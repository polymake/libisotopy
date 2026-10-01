#define CATCH_CONFIG_MAIN
#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <set>
#include <string>
#include <vector>

#include "catch.hpp"
#include "isotopy_batch.h"
#include "isotopy_graph.h"
#include "isotopy_viro.h"
#include "node.hpp"
#include "test_env.h"

namespace {

using IsotopyGPU::Result;
using IsotopyGPU::Tables;
using IsotopyGPU::TypeEntry;
using IsotopyGPU::TypeTable;

struct Packed {
  std::vector<uint8_t> edges;
  std::vector<uint64_t> signs;
};

Packed pack_batch(const Tables& t, const std::vector<std::vector<bool>>& signs,
                  const std::vector<std::vector<Isotopy::Edge>>& edges) {
  Packed out;
  const size_t n = signs.size();
  out.edges.assign(n * 2 * t.nedges, 0);
  out.signs.assign(n * IsotopyGPU::kQ1Words, 0ull);
  for (size_t i = 0; i < n; ++i) {
    for (int v = 0; v < (int)signs[i].size(); ++v) {
      if (signs[i][v]) out.signs[i * IsotopyGPU::kQ1Words + (v >> 6)] |= 1ull << (v & 63);
    }
    for (int e = 0; e < t.nedges; ++e) {
      out.edges[i * 2 * t.nedges + 2 * e] = (uint8_t)edges[i][e].first;
      out.edges[i * 2 * t.nedges + 2 * e + 1] = (uint8_t)edges[i][e].second;
    }
  }
  return out;
}

struct Instances {
  std::vector<std::vector<bool>> signs;
  std::vector<std::vector<Isotopy::Edge>> edges;
  std::vector<std::vector<Isotopy::Triangle>> triangles;
  std::vector<std::string> viro;
};

Instances random_instances(int delta, int seeds, int trials) {
  Instances out;
  const int nverts = Isotopy::num_vertices(delta);
  for (int seed = 0; seed < seeds; ++seed) {
    std::srand((unsigned)(seed * 15485863 + delta));
    const std::vector<Isotopy::Edge> edges = Utils::get_random_triangulation(delta);
    const std::vector<Isotopy::Triangle> triangles = Isotopy::edges_to_triangles(edges, delta);
    if ((int)triangles.size() != delta * delta) continue;
    for (int trial = 0; trial < trials; ++trial) {
      std::vector<bool> signs(nverts);
      for (int v = 0; v < nverts; ++v) signs[v] = (std::rand() & 1) != 0;
      Isotopy::Graph g(delta, signs, triangles);
      g.isotopy_type();
      out.signs.push_back(signs);
      out.edges.push_back(edges);
      out.triangles.push_back(triangles);
      out.viro.push_back(g.viro_notation());
    }
  }
  return out;
}

TypeTable table_from_instances(int delta, const Instances& in,
                               const std::string& omit = std::string()) {
  std::set<std::string> distinct(in.viro.begin(), in.viro.end());
  std::vector<TypeEntry> raw;
  for (const std::string& v : distinct) {
    if (!omit.empty() && v == omit) continue;
    const IsotopyGPU::Viro::Parsed p = IsotopyGPU::Viro::parse(v);
    TypeEntry e;
    e.viro = v;
    IsotopyGPU::Viro::region_pn(p.tree, p.has_j, &e.p, &e.n);
    e.num_ovals = IsotopyGPU::Viro::num_ovals(p.tree, p.has_j);
    raw.push_back(e);
  }
  return IsotopyGPU::build_type_table(delta, raw);
}

const Tables& tables_for(int delta) {
  static std::vector<Tables> cache(IsotopyGPU::kMaxDelta + 1);
  static std::vector<char> ready(IsotopyGPU::kMaxDelta + 1, 0);
  if (!ready[delta]) {
    cache[delta] = IsotopyGPU::build_tables(delta);
    ready[delta] = 1;
  }
  return cache[delta];
}

#ifdef ISOTOPY_GPU_HOST_WARP
uint16_t lookup_threaded(const std::vector<uint64_t>& keys, uint64_t code) {
  IsotopyGPU::HostWarpContext ctx;
  uint16_t lane_result[32] = {0};
  std::vector<std::thread> lanes;
  lanes.reserve(32);
  for (int i = 0; i < 32; ++i) {
    lanes.emplace_back([&, i]() {
      lane_result[i] = IsotopyGPU::lookup_type(keys.data(), (int)keys.size(), code,
                                               IsotopyGPU::Warp<32>{&ctx, i});
    });
  }
  for (auto& th : lanes) th.join();
  for (int i = 1; i < 32; ++i) {
    if (lane_result[i] != lane_result[0]) return IsotopyGPU::kTypeFallback;
  }
  return lane_result[0];
}
#endif

std::string json_string_field(const std::string& text, const std::string& key) {
  const std::string needle = "\"" + key + "\"";
  size_t at = text.find(needle);
  if (at == std::string::npos) return std::string();
  at = text.find('"', at + needle.size());
  if (at == std::string::npos) return std::string();
  const size_t end = text.find('"', at + 1);
  if (end == std::string::npos) return std::string();
  return text.substr(at + 1, end - at - 1);
}

}  // namespace

TEST_CASE("the Viro grammar parses every shape libisotopy emits", "[viro][grammar]") {
  struct Shape {
    const char* text;
    bool has_j;
    int ovals;
  };
  const Shape shapes[] = {
      {"<>", false, 0},        {"<J>", true, 1},          {"<1>", false, 1},
      {"<10>", false, 10},     {"<Jv1>", true, 2},        {"<Jv10v1<1>>", true, 13},
      {"<1<2>>", false, 3},    {"<9v1<1>v1<10>>", false, 22},
      {"<10v1<1>v2<2>>", false, 18}, {"<2v1<2v1<1>>>", false, 7},
  };

  for (const Shape& s : shapes) {
    INFO("shape " << s.text);
    const IsotopyGPU::Viro::Parsed p = IsotopyGPU::Viro::parse(s.text);
    CHECK(p.has_j == s.has_j);
    CHECK(IsotopyGPU::Viro::num_ovals(p.tree, p.has_j) == s.ovals);
    CHECK(IsotopyGPU::Viro::to_viro(p.tree, p.has_j, false) == s.text);
  }

  const IsotopyGPU::Viro::Parsed empty_set = IsotopyGPU::Viro::parse("∅");
  CHECK(empty_set.has_j == false);
  CHECK(IsotopyGPU::Viro::to_viro(empty_set.tree, false, false) == "<>");

  const std::string unicode = "⟨J⊔10⊔1⟨1⟩⟩";
  const IsotopyGPU::Viro::Parsed uni = IsotopyGPU::Viro::parse(unicode);
  CHECK(uni.has_j);
  CHECK(IsotopyGPU::Viro::to_viro(uni.tree, uni.has_j, true) == unicode);
  CHECK(IsotopyGPU::Viro::to_viro(uni.tree, uni.has_j, false) == "<Jv10v1<1>>");
  CHECK(IsotopyGPU::Viro::tree_code(uni.tree) ==
        IsotopyGPU::Viro::tree_code(IsotopyGPU::Viro::parse("<Jv10v1<1>>").tree));
}

TEST_CASE("malformed Viro notation is rejected", "[viro][grammar]") {
  const char* bogus[] = {"",      "<",     ">",      "<1",     "1>",      "<v1>",  "<1J>",
                         "<>x",   "<1<2>", "<<1>>",  "<1v>",   "<1v2>",   "<J",    "J<>",
                         "<-1>",  "<1>>",  "∅x"};
  for (const char* b : bogus) {
    INFO("input \"" << b << "\"");
    CHECK_THROWS_AS(IsotopyGPU::Viro::parse(b), std::invalid_argument);
  }
}

TEST_CASE("the tree code is blind to how children are ordered", "[viro][code]") {
  CHECK(IsotopyGPU::Viro::tree_code(IsotopyGPU::Viro::parse("<1<1>v1<2>>").tree) ==
        IsotopyGPU::Viro::tree_code(IsotopyGPU::Viro::parse("<1<2>v1<1>>").tree));
  CHECK(IsotopyGPU::Viro::tree_code(IsotopyGPU::Viro::parse("<2<1>>").tree) ==
        IsotopyGPU::Viro::tree_code(IsotopyGPU::Viro::parse("<1<1>v1<1>>").tree));
  CHECK(IsotopyGPU::Viro::tree_code(IsotopyGPU::Viro::parse("<1>").tree) !=
        IsotopyGPU::Viro::tree_code(IsotopyGPU::Viro::parse("<2>").tree));
}

TEST_CASE("the parsed corpus notation reproduces the kernel's tree code", "[viro][corpus]") {
  std::ifstream ifs("tests/isotopy_tests.yaml");
  if (!ifs.is_open()) {
    ISOTOPY_SKIP_UNLESS_REQUIRED("ISOTOPY_REQUIRE_CORPUS",
                                 "Could not open tests/isotopy_tests.yaml. "
                                 "Skipping corpus comparison.");
    return;
  }
  const std::string yaml((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
  auto doc = fkyaml::node::deserialize(yaml);

  int checked = 0;
  int skipped = 0;
  int bad_code = 0;
  int bad_roundtrip = 0;
  int bad_pn = 0;
  std::set<std::string> distinct;

  for (size_t i = 0; i < doc.size(); ++i) {
    auto& node = doc[i];
    const int delta = node.at("degree").get_value<int>();
    if (delta < 1 || delta > IsotopyGPU::kMaxDelta || !node.contains("viro")) {
      ++skipped;
      continue;
    }
    const std::string viro = node.at("viro").get_value<std::string>();
    const std::vector<bool> signs = node.at("polarisation").get_value<std::vector<bool>>();
    const std::vector<Isotopy::Triangle> triangles =
        node.at("triangulation").get_value<std::vector<Isotopy::Triangle>>();
    if (viro.empty() || (int)signs.size() != Isotopy::num_vertices(delta) ||
        (int)triangles.size() != delta * delta) {
      ++skipped;
      continue;
    }

    const Tables& t = tables_for(delta);
    const std::vector<Isotopy::Edge> edges = Isotopy::triangles_to_edges(triangles);
    if ((int)edges.size() != t.nedges) {
      ++skipped;
      continue;
    }

    Isotopy::Graph g(delta, signs, triangles);
    g.isotopy_type();
    const uint64_t want = IsotopyGPU::tree_code_from_adjacency(g.root_region, g.region_adjacency);

    const IsotopyGPU::Viro::Parsed parsed = IsotopyGPU::Viro::parse(viro);
    const uint64_t got = IsotopyGPU::Viro::tree_code(parsed.tree);

    const Packed packed = pack_batch(t, {signs}, {edges});
    const Result device = IsotopyGPU::classify_host(t, packed.edges.data(), packed.signs.data());

    if (got != want || (device.flags == IsotopyGPU::kOk && device.code != got)) {
      if (bad_code < 5) {
        INFO("case " << i << " degree " << delta << " viro " << viro << ": parsed " << got
                     << " reference " << want << " kernel " << device.code);
        CHECK(false);
      }
      ++bad_code;
    }
    if (IsotopyGPU::Viro::to_viro(parsed.tree, parsed.has_j, false) != viro) {
      if (bad_roundtrip < 5) {
        INFO("case " << i << " viro " << viro << " rendered "
                     << IsotopyGPU::Viro::to_viro(parsed.tree, parsed.has_j, false));
        CHECK(false);
      }
      ++bad_roundtrip;
    }
    int p = 0;
    int n = 0;
    IsotopyGPU::Viro::region_pn(parsed.tree, parsed.has_j, &p, &n);
    if (p != g.p_regions || n != g.n_regions) {
      if (bad_pn < 5) {
        INFO("case " << i << " viro " << viro << ": p/n " << p << "/" << n << " against "
                     << g.p_regions << "/" << g.n_regions);
        CHECK(false);
      }
      ++bad_pn;
    }
    distinct.insert(viro);
    ++checked;
  }

  WARN("corpus: checked " << checked << " (" << distinct.size() << " distinct notations), skipped "
                          << skipped << ", bad code " << bad_code << ", bad round-trip "
                          << bad_roundtrip << ", bad p/n " << bad_pn);
  CHECK(checked > 0);
  CHECK(bad_code == 0);
  CHECK(bad_roundtrip == 0);
  CHECK(bad_pn == 0);
}

TEST_CASE("lookup_type agrees with the host table search", "[viro][lookup]") {
  for (int delta = 5; delta <= IsotopyGPU::kMaxDelta; ++delta) {
    const Instances in = random_instances(delta, 12, 12);
    REQUIRE(!in.viro.empty());
    const TypeTable table = table_from_instances(delta, in);
    REQUIRE(!table.keys.empty());

    for (size_t k = 0; k < table.keys.size(); ++k) {
      INFO("degree " << delta << " key " << k << " of " << table.keys.size());
      REQUIRE(IsotopyGPU::lookup_type(table.keys.data(), (int)table.keys.size(), table.keys[k],
                                      IsotopyGPU::Warp<1>{}) == k);
    }

    const uint64_t absent[] = {0ull, table.keys.front() - 1, table.keys.back() + 1, ~0ull};
    for (uint64_t code : absent) {
      if (std::binary_search(table.keys.begin(), table.keys.end(), code)) continue;
      INFO("degree " << delta << " absent code " << code);
      REQUIRE(IsotopyGPU::lookup_type(table.keys.data(), (int)table.keys.size(), code,
                                      IsotopyGPU::Warp<1>{}) == IsotopyGPU::kTypeUnknown);
      REQUIRE(table.lookup(code) == IsotopyGPU::kTypeUnknown);
    }

    for (size_t k = 0; k + 1 < table.keys.size(); ++k) {
      const uint64_t between = table.keys[k] + 1;
      if (between >= table.keys[k + 1]) continue;
      INFO("degree " << delta << " gap after key " << k);
      REQUIRE(IsotopyGPU::lookup_type(table.keys.data(), (int)table.keys.size(), between,
                                      IsotopyGPU::Warp<1>{}) == IsotopyGPU::kTypeUnknown);
      REQUIRE(table.lookup(between) == IsotopyGPU::kTypeUnknown);
    }
  }
}

#ifdef ISOTOPY_GPU_HOST_WARP
TEST_CASE("the 32-lane table search agrees with the single-lane one", "[viro][threaded]") {
  for (int delta = 6; delta <= IsotopyGPU::kMaxDelta; ++delta) {
    const Instances in = random_instances(delta, 10, 10);
    REQUIRE(!in.viro.empty());
    const TypeTable table = table_from_instances(delta, in);
    REQUIRE(table.keys.size() > 1);

    for (size_t k = 0; k < table.keys.size(); ++k) {
      INFO("degree " << delta << " key " << k << " of " << table.keys.size());
      REQUIRE(lookup_threaded(table.keys, table.keys[k]) == k);
    }

    const uint64_t absent[] = {0ull, table.keys.front() - 1, table.keys.back() + 1, ~0ull};
    for (uint64_t code : absent) {
      if (std::binary_search(table.keys.begin(), table.keys.end(), code)) continue;
      INFO("degree " << delta << " absent code " << code);
      REQUIRE(lookup_threaded(table.keys, code) == IsotopyGPU::kTypeUnknown);
    }

    for (size_t k = 0; k + 1 < table.keys.size(); ++k) {
      const uint64_t between = table.keys[k] + 1;
      if (between >= table.keys[k + 1]) continue;
      INFO("degree " << delta << " gap after key " << k);
      REQUIRE(lookup_threaded(table.keys, between) == IsotopyGPU::kTypeUnknown);
    }
  }
}
#endif

TEST_CASE("every clean instance resolves to its own type", "[viro][batch]") {
  for (int delta = 4; delta <= IsotopyGPU::kMaxDelta; ++delta) {
    const Instances in = random_instances(delta, 10, 10);
    REQUIRE(!in.viro.empty());
    const Tables& t = tables_for(delta);
    const TypeTable table = table_from_instances(delta, in);
    const Packed packed = pack_batch(t, in.signs, in.edges);

    IsotopyGPU::BatchOutput out;
    IsotopyGPU::classify_batch_host(t, table, packed.edges.data(), packed.signs.data(),
                                    (int)in.signs.size(), &out);

    INFO("degree " << delta);
    CHECK(out.unknown == 0);
    CHECK(out.recompute_failed == 0);
    CHECK(out.unknown_overflow == false);
    CHECK(out.device_disagreements == 0);
    CHECK(out.table_fingerprint == table.fingerprint);
    for (size_t i = 0; i < in.signs.size(); ++i) {
      REQUIRE(out.type_id[i] != IsotopyGPU::kTypeUnknown);
      REQUIRE(out.type_id[i] != IsotopyGPU::kTypeFallback);
      const TypeEntry* e = table.entry(out.type_id[i]);
      REQUIRE(e != nullptr);
      REQUIRE(e->viro == in.viro[i]);
      REQUIRE(e->p == out.results[i].p);
      REQUIRE(e->n == out.results[i].n);
    }

    IsotopyGPU::BatchOutput reused = out;
    IsotopyGPU::classify_batch_host(t, table, packed.edges.data(), packed.signs.data(),
                                    (int)in.signs.size(), &reused);
    CHECK(reused.flagged == out.flagged);
    CHECK(reused.recomputed == out.recomputed);
    CHECK(reused.unknown == out.unknown);
    CHECK(reused.unknown_codes.size() == out.unknown_codes.size());
    CHECK(reused.type_id == out.type_id);
  }
}

TEST_CASE("a device id that contradicts the host table is counted", "[viro][batch]") {
  const int delta = 6;
  const Instances in = random_instances(delta, 6, 6);
  REQUIRE(!in.viro.empty());
  const int batch = (int)in.signs.size();
  const Tables& t = tables_for(delta);
  const TypeTable table = table_from_instances(delta, in);
  const Packed packed = pack_batch(t, in.signs, in.edges);

  IsotopyGPU::BatchOutput out;
  IsotopyGPU::classify_batch_host(t, table, packed.edges.data(), packed.signs.data(), batch, &out);
  REQUIRE(out.device_disagreements == 0);

  int clean = -1;
  for (int i = 0; i < batch; ++i) {
    if (out.results[i].flags == IsotopyGPU::kOk) {
      clean = i;
      break;
    }
  }
  REQUIRE(clean >= 0);

  const uint16_t honest = out.type_id[clean];
  REQUIRE(honest != IsotopyGPU::kTypeUnknown);
  REQUIRE(honest != IsotopyGPU::kTypeFallback);
  out.results[clean].type_id = static_cast<uint16_t>(honest == 0 ? 1 : honest - 1);

  IsotopyGPU::resolve_batch(t, table, packed.edges.data(), packed.signs.data(), batch, &out, 4096);

  CHECK(out.device_disagreements == 1);
  CHECK(out.type_id[clean] == honest);
  CHECK(out.results[clean].type_id == honest);
}

TEST_CASE("the CPU fallback is transparent when the kernel gives up", "[viro][fallback]") {
  struct Squeeze {
    const char* name;
    int max_rounds;
    int max_regions;
    int code_bits;
  };
  const Squeeze squeezes[] = {
      {"rounds", 1, IsotopyGPU::kMaxRegions, IsotopyGPU::kCodeBits},
      {"regions", IsotopyGPU::kMaxRounds, 3, IsotopyGPU::kCodeBits},
      {"codebits", IsotopyGPU::kMaxRounds, IsotopyGPU::kMaxRegions, 6},
  };

  for (int delta = 6; delta <= IsotopyGPU::kMaxDelta; ++delta) {
    const Instances in = random_instances(delta, 8, 8);
    REQUIRE(!in.viro.empty());
    const int batch = (int)in.signs.size();
    const Tables& clean = tables_for(delta);
    const TypeTable table = table_from_instances(delta, in);
    const Packed packed = pack_batch(clean, in.signs, in.edges);

    IsotopyGPU::BatchOutput reference;
    IsotopyGPU::classify_batch_host(clean, table, packed.edges.data(), packed.signs.data(), batch,
                                    &reference);

    for (const Squeeze& sq : squeezes) {
      Tables tight = clean;
      tight.max_rounds = sq.max_rounds;
      tight.max_regions = sq.max_regions;
      tight.code_bits = sq.code_bits;

      IsotopyGPU::BatchOutput out;
      IsotopyGPU::classify_batch_host(tight, table, packed.edges.data(), packed.signs.data(), batch,
                                      &out);

      INFO("degree " << delta << " squeeze " << sq.name << ": flagged " << out.flagged << " of "
                     << batch << ", recomputed " << out.recomputed);
      CHECK(out.flagged > 0);
      CHECK(out.recompute_failed == 0);
      for (int i = 0; i < batch; ++i) {
        REQUIRE(out.type_id[i] == reference.type_id[i]);
        REQUIRE(out.results[i].p == reference.results[i].p);
        REQUIRE(out.results[i].n == reference.results[i].n);
        REQUIRE(out.results[i].code == reference.results[i].code);
      }
    }
  }
}

TEST_CASE("a type missing from the table comes back as unknown with a witness", "[viro][unknown]") {
  const int delta = 8;
  const Instances in = random_instances(delta, 10, 10);
  REQUIRE(!in.viro.empty());
  const int batch = (int)in.signs.size();
  const Tables& t = tables_for(delta);
  const Packed packed = pack_batch(t, in.signs, in.edges);

  const std::string omitted = in.viro.front();
  const TypeTable full = table_from_instances(delta, in);
  const TypeTable holed = table_from_instances(delta, in, omitted);
  REQUIRE(holed.entries.size() + 1 == full.entries.size());

  const uint64_t omitted_code = IsotopyGPU::Viro::tree_code(IsotopyGPU::Viro::parse(omitted).tree);

  IsotopyGPU::BatchOutput out;
  IsotopyGPU::classify_batch_host(t, holed, packed.edges.data(), packed.signs.data(), batch, &out);

  int expected_unknown = 0;
  for (const std::string& v : in.viro) {
    if (v == omitted) ++expected_unknown;
  }
  REQUIRE(expected_unknown > 0);

  CHECK(out.unknown == expected_unknown);
  CHECK(out.unknown_codes.size() == 1);
  CHECK(out.unknown_codes.front() == omitted_code);
  CHECK(out.unknown_overflow == false);
  REQUIRE(out.unknown_witness.size() == 1);
  CHECK(in.viro[out.unknown_witness.front()] == omitted);
  for (int i = 0; i < batch; ++i) {
    REQUIRE((out.type_id[i] == IsotopyGPU::kTypeUnknown) == (in.viro[i] == omitted));
  }

  IsotopyGPU::BatchOutput empty_table_run;
  const TypeTable empty = IsotopyGPU::build_type_table(delta, std::vector<TypeEntry>());
  IsotopyGPU::classify_batch_host(t, empty, packed.edges.data(), packed.signs.data(), batch,
                                  &empty_table_run, 1);
  CHECK(empty_table_run.unknown == batch);
  CHECK(empty_table_run.unknown_codes.size() == 1);
  CHECK(empty_table_run.unknown_overflow == true);
}

TEST_CASE("two spellings of one type merge into a primary and an alias", "[viro][table]") {
  const std::vector<TypeEntry> raw = {TypeEntry{0ull, "∅", {}, 0, 0, 0},
                                      TypeEntry{0ull, "<>", {}, 0, 0, 0}};
  const TypeTable table = IsotopyGPU::build_type_table(8, raw);

  REQUIRE(table.rejected.empty());
  REQUIRE(table.entries.size() == 1);
  CHECK(table.entries[0].viro == "<>");
  CHECK(table.entries[0].aliases == std::vector<std::string>{"∅"});
  CHECK(table.keys.size() == 1);
  CHECK(table.lookup(table.entries[0].code) == 0);
}

TEST_CASE("wrong-parity and unparseable entries land in rejected", "[viro][table]") {
  const std::vector<TypeEntry> raw = {TypeEntry{0ull, "<J>", {}, 0, 0, 0},
                                      TypeEntry{0ull, "<1", {}, 0, 0, 0},
                                      TypeEntry{0ull, "<1>", {}, 0, 0, 0}};
  const TypeTable table = IsotopyGPU::build_type_table(8, raw);

  CHECK(table.rejected == std::vector<std::string>{"<J>", "<1"});
  REQUIRE(table.entries.size() == 1);
  CHECK(table.entries[0].viro == "<1>");

  const std::vector<TypeEntry> contradictory = {TypeEntry{0ull, "<1>", {}, 1, 0, 1},
                                                TypeEntry{0ull, "<1>", {}, 7, 0, 1}};
  CHECK_THROWS_AS(IsotopyGPU::build_type_table(8, contradictory), std::runtime_error);
}

TEST_CASE("the fingerprint moves when the table moves", "[viro][table]") {
  const std::vector<TypeEntry> raw = {TypeEntry{0ull, "<>", {}, 0, 0, 0},
                                      TypeEntry{0ull, "<1>", {}, 0, 0, 1},
                                      TypeEntry{0ull, "<2>", {}, 0, 0, 2}};

  const uint64_t full = IsotopyGPU::build_type_table(8, raw).fingerprint;
  const uint64_t same = IsotopyGPU::build_type_table(8, raw).fingerprint;
  const uint64_t shorter =
      IsotopyGPU::build_type_table(8, {raw.begin(), raw.end() - 1}).fingerprint;
  const uint64_t other_delta = IsotopyGPU::build_type_table(6, raw).fingerprint;

  CHECK(full == same);
  CHECK(full != shorter);
  CHECK(full != other_delta);
  CHECK(shorter != other_delta);
}

TEST_CASE("the published type catalog parses and agrees on p, n and ovals", "[viro][catalog]") {
  const std::string path = env_or_empty("ISOTOPY_TYPE_INFO");
  if (path.empty()) {
    ISOTOPY_SKIP_UNLESS_REQUIRED("ISOTOPY_REQUIRE_CATALOG",
                                 "ISOTOPY_TYPE_INFO is unset. "
                                 "Skipping the type_info.yaml cross-check.");
    return;
  }

  const std::map<int, std::vector<TypeEntry>> raw = IsotopyGPU::read_type_info(path);
  REQUIRE(!raw.empty());

  int total = 0;
  int aliases = 0;
  for (const auto& kv : raw) {
    const int delta = kv.first;
    const TypeTable table = IsotopyGPU::build_type_table(delta, kv.second);

    INFO("degree " << delta);
    CHECK(table.rejected.empty());

    for (const TypeEntry& e : table.entries) {
      const IsotopyGPU::Viro::Parsed p = IsotopyGPU::Viro::parse(e.viro);
      INFO("degree " << delta << " type " << e.viro);
      REQUIRE(p.has_j == ((delta % 2) != 0));
      REQUIRE(IsotopyGPU::Viro::tree_code(p.tree) == e.code);
      REQUIRE(IsotopyGPU::Viro::to_viro(p.tree, p.has_j, false) == e.viro);
      int pp = 0;
      int nn = 0;
      IsotopyGPU::Viro::region_pn(p.tree, p.has_j, &pp, &nn);
      REQUIRE(pp == e.p);
      REQUIRE(nn == e.n);
      REQUIRE(IsotopyGPU::Viro::num_ovals(p.tree, p.has_j) == e.num_ovals);
      REQUIRE(table.lookup(e.code) < table.entries.size());
      aliases += (int)e.aliases.size();
      ++total;
    }
  }
  WARN("catalog: " << total << " distinct types over " << raw.size() << " degrees, " << aliases
                   << " alias spellings merged");
  CHECK(total > 0);
}

TEST_CASE("realizing patchworks classify to their own catalog type", "[viro][witness]") {
  const std::string type_info = env_or_empty("ISOTOPY_TYPE_INFO");
  const std::string root = env_or_empty("ISOTOPY_REALIZABLE");
  if (type_info.empty() || root.empty()) {
    ISOTOPY_SKIP_UNLESS_REQUIRED("ISOTOPY_REQUIRE_CATALOG",
                                 "ISOTOPY_TYPE_INFO or ISOTOPY_REALIZABLE is unset. "
                                 "Skipping the witness cross-check.");
    return;
  }
  if (!std::filesystem::exists(root)) {
    ISOTOPY_SKIP_UNLESS_REQUIRED("ISOTOPY_REQUIRE_CATALOG",
                                 "ISOTOPY_REALIZABLE does not exist. "
                                 "Skipping the witness cross-check.");
    return;
  }

  const std::map<int, TypeTable> tables = IsotopyGPU::load_type_tables(type_info);

  int checked = 0;
  int skipped = 0;
  int mismatched = 0;
  for (const auto& entry : std::filesystem::recursive_directory_iterator(root)) {
    if (!entry.is_regular_file() || entry.path().extension() != ".pcom") continue;

    std::ifstream ifs(entry.path());
    const std::string text((std::istreambuf_iterator<char>(ifs)),
                           std::istreambuf_iterator<char>());
    const std::string want = json_string_field(text, "TYPE");
    if (want.empty()) {
      ++skipped;
      continue;
    }

    std::vector<bool> signs;
    std::vector<Isotopy::Triangle> triangles;
    try {
      const auto parsed = Utils::pcom_to_signs_and_triangles_vec(text);
      signs = parsed.first;
      triangles = parsed.second;
    } catch (const std::exception&) {
      ++skipped;
      continue;
    }

    int delta = 0;
    for (int d = 1; d <= IsotopyGPU::kMaxDelta; ++d) {
      if (Isotopy::num_vertices(d) == (int)signs.size() && (int)triangles.size() == d * d) delta = d;
    }
    if (delta == 0 || tables.find(delta) == tables.end()) {
      ++skipped;
      continue;
    }

    const Tables& t = tables_for(delta);
    const std::vector<Isotopy::Edge> edges = Isotopy::triangles_to_edges(triangles);
    if ((int)edges.size() != t.nedges) {
      ++skipped;
      continue;
    }

    const TypeTable& table = tables.at(delta);
    const Packed packed = pack_batch(t, {signs}, {edges});
    IsotopyGPU::BatchOutput out;
    IsotopyGPU::classify_batch_host(t, table, packed.edges.data(), packed.signs.data(), 1, &out);

    const TypeEntry* got = table.entry(out.type_id[0]);
    const bool alias_hit =
        got != nullptr &&
        std::find(got->aliases.begin(), got->aliases.end(), want) != got->aliases.end();
    if (got == nullptr || (got->viro != want && !alias_hit)) {
      if (mismatched < 5) {
        INFO(entry.path().string() << ": got "
                                   << (got == nullptr ? std::string("<unresolved>") : got->viro)
                                   << " want " << want);
        CHECK(false);
      }
      ++mismatched;
    }
    ++checked;
  }

  WARN("witnesses: checked " << checked << ", skipped " << skipped << ", mismatched "
                             << mismatched);
  CHECK(checked > 0);
  CHECK(mismatched == 0);
}
