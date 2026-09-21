#define CATCH_CONFIG_MAIN
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <random>
#include <string>
#include <thread>
#include <vector>

#include "catch.hpp"
#include "isotopy_gpu.cuh"
#include "isotopy_graph.h"
#include "node.hpp"

namespace {

struct Packed {
  std::vector<uint8_t> edges;
  std::vector<uint64_t> signs;
};

Packed pack_instance(const IsotopyGPU::Tables& t, const std::vector<bool>& signs,
                     const std::vector<Isotopy::Edge>& edges) {
  Packed out;
  out.signs.assign(IsotopyGPU::kQ1Words, 0ull);
  for (int i = 0; i < (int)signs.size(); ++i) {
    if (signs[i]) out.signs[i >> 6] |= 1ull << (i & 63);
  }
  out.edges.resize(2 * t.nedges);
  for (int e = 0; e < t.nedges; ++e) {
    out.edges[2 * e] = (uint8_t)edges[e].first;
    out.edges[2 * e + 1] = (uint8_t)edges[e].second;
  }
  return out;
}

struct Expectation {
  uint64_t code;
  int p;
  int n;
};

struct Case {
  int delta;
  std::vector<bool> signs;
  std::vector<Isotopy::Triangle> triangles;
  std::vector<Isotopy::Edge> edges;
};

std::vector<Case> random_cases(int delta, int seeds, int trials) {
  std::vector<Case> out;
  const int nverts = Isotopy::num_vertices(delta);
  for (int seed = 0; seed < seeds; ++seed) {
    std::srand((unsigned)(seed * 15485863 + delta));
    Case base;
    base.delta = delta;
    base.edges = Utils::get_random_triangulation(delta);
    base.triangles = Isotopy::edges_to_triangles(base.edges, delta);
    if ((int)base.triangles.size() != delta * delta) continue;
    for (int trial = 0; trial < trials; ++trial) {
      Case c = base;
      c.signs.resize(nverts);
      for (int i = 0; i < nverts; ++i) c.signs[i] = (std::rand() & 1) != 0;
      out.push_back(std::move(c));
    }
  }
  return out;
}

Expectation reference(int delta, const std::vector<bool>& signs,
                      const std::vector<Isotopy::Triangle>& triangles) {
  Isotopy::Graph g(delta, signs, triangles);
  g.isotopy_type();
  return {IsotopyGPU::tree_code_from_adjacency(g.root_region, g.region_adjacency), g.p_regions,
          g.n_regions};
}

const IsotopyGPU::Tables& tables_for(int delta) {
  static std::vector<IsotopyGPU::Tables> cache(IsotopyGPU::kMaxDelta + 1);
  static std::vector<char> ready(IsotopyGPU::kMaxDelta + 1, 0);
  if (!ready[delta]) {
    cache[delta] = IsotopyGPU::build_tables(delta);
    ready[delta] = 1;
  }
  return cache[delta];
}

bool compare(int delta, const std::vector<bool>& signs,
             const std::vector<Isotopy::Triangle>& triangles, std::string* why) {
  const IsotopyGPU::Tables& t = tables_for(delta);
  const std::vector<Isotopy::Edge> edges = Isotopy::triangles_to_edges(triangles);
  if ((int)edges.size() != t.nedges) {
    *why = "edge count " + std::to_string(edges.size()) + " != " + std::to_string(t.nedges);
    return false;
  }

  const Packed packed = pack_instance(t, signs, edges);
  const IsotopyGPU::Result got =
      IsotopyGPU::classify_host(t, packed.edges.data(), packed.signs.data());
  if (got.flags != IsotopyGPU::kOk) {
    *why = "gpu flags 0x" + std::to_string(got.flags) + " rounds " + std::to_string(got.rounds);
    return false;
  }

  const Expectation want = reference(delta, signs, triangles);
  if (got.code != want.code) {
    *why = "code " + std::to_string(got.code) + " != " + std::to_string(want.code);
    return false;
  }
  if (got.p != want.p || got.n != want.n) {
    *why = "p/n " + std::to_string(got.p) + "/" + std::to_string(got.n) +
           " != " + std::to_string(want.p) + "/" + std::to_string(want.n);
    return false;
  }
  return true;
}

}  // namespace

TEST_CASE("compile-time widths stay inside the types that carry them", "[gpu][limits]") {
  REQUIRE(Isotopy::num_vertices(IsotopyGPU::kMaxDelta) <= 256);
  REQUIRE(IsotopyGPU::kMaxRegions <= 32);
  REQUIRE(IsotopyGPU::kCodeBits < 64);

  const int worst_regions = (IsotopyGPU::kMaxDelta - 1) * (IsotopyGPU::kMaxDelta - 2) / 2 + 1;
  INFO("Harnack bound at delta=" << IsotopyGPU::kMaxDelta << " is " << worst_regions
                                 << " regions, needing " << 2 * (worst_regions - 1)
                                 << " code bits");
  REQUIRE(worst_regions <= IsotopyGPU::kMaxRegions);
  REQUIRE(2 * (worst_regions - 1) < IsotopyGPU::kCodeBits);
}

TEST_CASE("the tree code is invariant under region relabeling", "[gpu][relabel]") {
  const std::vector<std::vector<std::vector<int>>> trees = {
      {{}},
      {{1}, {0}},
      {{1, 2}, {0}, {0}},
      {{1}, {0, 2}, {1}},
      {{1, 2, 3}, {0}, {0, 4}, {0}, {2}},
      {{1, 2}, {0, 3, 4}, {0}, {1}, {1, 5}, {4}},
  };

  for (const auto& adjacency : trees) {
    const int n = (int)adjacency.size();
    const uint64_t base = IsotopyGPU::tree_code_from_adjacency(0, adjacency);

    std::vector<int> perm(n);
    for (int i = 0; i < n; ++i) perm[i] = i;

    std::mt19937 rng(9001u + n);
    for (int trial = 0; trial < 50; ++trial) {
      std::shuffle(perm.begin(), perm.end(), rng);

      std::vector<std::vector<int>> relabeled(n);
      for (int u = 0; u < n; ++u) {
        for (int v : adjacency[u]) relabeled[perm[u]].push_back(perm[v]);
      }
      for (auto& row : relabeled) std::sort(row.begin(), row.end());

      INFO("n=" << n << " trial=" << trial);
      REQUIRE(IsotopyGPU::tree_code_from_adjacency(perm[0], relabeled) == base);
    }
  }
}

TEST_CASE("classify_host is reentrant across threads", "[gpu][reentrant]") {
  for (int delta = 5; delta <= 8; ++delta) {
    const IsotopyGPU::Tables& t = tables_for(delta);
    const std::vector<Case> cases = random_cases(delta, 8, 8);
    REQUIRE(!cases.empty());

    std::vector<Packed> packed;
    packed.reserve(cases.size());
    for (const Case& c : cases) packed.push_back(pack_instance(t, c.signs, c.edges));

    std::vector<IsotopyGPU::Result> serial(cases.size());
    for (size_t i = 0; i < cases.size(); ++i) {
      serial[i] = IsotopyGPU::classify_host(t, packed[i].edges.data(), packed[i].signs.data());
    }

    std::vector<IsotopyGPU::Result> concurrent(cases.size());
    std::vector<std::thread> workers;
    const int nthreads = 4;
    for (int w = 0; w < nthreads; ++w) {
      workers.emplace_back([&, w]() {
        for (size_t i = w; i < cases.size(); i += nthreads) {
          concurrent[i] =
              IsotopyGPU::classify_host(t, packed[i].edges.data(), packed[i].signs.data());
        }
      });
    }
    for (auto& th : workers) th.join();

    for (size_t i = 0; i < cases.size(); ++i) {
      INFO("degree " << delta << " case " << i);
      REQUIRE(concurrent[i].flags == serial[i].flags);
      REQUIRE(concurrent[i].code == serial[i].code);
      REQUIRE(concurrent[i].p == serial[i].p);
      REQUIRE(concurrent[i].n == serial[i].n);
      REQUIRE(concurrent[i].regions == serial[i].regions);
    }
  }
}

TEST_CASE("squeezed limits flag instead of answering wrongly", "[gpu][fallback]") {
  struct Squeeze {
    const char* name;
    int max_rounds;
    int max_regions;
    int code_bits;
    uint32_t expect;
  };
  const Squeeze squeezes[] = {
      {"rounds", 1, IsotopyGPU::kMaxRegions, IsotopyGPU::kCodeBits, IsotopyGPU::kNoConverge},
      {"regions", IsotopyGPU::kMaxRounds, 3, IsotopyGPU::kCodeBits, IsotopyGPU::kTooManyRegions},
      {"codebits", IsotopyGPU::kMaxRounds, IsotopyGPU::kMaxRegions, 6,
       IsotopyGPU::kCodeOverflow},
  };

  for (const Squeeze& sq : squeezes) {
    int observed = 0;
    int wrong = 0;
    int total = 0;

    for (int delta = 6; delta <= 8; ++delta) {
      IsotopyGPU::Tables t = IsotopyGPU::build_tables(delta);
      t.max_rounds = sq.max_rounds;
      t.max_regions = sq.max_regions;
      t.code_bits = sq.code_bits;

      for (const Case& c : random_cases(delta, 6, 6)) {
        const Packed packed = pack_instance(t, c.signs, c.edges);
        const IsotopyGPU::Result got =
            IsotopyGPU::classify_host(t, packed.edges.data(), packed.signs.data());
        ++total;

        if (got.flags & sq.expect) ++observed;
        if (got.flags != IsotopyGPU::kOk) continue;

        const Expectation want = reference(delta, c.signs, c.triangles);
        if (got.code != want.code || got.p != want.p || got.n != want.n) ++wrong;
      }
    }

    WARN("squeeze " << sq.name << ": " << observed << " of " << total << " hit the target flag, "
                    << wrong << " wrong-and-unflagged");
    CHECK(total > 0);
    CHECK(observed > 0);
    CHECK(wrong == 0);
  }
}

TEST_CASE("sphere tables agree with libisotopy quadrant conventions", "[gpu][tables]") {
  for (int delta = 1; delta <= IsotopyGPU::kMaxDelta; ++delta) {
    const IsotopyGPU::Tables& t = IsotopyGPU::build_tables(delta);

    REQUIRE(t.nverts == Isotopy::num_vertices(delta));
    REQUIRE(t.nedges == Isotopy::num_edges(delta));
    REQUIRE(t.nsphere == 4 * delta * delta + 2);

    for (int v = 0; v < t.nsphere; ++v) {
      REQUIRE(t.anti[t.anti[v]] == v);
      REQUIRE(t.anti[v] != v);
    }

    for (int x = 0; x <= delta; ++x) {
      for (int y = 0; y + x <= delta; ++y) {
        const int v = IsotopyGPU::q1_index(delta, x, y);
        REQUIRE(t.flip[0 * t.nverts + v] == 0);
        REQUIRE(t.flip[1 * t.nverts + v] == (uint8_t)(x & 1));
        REQUIRE(t.flip[2 * t.nverts + v] == (uint8_t)(y & 1));
        REQUIRE(t.flip[3 * t.nverts + v] == (uint8_t)((x & 1) ^ (y & 1)));

        const int q1 = t.oct[0 * t.nverts + v];
        REQUIRE(t.oct[1 * t.nverts + v] == (x == 0 ? q1 : t.oct[1 * t.nverts + v]));
        REQUIRE(t.oct[2 * t.nverts + v] == (y == 0 ? q1 : t.oct[2 * t.nverts + v]));
      }
    }
  }
}

TEST_CASE("antipodal map flips the sign by (-1)^delta", "[gpu][tables]") {
  for (int delta = 1; delta <= IsotopyGPU::kMaxDelta; ++delta) {
    const IsotopyGPU::Tables& t = IsotopyGPU::build_tables(delta);
    for (int x = 0; x <= delta; ++x) {
      for (int y = 0; y + x <= delta; ++y) {
        const int z = delta - x - y;
        const int v = IsotopyGPU::q1_index(delta, x, y);
        for (int o = 0; o < 8; ++o) {
          const int opp = o ^ 7;
          if (x == 0 || y == 0 || z == 0) continue;
          REQUIRE(t.oct[opp * t.nverts + v] == t.anti[t.oct[o * t.nverts + v]]);
          const int d = t.flip[o * t.nverts + v] ^ t.flip[opp * t.nverts + v];
          REQUIRE(d == (delta & 1));
        }
      }
    }
  }
}

TEST_CASE("GPU pipeline matches libisotopy on the YAML corpus", "[gpu][corpus]") {
  std::ifstream ifs("tests/isotopy_tests.yaml");
  if (!ifs.is_open()) {
    WARN("Could not open tests/isotopy_tests.yaml. Skipping corpus comparison.");
    return;
  }
  const std::string yaml((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
  auto doc = fkyaml::node::deserialize(yaml);

  int checked = 0;
  int skipped = 0;
  int failed = 0;
  for (size_t i = 0; i < doc.size(); ++i) {
    auto& node = doc[i];
    const int delta = node.at("degree").get_value<int>();
    if (delta < 1 || delta > IsotopyGPU::kMaxDelta) {
      ++skipped;
      continue;
    }
    const std::vector<bool> signs = node.at("polarisation").get_value<std::vector<bool>>();
    const std::vector<Isotopy::Triangle> triangles =
        node.at("triangulation").get_value<std::vector<Isotopy::Triangle>>();
    if ((int)signs.size() != Isotopy::num_vertices(delta) ||
        (int)triangles.size() != delta * delta) {
      ++skipped;
      continue;
    }

    std::string why;
    if (!compare(delta, signs, triangles, &why)) {
      if (failed < 10) {
        INFO("case " << i << " degree " << delta << ": " << why);
        CHECK(false);
      }
      ++failed;
    }
    ++checked;
  }

  WARN("corpus: checked " << checked << ", skipped " << skipped << ", failed " << failed);
  CHECK(failed == 0);
  CHECK(checked > 0);
}

#ifdef ISOTOPY_GPU_HOST_WARP
TEST_CASE("32 host lanes agree with the serial lane", "[gpu][threaded]") {
  for (int delta = 4; delta <= 8; ++delta) {
    const IsotopyGPU::Tables& t = tables_for(delta);
    const int nverts = Isotopy::num_vertices(delta);

    for (int seed = 0; seed < 4; ++seed) {
      std::srand((unsigned)(seed * 104729 + delta));
      const std::vector<Isotopy::Edge> edges = Utils::get_random_triangulation(delta);
      const std::vector<Isotopy::Triangle> triangles = Isotopy::edges_to_triangles(edges, delta);
      if ((int)triangles.size() != delta * delta) continue;

      for (int trial = 0; trial < 4; ++trial) {
        std::vector<bool> signs(nverts);
        for (int i = 0; i < nverts; ++i) signs[i] = (std::rand() & 1) != 0;

        const Packed packed = pack_instance(t, signs, edges);
        const IsotopyGPU::Result serial =
            IsotopyGPU::classify_host(t, packed.edges.data(), packed.signs.data());
        const IsotopyGPU::Result threaded =
            IsotopyGPU::classify_host_threaded(t, packed.edges.data(), packed.signs.data());

        INFO("degree " << delta << " seed " << seed << " trial " << trial);
        REQUIRE(threaded.flags == serial.flags);
        REQUIRE(threaded.code == serial.code);
        REQUIRE(threaded.p == serial.p);
        REQUIRE(threaded.n == serial.n);
        REQUIRE(threaded.regions == serial.regions);
        REQUIRE(threaded.components == serial.components);
      }
    }
  }
}
#endif

TEST_CASE("GPU pipeline matches libisotopy on random triangulations", "[gpu][random]") {
  for (int delta = 3; delta <= IsotopyGPU::kMaxDelta; ++delta) {
    const int nverts = Isotopy::num_vertices(delta);
    int checked = 0;
    int failed = 0;

    for (int seed = 0; seed < 40; ++seed) {
      std::srand((unsigned)(seed * 7919 + delta));
      const std::vector<Isotopy::Edge> edges = Utils::get_random_triangulation(delta);
      const std::vector<Isotopy::Triangle> triangles = Isotopy::edges_to_triangles(edges, delta);
      if ((int)triangles.size() != delta * delta) continue;

      for (int trial = 0; trial < 8; ++trial) {
        std::vector<bool> signs(nverts);
        for (int i = 0; i < nverts; ++i) signs[i] = (std::rand() & 1) != 0;

        std::string why;
        if (!compare(delta, signs, triangles, &why)) {
          if (failed < 5) {
            INFO("degree " << delta << " seed " << seed << " trial " << trial << ": " << why);
            CHECK(false);
          }
          ++failed;
        }
        ++checked;
      }
    }

    WARN("degree " << delta << ": checked " << checked << ", failed " << failed);
    CHECK(failed == 0);
    CHECK(checked > 0);
  }
}
