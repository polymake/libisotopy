#pragma once

#include <cstdint>
#include <unordered_map>
#include <vector>

#include "isotopy_graph.h"
#include "isotopy_types.h"

namespace IsotopyGPU {

struct BatchOutput {
  std::vector<Result> results;
  std::vector<uint16_t> type_id;
  std::vector<uint64_t> unknown_codes;
  std::vector<int> unknown_witness;
  bool unknown_overflow = false;
  uint64_t table_fingerprint = 0ull;
  int flagged = 0;
  int recomputed = 0;
  int recompute_failed = 0;
  int unknown = 0;
  int device_disagreements = 0;
};

inline std::vector<bool> unpack_signs(const Tables& t, const uint64_t* q1signs) {
  std::vector<bool> signs(t.nverts);
  for (int v = 0; v < t.nverts; ++v) signs[v] = get_bit(q1signs, v);
  return signs;
}

inline std::vector<Isotopy::Edge> unpack_edges(const Tables& t, const uint8_t* edges) {
  std::vector<Isotopy::Edge> out;
  out.reserve(t.nedges);
  for (int e = 0; e < t.nedges; ++e) {
    out.push_back(Isotopy::Edge(edges[2 * e], edges[2 * e + 1]));
  }
  return out;
}

inline uint16_t resolve_one_cpu(const Tables& t, const TypeTable& table, const uint8_t* edges,
                                const uint64_t* q1signs, Result* r) {
  try {
    const std::vector<Isotopy::Triangle> triangles =
        Isotopy::edges_to_triangles(unpack_edges(t, edges), t.delta);
    Isotopy::Graph g(t.delta, unpack_signs(t, q1signs), triangles);
    g.isotopy_type();

    const uint64_t code = tree_code_from_adjacency(g.root_region, g.region_adjacency);
    r->code = code;
    r->p = g.p_regions;
    r->n = g.n_regions;
    r->regions = g.region_count;
    r->flags = kOk;
    return table.lookup(code);
  } catch (const std::exception&) {
    return kTypeFallback;
  }
}

inline void resolve_batch(const Tables& t, const TypeTable& table, const uint8_t* edges,
                          const uint64_t* q1signs, int batch, BatchOutput* out, int unknown_max) {
  const size_t edge_stride = 2 * (size_t)t.nedges;
  const size_t sign_stride = kQ1Words;

  out->type_id.assign(batch, kTypeFallback);
  out->unknown_codes.clear();
  out->unknown_witness.clear();
  out->unknown_overflow = false;
  out->table_fingerprint = table.fingerprint;
  out->flagged = 0;
  out->recomputed = 0;
  out->recompute_failed = 0;
  out->unknown = 0;
  out->device_disagreements = 0;

  std::unordered_map<uint64_t, int> seen;
  for (int i = 0; i < batch; ++i) {
    Result& r = out->results[i];
    const bool was_flagged = r.flags != kOk;
    uint16_t id;

    if (was_flagged) {
      ++out->flagged;
      id = resolve_one_cpu(t, table, edges + i * edge_stride, q1signs + i * sign_stride, &r);
      if (id == kTypeFallback) {
        ++out->recompute_failed;
      } else {
        ++out->recomputed;
      }
    } else {
      id = table.lookup(r.code);
      if (r.type_id != kTypeFallback && r.type_id != id) ++out->device_disagreements;
    }

    r.type_id = id;
    out->type_id[i] = id;
    if (id != kTypeUnknown) continue;

    ++out->unknown;
    if (seen.find(r.code) != seen.end()) continue;
    if ((int)out->unknown_codes.size() >= unknown_max) {
      out->unknown_overflow = true;
      continue;
    }
    seen.emplace(r.code, (int)out->unknown_codes.size());
    out->unknown_codes.push_back(r.code);
    out->unknown_witness.push_back(i);
  }
}

void classify_batch(const Tables& t, const TypeTable& table, const uint8_t* edges,
                    const uint64_t* q1signs, int batch, BatchOutput* out, int unknown_max = 4096);

inline void classify_batch_host(const Tables& t, const TypeTable& table, const uint8_t* edges,
                                const uint64_t* q1signs, int batch, BatchOutput* out,
                                int unknown_max = 4096) {
  const size_t edge_stride = 2 * (size_t)t.nedges;
  const size_t sign_stride = kQ1Words;

  out->results.assign(batch, Result{});
  for (int i = 0; i < batch; ++i) {
    out->results[i] = classify_host(t, edges + i * edge_stride, q1signs + i * sign_stride);
  }
  resolve_batch(t, table, edges, q1signs, batch, out, unknown_max);
}

}  // namespace IsotopyGPU
