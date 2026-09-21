#pragma once

#include <cstdint>

#if !defined(__CUDACC__) && defined(ISOTOPY_GPU_HOST_WARP)
#include <atomic>
#include <barrier>
#include <thread>
#endif

#ifndef ISOTOPY_GPU_MAX_DELTA
#define ISOTOPY_GPU_MAX_DELTA 9
#endif

// Stage-truncation profiling: build with -DISOTOPY_GPU_STOP_AFTER=k to return after stage k.
// Each truncation point folds the stage's own output into out.code, so nvcc cannot eliminate
// the work being measured. Differencing the wall times of successive k gives a per-stage cost
// without needing Nsight or profiling permissions.
#ifndef ISOTOPY_GPU_STOP_AFTER
#define ISOTOPY_GPU_STOP_AFTER 6
#endif

#ifdef __CUDACC__
#define ISO_HD __host__ __device__
#else
#define ISO_HD
#endif

namespace IsotopyGPU {

inline constexpr int kMaxDelta = ISOTOPY_GPU_MAX_DELTA;
inline constexpr int kMaxQ1Verts = (kMaxDelta + 1) * (kMaxDelta + 2) / 2;
inline constexpr int kMaxSphereVerts = 4 * kMaxDelta * kMaxDelta + 2;
inline constexpr int kMaxQ1Edges = 3 * kMaxDelta * (kMaxDelta + 1) / 2;
inline constexpr int kMaxRegions = 32;
inline constexpr int kMaxRounds = 64;
inline constexpr int kSphereWords = (kMaxSphereVerts + 63) / 64;
inline constexpr int kQ1Words = (kMaxQ1Verts + 63) / 64;
inline constexpr int kCodeBits = 63;

static_assert(kMaxQ1Verts <= 256, "Q1 vertex indices cross the API as uint8_t");
static_assert(kMaxRegions <= 32, "region adjacency rows are uint32_t bitsets");
static_assert(kCodeBits < 64, "the tree code plus its length tag must fit one uint64_t");
static_assert(kMaxSphereVerts <= 32768, "sphere vertex ids are stored as uint16_t");

enum Flag : uint32_t {
  kOk = 0u,
  kNoConverge = 1u << 0,
  kLabelsInconsistent = 1u << 1,
  kComponentCountMismatch = 1u << 2,
  kNoRoot = 1u << 3,
  kTooManyRegions = 1u << 4,
  kBfsIncomplete = 1u << 5,
  kCodeOverflow = 1u << 6,
  kDeltaOutOfRange = 1u << 7,
  kAmbiguousRoot = 1u << 8,
};

ISO_HD inline int popcount32(uint32_t x) {
#ifdef __CUDA_ARCH__
  return __popc(x);
#else
  int c = 0;
  while (x) {
    x &= x - 1u;
    ++c;
  }
  return c;
#endif
}

ISO_HD inline int lowest_bit(uint32_t x) {
  int p = 0;
  while (!((x >> p) & 1u)) ++p;
  return p;
}

struct Result {
  uint64_t code;
  int32_t p;
  int32_t n;
  int32_t regions;
  int32_t components;
  int32_t rounds;
  uint32_t flags;
};

struct Tables {
  int32_t delta;
  int32_t nverts;
  int32_t nsphere;
  int32_t nedges;
  uint32_t edge_magic;
  int32_t max_rounds;
  int32_t max_regions;
  int32_t code_bits;
  uint16_t oct[8 * kMaxQ1Verts];
  uint8_t flip[8 * kMaxQ1Verts];
  uint16_t anti[kMaxSphereVerts];
};

struct Scratch {
  int32_t label[kMaxSphereVerts];
  uint64_t sign[kSphereWords];
  uint8_t rid[kMaxSphereVerts];
  uint32_t adj[kMaxRegions];
  uint8_t level[kMaxRegions];
  uint8_t parent[kMaxRegions];
  uint64_t code[kMaxRegions];
  uint8_t codelen[kMaxRegions];
  uint64_t result_code;
};

// s / nedges for s < 8*nedges, without a runtime integer divide (sm_75 has no IDIV).
ISO_HD inline int slot_octant(int s, uint32_t magic) {
#ifdef __CUDA_ARCH__
  return (int)__umulhi((unsigned)s, magic);
#else
  return (int)(((uint64_t)(unsigned)s * (uint64_t)magic) >> 32);
#endif
}

ISO_HD inline bool get_bit(const uint64_t* words, int i) {
  return (words[i >> 6] >> (i & 63)) & 1ull;
}

template <int LANES>
struct Warp;

template <>
struct Warp<1> {
  static constexpr int lanes = 1;
  ISO_HD int lane() const { return 0; }
  ISO_HD bool leader() const { return true; }
  ISO_HD void sync() const {}
  ISO_HD bool any(bool p) const { return p; }
  ISO_HD uint32_t reduce_or(uint32_t v) const { return v; }
  ISO_HD int reduce_min(int v) const { return v; }
  ISO_HD int reduce_add(int v) const { return v; }
  ISO_HD int exclusive_scan(int v, int* total) const {
    *total = v;
    return 0;
  }
  ISO_HD int32_t load(const int32_t* addr) const { return *addr; }
  ISO_HD void atomic_min(int32_t* addr, int32_t v) const {
    if (v < *addr) *addr = v;
  }
  ISO_HD void atomic_or(uint32_t* addr, uint32_t v) const { *addr |= v; }
  ISO_HD void atomic_or64(uint64_t* addr, uint64_t v) const { *addr |= v; }
};

#ifdef __CUDACC__
template <>
struct Warp<32> {
  static constexpr int lanes = 32;
  static constexpr uint32_t kAll = 0xffffffffu;
  __device__ int lane() const { return threadIdx.x & 31; }
  __device__ bool leader() const { return lane() == 0; }
  __device__ void sync() const { __syncwarp(kAll); }
  __device__ bool any(bool p) const { return __any_sync(kAll, p); }
  __device__ uint32_t reduce_or(uint32_t v) const {
    for (int off = 16; off > 0; off >>= 1) v |= __shfl_xor_sync(kAll, v, off);
    return v;
  }
  __device__ int reduce_min(int v) const {
    for (int off = 16; off > 0; off >>= 1) {
      int o = __shfl_xor_sync(kAll, v, off);
      if (o < v) v = o;
    }
    return v;
  }
  __device__ int reduce_add(int v) const {
    for (int off = 16; off > 0; off >>= 1) v += __shfl_xor_sync(kAll, v, off);
    return v;
  }
  __device__ int exclusive_scan(int v, int* total) const {
    int x = v;
    for (int off = 1; off < 32; off <<= 1) {
      int y = __shfl_up_sync(kAll, x, off);
      if (lane() >= off) x += y;
    }
    *total = __shfl_sync(kAll, x, 31);
    return x - v;
  }
  __device__ int32_t load(const int32_t* addr) const {
    return *static_cast<const volatile int32_t*>(addr);
  }
  __device__ void atomic_min(int32_t* addr, int32_t v) const { atomicMin(addr, v); }
  __device__ void atomic_or(uint32_t* addr, uint32_t v) const { atomicOr(addr, v); }
  __device__ void atomic_or64(uint64_t* addr, uint64_t v) const {
    atomicOr(reinterpret_cast<unsigned long long*>(addr), (unsigned long long)v);
  }
};
#endif

#if !defined(__CUDACC__) && defined(ISOTOPY_GPU_HOST_WARP)

struct HostWarpContext {
  std::barrier<> gate{32};
  int64_t slot[32];
};

template <>
struct Warp<32> {
  static constexpr int lanes = 32;
  HostWarpContext* ctx;
  int id;

  int lane() const { return id; }
  bool leader() const { return id == 0; }
  void sync() const { ctx->gate.arrive_and_wait(); }

  int64_t exchange(int64_t v) const {
    ctx->slot[id] = v;
    ctx->gate.arrive_and_wait();
    return 0;
  }
  void release() const { ctx->gate.arrive_and_wait(); }

  bool any(bool p) const {
    exchange(p ? 1 : 0);
    bool r = false;
    for (int i = 0; i < 32; ++i) r = r || (ctx->slot[i] != 0);
    release();
    return r;
  }
  uint32_t reduce_or(uint32_t v) const {
    exchange((int64_t)v);
    uint32_t r = 0u;
    for (int i = 0; i < 32; ++i) r |= (uint32_t)ctx->slot[i];
    release();
    return r;
  }
  int reduce_min(int v) const {
    exchange(v);
    int64_t r = ctx->slot[0];
    for (int i = 1; i < 32; ++i) r = ctx->slot[i] < r ? ctx->slot[i] : r;
    release();
    return (int)r;
  }
  int reduce_add(int v) const {
    exchange(v);
    int64_t r = 0;
    for (int i = 0; i < 32; ++i) r += ctx->slot[i];
    release();
    return (int)r;
  }
  int exclusive_scan(int v, int* total) const {
    exchange(v);
    int64_t before = 0;
    int64_t sum = 0;
    for (int i = 0; i < 32; ++i) {
      if (i < id) before += ctx->slot[i];
      sum += ctx->slot[i];
    }
    release();
    *total = (int)sum;
    return (int)before;
  }

  int32_t load(const int32_t* addr) const {
    return std::atomic_ref<int32_t>(*const_cast<int32_t*>(addr)).load(std::memory_order_relaxed);
  }
  void atomic_min(int32_t* addr, int32_t v) const {
    std::atomic_ref<int32_t> a(*addr);
    int32_t cur = a.load(std::memory_order_relaxed);
    while (v < cur && !a.compare_exchange_weak(cur, v, std::memory_order_relaxed)) {
    }
  }
  void atomic_or(uint32_t* addr, uint32_t v) const {
    std::atomic_ref<uint32_t>(*addr).fetch_or(v, std::memory_order_relaxed);
  }
  void atomic_or64(uint64_t* addr, uint64_t v) const {
    std::atomic_ref<uint64_t>(*addr).fetch_or(v, std::memory_order_relaxed);
  }
};

#endif

template <int LANES>
ISO_HD void sphere_signs(const Tables& t, const uint64_t* q1signs, Scratch& sh, Warp<LANES> w) {
  for (int i = w.lane(); i < kSphereWords; i += LANES) sh.sign[i] = 0ull;
  w.sync();
  for (int v = w.lane(); v < t.nverts; v += LANES) {
    const bool s = get_bit(q1signs, v);
    for (int o = 0; o < 8; ++o) {
      const int slot = o * t.nverts + v;
      if (!(s ^ (bool)t.flip[slot])) continue;
      const int id = t.oct[slot];
      w.atomic_or64(&sh.sign[id >> 6], 1ull << (id & 63));
    }
  }
  w.sync();
}

template <int LANES>
ISO_HD int propagate_labels(const Tables& t, const uint8_t* edges, Scratch& sh, Warp<LANES> w) {
  for (int v = w.lane(); v < t.nsphere; v += LANES) sh.label[v] = v;
  w.sync();

  int round = 0;
  for (; round < t.max_rounds; ++round) {
    bool changed = false;

    const int slots = 8 * t.nedges;
    for (int sidx = w.lane(); sidx < slots; sidx += LANES) {
      const int o = slot_octant(sidx, t.edge_magic);
      const int e = sidx - o * t.nedges;
      const uint16_t* oct_o = t.oct + o * t.nverts;
        const int u = oct_o[edges[2 * e]];
        const int v = oct_o[edges[2 * e + 1]];
        if (get_bit(sh.sign, u) != get_bit(sh.sign, v)) continue;
        const int lu = w.load(&sh.label[u]);
        const int lv = w.load(&sh.label[v]);
        if (lu == lv) continue;
        const int m = lu < lv ? lu : lv;
        if (m < lu) {
          w.atomic_min(&sh.label[u], m);
          changed = true;
        }
        if (m < lv) {
          w.atomic_min(&sh.label[v], m);
          changed = true;
        }
    }
    w.sync();

    for (int v = w.lane(); v < t.nsphere; v += LANES) {
      const int old = w.load(&sh.label[v]);
      const int jumped = w.load(&sh.label[old]);
      if (jumped != old) {
        w.atomic_min(&sh.label[v], jumped);
        changed = true;
      }
    }
    w.sync();

    if (!w.any(changed)) break;
  }
  return round;
}

template <int LANES>
ISO_HD bool labels_consistent(const Tables& t, const uint8_t* edges, Scratch& sh, Warp<LANES> w) {
  bool bad = false;
  const int slots = 8 * t.nedges;
  for (int sidx = w.lane(); sidx < slots; sidx += LANES) {
    const int o = slot_octant(sidx, t.edge_magic);
    const int e = sidx - o * t.nedges;
    const uint16_t* oct_o = t.oct + o * t.nverts;
    const int u = oct_o[edges[2 * e]];
    const int v = oct_o[edges[2 * e + 1]];
    if (get_bit(sh.sign, u) == get_bit(sh.sign, v) && sh.label[u] != sh.label[v]) bad = true;
  }
  for (int v = w.lane(); v < t.nsphere; v += LANES) {
    if (sh.label[sh.label[v]] != sh.label[v]) bad = true;
  }
  return !w.any(bad);
}

ISO_HD inline int region_key(const Scratch& sh, const Tables& t, int v) {
  const int a = sh.label[t.anti[v]];
  const int b = sh.label[v];
  return a < b ? a : b;
}

template <int LANES>
ISO_HD void assign_region_ids(const Tables& t, Scratch& sh, Warp<LANES> w, int* n_reg) {
  constexpr int kPerLane = (kMaxSphereVerts + LANES - 1) / LANES;

  int local = 0;
  for (int v = w.lane(); v < t.nsphere; v += LANES) {
    if (region_key(sh, t, v) == v) ++local;
  }
  int total = 0;
  int base = w.exclusive_scan(local, &total);
  *n_reg = total;
  if (total > t.max_regions) return;

  int next = base;
  for (int v = w.lane(); v < t.nsphere; v += LANES) {
    if (region_key(sh, t, v) == v) sh.rid[v] = (uint8_t)next++;
  }
  w.sync();

  uint8_t gathered[kPerLane];
  int k = 0;
  for (int v = w.lane(); v < t.nsphere; v += LANES) gathered[k++] = sh.rid[region_key(sh, t, v)];
  w.sync();

  k = 0;
  for (int v = w.lane(); v < t.nsphere; v += LANES) sh.rid[v] = gathered[k++];
  w.sync();
}

template <int LANES>
ISO_HD void build_adjacency(const Tables& t, const uint8_t* edges, Scratch& sh, Warp<LANES> w,
                            int n_reg) {
  for (int r = w.lane(); r < kMaxRegions; r += LANES) sh.adj[r] = 0u;
  w.sync();

  const int slots = 8 * t.nedges;
  for (int sidx = w.lane(); sidx < slots; sidx += LANES) {
    const int o = slot_octant(sidx, t.edge_magic);
    const int e = sidx - o * t.nedges;
    const uint16_t* oct_o = t.oct + o * t.nverts;
    const int u = oct_o[edges[2 * e]];
    const int v = oct_o[edges[2 * e + 1]];
    if (get_bit(sh.sign, u) == get_bit(sh.sign, v)) continue;
    const int ru = sh.rid[u];
    const int rv = sh.rid[v];
    w.atomic_or(&sh.adj[ru], 1u << rv);
    w.atomic_or(&sh.adj[rv], 1u << ru);
  }
  w.sync();

  for (int r = w.lane(); r < n_reg; r += LANES) sh.adj[r] &= ~(1u << r);
  w.sync();
}

template <int LANES>
ISO_HD bool breadth_first(Scratch& sh, Warp<LANES> w, int n_reg, int root) {
  for (int r = w.lane(); r < kMaxRegions; r += LANES) {
    sh.level[r] = 0;
    sh.parent[r] = 0;
  }
  w.sync();

  uint32_t visited = 1u << root;
  uint32_t frontier = visited;

  for (int depth = 1; depth <= n_reg; ++depth) {
    uint32_t contrib = 0u;
    for (int r = w.lane(); r < n_reg; r += LANES) {
      if ((frontier >> r) & 1u) contrib |= sh.adj[r];
    }
    const uint32_t next = w.reduce_or(contrib) & ~visited;
    if (next == 0u) break;

    for (int r = w.lane(); r < n_reg; r += LANES) {
      if (!((next >> r) & 1u)) continue;
      const uint32_t from = sh.adj[r] & frontier;
      if (from == 0u) continue;
      sh.parent[r] = (uint8_t)lowest_bit(from);
      sh.level[r] = (uint8_t)depth;
    }
    w.sync();

    visited |= next;
    frontier = next;
  }

  const uint32_t full = n_reg >= 32 ? 0xffffffffu : ((1u << n_reg) - 1u);
  return visited == full;
}

ISO_HD inline bool code_less(uint64_t ca, uint8_t la, uint64_t cb, uint8_t lb) {
  if (la != lb) return la < lb;
  return ca < cb;
}

ISO_HD inline bool encode_tree(Scratch& sh, int n_reg, int root, int code_bits,
                               uint64_t* out_code) {
  uint8_t first_child[kMaxRegions];
  uint8_t next_sibling[kMaxRegions];
  uint8_t order[kMaxRegions];

  int max_level = 0;
  for (int r = 0; r < n_reg; ++r) {
    sh.code[r] = 0ull;
    sh.codelen[r] = 0;
    first_child[r] = (uint8_t)kMaxRegions;
    next_sibling[r] = (uint8_t)kMaxRegions;
    if (sh.level[r] > max_level) max_level = sh.level[r];
  }

  for (int c = 0; c < n_reg; ++c) {
    if (c == root) continue;
    const int p = sh.parent[c];
    next_sibling[c] = first_child[p];
    first_child[p] = (uint8_t)c;
  }

  for (int depth = max_level; depth >= 1; --depth) {
    for (int p = 0; p < n_reg; ++p) {
      if (sh.level[p] != depth - 1) continue;

      int m = 0;
      for (int c = first_child[p]; c != kMaxRegions; c = next_sibling[c]) {
        int j = m - 1;
        while (j >= 0 && code_less(sh.code[c], sh.codelen[c], sh.code[order[j]],
                                   sh.codelen[order[j]])) {
          order[j + 1] = order[j];
          --j;
        }
        order[j + 1] = (uint8_t)c;
        ++m;
      }
      if (m == 0) continue;

      int off = 0;
      for (int i = 0; i < m; ++i) {
        const int c = order[i];
        const int span = sh.codelen[c] + 2;
        if (off + span > code_bits) return false;
        sh.code[p] |= (1ull | (sh.code[c] << 1)) << off;
        off += span;
      }
      sh.codelen[p] = (uint8_t)off;
    }
  }

  if (sh.codelen[root] >= code_bits) return false;
  *out_code = (1ull << sh.codelen[root]) | sh.code[root];
  return true;
}

template <int LANES>
ISO_HD void classify_one(const Tables& t, const uint8_t* edges, const uint64_t* q1signs,
                         Scratch& sh, Warp<LANES> w, Result& out) {
  out.code = 0ull;
  out.p = 0;
  out.n = 0;
  out.regions = 0;
  out.components = 0;
  out.rounds = 0;
  out.flags = kOk;

  if (t.delta < 1 || t.delta > kMaxDelta) {
    out.flags |= kDeltaOutOfRange;
    return;
  }

  sphere_signs(t, q1signs, sh, w);
#if ISOTOPY_GPU_STOP_AFTER <= 1
  { uint32_t dig = 0u;
    for (int i = w.lane(); i < kSphereWords; i += LANES) dig ^= (uint32_t)sh.sign[i];
    out.code = w.reduce_or(dig); return; }
#endif
  out.rounds = propagate_labels(t, edges, sh, w);
  if (out.rounds >= t.max_rounds) out.flags |= kNoConverge;
  if (!labels_consistent(t, edges, sh, w)) out.flags |= kLabelsInconsistent;
  if (out.flags) return;
#if ISOTOPY_GPU_STOP_AFTER <= 2
  { int dig = 0;
    for (int v = w.lane(); v < t.nsphere; v += LANES) dig += sh.label[v];
    out.code = (uint64_t)w.reduce_add(dig); return; }
#endif

  int local_components = 0;
  for (int v = w.lane(); v < t.nsphere; v += LANES) {
    if (sh.label[v] == v) ++local_components;
  }
  out.components = w.reduce_add(local_components);

  int n_reg = 0;
  assign_region_ids(t, sh, w, &n_reg);
  out.regions = n_reg;
  if (n_reg > t.max_regions) {
    out.flags |= kTooManyRegions;
    return;
  }

  const bool even = (t.delta % 2) == 0;
  const int expected = even ? 2 * n_reg - 1 : 2 * n_reg;
  if (out.components != expected) {
    out.flags |= kComponentCountMismatch;
    return;
  }

#if ISOTOPY_GPU_STOP_AFTER <= 3
  { int dig = out.components * 131 + n_reg;
    for (int v = w.lane(); v < t.nsphere; v += LANES) dig += sh.rid[v];
    out.code = (uint64_t)w.reduce_add(dig); return; }
#endif

  build_adjacency(t, edges, sh, w, n_reg);
#if ISOTOPY_GPU_STOP_AFTER <= 4
  { uint32_t dig = 0u;
    for (int r = w.lane(); r < n_reg; r += LANES) dig ^= sh.adj[r];
    out.code = w.reduce_or(dig); return; }
#endif

  uint32_t witness = 0u;
  if (even) {
    for (int v = w.lane(); v < t.nsphere; v += LANES) {
      if (sh.label[v] != sh.label[t.anti[v]]) continue;
      witness |= 1u << sh.rid[v];
    }
  } else {
    const int slots = 8 * t.nedges;
    for (int sidx = w.lane(); sidx < slots; sidx += LANES) {
      const int o = slot_octant(sidx, t.edge_magic);
      const int e = sidx - o * t.nedges;
      const uint16_t* oct_o = t.oct + o * t.nverts;
      const int u = oct_o[edges[2 * e]];
      const int v = oct_o[edges[2 * e + 1]];
      if (get_bit(sh.sign, u) == get_bit(sh.sign, v)) continue;
      if (sh.label[u] != sh.label[t.anti[v]]) continue;
      witness |= 1u << sh.rid[u];
    }
  }
  witness = w.reduce_or(witness);
  if (witness == 0u) {
    out.flags |= kNoRoot;
    return;
  }
  if (popcount32(witness) != 1) {
    out.flags |= kAmbiguousRoot;
    return;
  }
  const int root = lowest_bit(witness);

  if (!breadth_first(sh, w, n_reg, root)) {
    out.flags |= kBfsIncomplete;
    return;
  }

  int local_even = 0;
  for (int r = w.lane(); r < n_reg; r += LANES) {
    if ((sh.level[r] & 1) == 0) ++local_even;
  }
  const int even_depth = w.reduce_add(local_even);
  const int odd_depth = n_reg - even_depth;
  out.p = even ? odd_depth : even_depth;
  out.n = even ? even_depth - 1 : odd_depth;
#if ISOTOPY_GPU_STOP_AFTER <= 5
  { int dig = root;
    for (int r = w.lane(); r < n_reg; r += LANES) dig += sh.level[r] * 7 + sh.parent[r];
    out.code = (uint64_t)w.reduce_add(dig); return; }
#endif

  if (w.leader()) {
    uint64_t code = 0ull;
    if (!encode_tree(sh, n_reg, root, t.code_bits, &code)) {
      out.flags |= kCodeOverflow;
    } else {
      sh.result_code = code;
    }
  }
  w.sync();
  out.flags = w.reduce_or(out.flags);
  if (!out.flags) out.code = sh.result_code;
}

}  // namespace IsotopyGPU

#ifndef __CUDA_ARCH__

#include <algorithm>
#include <array>
#include <map>
#include <stdexcept>
#include <vector>

namespace IsotopyGPU {

inline int q1_index(int delta, int x, int y) {
  return x * (delta + 1) - x * (x - 1) / 2 + y;
}

inline std::vector<std::array<int, 3>> sphere_points(int delta) {
  std::vector<std::array<int, 3>> pts;
  for (int a = -delta; a <= delta; ++a) {
    const int ra = delta - (a < 0 ? -a : a);
    for (int b = -ra; b <= ra; ++b) {
      const int rb = ra - (b < 0 ? -b : b);
      pts.push_back({a, b, rb});
      if (rb != 0) pts.push_back({a, b, -rb});
    }
  }
  std::sort(pts.begin(), pts.end());
  return pts;
}

inline Tables build_tables(int delta) {
  if (delta < 1 || delta > kMaxDelta) {
    throw std::invalid_argument("IsotopyGPU::build_tables: delta out of range");
  }

  Tables t{};
  t.delta = delta;
  t.nverts = (delta + 1) * (delta + 2) / 2;
  t.nsphere = 4 * delta * delta + 2;
  t.nedges = 3 * delta * (delta + 1) / 2;
  {
    const uint32_t d = (uint32_t)t.nedges;
    t.edge_magic = (uint32_t)((0x100000000ULL + d - 1) / d);
    for (int sidx = 0; sidx < 8 * t.nedges; ++sidx) {
      if (slot_octant(sidx, t.edge_magic) != sidx / t.nedges) {
        throw std::logic_error("IsotopyGPU::build_tables: edge_magic is wrong");
      }
    }
  }
  t.max_rounds = kMaxRounds;
  t.max_regions = kMaxRegions;
  t.code_bits = kCodeBits;

  const auto pts = sphere_points(delta);
  if ((int)pts.size() != t.nsphere) {
    throw std::logic_error("IsotopyGPU::build_tables: sphere vertex count mismatch");
  }

  std::map<std::array<int, 3>, int> id;
  for (int i = 0; i < (int)pts.size(); ++i) id[pts[i]] = i;

  for (int i = 0; i < (int)pts.size(); ++i) {
    t.anti[i] = (uint16_t)id.at({-pts[i][0], -pts[i][1], -pts[i][2]});
  }

  for (int x = 0; x <= delta; ++x) {
    for (int y = 0; y + x <= delta; ++y) {
      const int z = delta - x - y;
      const int v = q1_index(delta, x, y);
      for (int o = 0; o < 8; ++o) {
        const int sx = o & 1;
        const int sy = (o >> 1) & 1;
        const int sz = (o >> 2) & 1;
        const std::array<int, 3> p{sx ? -x : x, sy ? -y : y, sz ? -z : z};
        const int slot = o * t.nverts + v;
        t.oct[slot] = (uint16_t)id.at(p);
        t.flip[slot] = (uint8_t)((sx & (x & 1)) ^ (sy & (y & 1)) ^ (sz & (z & 1)));
      }
    }
  }
  return t;
}

inline Result classify_host(const Tables& t, const uint8_t* edges, const uint64_t* q1signs) {
  Scratch sh;
  Result out{};
  classify_one(t, edges, q1signs, sh, Warp<1>{}, out);
  return out;
}

#ifdef ISOTOPY_GPU_HOST_WARP
inline Result classify_host_threaded(const Tables& t, const uint8_t* edges,
                                     const uint64_t* q1signs) {
  Scratch sh;
  HostWarpContext ctx;
  Result out[32]{};

  std::vector<std::thread> lanes;
  lanes.reserve(32);
  for (int i = 0; i < 32; ++i) {
    lanes.emplace_back([&, i]() {
      classify_one(t, edges, q1signs, sh, Warp<32>{&ctx, i}, out[i]);
    });
  }
  for (auto& th : lanes) th.join();
  return out[0];
}
#endif

inline uint64_t tree_code_from_adjacency(int root, const std::vector<std::vector<int>>& adjacency) {
  const int n = (int)adjacency.size();
  if (root < 0 || root >= n) throw std::invalid_argument("tree_code_from_adjacency: bad root");

  std::vector<int> level(n, -1);
  std::vector<int> parent(n, -1);
  std::vector<int> order;
  order.reserve(n);
  level[root] = 0;
  order.push_back(root);
  for (size_t i = 0; i < order.size(); ++i) {
    const int cur = order[i];
    for (int nb : adjacency[cur]) {
      if (nb == cur || level[nb] != -1) continue;
      level[nb] = level[cur] + 1;
      parent[nb] = cur;
      order.push_back(nb);
    }
  }

  std::vector<uint64_t> code(n, 0ull);
  std::vector<int> len(n, 0);
  for (auto it = order.rbegin(); it != order.rend(); ++it) {
    const int cur = *it;
    std::vector<int> kids;
    for (int c = 0; c < n; ++c) {
      if (parent[c] == cur) kids.push_back(c);
    }
    std::sort(kids.begin(), kids.end(), [&](int a, int b) {
      if (len[a] != len[b]) return len[a] < len[b];
      return code[a] < code[b];
    });
    int off = 0;
    for (int c : kids) {
      if (off + len[c] + 2 > kCodeBits) throw std::overflow_error("tree_code_from_adjacency");
      code[cur] |= (1ull | (code[c] << 1)) << off;
      off += len[c] + 2;
    }
    len[cur] = off;
  }
  return (1ull << len[root]) | code[root];
}

}  // namespace IsotopyGPU

#endif
