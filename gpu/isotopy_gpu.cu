#include <cstdio>
#include <cstdlib>

#include "isotopy_gpu.cuh"

#ifndef __CUDA_ARCH__
#include "isotopy_batch.h"
#endif

namespace IsotopyGPU {

inline constexpr int kWarpsPerBlock = 8;
inline constexpr int kThreadsPerBlock = 32 * kWarpsPerBlock;

static_assert(sizeof(Tables) % sizeof(uint32_t) == 0,
              "Tables is staged into shared memory a uint32_t at a time");
static_assert(sizeof(Scratch) * kWarpsPerBlock + sizeof(Tables) <= 48 * 1024,
              "static __shared__ per block exceeds the 48 KB limit; lower "
              "ISOTOPY_GPU_MAX_DELTA or kWarpsPerBlock");

__constant__ Tables g_tables;

__global__ void classify_kernel(const uint8_t* edges, const uint64_t* signs, Result* out,
                                int batch, int edge_stride, int sign_stride, const uint64_t* keys,
                                int nkeys) {
  __shared__ Tables tables;
  __shared__ Scratch scratch[kWarpsPerBlock];

  uint32_t* dst = reinterpret_cast<uint32_t*>(&tables);
  const uint32_t* src = reinterpret_cast<const uint32_t*>(&g_tables);
  for (int i = threadIdx.x; i < (int)(sizeof(Tables) / sizeof(uint32_t)); i += blockDim.x) {
    dst[i] = src[i];
  }
  __syncthreads();

  const int slot = threadIdx.x >> 5;
  const int instance = (int)blockIdx.x * kWarpsPerBlock + slot;
  if (instance >= batch) return;

  Result r;
  classify_one(tables, edges + (size_t)instance * edge_stride,
               signs + (size_t)instance * sign_stride, scratch[slot], Warp<32>{}, r);

#ifdef ISOTOPY_GPU_DEVICE_TABLE
  if (r.flags == kOk && nkeys > 0) r.type_id = lookup_type(keys, nkeys, r.code, Warp<32>{});
#else
  (void)keys;
  (void)nkeys;
#endif

  if ((threadIdx.x & 31) == 0) out[instance] = r;
}

#define ISO_CUDA_CHECK(expr)                                                              \
  do {                                                                                    \
    cudaError_t err_ = (expr);                                                             \
    if (err_ != cudaSuccess) {                                                             \
      std::fprintf(stderr, "%s:%d %s -> %s\n", __FILE__, __LINE__, #expr,                  \
                   cudaGetErrorString(err_));                                              \
      std::abort();                                                                        \
    }                                                                                      \
  } while (0)

struct DeviceBatch {
  uint8_t* edges = nullptr;
  uint64_t* signs = nullptr;
  Result* out = nullptr;
  uint64_t* keys = nullptr;
  int nkeys = 0;
  int capacity = 0;
  int edge_stride = 0;
  int sign_stride = 0;
};

void device_batch_alloc(DeviceBatch& d, const Tables& t, int capacity, const uint64_t* keys,
                        int nkeys) {
  d.capacity = capacity;
  d.edge_stride = 2 * t.nedges;
  d.sign_stride = kQ1Words;
  ISO_CUDA_CHECK(cudaMemcpyToSymbol(g_tables, &t, sizeof(Tables)));
  ISO_CUDA_CHECK(cudaMalloc(&d.edges, (size_t)capacity * d.edge_stride * sizeof(uint8_t)));
  ISO_CUDA_CHECK(cudaMalloc(&d.signs, (size_t)capacity * d.sign_stride * sizeof(uint64_t)));
  ISO_CUDA_CHECK(cudaMalloc(&d.out, (size_t)capacity * sizeof(Result)));
  if (keys != nullptr && nkeys > 0) {
    d.nkeys = nkeys;
    ISO_CUDA_CHECK(cudaMalloc(&d.keys, (size_t)nkeys * sizeof(uint64_t)));
    ISO_CUDA_CHECK(
        cudaMemcpy(d.keys, keys, (size_t)nkeys * sizeof(uint64_t), cudaMemcpyHostToDevice));
  }
}

void device_batch_alloc(DeviceBatch& d, const Tables& t, int capacity) {
  device_batch_alloc(d, t, capacity, nullptr, 0);
}

void device_batch_free(DeviceBatch& d) {
  ISO_CUDA_CHECK(cudaFree(d.edges));
  ISO_CUDA_CHECK(cudaFree(d.signs));
  ISO_CUDA_CHECK(cudaFree(d.out));
  if (d.keys != nullptr) ISO_CUDA_CHECK(cudaFree(d.keys));
  d = DeviceBatch{};
}

void device_batch_upload(DeviceBatch& d, const uint8_t* edges, const uint64_t* signs, int batch) {
  ISO_CUDA_CHECK(cudaMemcpy(d.edges, edges, (size_t)batch * d.edge_stride * sizeof(uint8_t),
                            cudaMemcpyHostToDevice));
  ISO_CUDA_CHECK(cudaMemcpy(d.signs, signs, (size_t)batch * d.sign_stride * sizeof(uint64_t),
                            cudaMemcpyHostToDevice));
}

void device_batch_launch(DeviceBatch& d, int batch) {
  const int blocks = (batch + kWarpsPerBlock - 1) / kWarpsPerBlock;
  classify_kernel<<<blocks, kThreadsPerBlock>>>(d.edges, d.signs, d.out, batch, d.edge_stride,
                                                d.sign_stride, d.keys, d.nkeys);
  ISO_CUDA_CHECK(cudaGetLastError());
  ISO_CUDA_CHECK(cudaDeviceSynchronize());
}

void device_batch_download(DeviceBatch& d, Result* out, int batch) {
  ISO_CUDA_CHECK(cudaMemcpy(out, d.out, (size_t)batch * sizeof(Result), cudaMemcpyDeviceToHost));
}

#ifndef __CUDA_ARCH__
void classify_batch(const Tables& t, const TypeTable& table, const uint8_t* edges,
                    const uint64_t* signs, int batch, BatchOutput* out, int unknown_max) {
  DeviceBatch d;
  device_batch_alloc(d, t, batch, table.keys.data(), (int)table.keys.size());
  device_batch_upload(d, edges, signs, batch);
  device_batch_launch(d, batch);
  out->results.assign(batch, Result{});
  device_batch_download(d, out->results.data(), batch);
  device_batch_free(d);
  resolve_batch(t, table, edges, signs, batch, out, unknown_max);
}
#endif

void classify_batch_device(const Tables& t, const uint8_t* edges, const uint64_t* signs,
                           Result* out, int batch) {
  DeviceBatch d;
  device_batch_alloc(d, t, batch);
  device_batch_upload(d, edges, signs, batch);
  device_batch_launch(d, batch);
  device_batch_download(d, out, batch);
  device_batch_free(d);
}

}  // namespace IsotopyGPU


#if defined(ISOTOPY_GPU_MAIN) && !defined(__CUDA_ARCH__)

#include <algorithm>
#include <chrono>
#include <cstring>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "isotopy_graph.h"

namespace {

using Clock = std::chrono::steady_clock;

double seconds_since(const Clock::time_point& t0) {
  return std::chrono::duration<double>(Clock::now() - t0).count();
}

struct Triangulation {
  std::vector<Isotopy::Edge> edges;
  std::vector<Isotopy::Triangle> triangles;
};

std::vector<Triangulation> make_pool(int delta, int pool_size, unsigned seed) {
  std::vector<Triangulation> pool;
  pool.reserve(pool_size);
  std::srand(seed);
  int attempts = 0;
  while ((int)pool.size() < pool_size && attempts < pool_size * 50) {
    ++attempts;
    Triangulation t;
    t.edges = Utils::get_random_triangulation(delta);
    t.triangles = Isotopy::edges_to_triangles(t.edges, delta);
    if ((int)t.triangles.size() != delta * delta) continue;
    pool.push_back(std::move(t));
  }
  if (pool.empty()) {
    std::fprintf(stderr, "could not build any triangulation at delta %d\n", delta);
    std::abort();
  }
  return pool;
}

struct Batch {
  std::vector<uint8_t> edges;
  std::vector<uint64_t> signs;
  std::vector<std::vector<bool>> sign_vectors;
  std::vector<int> pool_index;
};

Batch make_batch(const IsotopyGPU::Tables& t, const std::vector<Triangulation>& pool, int count,
                 unsigned seed) {
  Batch b;
  std::srand(seed);
  b.edges.resize((size_t)count * 2 * t.nedges);
  b.signs.assign((size_t)count * IsotopyGPU::kQ1Words, 0ull);
  b.sign_vectors.resize(count);
  b.pool_index.resize(count);

  for (int i = 0; i < count; ++i) {
    const int pi = i % (int)pool.size();
    b.pool_index[i] = pi;
    const std::vector<Isotopy::Edge>& e = pool[pi].edges;
    for (int k = 0; k < t.nedges; ++k) {
      b.edges[(size_t)i * 2 * t.nedges + 2 * k] = (uint8_t)e[k].first;
      b.edges[(size_t)i * 2 * t.nedges + 2 * k + 1] = (uint8_t)e[k].second;
    }
    b.sign_vectors[i].resize(t.nverts);
    for (int v = 0; v < t.nverts; ++v) {
      const bool s = ((std::rand() >> 8) & 1) != 0;
      b.sign_vectors[i][v] = s;
      if (s) b.signs[(size_t)i * IsotopyGPU::kQ1Words + (v >> 6)] |= 1ull << (v & 63);
    }
  }
  return b;
}

IsotopyGPU::TypeTable table_for(int delta, const std::vector<Triangulation>& pool, const Batch& b,
                                const std::string& yaml_path, std::string* source) {
#ifdef ISOTOPY_HAVE_FKYAML
  if (!yaml_path.empty()) {
    *source = yaml_path;
    return IsotopyGPU::load_type_table(yaml_path, delta);
  }
#else
  if (!yaml_path.empty()) {
    std::fprintf(stderr, "this binary was built without ISOTOPY_HAVE_FKYAML; ignoring %s\n",
                 yaml_path.c_str());
  }
#endif
  *source = "libisotopy over this batch";
  std::set<std::string> distinct;
  for (size_t i = 0; i < b.sign_vectors.size(); ++i) {
    Isotopy::Graph g(delta, b.sign_vectors[i], pool[b.pool_index[i]].triangles);
    g.isotopy_type();
    distinct.insert(g.viro_notation());
  }
  std::vector<IsotopyGPU::TypeEntry> raw;
  for (const std::string& v : distinct) {
    IsotopyGPU::TypeEntry e;
    e.viro = v;
    const IsotopyGPU::Viro::Parsed parsed = IsotopyGPU::Viro::parse(v);
    IsotopyGPU::Viro::region_pn(parsed.tree, parsed.has_j, &e.p, &e.n);
    e.num_ovals = IsotopyGPU::Viro::num_ovals(parsed.tree, parsed.has_j);
    raw.push_back(e);
  }
  return IsotopyGPU::build_type_table(delta, raw);
}

struct Stats {
  double best = 0.0;
  double mean = 0.0;
};

Stats rate_stats(const std::vector<double>& secs, int count) {
  Stats s;
  double total = 0.0;
  for (double x : secs) {
    const double r = count / x;
    if (r > s.best) s.best = r;
    total += r;
  }
  s.mean = total / (double)secs.size();
  return s;
}

int verify_sample(const IsotopyGPU::Tables& t, const Batch& b,
                  const std::vector<Triangulation>& pool,
                  const std::vector<IsotopyGPU::Result>& out, int sample, int* flagged) {
  int mismatched = 0;
  *flagged = 0;
  for (int i = 0; i < sample; ++i) {
    if (out[i].flags != IsotopyGPU::kOk) {
      ++(*flagged);
      continue;
    }
    Isotopy::Graph g(t.delta, b.sign_vectors[i], pool[b.pool_index[i]].triangles);
    g.isotopy_type();
    const uint64_t want =
        IsotopyGPU::tree_code_from_adjacency(g.root_region, g.region_adjacency);
    if (out[i].code != want || out[i].p != g.p_regions || out[i].n != g.n_regions) ++mismatched;
  }
  return mismatched;
}

void report_rounds(const std::vector<IsotopyGPU::Result>& out) {
  int maxr = 0;
  long long total = 0;
  int hist[IsotopyGPU::kMaxRounds + 1] = {0};
  for (const auto& r : out) {
    const int k = r.rounds < 0 ? 0 : (r.rounds > IsotopyGPU::kMaxRounds ? IsotopyGPU::kMaxRounds
                                                                       : r.rounds);
    ++hist[k];
    total += r.rounds;
    if (r.rounds > maxr) maxr = r.rounds;
  }
  std::printf("  rounds mean %.2f  max %d  histogram", (double)total / (double)out.size(), maxr);
  for (int k = 0; k <= maxr; ++k) {
    if (hist[k]) std::printf(" %d:%d", k, hist[k]);
  }
  std::printf("\n");
}

double cpu_baseline(int delta, const std::vector<Triangulation>& pool, const Batch& b,
                    double budget_seconds) {
  const auto t0 = Clock::now();
  long long done = 0;
  const int n = (int)b.sign_vectors.size();
  while (seconds_since(t0) < budget_seconds) {
    const int i = (int)(done % n);
    Isotopy::Graph g(delta, b.sign_vectors[i], pool[b.pool_index[i]].triangles);
    g.isotopy_type();
    volatile uint64_t code =
        IsotopyGPU::tree_code_from_adjacency(g.root_region, g.region_adjacency);
    (void)code;
    ++done;
  }
  return (double)done / seconds_since(t0);
}

}  // namespace

int main(int argc, char** argv) {
  const std::string mode = argc > 1 ? argv[1] : "bench";
  const int delta = argc > 2 ? std::atoi(argv[2]) : 8;
  const int count = argc > 3 ? std::atoi(argv[3]) : 65536;
  const int pool_size = argc > 4 ? std::atoi(argv[4]) : 1024;
  const int reps = argc > 5 ? std::atoi(argv[5]) : 7;
  const std::string yaml = argc > 6 ? argv[6] : "";

  if (delta < 1 || delta > IsotopyGPU::kMaxDelta) {
    std::fprintf(stderr,
                 "delta %d is outside this build's range (ISOTOPY_GPU_MAX_DELTA=%d).\n"
                 "Rebuild with GPU_MAX_DELTA=%d; note that raising it costs occupancy.\n",
                 delta, IsotopyGPU::kMaxDelta, delta);
    return 2;
  }
  const IsotopyGPU::Tables t = IsotopyGPU::build_tables(delta);

  if (mode == "types") {
    const std::vector<Triangulation> pool = make_pool(delta, pool_size, 4242u);
    const Batch b = make_batch(t, pool, count, 12345u);

    std::string source;
    const IsotopyGPU::TypeTable table = table_for(delta, pool, b, yaml, &source);
    std::printf("delta %d  batch %d  pool %zu\n", delta, count, pool.size());
    std::printf("  table       %zu types from %s (fingerprint %016llx, %zu unparsed)\n",
                table.entries.size(), source.c_str(),
                (unsigned long long)table.fingerprint, table.rejected.size());

    IsotopyGPU::BatchOutput out;
    IsotopyGPU::classify_batch(t, table, b.edges.data(), b.signs.data(), count, &out);

    std::map<uint16_t, int> histogram;
    for (uint16_t id : out.type_id) ++histogram[id];
    std::printf("  resolved    %zu distinct ids over %d instances\n", histogram.size(), count);
    std::printf("  flagged     %d  recomputed %d  recompute_failed %d\n", out.flagged,
                out.recomputed, out.recompute_failed);
    std::printf("  unknown     %d instances, %zu distinct codes, overflow %d\n", out.unknown,
                out.unknown_codes.size(), (int)out.unknown_overflow);
    std::printf("  device ids  %d disagreements with the host search\n", out.device_disagreements);

    std::vector<std::pair<int, uint16_t>> ranked;
    for (const auto& kv : histogram) ranked.push_back(std::make_pair(kv.second, kv.first));
    std::sort(ranked.rbegin(), ranked.rend());
    for (size_t k = 0; k < ranked.size() && k < 8; ++k) {
      const IsotopyGPU::TypeEntry* e = table.entry(ranked[k].second);
      std::printf("    %7d x id %5u  %s\n", ranked[k].first, (unsigned)ranked[k].second,
                  e == nullptr ? "<unresolved>" : e->viro.c_str());
    }
    for (size_t k = 0; k < out.unknown_codes.size() && k < 8; ++k) {
      std::printf("    unknown code %020llu  witness instance %d\n",
                  (unsigned long long)out.unknown_codes[k], out.unknown_witness[k]);
    }

    std::printf("METRIC gpu_types_delta%d_unknown %d\n", delta, out.unknown);
    std::printf("METRIC gpu_types_delta%d_flagged %d\n", delta, out.flagged);
    std::printf("METRIC gpu_types_delta%d_distinct %zu\n", delta, histogram.size());
    return (out.recompute_failed == 0 && out.device_disagreements == 0) ? 0 : 1;
  }

  if (mode == "sweep") {
    const std::vector<Triangulation> pool = make_pool(delta, pool_size, 4242u);
    std::printf("delta %d  pool %zu  batch-size sweep (kernel-only, best of %d)\n", delta,
                pool.size(), reps);
    for (int n = 1024; n <= count; n *= 4) {
      const Batch b = make_batch(t, pool, n, 12345u);
      std::vector<IsotopyGPU::Result> out(n);
      IsotopyGPU::DeviceBatch d;
      IsotopyGPU::device_batch_alloc(d, t, n);
      IsotopyGPU::device_batch_upload(d, b.edges.data(), b.signs.data(), n);
      IsotopyGPU::device_batch_launch(d, n);
      std::vector<double> secs;
      for (int r = 0; r < reps; ++r) {
        const auto t0 = Clock::now();
        IsotopyGPU::device_batch_launch(d, n);
        secs.push_back(seconds_since(t0));
      }
      IsotopyGPU::device_batch_download(d, out.data(), n);
      IsotopyGPU::device_batch_free(d);
      const Stats s = rate_stats(secs, n);
      std::printf("  batch %7d   kernel %10.0f /s (mean %10.0f)\n", n, s.best, s.mean);
    }
    return 0;
  }

  const std::vector<Triangulation> pool = make_pool(delta, pool_size, 4242u);
  const Batch b = make_batch(t, pool, count, 12345u);
  std::vector<IsotopyGPU::Result> out(count);

  IsotopyGPU::DeviceBatch d;
  IsotopyGPU::device_batch_alloc(d, t, count);
  IsotopyGPU::device_batch_upload(d, b.edges.data(), b.signs.data(), count);
  IsotopyGPU::device_batch_launch(d, count);
  IsotopyGPU::device_batch_download(d, out.data(), count);

  if (mode == "time") {
    std::vector<double> secs;
    for (int r = 0; r < reps; ++r) {
      const auto t0 = Clock::now();
      IsotopyGPU::device_batch_launch(d, count);
      secs.push_back(seconds_since(t0));
    }
    const Stats s = rate_stats(secs, count);
    IsotopyGPU::device_batch_free(d);
    std::printf("delta %d  batch %d  kernel best %.0f /s  mean %.0f /s\n", delta, count, s.best,
                s.mean);
    return 0;
  }

  int flagged = 0;
  const int sample = count < 4096 ? count : 4096;
  const int mismatched = verify_sample(t, b, pool, out, sample, &flagged);
  std::printf("delta %d  batch %d  pool %zu  verified %d  mismatched %d  flagged %d\n", delta,
              count, pool.size(), sample, mismatched, flagged);
  report_rounds(out);

  if (mode == "verify") return (mismatched == 0 && flagged == 0) ? 0 : 1;

  std::vector<double> kernel_secs;
  for (int r = 0; r < reps; ++r) {
    const auto t0 = Clock::now();
    IsotopyGPU::device_batch_launch(d, count);
    kernel_secs.push_back(seconds_since(t0));
  }
  const Stats ks = rate_stats(kernel_secs, count);

  std::vector<double> e2e_secs;
  for (int r = 0; r < reps; ++r) {
    const auto t0 = Clock::now();
    IsotopyGPU::device_batch_upload(d, b.edges.data(), b.signs.data(), count);
    IsotopyGPU::device_batch_launch(d, count);
    IsotopyGPU::device_batch_download(d, out.data(), count);
    e2e_secs.push_back(seconds_since(t0));
  }
  const Stats es = rate_stats(e2e_secs, count);

  IsotopyGPU::device_batch_free(d);

  const double cpu = cpu_baseline(delta, pool, b, 5.0);

  std::printf("  gpu kernel   best %10.0f /s   mean %10.0f /s\n", ks.best, ks.mean);
  std::printf("  gpu e2e      best %10.0f /s   mean %10.0f /s\n", es.best, es.mean);
  std::printf("  cpu 1 core        %10.0f /s   (same unit: Graph + isotopy_type + tree code)\n",
              cpu);
  std::printf("  speedup vs 1 core %10.1f x   vs 16 cores %8.1f x  (kernel, best)\n",
              ks.best / cpu, ks.best / (cpu * 16.0));
  std::printf("METRIC gpu_kernel_delta%d %.0f\n", delta, ks.best);
  std::printf("METRIC gpu_e2e_delta%d %.0f\n", delta, es.best);
  std::printf("METRIC cpu_code_delta%d %.0f\n", delta, cpu);

  return (mismatched == 0 && flagged == 0) ? 0 : 1;
}

#endif
