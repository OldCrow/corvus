// Standalone probe (K, v1.2.1): no corvus code. Does an 8-byte load cost more
// when it FOLLOWS an AVX2 gather and its address has a particular offset
// within its 4 KB page? Each iteration: one vgatherqpd from a small table,
// then one 8-byte load from the probe address; both feed the next
// iteration's gather index, so the loop runs at the latency of that path.
// The probe address is swept over two pages in 8-byte steps.
//   c++ -O2 -mavx2 -mfma gather_alias.cpp -o gather_alias
// usage: gather_alias [nogather]   ("nogather" replaces the gather by a
//        plain 32-byte load of the same table, as a control)
// Output: page_offset ns_per_iteration (min of 7 repetitions)
#include <immintrin.h>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>

alignas(4096) static double g_probe[1024];
alignas(64) static double g_table[64];

int main(int argc, char** argv) {
  const bool gather = !(argc > 1 && std::strcmp(argv[1], "nogather") == 0);
  for (double& v : g_probe) v = 0.0;
  for (int i = 0; i < 64; ++i) g_table[i] = 1.0;
  using Clock = std::chrono::steady_clock;
  constexpr int kIters = 200000, kReps = 7;
  const __m256i lane = _mm256_set_epi64x(3, 2, 1, 0);
  for (int off = 0; off < 1024; ++off) {
    const volatile double* p = &g_probe[off];
    double best = 1e30;
    for (int r = 0; r < kReps; ++r) {
      __m256i idx = lane;
      double acc = 0.0;
      const auto t0 = Clock::now();
      for (int i = 0; i < kIters; ++i) {
        const __m256d g = gather ? _mm256_i64gather_pd(g_table, idx, 8)
                                 : _mm256_loadu_pd(g_table);
        const double c = *p;  // the probed load, after the gather
        const double s = _mm256_cvtsd_f64(g) + c;  // 1.0 + 0.0
        acc += s;
        // next index depends on both the gather and the probed load
        const long long z = static_cast<long long>(s) - 1;  // 0
        idx = _mm256_add_epi64(lane, _mm256_set1_epi64x(z));
      }
      const auto t1 = Clock::now();
      if (acc != kIters) std::abort();
      best = std::min(best, std::chrono::duration<double, std::nano>(t1 - t0)
                                    .count() / kIters);
    }
    std::printf("%4d %7.2f\n", (off * 8) & 4095, best);
  }
  return 0;
}
