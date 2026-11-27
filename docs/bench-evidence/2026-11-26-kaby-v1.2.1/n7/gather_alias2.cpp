// Standalone probe 2 (K, v1.2.1): the window of gather_alias.cpp at 1-byte
// resolution, for several positions of the gathered table within its page
// and for 32- and 64-bit-index gathers of doubles. No corvus code.
//   c++ -O2 -mavx2 -mfma gather_alias2.cpp -o gather_alias2
// Output: kind table_page_offset probe_page_offset ns (min of 7 reps)
#include <immintrin.h>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>

alignas(4096) static unsigned char g_probe[8192];
alignas(4096) static unsigned char g_tab[8192];

template <int kKind>
static double Run(const double* table, const unsigned char* pp) {
  using Clock = std::chrono::steady_clock;
  constexpr int kIters = 200000, kReps = 7;
  const __m256i lane = _mm256_set_epi64x(3, 2, 1, 0);
  const __m128i lane32 = _mm_set_epi32(3, 2, 1, 0);
  double best = 1e30;
  for (int r = 0; r < kReps; ++r) {
    __m256i idx = lane;
    __m128i idx32 = lane32;
    double acc = 0.0;
    const auto t0 = Clock::now();
    for (int i = 0; i < kIters; ++i) {
      __m256d g;
      if (kKind == 0) g = _mm256_i64gather_pd(table, idx, 8);
      else if (kKind == 1) g = _mm256_i32gather_pd(table, idx32, 8);
      else g = _mm256_loadu_pd(table);
      double c;
      const unsigned char* q = pp;
      __asm__ volatile("" : "+r"(q));          // the address is opaque
      std::memcpy(&c, q, 8);                   // the probed 8-byte load
      const double s = _mm256_cvtsd_f64(g) + c;
      acc += s;
      const long long z = static_cast<long long>(s) - 1;
      idx = _mm256_add_epi64(lane, _mm256_set1_epi64x(z));
      idx32 = _mm_add_epi32(lane32, _mm_set1_epi32(static_cast<int>(z)));
    }
    const auto t1 = Clock::now();
    if (acc != kIters) std::abort();
    best = std::min(best, std::chrono::duration<double, std::nano>(t1 - t0)
                                  .count() / kIters);
  }
  return best;
}

int main() {
  std::memset(g_probe, 0, sizeof g_probe);
  const char* names[] = {"gatherq", "gatherd", "load"};
  for (int toff : {0x000, 0x340, 0x7c0, 0xa80, 0xff0}) {
    double* table = reinterpret_cast<double*>(g_tab + toff);
    for (int i = 0; i < 8; ++i) table[i] = 1.0;
    for (int kind = 0; kind < 3; ++kind) {
      for (int po = 0x780; po < 0x820; ++po) {
        const double ns = kind == 0 ? Run<0>(table, g_probe + po)
                        : kind == 1 ? Run<1>(table, g_probe + po)
                                    : Run<2>(table, g_probe + po);
        std::printf("%s %#05x %#05x %7.2f\n", names[kind], toff, po, ns);
      }
    }
    for (int i = 0; i < 8; ++i) table[i] = 0.0;
  }
  return 0;
}
