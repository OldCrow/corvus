// Benchmark: corvus::Erf (SIMD batch) vs scalar std::erf, ns per element.
// Not a correctness test; run manually on a quiet machine, Release build only.
//
//   ./build/tests/bench_erf
//
// Repeat under tier caps (CORVUS_DISABLED_TARGETS) for per-tier numbers.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include "corvus/corvus.h"

namespace {

using Clock = std::chrono::steady_clock;

volatile double g_sink;  // defeat dead-code elimination

template <class F>
double NsPerElement(F&& fn, size_t n, int reps) {
  // Median of reps, one untimed warmup.
  fn();
  std::vector<double> times(static_cast<size_t>(reps));
  for (auto& t : times) {
    const auto t0 = Clock::now();
    fn();
    const auto t1 = Clock::now();
    t = std::chrono::duration<double, std::nano>(t1 - t0).count() /
        static_cast<double>(n);
  }
  std::nth_element(times.begin(), times.begin() + reps / 2, times.end());
  return times[static_cast<size_t>(reps) / 2];
}

}  // namespace

int main() {
  std::printf("corvus active SIMD target: %s\n", corvus::ActiveTarget());
  std::printf("%10s %14s %14s %10s\n", "n", "corvus ns/el", "libm ns/el", "speedup");

  std::mt19937_64 rng(20260720);
  std::uniform_real_distribution<double> dist(-6.0, 6.0);

  for (size_t n : {1000UL, 10000UL, 100000UL, 1000000UL}) {
    std::vector<double> in(n), out(n);
    for (auto& v : in) {
      v = dist(rng);
    }
    const int reps = n >= 1000000 ? 11 : 51;

    const double simd = NsPerElement(
        [&] {
          corvus::Erf(in, out);
          g_sink = out[n / 2];
        },
        n, reps);

    const double scalar = NsPerElement(
        [&] {
          for (size_t i = 0; i < n; ++i) {
            out[i] = std::erf(in[i]);
          }
          g_sink = out[n / 2];
        },
        n, reps);

    std::printf("%10zu %14.2f %14.2f %9.2fx\n", n, simd, scalar, scalar / simd);
  }
  return 0;
}
