// Benchmark: corvus::erfc vs scalar std::erfc, ns per element. Same rules
// as bench_erf: manual, quiet machine, Release build only.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include "corvus/corvus.h"

namespace {

using Clock = std::chrono::steady_clock;

volatile double g_sink;

template <class F>
double NsPerElement(F&& fn, size_t n, int reps) {
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

void Run(const char* label, double lo, double hi) {
  std::printf("input range [%g, %g] (%s)\n", lo, hi, label);
  std::printf("%10s %14s %14s %10s\n", "n", "corvus ns/el", "libm ns/el", "speedup");
  std::mt19937_64 rng(20260721);
  std::uniform_real_distribution<double> dist(lo, hi);
  for (size_t n : {10000UL, 1000000UL}) {
    std::vector<double> in(n), out(n);
    for (auto& v : in) {
      v = dist(rng);
    }
    const int reps = n >= 1000000 ? 11 : 51;
    const double simd = NsPerElement(
        [&] {
          corvus::erfc(in, out);
          g_sink = out[n / 2];
        },
        n, reps);
    const double scalar = NsPerElement(
        [&] {
          for (size_t i = 0; i < n; ++i) {
            out[i] = std::erfc(in[i]);
          }
          g_sink = out[n / 2];
        },
        n, reps);
    std::printf("%10zu %14.2f %14.2f %9.2fx\n", n, simd, scalar, scalar / simd);
  }
}

}  // namespace

int main() {
  std::printf("corvus active SIMD target: %s\n", corvus::active_target());
  Run("core-dominated", -6.0, 6.0);
  Run("mixed", -6.5, 28.0);
  Run("tail-only", 6.0, 28.0);
  return 0;
}
