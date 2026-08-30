// Benchmark: corvus::cos / corvus::sin (SIMD batch) vs scalar std::cos /
// std::sin, ns per element. Not a correctness test; run manually on a quiet
// machine, Release build only.
//
//   ./build/tests/bench_trig
//
// Inputs are uniform in [-1e6, 1e6] — the small-region hot path, not the
// huge-argument reduction path.
// Repeat under tier caps (CORVUS_DISABLED_TARGETS) for per-tier numbers.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include "corvus/corvus.h"
#include "expect_target.h"

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
  // Numbers attributed to the wrong tier are worse than no numbers.
  if (!corvus_test::ReportAndCheckTarget()) return 2;
  std::printf("%10s %4s %14s %14s %10s\n", "n", "fn", "corvus ns/el",
              "libm ns/el", "speedup");

  std::mt19937_64 rng(20260720);
  std::uniform_real_distribution<double> dist(-1e6, 1e6);

  for (size_t n : {1000UL, 10000UL, 100000UL, 1000000UL}) {
    std::vector<double> in(n), out(n);
    for (auto& v : in) {
      v = dist(rng);
    }
    const int reps = n >= 1000000 ? 11 : 51;

    const double simd_cos = NsPerElement(
        [&] {
          corvus::cos(in, out);
          g_sink = out[n / 2];
        },
        n, reps);

    const double scalar_cos = NsPerElement(
        [&] {
          for (size_t i = 0; i < n; ++i) {
            out[i] = std::cos(in[i]);
          }
          g_sink = out[n / 2];
        },
        n, reps);

    std::printf("%10zu %4s %14.2f %14.2f %9.2fx\n", n, "cos", simd_cos,
                scalar_cos, scalar_cos / simd_cos);

    const double simd_sin = NsPerElement(
        [&] {
          corvus::sin(in, out);
          g_sink = out[n / 2];
        },
        n, reps);

    const double scalar_sin = NsPerElement(
        [&] {
          for (size_t i = 0; i < n; ++i) {
            out[i] = std::sin(in[i]);
          }
          g_sink = out[n / 2];
        },
        n, reps);

    std::printf("%10zu %4s %14.2f %14.2f %9.2fx\n", n, "sin", simd_sin,
                scalar_sin, scalar_sin / simd_sin);
  }
  return 0;
}
