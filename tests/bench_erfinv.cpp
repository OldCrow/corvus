// Benchmark: corvus::erfinv / corvus::erfcinv vs a SCALAR baseline.
//
// libm has no erfinv/erfcinv (neither C99 nor POSIX defines one), so unlike
// bench_erf/bench_erfc there is no libm scalar to compare against. The
// baseline here is a plain scalar walk of corvus's own SIMD kernel
// (corvus::erfinv/erfcinv called with a length-1 span per element), which
// isolates the SIMD width's contribution from everything else the kernel
// does (table lookups, the dd Halley step, region branches) -- it is NOT an
// independent implementation, just the same code path minus vectorization.
// Same rules as the other bench_* files: manual, quiet machine, Release
// build only, not ctest-registered.
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

void ScalarErfinv(const std::vector<double>& in, std::vector<double>& out) {
  for (size_t i = 0; i < in.size(); ++i) {
    corvus::erfinv(std::span<const double>(&in[i], 1),
                  std::span<double>(&out[i], 1));
  }
}

void ScalarErfcinv(const std::vector<double>& in, std::vector<double>& out) {
  for (size_t i = 0; i < in.size(); ++i) {
    corvus::erfcinv(std::span<const double>(&in[i], 1),
                    std::span<double>(&out[i], 1));
  }
}

void RunErfinv(const char* label, double lo, double hi) {
  std::printf("erfinv, input range [%g, %g] (%s)\n", lo, hi, label);
  std::printf("%10s %14s %14s %10s\n", "n", "simd ns/el", "scalar ns/el",
              "speedup");
  std::mt19937_64 rng(20260725);
  std::uniform_real_distribution<double> dist(lo, hi);
  for (size_t n : {10000UL, 1000000UL}) {
    std::vector<double> in(n), out(n);
    for (auto& v : in) v = dist(rng);
    const int reps = n >= 1000000 ? 11 : 51;
    const double simd = NsPerElement(
        [&] {
          corvus::erfinv(in, out);
          g_sink = out[n / 2];
        },
        n, reps);
    const double scalar = NsPerElement(
        [&] {
          ScalarErfinv(in, out);
          g_sink = out[n / 2];
        },
        n, n >= 1000000 ? 3 : 11);  // scalar path is span-call-per-element, slow
    std::printf("%10zu %14.2f %14.2f %9.2fx\n", n, simd, scalar,
                scalar / simd);
  }
}

void RunErfcinv(const char* label, double lo, double hi, bool log_spaced = false) {
  std::printf("erfcinv, input range [%g, %g] (%s)\n", lo, hi, label);
  std::printf("%10s %14s %14s %10s\n", "n", "simd ns/el", "scalar ns/el",
              "speedup");
  std::mt19937_64 rng(20260725);
  // Uniform-in-exponent sampling for the far-tail bench: a linear
  // uniform_real_distribution over a range spanning hundreds of orders of
  // magnitude would draw almost every sample from the top decade.
  std::uniform_real_distribution<double> dist(lo, hi);
  std::uniform_real_distribution<double> exp_dist(std::log2(lo), std::log2(hi));
  for (size_t n : {10000UL, 1000000UL}) {
    std::vector<double> in(n), out(n);
    for (auto& v : in) {
      v = log_spaced ? std::ldexp(1.0, static_cast<int>(exp_dist(rng)))
                     : dist(rng);
    }
    const int reps = n >= 1000000 ? 11 : 51;
    const double simd = NsPerElement(
        [&] {
          corvus::erfcinv(in, out);
          g_sink = out[n / 2];
        },
        n, reps);
    const double scalar = NsPerElement(
        [&] {
          ScalarErfcinv(in, out);
          g_sink = out[n / 2];
        },
        n, n >= 1000000 ? 3 : 11);
    std::printf("%10zu %14.2f %14.2f %9.2fx\n", n, simd, scalar,
                scalar / simd);
  }
}

}  // namespace

int main() {
  if (!corvus_test::ReportAndCheckTarget()) return 2;
  RunErfinv("central C", -0.5, 0.5);
  RunErfinv("mixed", -0.999, 0.999);
  RunErfcinv("central C", 0.5, 1.5);
  RunErfcinv("T-mid", 1e-6, 0.5);
  RunErfcinv("T-far, subnormal z", 5e-324, 2.2e-308, /*log_spaced=*/true);
  return 0;
}
