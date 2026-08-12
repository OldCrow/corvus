// Benchmark: corvus::lbeta vs a SCALAR baseline (a plain scalar walk of the
// same SIMD kernel, length-1 spans -- see bench_bessel.cpp's header for why
// the ratio is an UPPER BOUND on any claim against a vendor scalar). Two
// samplers: the main band (the LgammaPosDd/LgammaDiffDd assembly, by far
// the common case) and the big band (min > 2^990, the Stirling-direct
// path). Span length 4001 (not a lane multiple; the masked tail runs).
// Manual, quiet machine, Release only, not ctest-registered.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <random>
#include <span>
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

void FillLog(std::mt19937_64& rng, std::vector<double>& v, double lo,
             double hi) {
  std::uniform_real_distribution<double> le(std::log(lo), std::log(hi));
  for (auto& x : v) x = std::exp(le(rng));
}

}  // namespace

int main() {
  if (!corvus_test::ReportAndCheckTarget()) return 2;

  constexpr size_t kN = 4001;
  constexpr int kReps = 41;
  std::mt19937_64 rng(20260811);

  std::vector<double> a(kN), b(kN), out(kN);

  struct Band {
    const char* name;
    double lo, hi;
  };
  for (const Band& band : {Band{"main band", 1e-3, 1e6},
                           Band{"big band ", 0x1.1p+990, 0x1.0p+1000}}) {
    FillLog(rng, a, band.lo, band.hi);
    FillLog(rng, b, band.lo, band.hi);

    const double simd = NsPerElement(
        [&] {
          corvus::lbeta(a, b, out);
          g_sink = out[0];
        },
        kN, kReps);
    const double scalar = NsPerElement(
        [&] {
          for (size_t i = 0; i < kN; ++i) {
            corvus::lbeta(std::span<const double>(&a[i], 1),
                          std::span<const double>(&b[i], 1),
                          std::span<double>(&out[i], 1));
          }
          g_sink = out[0];
        },
        kN, kReps);
    std::printf("lbeta %s: %7.2f ns/el SIMD, %7.2f ns/el scalar-walk, "
                "%.2fx\n",
                band.name, simd, scalar, scalar / simd);
  }
  return 0;
}
