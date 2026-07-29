// Benchmark: corvus::gamma_p / gamma_q vs a SCALAR baseline.
//
// libm has no incomplete gamma (neither C99 nor POSIX defines one), so as in
// bench_erfinv there is no independent scalar to compare against. The
// baseline is a plain scalar walk of corvus's own SIMD kernel (called with
// length-1 spans), which isolates the SIMD width's contribution from
// everything else the kernel does -- region routing, the dd cores, the
// Temme table's Clenshaw passes. It is NOT an independent implementation.
//
// Ranges are chosen so each run exercises ONE region: the four cores differ
// by more than an order of magnitude in cost (R3 evaluates eleven Chebyshev
// rows; R1 can exit after a handful of terms), so a mixed range reports a
// number that describes no real workload. Same rules as the other bench_*
// files: manual, quiet machine, Release build only, not ctest-registered.
#include <algorithm>
#include <chrono>
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

using Fn = void (*)(std::span<const double>, std::span<const double>,
                    std::span<double>);

void ScalarWalk(Fn fn, const std::vector<double>& a,
                const std::vector<double>& x, std::vector<double>& out) {
  for (size_t i = 0; i < a.size(); ++i) {
    fn(std::span<const double>(&a[i], 1), std::span<const double>(&x[i], 1),
       std::span<double>(&out[i], 1));
  }
}

// lam_lo/lam_hi bracket x/a; a is drawn log-uniformly in [a_lo, a_hi].
void Run(const char* label, Fn fn, const char* fname, double a_lo, double a_hi,
         double lam_lo, double lam_hi) {
  std::printf("%s, %s: a in [%g, %g], x/a in [%g, %g]\n", fname, label, a_lo,
              a_hi, lam_lo, lam_hi);
  std::printf("%10s %14s %14s %10s\n", "n", "simd ns/el", "scalar ns/el",
              "speedup");
  std::mt19937_64 rng(20260727);
  std::uniform_real_distribution<double> la(std::log(a_lo), std::log(a_hi));
  std::uniform_real_distribution<double> ll(lam_lo, lam_hi);
  for (size_t n : {10000UL, 1000000UL}) {
    std::vector<double> a(n), x(n), out(n);
    for (size_t i = 0; i < n; ++i) {
      a[i] = std::exp(la(rng));
      x[i] = a[i] * ll(rng);
    }
    const int reps = n >= 1000000 ? 11 : 51;
    const double simd = NsPerElement(
        [&] {
          fn(a, x, out);
          g_sink = out[n / 2];
        },
        n, reps);
    const double scalar = NsPerElement(
        [&] {
          ScalarWalk(fn, a, x, out);
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
  // R1: small a with x <= a+1 (lambda <= 1 keeps it inside for a >= 1).
  Run("R1 series", corvus::gamma_p, "gamma_p", 2.0, 19.0, 0.2, 1.0);
  // R2: small a well past x = a+1.
  Run("R2 backward CF", corvus::gamma_q, "gamma_q", 2.0, 19.0, 3.0, 20.0);
  // R3: the ridge, a >= 20 and lambda strictly inside (1/2, 2).
  Run("R3 Temme", corvus::gamma_p, "gamma_p", 50.0, 1e6, 0.7, 1.4);
  // R4: the small-a box, where gamma_q owns the whole thing.
  Run("R4 small-a", corvus::gamma_q, "gamma_q", 0.05, 1.5, 0.5, 2.5);
  return 0;
}
