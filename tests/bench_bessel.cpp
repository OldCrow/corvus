// Benchmark: corvus::i0/i1/i0e/i1e vs a SCALAR baseline.
//
// libm's bessel_i0/i1 (C++17 std::cyl_bessel_i-family aside, which is
// unscaled-only and not implemented everywhere) is not a drop-in comparison
// for a SIMD batch API, so as in bench_trigamma/bench_digamma the baseline
// is a plain scalar walk of corvus's own SIMD kernel (length-1 spans). It
// isolates the SIMD width's contribution from everything else the kernel
// does -- region routing, the dd series/tail cores, the exp_dd folds -- and
// is NOT an independent implementation, so the ratio is an UPPER BOUND on
// the speedup a vendor scalar Bessel would show, not a claim against one.
//
// Each region is sampled separately: series and tail have materially
// different cost (the series pays for the exact-q residual and a dd-lead
// Horner; the tail pays for a dd sqrt+reciprocal on top of its own Horner),
// and a mixed range would describe no real workload. Span lengths are 4001,
// deliberately not a multiple of AVX3_ZEN4's 8-lane width, so every run
// also exercises the masked-tail path. Same rules as the other bench_*
// files: manual, quiet machine, Release build only, not ctest-registered.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <random>
#include <span>
#include <vector>

#include "corvus/corvus.h"
#include "expect_target.h"
#include "src/bessel_data.h"

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

using Fn = void (*)(std::span<const double>, std::span<double>);

void ScalarWalk(Fn fn, const std::vector<double>& in,
                std::vector<double>& out) {
  for (size_t i = 0; i < in.size(); ++i) {
    fn(std::span<const double>(&in[i], 1), std::span<double>(&out[i], 1));
  }
}

// Each sampler fills `in` (already sized) from its region's own bounds.
using Sampler = void (*)(std::mt19937_64&, std::vector<double>&);

void SampleSeries(std::mt19937_64& rng, std::vector<double>& in) {
  // (0, kBesselSplit]: log-spaced so the exact-q residual and the dd-lead
  // Horner are exercised across many decades, not just near the split.
  std::uniform_real_distribution<double> le(std::log(1e-6),
                                            std::log(corvus::detail::kBesselSplit));
  for (auto& v : in) v = std::exp(le(rng));
}

void SampleTail(std::mt19937_64& rng, std::vector<double>& in) {
  // (kBesselSplit, 700]: well inside the tail, away from the overflow
  // boundary (~714), log-spaced.
  std::uniform_real_distribution<double> le(
      std::log(corvus::detail::kBesselSplit), std::log(700.0));
  for (auto& v : in) v = std::exp(le(rng));
}

void SampleMixedSign(std::mt19937_64& rng, std::vector<double>& in) {
  // Both regions, both signs, uniform in log-magnitude -- the realistic
  // case, and the one the routing selects cost the most on.
  std::uniform_real_distribution<double> le(std::log(1e-3), std::log(700.0));
  std::uniform_real_distribution<double> sign(0.0, 1.0);
  for (auto& v : in) {
    v = std::exp(le(rng));
    if (sign(rng) < 0.5) v = -v;
  }
}

void Run(const char* label, size_t n, Sampler sampler) {
  std::mt19937_64 rng(20260811);
  std::vector<double> in(n), out(n);
  sampler(rng, in);

  std::printf("bessel, %s: n=%zu\n", label, n);
  std::printf("%8s %14s %14s %10s\n", "fn", "simd ns/el", "scalar ns/el",
              "upper bnd");
  const struct {
    const char* name;
    Fn fn;
  } kFns[] = {
      {"i0", corvus::i0},
      {"i1", corvus::i1},
      {"i0e", corvus::i0e},
      {"i1e", corvus::i1e},
  };
  for (const auto& f : kFns) {
    const double simd = NsPerElement(
        [&] {
          f.fn(in, out);
          g_sink = out[n / 2];
        },
        n, 21);
    const double scalar = NsPerElement(
        [&] {
          ScalarWalk(f.fn, in, out);
          g_sink = out[n / 2];
        },
        n, 5);
    std::printf("%8s %14.2f %14.2f %9.2fx\n", f.name, simd, scalar,
               scalar / simd);
  }
}

}  // namespace

int main() {
  // Numbers attributed to the wrong tier are worse than no numbers.
  if (!corvus_test::ReportAndCheckTarget()) return 2;
  constexpr size_t kN = 4001;  // not a multiple of AVX3_ZEN4's 8-lane width
  Run("series (0,8]", kN, SampleSeries);
  Run("tail (8,700]", kN, SampleTail);
  Run("mixed sign, both regions", kN, SampleMixedSign);
  return 0;
}
