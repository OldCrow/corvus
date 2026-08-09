// Benchmark: corvus::gamma_p_inv / gamma_q_inv vs a SCALAR baseline.
//
// libm has no inverse incomplete gamma, and neither does anything else outside
// SciPy, so as in bench_gamma/bench_beta the baseline is a plain scalar walk of
// corvus's own SIMD kernel (length-1 spans). It isolates the SIMD width's
// contribution from everything else the kernel does -- three candidate seeds,
// the residual comparison, five forward evaluations -- and is NOT an
// independent implementation, so the ratio is an UPPER BOUND on the speedup a
// vendor scalar inverse would show, not a claim against one.
//
// The regimes have wildly different costs and a mixed range would describe no
// real workload: the deep-small closed form is one dd division and one
// exponential and skips the entire pipeline (there is an explicit all-lanes
// fast path for it), while a mid-band point at a < a_T pays for all three
// candidate seeds -- including erfcinv, six Picard corrections each with their
// own series, and a fixed-point iteration -- before its three Newton steps
// start. Each point set is built from the kernel's own pinned constants in
// src/gammainv_data.h with margin from the routing walls, so membership is
// self-evident and needs no router replica. Span lengths are 4001,
// deliberately not a multiple of AVX3_ZEN4's 8-lane width, so every run also
// exercises the masked-tail path. Same rules as the other bench_* files:
// manual, quiet machine, Release build only, not ctest-registered.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <random>
#include <span>
#include <vector>

#include "corvus/corvus.h"
#include "expect_target.h"
#include "src/gamma_data.h"
#include "src/gammainv_data.h"

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

void ScalarWalkP(const std::vector<double>& a, const std::vector<double>& s,
                 std::vector<double>& out) {
  for (size_t i = 0; i < a.size(); ++i) {
    corvus::gamma_p_inv(std::span<const double>(&a[i], 1),
                        std::span<const double>(&s[i], 1),
                        std::span<double>(&out[i], 1));
  }
}

using Sampler = void (*)(std::mt19937_64&, std::vector<double>&,
                         std::vector<double>&);

void SampleDeepSmall(std::mt19937_64& rng, std::vector<double>& a,
                     std::vector<double>& s) {
  // The closed form's own territory: a tiny enough that x collapses far below
  // the cut for any p in range.
  std::uniform_real_distribution<double> la(std::log(1e-40), std::log(1e-6));
  std::uniform_real_distribution<double> ls(std::log(1e-300), std::log(0.4));
  for (size_t i = 0; i < a.size(); ++i) {
    a[i] = std::exp(la(rng));
    s[i] = std::exp(ls(rng));
  }
}

void SampleSmallMid(std::mt19937_64& rng, std::vector<double>& a,
                    std::vector<double>& s) {
  // a < a_T with a target near the median: the weak-seed middle band, where
  // all three candidates are computed and compared.
  std::uniform_real_distribution<double> ua(0.5,
                                            corvus::detail::kGammaAT - 1.0);
  std::uniform_real_distribution<double> us(0.05, 0.5);
  for (size_t i = 0; i < a.size(); ++i) {
    a[i] = ua(rng);
    s[i] = us(rng);
  }
}

void SampleSmallTail(std::mt19937_64& rng, std::vector<double>& a,
                     std::vector<double>& s) {
  // a < a_T, target well below the shallow threshold: S2/S3's genuine tail.
  std::uniform_real_distribution<double> ua(0.5,
                                            corvus::detail::kGammaAT - 1.0);
  std::uniform_real_distribution<double> ls(
      std::log(1e-200),
      std::log(corvus::detail::kGammaInvShallowThreshold * 0.5));
  for (size_t i = 0; i < a.size(); ++i) {
    a[i] = ua(rng);
    s[i] = std::exp(ls(rng));
  }
}

void SampleRidge(std::mt19937_64& rng, std::vector<double>& a,
                 std::vector<double>& s) {
  // a >= a_T near the median: S1 plus the forward's Temme core on every step.
  std::uniform_real_distribution<double> la(std::log(corvus::detail::kGammaAT),
                                            std::log(1e8));
  std::uniform_real_distribution<double> us(0.05, 0.5);
  for (size_t i = 0; i < a.size(); ++i) {
    a[i] = std::exp(la(rng));
    s[i] = us(rng);
  }
}

void SampleBigTail(std::mt19937_64& rng, std::vector<double>& a,
                   std::vector<double>& s) {
  // a >= a_T with a deep target: the forward lands in R1/R2 rather than on the
  // ridge, and the seed comparison has real work to do.
  std::uniform_real_distribution<double> la(std::log(corvus::detail::kGammaAT),
                                            std::log(1e6));
  std::uniform_real_distribution<double> ls(std::log(1e-250), std::log(1e-6));
  for (size_t i = 0; i < a.size(); ++i) {
    a[i] = std::exp(la(rng));
    s[i] = std::exp(ls(rng));
  }
}

void SampleBeyond(std::mt19937_64& rng, std::vector<double>& a,
                  std::vector<double>& s) {
  // Beyond resolution: the seed is the answer and every step self-freezes.
  std::uniform_real_distribution<double> la(std::log(1e40), std::log(1e300));
  std::uniform_real_distribution<double> us(0.05, 0.5);
  for (size_t i = 0; i < a.size(); ++i) {
    a[i] = std::exp(la(rng));
    s[i] = us(rng);
  }
}

void SampleLargeSide(std::mt19937_64& rng, std::vector<double>& a,
                     std::vector<double>& s) {
  // s > 1/2: the exact input flip, otherwise the same work as small-mid.
  std::uniform_real_distribution<double> ua(0.5,
                                            corvus::detail::kGammaAT - 1.0);
  std::uniform_real_distribution<double> us(0.5, 0.9999);
  for (size_t i = 0; i < a.size(); ++i) {
    a[i] = ua(rng);
    s[i] = us(rng);
  }
}

void Run(const char* label, size_t n, Sampler sampler) {
  std::mt19937_64 rng(20260809);
  std::vector<double> a(n), s(n), out(n);
  sampler(rng, a, s);

  std::printf("gammainv, %s: n=%zu\n", label, n);
  std::printf("%10s %14s %14s %14s %10s\n", "n", "p_inv ns/el", "q_inv ns/el",
              "scalar ns/el", "upper bnd");
  const double simd_p = NsPerElement(
      [&] {
        corvus::gamma_p_inv(a, s, out);
        g_sink = out[n / 2];
      },
      n, 21);
  const double simd_q = NsPerElement(
      [&] {
        corvus::gamma_q_inv(a, s, out);
        g_sink = out[n / 2];
      },
      n, 21);
  const double scalar = NsPerElement(
      [&] {
        ScalarWalkP(a, s, out);
        g_sink = out[n / 2];
      },
      n, 5);
  std::printf("%10zu %14.2f %14.2f %14.2f %9.2fx\n", n, simd_p, simd_q, scalar,
              scalar / simd_p);
}

}  // namespace

int main() {
  // Numbers attributed to the wrong tier are worse than no numbers.
  if (!corvus_test::ReportAndCheckTarget()) return 2;
  constexpr size_t kN = 4001;  // not a multiple of AVX3_ZEN4's 8-lane width
  Run("deep-small closed form", kN, SampleDeepSmall);
  Run("small-a mid band", kN, SampleSmallMid);
  Run("small-a deep tail", kN, SampleSmallTail);
  Run("a >= a_T ridge", kN, SampleRidge);
  Run("a >= a_T deep tail", kN, SampleBigTail);
  Run("beyond resolution", kN, SampleBeyond);
  Run("large-side input (flip)", kN, SampleLargeSide);
  return 0;
}
