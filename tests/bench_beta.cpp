// Benchmark: corvus::beta_p / corvus::beta_q vs a SCALAR baseline.
//
// libm has no incomplete beta (neither C99 nor POSIX defines one), so as in
// bench_gamma/bench_erfinv there is no independent scalar to compare
// against. The baseline is a plain scalar walk of corvus's own SIMD kernel
// (called with length-1 spans), which isolates the SIMD width's
// contribution from everything else the kernel does -- region routing, the
// dd cores, the R3 Clenshaw table. It is NOT an independent implementation.
//
// Ranges are chosen so each run exercises ONE region: the beta router (see
// the Route replica below, mirroring tests/test_beta_ulp.cpp's) has SIX
// destinations with materially different cost -- R1's power series can exit
// in a handful of terms, R3 walks a ten-row Chebyshev table, and the two
// "postroute" flavors (Pr, Gl) fold extra assembly on top of R1/R2 -- so a
// mixed range reports a number that describes no real workload. Each point
// set is built directly from the region's defining inequalities (same
// src/beta_data.h constants the kernel and test_beta_ulp read), with enough
// margin from the routing walls that the construction is self-evidently in
// region; the Route replica then re-checks membership on the actual (a, b,
// x) triples and the hit rate is printed as a diagnostic, never a gate.
// Span lengths are 4001, deliberately not a multiple of AVX3_ZEN4's 8-lane
// width, so every run also exercises the masked-tail path. Same rules as
// the other bench_* files: manual, quiet machine, Release build only, not
// ctest-registered.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include "corvus/corvus.h"
#include "expect_target.h"
#include "src/beta_data.h"

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
                    std::span<const double>, std::span<double>);

void ScalarWalk(Fn fn, const std::vector<double>& a,
                const std::vector<double>& b, const std::vector<double>& x,
                std::vector<double>& out) {
  for (size_t i = 0; i < a.size(); ++i) {
    fn(std::span<const double>(&a[i], 1), std::span<const double>(&b[i], 1),
       std::span<const double>(&x[i], 1), std::span<double>(&out[i], 1));
  }
}

// ---- Router replica -------------------------------------------------------
// Verbatim re-derivation of BetaVec's router, mirroring the one in
// tests/test_beta_ulp.cpp (kept identical on purpose: two independent
// transcriptions of src/beta_data.h drifting apart is exactly the failure
// mode both files exist to catch). Used here ONLY to verify each point set's
// membership as a printed diagnostic -- never to gate the benchmark.
enum : int { kR1 = 0, kR2 = 1, kR3 = 2, kR4 = 3, kSp = 4, kPr = 5, kGl = 6 };

int Route(double a, double b, double x, double pref, double qref,
         bool* direct_is_p) {
  const bool in_domain = a > 0.0 && b > 0.0 && std::isfinite(a) &&
                         std::isfinite(b) && x > 0.0 && x < 1.0;
  if (!in_domain) {
    *direct_is_p = true;
    return kSp;
  }
  using corvus::detail::kBetaB1;
  using corvus::detail::kBetaEpsR4;
  using corvus::detail::kBetaLn2;
  using corvus::detail::kBetaTRidge;
  using corvus::detail::kBetaXi1;
  using corvus::detail::kBetaXiRatioHi;
  using corvus::detail::kBetaXiRatioLo;

  const double y = 1.0 - x;
  const double tau = std::min(a, b);
  const double bmax = std::max(a, b);
  const bool sw4 = b < a;
  const double xt = sw4 ? y : x;

  const double thr_t = 1.0 / (1.0 + (bmax + 1.0) / (tau + 1.0));
  const bool r4 = tau <= kBetaEpsR4 &&
                  tau * std::fabs(std::log(xt)) <= kBetaLn2 &&
                  (xt <= kBetaXi1 || xt < thr_t) && bmax * xt <= kBetaB1;
  if (r4) {
    *direct_is_p = sw4;
    return kR4;
  }
  const bool r1n = x <= kBetaXi1 && b * x <= kBetaB1;
  const bool r1s = y <= kBetaXi1 && a * y <= kBetaB1;
  if (r1n || r1s) {
    const bool sw = !r1n && r1s;
    const double eval_v = sw ? qref : pref;
    const double fired_alpha = sw ? b : a;
    if (eval_v > corvus::detail::kBetaNearOne &&
        fired_alpha <= corvus::detail::kBetaPrTauMax) {
      *direct_is_p = sw;
      return kPr;
    }
    *direct_is_p = !sw;
    return kR1;
  }
  const double nu = 1.0 / (1.0 / a + 1.0 / b);
  const double rat1 = x * (1.0 + b / a);
  const double rat2 = y * (1.0 + a / b);
  const bool band = rat1 >= kBetaXiRatioLo && rat1 <= kBetaXiRatioHi &&
                    rat2 >= kBetaXiRatioLo && rat2 <= kBetaXiRatioHi;
  const bool gl_hi = bmax >= corvus::detail::kBetaGammaLim;
  if (band && (nu >= kBetaTRidge ||
               (gl_hi && nu >= corvus::detail::kBetaGlRidgeMin))) {
    const bool sw = rat1 > 1.0;
    const double lam = a * y - b * x;
    *direct_is_p = sw ? (lam > 0.0) : (lam >= 0.0);
    return kR3;
  }
  const double thr = 1.0 / (1.0 + (b + 1.0) / (a + 1.0));
  const bool sw = !(x < thr);
  if (gl_hi && std::min(a, b) < corvus::detail::kBetaGammaLim) {
    const double ra = sw ? b : a;
    const double rxi = sw ? y : x;
    const bool hf = ra >= corvus::detail::kBetaGammaLim;
    const double s = hf ? (sw ? a : b) : ra;
    const double huge = hf ? ra : (sw ? a : b);
    const double t = -huge * std::log(hf ? rxi : 1.0 - rxi);
    const bool ser = (s < 20.0 && t <= s + 1.0) || (s >= 20.0 && s >= 2.0 * t);
    const bool agree = (ser == !hf);
    *direct_is_p = agree ? !sw : sw;
    return kGl;
  }
  *direct_is_p = !sw;
  return kR2;
}

const char* RegionName(int code) {
  static const char* kNames[] = {"R1", "R2", "R3", "R4", "Sp", "Pr", "Gl"};
  return kNames[code];
}

// ---- per-region point-set samplers ----------------------------------------
// Each fills a, b, x (already sized to n) from the defining inequalities of
// its target region, with margin from the routing walls -- see the header
// comment. `Run` re-checks membership with the Route replica.

void SampleR1(std::mt19937_64& rng, std::vector<double>& a,
              std::vector<double>& b, std::vector<double>& x) {
  // R1 power series, "near" orientation: x well under kBetaXi1 (~0.45) and
  // b*x well under kBetaB1 (8) for a, b in [2, 19].
  std::uniform_real_distribution<double> la(std::log(2.0), std::log(19.0));
  std::uniform_real_distribution<double> xu(0.01, 0.30);
  for (size_t i = 0; i < a.size(); ++i) {
    a[i] = std::exp(la(rng));
    b[i] = std::exp(la(rng));
    x[i] = xu(rng);
  }
}

void SampleR2(std::mt19937_64& rng, std::vector<double>& a,
              std::vector<double>& b, std::vector<double>& x) {
  // R2 continued fraction: a, b in [2, 19] never reach nu >= kBetaTRidge
  // (32) or bmax >= kBetaGammaLim, so the only wall to clear is R1's -- x
  // held in a band straddling 0.5 keeps both x and y = 1-x above kBetaXi1.
  std::uniform_real_distribution<double> la(std::log(2.0), std::log(19.0));
  std::uniform_real_distribution<double> xu(0.46, 0.54);
  for (size_t i = 0; i < a.size(); ++i) {
    a[i] = std::exp(la(rng));
    b[i] = std::exp(la(rng));
    x[i] = xu(rng);
  }
}

void SampleR3(std::mt19937_64& rng, std::vector<double>& a,
              std::vector<double>& b, std::vector<double>& x) {
  // R3 Temme ridge: a == b large (nu = a/2 >= 40, well past kBetaTRidge)
  // with x within +-0.03 of the mean 0.5 -- safely inside the ratio band
  // and safely away from kBetaXi1/1-kBetaXi1 (~0.45/0.55).
  std::uniform_real_distribution<double> lt(std::log(80.0), std::log(5000.0));
  std::uniform_real_distribution<double> du(-0.03, 0.03);
  for (size_t i = 0; i < a.size(); ++i) {
    const double t = std::exp(lt(rng));
    a[i] = t;
    b[i] = t;
    x[i] = 0.5 + du(rng);
  }
}

void SampleR4(std::mt19937_64& rng, std::vector<double>& a,
              std::vector<double>& b, std::vector<double>& x) {
  // R4 tiny-min: tau = a in [1e-4, 1e-2], comfortably under kBetaEpsR4
  // (2^-6 ~ 0.0156); bmax = b in [1, 20] and x in [0.001, 0.3] keep
  // bmax*x <= 6 < kBetaB1 and x <= kBetaXi1 throughout.
  std::uniform_real_distribution<double> lt(std::log(1e-4), std::log(1e-2));
  std::uniform_real_distribution<double> lbmax(std::log(1.0), std::log(20.0));
  std::uniform_real_distribution<double> xu(0.001, 0.3);
  for (size_t i = 0; i < a.size(); ++i) {
    a[i] = std::exp(lt(rng));
    b[i] = std::exp(lbmax(rng));
    x[i] = xu(rng);
  }
}

void SamplePr(std::mt19937_64& rng, std::vector<double>& a,
              std::vector<double>& b, std::vector<double>& x) {
  // R4-postroute: a fixed small (0.5, under kBetaPrTauMax = 2.5) with b
  // large (80-150) concentrates almost all of Beta(a, b)'s mass near 0, so
  // P(a, b, x) is already within a part in 1e3 of 1 by x = kBetaB1/b -- the
  // R1 power-series box's own edge. x is drawn at 85-98% of that edge to
  // land past kBetaNearOne while keeping b*x <= kBetaB1 for r1n to fire.
  std::uniform_real_distribution<double> lb(std::log(80.0), std::log(150.0));
  std::uniform_real_distribution<double> fu(0.85, 0.98);
  for (size_t i = 0; i < a.size(); ++i) {
    a[i] = 0.5;
    b[i] = std::exp(lb(rng));
    x[i] = (8.0 / b[i]) * fu(rng);
  }
}

void SampleGl(std::mt19937_64& rng, std::vector<double>& a,
              std::vector<double>& b, std::vector<double>& x) {
  // R2 gamma-limit: a small-moderate (2-19), b well past kBetaGammaLim
  // (2^59). b*x is then astronomically past kBetaB1, so R1 cannot fire for
  // this orientation; y stays above kBetaXi1 (x <= 0.45), so R1's swapped
  // orientation cannot fire either; b/a is astronomically past the ratio
  // band, so R3 cannot fire. What's left is the (C) gamma-limit CF slice.
  std::uniform_real_distribution<double> la(std::log(2.0), std::log(19.0));
  std::uniform_real_distribution<double> lb(std::log(0x1p59 * 1.1),
                                            std::log(0x1p62));
  std::uniform_real_distribution<double> xu(0.05, 0.45);
  for (size_t i = 0; i < a.size(); ++i) {
    a[i] = std::exp(la(rng));
    b[i] = std::exp(lb(rng));
    x[i] = xu(rng);
  }
}

using Sampler = void (*)(std::mt19937_64&, std::vector<double>&,
                         std::vector<double>&, std::vector<double>&);

void Run(const char* fname, Fn fn, const char* label, int want_region,
         size_t n, Sampler sampler) {
  std::mt19937_64 rng(20260806);
  std::vector<double> a(n), b(n), x(n), out(n);
  sampler(rng, a, b, x);

  // Diagnostic-only membership check against the router replica -- printed,
  // never gated (this is a benchmark, not a correctness test).
  std::vector<double> pref(n), qref(n);
  corvus::beta_p(a, b, x, pref);
  corvus::beta_q(a, b, x, qref);
  size_t hit = 0;
  for (size_t i = 0; i < n; ++i) {
    bool direct_is_p = false;
    if (Route(a[i], b[i], x[i], pref[i], qref[i], &direct_is_p) ==
        want_region) {
      ++hit;
    }
  }
  std::printf("%s, %s: n=%zu, region=%s membership %zu/%zu (%.1f%%)\n", fname,
              label, n, RegionName(want_region), hit, n,
              100.0 * static_cast<double>(hit) / static_cast<double>(n));
  std::printf("%10s %14s %14s %10s\n", "n", "simd ns/el", "scalar ns/el",
              "speedup");
  const int simd_reps = 51;
  const int scalar_reps = 11;
  const double simd = NsPerElement(
      [&] {
        fn(a, b, x, out);
        g_sink = out[n / 2];
      },
      n, simd_reps);
  const double scalar = NsPerElement(
      [&] {
        ScalarWalk(fn, a, b, x, out);
        g_sink = out[n / 2];
      },
      n, scalar_reps);
  std::printf("%10zu %14.2f %14.2f %9.2fx\n", n, simd, scalar, scalar / simd);
}

}  // namespace

int main() {
  if (!corvus_test::ReportAndCheckTarget()) return 2;
  constexpr size_t kN = 4001;  // not a multiple of AVX3_ZEN4's 8-lane width
  Run("beta_p", corvus::beta_p, "R1 series", kR1, kN, SampleR1);
  Run("beta_q", corvus::beta_q, "R2 backward CF", kR2, kN, SampleR2);
  Run("beta_p", corvus::beta_p, "R3 Temme ridge", kR3, kN, SampleR3);
  Run("beta_q", corvus::beta_q, "R4 tiny-min", kR4, kN, SampleR4);
  Run("beta_p", corvus::beta_p, "R4 postroute near-one", kPr, kN, SamplePr);
  Run("beta_q", corvus::beta_q, "R2 gamma-limit", kGl, kN, SampleGl);
  return 0;
}
