// Benchmark: corvus::beta_p_inv / beta_q_inv vs a SCALAR baseline.
//
// libm has no inverse incomplete beta, and neither does anything outside
// SciPy, so as in bench_gamma/bench_beta/bench_gammainv the baseline is a
// plain scalar walk of corvus's own SIMD kernel (length-1 spans). It isolates
// the SIMD width's contribution from everything else the kernel does -- up to
// seven candidate seeds, the residual comparison and four Newton steps, each
// of which is a full region-routed forward -- and is NOT an independent
// implementation, so the ratio is an UPPER BOUND on the speedup a vendor
// scalar inverse would show, not a claim against one.
//
// THE REGIMES HAVE WILDLY DIFFERENT COSTS and a mixed range would describe no
// real workload. The deep-small closed form is one dd division and one
// exponential and skips the entire pipeline (there is an explicit all-lanes
// fast path for it); a moderate interior point pays for five seed families and
// four forwards; a gamma-limit point additionally instantiates gammainv's
// three seeds and beta's gamma-limit slice. Each point set is drawn with
// margin from the routing walls, and because an inverse's regime depends on
// its ANSWER, membership is not self-evident from (a, b, sigma) the way a
// forward's is -- so this bench carries a ROUTER REPLICA: it evaluates the
// kernel once per set and reports the measured regime histogram, which is what
// makes an unexpected timing interpretable rather than mysterious. Span
// lengths are 4001, deliberately not a multiple of AVX3_ZEN4's 8-lane width,
// so every run also exercises the masked-tail path. Same rules as the other
// bench_* files: manual, quiet machine, Release build only, not
// ctest-registered.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <random>
#include <span>
#include <vector>

#include "corvus/corvus.h"
#include "expect_target.h"
#include "src/beta_data.h"
#include "src/betainv_data.h"

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

void ScalarWalkP(const std::vector<double>& a, const std::vector<double>& b,
                 const std::vector<double>& s, std::vector<double>& out) {
  for (size_t i = 0; i < a.size(); ++i) {
    corvus::beta_p_inv(std::span<const double>(&a[i], 1),
                       std::span<const double>(&b[i], 1),
                       std::span<const double>(&s[i], 1),
                       std::span<double>(&out[i], 1));
  }
}

// --- router replica -------------------------------------------------------
// The kernel's own regime predicates, evaluated on the answer the kernel
// actually returned -- including its internal frame (input flip + orientation
// swap), because the deep-small cut is only meaningful there.
enum : int { kDeep = 0, kHuge = 1, kPlateau = 2, kRidge = 3, kGl = 4, kRest = 5,
             kNReg = 6 };
const char* kRegName[kNReg] = {"deep", "hugeNu", "plateau",
                               "ridge", "gammaLim", "rest"};
constexpr double kHugeNu = 1e31;
constexpr double kPlateauMin = 1.1e-16;
constexpr double kGammaLimit = 0x1.0p+20;

int Regime(double a, double b, double s, double x) {
  const bool swap = s > 0.5;  // p side: swap == flip
  const double alpha = swap ? b : a;
  const double beta = swap ? a : b;
  const double y = swap ? 1.0 - x : x;
  if (y == 0.0) return kDeep;
  if (y > 0.0 && y < 1.0) {
    const double corr = y < 1e-8 ? 1.0 : -std::log1p(-y) / y;
    if (std::fabs(1.0 - beta) * y / (1.0 + alpha) * corr <
        corvus::detail::kBetaInvDeepSmallCut) {
      return kDeep;
    }
  }
  const double nu = a * (b / (a + b));
  if (nu >= kHugeNu) return kHuge;
  if (std::fmin(a, b) <= kPlateauMin) return kPlateau;
  if (std::fmax(a, b) >= kGammaLimit) return kGl;
  if (nu >= corvus::detail::kBetaInvS1NuMin) return kRidge;
  return kRest;
}

using Sampler = void (*)(std::mt19937_64&, std::vector<double>&,
                         std::vector<double>&, std::vector<double>&);

double LogU(std::mt19937_64& rng, double lo, double hi) {
  std::uniform_real_distribution<double> u(std::log(lo), std::log(hi));
  return std::exp(u(rng));
}

void SampleDeepSmall(std::mt19937_64& rng, std::vector<double>& a,
                     std::vector<double>& b, std::vector<double>& s) {
  // alpha tiny enough that y collapses far below the cut for any sigma in
  // range: the closed form's own territory, and the one path with an
  // all-lanes fast exit.
  for (size_t i = 0; i < a.size(); ++i) {
    a[i] = LogU(rng, 1e-40, 1e-6);
    b[i] = LogU(rng, 0.5, 20.0);
    s[i] = LogU(rng, 1e-300, 0.4);
  }
}

void SampleR1Tiny(std::mt19937_64& rng, std::vector<double>& a,
                  std::vector<double>& b, std::vector<double>& s) {
  // Small alpha, deep target: S2's Picard series does real work and the
  // forward lands in R1 every step.
  for (size_t i = 0; i < a.size(); ++i) {
    a[i] = LogU(rng, 1e-2, 2.0);
    b[i] = LogU(rng, 0.5, 50.0);
    s[i] = LogU(rng, 1e-200, 1e-4);
  }
}

void SampleModerate(std::mt19937_64& rng, std::vector<double>& a,
                    std::vector<double>& b, std::vector<double>& s) {
  // The weak-seed middle band PLAN's SECOND correction is about: no candidate
  // exceeds a few bits, so every seed is computed and all four steps run.
  std::uniform_real_distribution<double> us(0.05, 0.5);
  for (size_t i = 0; i < a.size(); ++i) {
    a[i] = LogU(rng, 0.02, 0.5);
    b[i] = a[i] * LogU(rng, 3.0, 10.0);
    s[i] = us(rng);
  }
}

void SampleRidge(std::mt19937_64& rng, std::vector<double>& a,
                 std::vector<double>& b, std::vector<double>& s) {
  // nu well above kBetaTRidge with a balanced mean: S1 plus the forward's
  // Temme core and its 2D tensor on every evaluation.
  std::uniform_real_distribution<double> us(0.05, 0.5);
  for (size_t i = 0; i < a.size(); ++i) {
    a[i] = LogU(rng, 4.0 * corvus::detail::kBetaTRidge, 1e7);
    b[i] = a[i] * LogU(rng, 0.5, 2.0);
    s[i] = us(rng);
  }
}

void SampleSkewed(std::mt19937_64& rng, std::vector<double>& a,
                  std::vector<double>& b, std::vector<double>& s) {
  // Strong skew off the ridge band: R2's continued fraction on every step.
  std::uniform_real_distribution<double> us(0.05, 0.5);
  for (size_t i = 0; i < a.size(); ++i) {
    a[i] = LogU(rng, 1.0, 100.0);
    b[i] = a[i] * LogU(rng, 1e3, 1e6);
    s[i] = us(rng);
  }
}

void SampleGammaLimit(std::mt19937_64& rng, std::vector<double>& a,
                      std::vector<double>& b, std::vector<double>& s) {
  // One parameter past B_GL: gammainv's three seeds join the comparison and
  // beta's gamma-limit slice answers the forward.
  std::uniform_real_distribution<double> us(0.05, 0.5);
  for (size_t i = 0; i < a.size(); ++i) {
    a[i] = LogU(rng, 0.5, 30.0);
    b[i] = LogU(rng, 1e18, 1e120);
    s[i] = us(rng);
  }
}

void SampleJointTiny(std::mt19937_64& rng, std::vector<double>& a,
                     std::vector<double>& b, std::vector<double>& s) {
  // Both parameters tiny with an interior target: the plateau band, where the
  // contract is backward error rather than a y-ULP.
  std::uniform_real_distribution<double> us(0.05, 0.5);
  for (size_t i = 0; i < a.size(); ++i) {
    a[i] = LogU(rng, 1e-30, 1e-18);
    b[i] = LogU(rng, 1e-30, 1e-18);
    s[i] = us(rng);
  }
}

void SampleBeyond(std::mt19937_64& rng, std::vector<double>& a,
                  std::vector<double>& b, std::vector<double>& s) {
  // Beyond resolution: the seed is the answer and every step self-freezes.
  std::uniform_real_distribution<double> us(0.05, 0.5);
  for (size_t i = 0; i < a.size(); ++i) {
    a[i] = LogU(rng, 1e36, 1e300);
    b[i] = a[i] * LogU(rng, 0.5, 2.0);
    s[i] = us(rng);
  }
}

void SampleLargeSide(std::mt19937_64& rng, std::vector<double>& a,
                     std::vector<double>& b, std::vector<double>& s) {
  // sigma > 1/2: the exact input flip AND the orientation swap, otherwise the
  // same work as the moderate band.
  std::uniform_real_distribution<double> us(0.5, 0.9999);
  for (size_t i = 0; i < a.size(); ++i) {
    a[i] = LogU(rng, 0.5, 20.0);
    b[i] = LogU(rng, 0.5, 20.0);
    s[i] = us(rng);
  }
}

void Run(const char* label, size_t n, Sampler sampler) {
  std::mt19937_64 rng(20260810);
  std::vector<double> a(n), b(n), s(n), out(n);
  sampler(rng, a, b, s);

  corvus::beta_p_inv(a, b, s, out);
  size_t hist[kNReg] = {0, 0, 0, 0, 0, 0};
  for (size_t i = 0; i < n; ++i) ++hist[Regime(a[i], b[i], s[i], out[i])];

  std::printf("betainv, %s: n=%zu  membership:", label, n);
  for (int k = 0; k < kNReg; ++k) {
    if (hist[k]) {
      std::printf(" %s=%.0f%%", kRegName[k],
                  100.0 * static_cast<double>(hist[k]) /
                      static_cast<double>(n));
    }
  }
  std::printf("\n%10s %14s %14s %14s %10s\n", "n", "p_inv ns/el", "q_inv ns/el",
              "scalar ns/el", "upper bnd");
  const double simd_p = NsPerElement(
      [&] {
        corvus::beta_p_inv(a, b, s, out);
        g_sink = out[n / 2];
      },
      n, 15);
  const double simd_q = NsPerElement(
      [&] {
        corvus::beta_q_inv(a, b, s, out);
        g_sink = out[n / 2];
      },
      n, 15);
  const double scalar = NsPerElement(
      [&] {
        ScalarWalkP(a, b, s, out);
        g_sink = out[n / 2];
      },
      n, 3);
  std::printf("%10zu %14.2f %14.2f %14.2f %9.2fx\n", n, simd_p, simd_q, scalar,
              scalar / simd_p);
}

}  // namespace

int main() {
  // Numbers attributed to the wrong tier are worse than no numbers.
  if (!corvus_test::ReportAndCheckTarget()) return 2;
  constexpr size_t kN = 4001;  // not a multiple of AVX3_ZEN4's 8-lane width
  Run("deep-small closed form", kN, SampleDeepSmall);
  Run("R1-tiny deep target", kN, SampleR1Tiny);
  Run("moderate weak-seed band", kN, SampleModerate);
  Run("ridge (Temme)", kN, SampleRidge);
  Run("skewed (CF)", kN, SampleSkewed);
  Run("gamma-limit transfer", kN, SampleGammaLimit);
  Run("joint-tiny plateau", kN, SampleJointTiny);
  Run("beyond resolution", kN, SampleBeyond);
  Run("large-side input (flip)", kN, SampleLargeSide);
  return 0;
}
