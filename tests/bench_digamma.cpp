// Benchmark: corvus::digamma vs a SCALAR baseline.
//
// libm has no digamma (neither C99 nor POSIX defines one), so as in
// bench_gamma/bench_beta there is no independent scalar to compare against.
// The baseline is a plain scalar walk of corvus's own SIMD kernel (called
// with length-1 spans), which isolates the SIMD width's contribution from
// everything else the kernel does -- region routing, the dd cores, log_dd.
// It is NOT an independent implementation, so the ratio is an UPPER BOUND on
// the speedup a vendor scalar digamma would show, not a claim against one.
//
// Ranges are chosen so each run exercises ONE region: the four positive
// branches have materially different cost -- the zone is one Horner, the
// walk adds up to six dd reciprocals, the asymptotic pays for log_dd, and the
// negative axis pays for the whole positive pipeline PLUS the sinc/cos pair
// and a dd division -- so a mixed range reports a number that describes no
// real workload. Each point set is built directly from the region's defining
// inequalities (the same src/digamma_data.h constants the kernel reads) with
// margin from the routing walls, so membership is self-evident and needs no
// router replica. Span lengths are 4001, deliberately not a multiple of
// AVX3_ZEN4's 8-lane width, so every run also exercises the masked-tail path.
// Same rules as the other bench_* files: manual, quiet machine, Release build
// only, not ctest-registered.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include "corvus/corvus.h"
#include "expect_target.h"
#include "src/digamma_data.h"

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

void ScalarWalk(const std::vector<double>& in, std::vector<double>& out) {
  for (size_t i = 0; i < in.size(); ++i) {
    corvus::digamma(std::span<const double>(&in[i], 1),
                    std::span<double>(&out[i], 1));
  }
}

// Each sampler fills `in` (already sized) from its region's own bounds.
using Sampler = void (*)(std::mt19937_64&, std::vector<double>&);

void SampleUpStep(std::mt19937_64& rng, std::vector<double>& in) {
  // (0, kDigammaZoneLo): the one up-step, log-spaced so the -1/x limb is
  // exercised across many decades rather than only near 1.
  std::uniform_real_distribution<double> le(std::log(1e-6), std::log(0.99));
  for (auto& v : in) v = std::exp(le(rng));
}

void SampleZone(std::mt19937_64& rng, std::vector<double>& in) {
  // [kDigammaZoneLo, kDigammaZoneHi): the product-form zone, no walk steps.
  std::uniform_real_distribution<double> u(corvus::detail::kDigammaZoneLo,
                                           corvus::detail::kDigammaZoneHi);
  for (auto& v : in) v = u(rng);
}

void SampleWalk(std::mt19937_64& rng, std::vector<double>& in) {
  // [kDigammaZoneHi, kDigammaX0): the down-walk, uniform so all six depths
  // appear and every vector is mixed-depth (which is the realistic case, and
  // the one the fire masks cost the most on).
  std::uniform_real_distribution<double> u(corvus::detail::kDigammaZoneHi,
                                           corvus::detail::kDigammaX0);
  for (auto& v : in) v = u(rng);
}

void SampleAsym(std::mt19937_64& rng, std::vector<double>& in) {
  // [kDigammaX0, inf): log-spaced well inside the correction-carrying part of
  // the branch (far below the log-only cut at 2^85, which is a cheaper path
  // and would flatter the number).
  std::uniform_real_distribution<double> le(std::log(corvus::detail::kDigammaX0),
                                            std::log(1e6));
  for (auto& v : in) v = std::exp(le(rng));
}

void SampleNeg(std::mt19937_64& rng, std::vector<double>& in) {
  // The reflection: uniform on (-30, 0), rejecting the immediate
  // neighbourhood of each pole so no lane takes the driver's direct -1/u
  // shortcut (which would be a different, much cheaper measurement).
  std::uniform_real_distribution<double> u(-30.0, -0.01);
  for (auto& v : in) {
    double x;
    do {
      x = u(rng);
    } while (std::fabs(x - std::nearbyint(x)) < 0.02);
    v = x;
  }
}

void Run(const char* label, size_t n, Sampler sampler) {
  std::mt19937_64 rng(20260807);
  std::vector<double> in(n), out(n);
  sampler(rng, in);

  std::printf("digamma, %s: n=%zu\n", label, n);
  std::printf("%10s %14s %14s %10s\n", "n", "simd ns/el", "scalar ns/el",
              "upper bnd");
  const double simd = NsPerElement(
      [&] {
        corvus::digamma(in, out);
        g_sink = out[n / 2];
      },
      n, 51);
  const double scalar = NsPerElement(
      [&] {
        ScalarWalk(in, out);
        g_sink = out[n / 2];
      },
      n, 11);
  std::printf("%10zu %14.2f %14.2f %9.2fx\n", n, simd, scalar, scalar / simd);
}

}  // namespace

int main() {
  // Numbers attributed to the wrong tier are worse than no numbers.
  if (!corvus_test::ReportAndCheckTarget()) return 2;
  constexpr size_t kN = 4001;  // not a multiple of AVX3_ZEN4's 8-lane width
  Run("(0,1) up-step", kN, SampleUpStep);
  Run("[1,2) zone", kN, SampleZone);
  Run("[2,8) down-walk", kN, SampleWalk);
  Run("[8,inf) asymptotic", kN, SampleAsym);
  Run("negative axis (reflection)", kN, SampleNeg);
  return 0;
}
