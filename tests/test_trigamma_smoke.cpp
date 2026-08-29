// Smoke test for corvus::trigamma: the full specials table, a handful of
// closed-form values, and lane-mix determinism.
//
// THE SPECIALS TABLE IS NOT DIGAMMA'S, AND THE DIFFERENCE IS THE POINT. Every
// pole of psi_1 is a DOUBLE pole -- psi_1(x) ~ 1/u^2 near it -- so it has no
// sign ambiguity and the answer is +inf on both sides. digamma returns NaN at
// its negative-integer poles precisely because a simple pole approaches -inf
// from one side and +inf from the other; that reasoning does not transfer, and
// scipy's polygamma(1, .) agrees. Likewise psi_1(+-0) = +inf on both signs of
// zero (digamma's signed-zero convention has nothing to distinguish here), and
// psi_1(+inf) = +0 rather than digamma's +inf.
//
// Lane-mix determinism is why this file exists separately from the ULP gate.
// The [2, 8) down-walk stops per lane on a fire mask and lets the whole vector
// break out once no lane fires, so a lane's iteration count is decided by its
// NEIGHBOURS; the region cores additionally run on all lanes with the inactive
// ones scrubbed or clamped. If a frozen accumulator were merely having zeros
// added to it rather than being left alone, the extra dd renormalizations
// would show in the last bits -- and only ever in mixed vectors, which a
// single-point test never builds. So every probe here is evaluated alone, at
// several lane offsets, and interleaved with points from other regions, and
// all of it must be bit-identical.
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <vector>

#include "corvus/corvus.h"
#include "expect_target.h"
#include "ulp_utils.h"

namespace {

using corvus_test::SameBits;
using corvus_test::UlpDiff;

const double kInf = std::numeric_limits<double>::infinity();

int g_fail = 0;

double One(double x) {
  double out = 0.0;
  corvus::trigamma(std::span<const double>(&x, 1), std::span<double>(&out, 1));
  return out;
}

// A NaN expectation means "any NaN". Everything else is matched WITH its sign,
// which is the whole content of the +0 at +inf and the +inf at both zeros.
void CheckSpecial(const char* what, double x, double want) {
  const double got = One(x);
  bool ok;
  if (std::isnan(want)) {
    ok = std::isnan(got);
  } else {
    ok = got == want && std::signbit(got) == std::signbit(want);
  }
  if (!ok) {
    std::fprintf(stderr, "FAIL: %s: trigamma(%.17g) = %.17g, want %.17g\n",
                 what, x, got, want);
    g_fail = 1;
  }
}

void Specials() {
  // The double pole at zero: +inf from BOTH signs of zero.
  CheckSpecial("psi_1(+0)", 0.0, kInf);
  CheckSpecial("psi_1(-0)", -0.0, kInf);

  // Every negative integer is a double pole -> +inf (NOT digamma's NaN). The
  // last four carry the rule that every double <= -2^53 is an integer, so the
  // whole far negative axis is poles; -0x1p53 itself is the boundary case.
  const double kNegInts[] = {-1.0,     -2.0,    -3.0,          -17.0,
                             -1.0e6,   -1.0e15, -0x1p53,       -0x1p53 - 2.0,
                             -1.0e300, -1.7976931348623157e308};
  for (double v : kNegInts) CheckSpecial("negative integer", v, kInf);

  CheckSpecial("psi_1(+inf)", kInf, 0.0);   // +0, not merely zero
  CheckSpecial("psi_1(-inf)", -kInf, kInf);  // scipy convention

  // Arguments small enough that 1/x^2 overflows, on BOTH sides: near zero the
  // reciprocal-square limb is the whole answer, and on the negative side the
  // pole distance |u| is the argument that overflows.
  CheckSpecial("psi_1(smallest +subnormal)", 5e-324, kInf);
  CheckSpecial("psi_1(+subnormal)", 1e-320, kInf);
  CheckSpecial("psi_1(smallest -subnormal)", -5e-324, kInf);
  CheckSpecial("psi_1(-subnormal)", -1e-320, kInf);
  CheckSpecial("psi_1(1e-300)", 1e-300, kInf);
  CheckSpecial("psi_1(-1e-300)", -1e-300, kInf);
  // 2^-512 = 1/sqrt(DBL_MAX) is exactly where 1/x^2 reaches 2^1024.
  CheckSpecial("psi_1(2^-512)", 0x1p-512, kInf);

  // NaN propagates WITH its payload, matching every other corvus kernel.
  uint64_t bits = 0x7FF8000000ABCDEFULL;
  double payload;
  std::memcpy(&payload, &bits, sizeof(payload));
  const double got = One(payload);
  if (!SameBits(got, payload)) {
    std::fprintf(stderr, "FAIL: trigamma did not propagate the NaN payload\n");
    g_fail = 1;
  }
}

// Closed-form and reference values, one or more per region. Not an accuracy
// gate -- test_trigamma_ulp owns that -- just a check that each branch is
// wired to the right formula. Values are mpmath polygamma(1, .) at dps 60,
// rounded to nearest double.
void KnownValues() {
  struct Point {
    double x;
    double want;
    const char* region;
  };
  const Point kPoints[] = {
      // The first double above 2^-512: the largest FINITE psi_1 there is,
      // and the point where the scaled reciprocal-square must not overflow.
      {0x1.0000000000001p-512, 0x1.ffffffffffffcp+1023, "(0,1) deep tiny"},
      {0x1p-480, 0x1.0000000000000p+960, "(0,1) deep-tiny guard"},
      {1e-30, 0x1.3e9e4e4c2f343p+199, "(0,1) reciprocal-dominated"},
      {0.25, 0x1.1328429d927c6p+4, "(0,1)"},
      {0.5, 0x1.3bd3cc9be45dep+2, "(0,1) (= pi^2/2)"},
      {0.75, 0x1.455c4ff28f0bfp+1, "(0,1)"},
      {1.0, 0x1.a51a6625307d3p+0, "zone (= pi^2/6)"},
      {1.5, 0x1.de9e64df22ef3p-1, "zone (= pi^2/2 - 4, fit centre)"},
      {2.0, 0x1.4a34cc4a60fa6p-1, "walk (= pi^2/6 - 1)"},
      {3.0, 0x1.94699894c1f4dp-2, "walk"},
      {7.5, 0x1.2413cda19dd03p-3, "walk (deepest, 6 steps)"},
      {8.0, 0x1.10aa239ffbc61p-3, "asymptotic (X0)"},
      {10.0, 0x1.aec2e54649b87p-4, "asymptotic"},
      {100.0, 0x1.4952e891b603ap-7, "asymptotic"},
      {1.0e10, 0x1.b7cdfd9dda4e3p-34, "asymptotic"},
      {-0.25, 0x1.28ab89fe51e18p+4, "reflection"},
      {-0.5, 0x1.1de9e64df22efp+3, "reflection"},
      {-0.4957, 0x1.1ddb5aa1b2939p+3, "reflection (global min 8.933)"},
      {-1.5, 0x1.2c22c9dc2b128p+3, "reflection"},
      {-2.5, 0x1.3141822e1697ap+3, "reflection"},
      {-1000000.5, 0x1.3bd3ca83058d0p+3, "reflection (far)"},
      {-4222343377343339.5, 0x1.3bd3cc9be45dep+3, "reflection (near 2^52)"},
  };
  for (const Point& p : kPoints) {
    const double got = One(p.x);
    const uint64_t u = UlpDiff(got, p.want);
    if (u > 2) {
      std::fprintf(
          stderr,
          "FAIL: trigamma(%.17g) [%s] = %.17g, want %.17g (%llu ULP)\n", p.x,
          p.region, got, p.want, static_cast<unsigned long long>(u));
      g_fail = 1;
    }
  }
}

// A spread of in-domain points covering every region, both sides of every
// boundary the driver tests, and the walk at every depth 1..6.
std::vector<double> Probes() {
  return {
      // (0, 1), from the deep-tiny pole lane up to the zone wall.
      1e-320, 0x1p-512, 0x1.0000000000001p-512, 0x1p-480,
      0x1.0000000000001p-480, 1e-300, 1e-100, 1e-30, 1e-8, 0x1p-28, 0.001,
      0.1, 0.25, 0.5, 0.75, 0.9, 0.99999999999,
      // zone [1, 2), including the fit centre and both walls.
      1.0, 1.0000000001, 1.2, 1.5, 1.75, 1.9999999999,
      // walk [2, 8) at each depth 1..6, plus both sides of every step wall.
      2.0, 2.5, 2.9999999999, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0,
      7.5, 7.9999999999,
      // asymptotic [8, inf), both sides of X0 and out past the 1/x-only cut.
      7.9999999999999991, 8.0, 8.0000000000000018, 9.0, 12.0, 100.0, 1e6,
      1e17, 1e30, 0x1p88, 0x1p89, 0x1p90, 1e150, 1e300,
      // negative axis: between poles, at the global minimum, near poles
      // (including a lane deep enough to take the 1/u^2 shortcut), and large.
      -0.4957, -0.5, -0.75, -1.5, -1.9, -2.1, -2.5, -3.5, -4.5, -10.5, -20.5,
      -1e6 - 0.5, -1e15 - 0.5, -1.0000000000000002, -0.9999999999999998,
      -2.0000000000000004, -1e-320, -1e-300,
  };
}

void LaneMix() {
  const auto probes = Probes();
  const size_t n = probes.size();

  std::vector<double> alone(n);
  for (size_t i = 0; i < n; ++i) alone[i] = One(probes[i]);

  // Offsets 0..8 cover every lane position on every tier up to AVX-512, and
  // the leading filler is drawn from the far end of the probe list so each
  // probe sits next to a different region's point every time.
  for (size_t off = 0; off < 9; ++off) {
    std::vector<double> in;
    for (size_t k = 0; k < off; ++k) in.push_back(probes[(n - 1 - k) % n]);
    for (double v : probes) in.push_back(v);
    std::vector<double> got(in.size());
    corvus::trigamma(in, got);
    for (size_t i = 0; i < n; ++i) {
      if (!SameBits(got[off + i], alone[i])) {
        std::fprintf(stderr,
                     "FAIL: lane-mix drift at x=%.17g (offset %zu): batched "
                     "%.17g vs alone %.17g\n",
                     probes[i], off, got[off + i], alone[i]);
        g_fail = 1;
      }
    }
  }
}

}  // namespace

int main() {
  if (!corvus_test::ReportAndCheckTarget()) return 2;
  Specials();
  KnownValues();
  LaneMix();
  if (g_fail == 0) std::printf("PASS: trigamma smoke\n");
  return g_fail;
}
