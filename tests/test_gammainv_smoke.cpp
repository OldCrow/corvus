// Smoke test for corvus::gamma_p_inv / corvus::gamma_q_inv: the full specials
// table, the closed-form rows, a round trip through the forward pair, and
// lane-mix determinism. The ULP gate lives in test_gammainv_ulp.
//
// THE SPECIALS ARE SciPy'S, AND THE TWO EXPORTS MIRROR EACH OTHER: p = 0 puts
// the answer at the bottom of the support and p = 1 at +inf, while q = 0 is
// the far tail (+inf) and q = 1 the bottom. Everything outside [0, 1] is NaN,
// as is a <= 0 and a = +inf -- the last because the whole distribution has
// escaped to infinity and no finite quantile exists.
//
// THE CLOSED-FORM ROWS ARE THE ONLY INDEPENDENT ARITHMETIC IN THIS FILE.
// a = 1 makes the incomplete gamma an exponential, P(1,x) = 1 - e^-x, so the
// inverse is -log1p(-p) and -log(q) from libm, computed without touching any
// corvus code. Everything else here is a consistency check (round trip
// through corvus::gamma_p / gamma_q) or a structural one (lane mixing).
//
// LANE-MIX DETERMINISM IS WHY THIS FILE EXISTS SEPARATELY FROM THE GATE. The
// kernel is full of vector-wide decisions: the seed series and lambda's Newton
// break out when NO lane still needs them, the forward's four region cores and
// its two E formulas are each skipped when no lane is in them, and every core
// runs on all lanes with the inactive ones clamped or scrubbed. A point's
// answer must not depend on which other points shared its vector, so each
// probe is evaluated alone, at several lane offsets, and interleaved with
// points from other regimes -- all of it bit-identical.
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <span>
#include <vector>

#include "corvus/corvus.h"
#include "expect_target.h"
#include "ulp_utils.h"

namespace {

using corvus_test::SameBits;
using corvus_test::UlpDiff;

const double kInf = std::numeric_limits<double>::infinity();

int g_fail = 0;

double OneP(double a, double s) {
  double out = 0.0;
  corvus::gamma_p_inv(std::span<const double>(&a, 1),
                      std::span<const double>(&s, 1),
                      std::span<double>(&out, 1));
  return out;
}
double OneQ(double a, double s) {
  double out = 0.0;
  corvus::gamma_q_inv(std::span<const double>(&a, 1),
                      std::span<const double>(&s, 1),
                      std::span<double>(&out, 1));
  return out;
}

// A NaN expectation means "any NaN"; everything else is matched WITH its sign,
// which is the whole content of the +0 at the bottom of the support.
void CheckSpecial(const char* what, bool want_q, double a, double s,
                  double want) {
  const double got = want_q ? OneQ(a, s) : OneP(a, s);
  bool ok;
  if (std::isnan(want)) {
    ok = std::isnan(got);
  } else {
    ok = got == want && std::signbit(got) == std::signbit(want);
  }
  if (!ok) {
    std::fprintf(stderr, "FAIL: %s: gamma_%c_inv(%.17g, %.17g) = %.17g, want %.17g\n",
                 what, want_q ? 'q' : 'p', a, s, got, want);
    g_fail = 1;
  }
}

void Specials() {
  const double kAs[] = {1e-300, 1e-3, 0.5, 1.0, 3.0, 20.0, 1e6, 1e300};
  for (double a : kAs) {
    CheckSpecial("p = 0", false, a, 0.0, 0.0);
    CheckSpecial("p = -0", false, a, -0.0, 0.0);
    CheckSpecial("p = 1", false, a, 1.0, kInf);
    CheckSpecial("q = 0", true, a, 0.0, kInf);
    CheckSpecial("q = 1", true, a, 1.0, 0.0);
    CheckSpecial("q = -0", true, a, -0.0, kInf);

    // Outside [0, 1], both directions, including the subnormal just below 0.
    CheckSpecial("s < 0", false, a, -1e-320, NAN);
    CheckSpecial("s < 0", true, a, -0.5, NAN);
    CheckSpecial("s > 1", false, a, 1.0000000000000002, NAN);
    CheckSpecial("s > 1", true, a, 1e300, NAN);
    CheckSpecial("s = inf", false, a, kInf, NAN);
  }

  // Degenerate a, on both exports and at an ordinary s.
  for (double s : {0.25, 0.5, 0.75}) {
    CheckSpecial("a = 0", false, 0.0, s, NAN);
    CheckSpecial("a = -0", false, -0.0, s, NAN);
    CheckSpecial("a < 0", false, -1.0, s, NAN);
    CheckSpecial("a = +inf", false, kInf, s, NAN);
    CheckSpecial("a = -inf", false, -kInf, s, NAN);
    CheckSpecial("a = 0", true, 0.0, s, NAN);
    CheckSpecial("a < 0", true, -1e-320, s, NAN);
    CheckSpecial("a = +inf", true, kInf, s, NAN);
  }

  // NaN propagates WITH its payload, in either argument, matching every other
  // corvus kernel.
  uint64_t bits = 0x7FF8000000ABCDEFULL;
  double payload;
  std::memcpy(&payload, &bits, sizeof(payload));
  for (bool q : {false, true}) {
    const double g1 = q ? OneQ(payload, 0.25) : OneP(payload, 0.25);
    const double g2 = q ? OneQ(2.0, payload) : OneP(2.0, payload);
    if (!SameBits(g1, payload) || !SameBits(g2, payload)) {
      std::fprintf(stderr, "FAIL: gamma_%c_inv did not propagate the NaN payload\n",
                   q ? 'q' : 'p');
      g_fail = 1;
    }
  }
}

// a = 1: P(1,x) = 1 - e^-x and Q(1,x) = e^-x, so both inverses are libm one
// liners and neither touches corvus. The tolerance is loose on purpose --
// this checks the wiring, not the last bit, and log1p carries its own error.
void ClosedForm() {
  const double kPs[] = {1e-320, 1e-300, 1e-100, 1e-16, 1e-8, 0.001, 0.1,
                        0.25,   0.5,    0.75,   0.9,   0.99, 1.0 - 1e-12};
  for (double p : kPs) {
    const double want = -std::log1p(-p);
    const double got = OneP(1.0, p);
    if (UlpDiff(got, want) > 8) {
      std::fprintf(stderr,
                   "FAIL: gamma_p_inv(1, %.17g) = %.17g, want ~%.17g (%llu ulp)\n",
                   p, got, want, static_cast<unsigned long long>(UlpDiff(got, want)));
      g_fail = 1;
    }
  }
  const double kQs[] = {1e-320, 1e-300, 1e-100, 1e-16, 0.001, 0.1,
                        0.5,    0.9,    0.999999};
  for (double q : kQs) {
    const double want = -std::log(q);
    const double got = OneQ(1.0, q);
    if (UlpDiff(got, want) > 8) {
      std::fprintf(stderr,
                   "FAIL: gamma_q_inv(1, %.17g) = %.17g, want ~%.17g (%llu ulp)\n",
                   q, got, want, static_cast<unsigned long long>(UlpDiff(got, want)));
      g_fail = 1;
    }
  }
}

// Round trip through the forward pair, as a wiring check -- the ULP gate is
// the accuracy claim.
//
// TWO THINGS MAKE THE TOLERANCE WHAT IT IS. First, a round trip in p-space is
// a measurement of the FORWARD's condition number: one ulp of x maps to
// kappa = x*g/P ulps of p, and kappa grows like sqrt(a) (it is ~1e8 at
// a = 1e16), so a fixed bound would report on the conditioning rather than on
// either kernel. Second, x = 0 is a legitimate answer over whole slabs of the
// domain -- P(1e-3, 5e-324) is already 0.475, so every p below that has its
// true root far under the smallest subnormal -- and those rows have no round
// trip to check.
void RoundTrip() {
  std::vector<double> as, ss;
  for (double a : {1e-3, 0.1, 0.5, 1.0, 2.5, 7.0, 19.0, 20.0, 21.0, 100.0,
                   1e4, 1e8, 1e16}) {
    for (double s : {1e-30, 1e-8, 1e-3, 0.1, 0.3, 0.5, 0.7, 0.9, 0.999}) {
      as.push_back(a);
      ss.push_back(s);
    }
  }
  const size_t n = as.size();
  std::vector<double> x(n), back(n);

  corvus::gamma_p_inv(as, ss, x);
  corvus::gamma_p(as, x, back);
  for (size_t i = 0; i < n; ++i) {
    if (x[i] == 0.0) continue;  // underflowed answer: nothing to round trip
    const double tol = 1e-11 + 64.0 * std::sqrt(as[i]) *
                                   std::numeric_limits<double>::epsilon();
    if (!(x[i] > 0.0) || !std::isfinite(x[i])) {
      std::fprintf(stderr, "FAIL: gamma_p_inv(%.17g, %.17g) = %.17g not in range\n",
                   as[i], ss[i], x[i]);
      g_fail = 1;
      continue;
    }
    const double rel = std::fabs(back[i] - ss[i]) / ss[i];
    if (!(rel < tol)) {
      std::fprintf(stderr,
                   "FAIL: p round trip a=%.17g p=%.17g -> x=%.17g -> %.17g "
                   "(rel %.3e)\n",
                   as[i], ss[i], x[i], back[i], rel);
      g_fail = 1;
    }
  }

  corvus::gamma_q_inv(as, ss, x);
  corvus::gamma_q(as, x, back);
  for (size_t i = 0; i < n; ++i) {
    if (x[i] == 0.0) continue;  // underflowed answer: nothing to round trip
    const double tol = 1e-11 + 64.0 * std::sqrt(as[i]) *
                                   std::numeric_limits<double>::epsilon();
    if (!(x[i] > 0.0) || !std::isfinite(x[i])) {
      std::fprintf(stderr, "FAIL: gamma_q_inv(%.17g, %.17g) = %.17g not in range\n",
                   as[i], ss[i], x[i]);
      g_fail = 1;
      continue;
    }
    const double rel = std::fabs(back[i] - ss[i]) / ss[i];
    if (!(rel < tol)) {
      std::fprintf(stderr,
                   "FAIL: q round trip a=%.17g q=%.17g -> x=%.17g -> %.17g "
                   "(rel %.3e)\n",
                   as[i], ss[i], x[i], back[i], rel);
      g_fail = 1;
    }
  }
}

// Beyond resolution: for a past ~3e34 the entire transition happens inside one
// ulp of x, so every interior target has the same answer, x = a. This is the
// no-branch case -- the Newton steps freeze themselves because the inverse's
// condition number is 2^-500-class there -- so it is worth pinning as a
// behaviour, not just as a bucket in the ULP table.
void BeyondResolution() {
  for (double a : {1e35, 1e100, 1e300, 1.7976931348623157e308}) {
    for (double s : {0.25, 0.5, 0.75}) {
      for (bool q : {false, true}) {
        const double got = q ? OneQ(a, s) : OneP(a, s);
        if (got != a) {
          std::fprintf(stderr,
                       "FAIL: beyond-resolution gamma_%c_inv(%.17g, %.17g) = "
                       "%.17g, want a\n",
                       q ? 'q' : 'p', a, s, got);
          g_fail = 1;
        }
      }
    }
  }
}

// Every probe alone, then at each lane offset inside a vector of unrelated
// points, then interleaved with points from other regimes. All bit-identical.
void LaneMix() {
  struct Pt {
    double a, s;
  };
  const Pt pts[] = {
      {1e-4, 0.3},   {1e-4, 1e-40},  {0.1, 0.5},    {0.5, 1e-8},
      {1.0, 0.25},   {3.0, 0.75},    {19.0, 1e-12}, {20.0, 0.5},
      {50.0, 0.999}, {1e4, 0.5},     {1e4, 1e-100}, {1e8, 0.3},
      {1e20, 0.5},   {1e100, 0.5},   {1e-300, 0.7}, {2.0, 1e-300},
  };
  constexpr size_t kNPts = sizeof(pts) / sizeof(pts[0]);
  // Filler points drawn from other regimes, so a mixed vector really does mix
  // branches rather than just lane positions.
  const Pt filler[] = {{7.0, 0.4}, {1e-2, 1e-9}, {1e6, 0.6}, {30.0, 1e-30}};

  double alone_p[kNPts], alone_q[kNPts];
  for (size_t i = 0; i < kNPts; ++i) {
    alone_p[i] = OneP(pts[i].a, pts[i].s);
    alone_q[i] = OneQ(pts[i].a, pts[i].s);
  }

  constexpr size_t kN = 13;  // not a multiple of any lane count in the fleet
  for (size_t i = 0; i < kNPts; ++i) {
    for (size_t off = 0; off < kN; ++off) {
      std::vector<double> a(kN), s(kN), gp(kN), gq(kN);
      for (size_t j = 0; j < kN; ++j) {
        const Pt& f = filler[(j + i) % 4];
        a[j] = f.a;
        s[j] = f.s;
      }
      a[off] = pts[i].a;
      s[off] = pts[i].s;
      corvus::gamma_p_inv(a, s, gp);
      corvus::gamma_q_inv(a, s, gq);
      if (!SameBits(gp[off], alone_p[i]) || !SameBits(gq[off], alone_q[i])) {
        std::fprintf(stderr,
                     "FAIL: lane mix changed the answer at a=%.17g s=%.17g "
                     "(offset %zu): p %.17g vs %.17g, q %.17g vs %.17g\n",
                     pts[i].a, pts[i].s, off, gp[off], alone_p[i], gq[off],
                     alone_q[i]);
        g_fail = 1;
      }
    }
  }

  // One vector holding all the probes at once (plus a non-lane-multiple tail).
  std::vector<double> a(kNPts), s(kNPts), gp(kNPts), gq(kNPts);
  for (size_t i = 0; i < kNPts; ++i) {
    a[i] = pts[i].a;
    s[i] = pts[i].s;
  }
  corvus::gamma_p_inv(a, s, gp);
  corvus::gamma_q_inv(a, s, gq);
  for (size_t i = 0; i < kNPts; ++i) {
    if (!SameBits(gp[i], alone_p[i]) || !SameBits(gq[i], alone_q[i])) {
      std::fprintf(stderr,
                   "FAIL: all-probe vector changed the answer at a=%.17g "
                   "s=%.17g\n",
                   pts[i].a, pts[i].s);
      g_fail = 1;
    }
  }
}

}  // namespace

int main() {
  if (!corvus_test::ReportAndCheckTarget()) return 2;
  Specials();
  ClosedForm();
  RoundTrip();
  BeyondResolution();
  LaneMix();
  if (g_fail == 0) std::printf("PASS: gammainv smoke\n");
  return g_fail;
}
