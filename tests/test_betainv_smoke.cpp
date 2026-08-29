// Smoke test for corvus::beta_p_inv / corvus::beta_q_inv: the full specials
// table, the closed-form rows, the swap identity, a round trip through the
// forward pair, and lane-mix determinism. The ULP gate lives in
// test_betainv_ulp.
//
// THE SPECIALS ARE beta's OWN TABLE READ AS QUANTILES. p = 0 puts the answer
// at the bottom of the support and p = 1 at the top, while q = 0 is the top
// and q = 1 the bottom. A single degenerate parameter puts all the mass at one
// endpoint and every quantile is that endpoint (a = 0 or b = +inf -> 0;
// b = 0 or a = +inf -> 1); two degeneracies give NaN, as do a negative
// parameter and an out-of-range probability.
//
// THE CLOSED-FORM ROWS ARE THE ONLY INDEPENDENT ARITHMETIC IN THIS FILE.
// b = 1 makes the incomplete beta a pure power, I_x(a,1) = x^a, and a = 1
// makes it 1 - (1-x)^b; both inverses are libm one-liners computed without
// touching any corvus code. Everything else is a consistency check (round trip
// through corvus::beta_p / beta_q, or the swap identity against the kernel's
// other orientation) or a structural one (lane mixing).
//
// LANE-MIX DETERMINISM IS WHY THIS FILE EXISTS SEPARATELY FROM THE GATE. The
// kernel is full of vector-wide decisions: the seed series and the lambda
// Newton break out when NO lane still needs them, each of the seven seed
// candidates and each of beta's region cores is skipped when no lane is in it,
// and every core runs on all lanes with the inactive ones scrubbed. A point's
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

using corvus_test::OrderedBits;
using corvus_test::SameBits;
using corvus_test::UlpDiff;

const double kInf = std::numeric_limits<double>::infinity();

int g_fail = 0;

double OneP(double a, double b, double s) {
  double out = 0.0;
  corvus::beta_p_inv(std::span<const double>(&a, 1),
                     std::span<const double>(&b, 1),
                     std::span<const double>(&s, 1),
                     std::span<double>(&out, 1));
  return out;
}
double OneQ(double a, double b, double s) {
  double out = 0.0;
  corvus::beta_q_inv(std::span<const double>(&a, 1),
                     std::span<const double>(&b, 1),
                     std::span<const double>(&s, 1),
                     std::span<double>(&out, 1));
  return out;
}

// A NaN expectation means "any NaN"; everything else is matched WITH its sign,
// which is the whole content of the +0 at the bottom of the support.
void CheckSpecial(const char* what, bool want_q, double a, double b, double s,
                  double want) {
  const double got = want_q ? OneQ(a, b, s) : OneP(a, b, s);
  bool ok;
  if (std::isnan(want)) {
    ok = std::isnan(got);
  } else {
    ok = got == want && std::signbit(got) == std::signbit(want);
  }
  if (!ok) {
    std::fprintf(stderr,
                 "FAIL: %s: beta_%c_inv(%.17g, %.17g, %.17g) = %.17g, "
                 "want %.17g\n",
                 what, want_q ? 'q' : 'p', a, b, s, got, want);
    g_fail = 1;
  }
}

void Specials() {
  const double kPars[][2] = {{1e-300, 1e-300}, {1e-3, 5.0},  {0.5, 0.5},
                             {1.0, 1.0},       {3.0, 7.0},   {20.0, 20.0},
                             {1e6, 1e-6},      {1e300, 2.0}, {2.0, 1e300}};
  for (const auto& ab : kPars) {
    const double a = ab[0], b = ab[1];
    CheckSpecial("p = 0", false, a, b, 0.0, 0.0);
    CheckSpecial("p = -0", false, a, b, -0.0, 0.0);
    CheckSpecial("p = 1", false, a, b, 1.0, 1.0);
    CheckSpecial("q = 0", true, a, b, 0.0, 1.0);
    CheckSpecial("q = -0", true, a, b, -0.0, 1.0);
    CheckSpecial("q = 1", true, a, b, 1.0, 0.0);

    // Outside [0, 1], both directions, including the subnormal just below 0.
    CheckSpecial("s < 0", false, a, b, -1e-320, NAN);
    CheckSpecial("s < 0", true, a, b, -0.5, NAN);
    CheckSpecial("s > 1", false, a, b, 1.0000000000000002, NAN);
    CheckSpecial("s > 1", true, a, b, 1e300, NAN);
    CheckSpecial("s = inf", false, a, b, kInf, NAN);
  }

  // One degenerate parameter: all the mass sits at one endpoint.
  for (double s : {0.25, 0.5, 0.75}) {
    for (bool q : {false, true}) {
      CheckSpecial("a = 0 (mass at 0)", q, 0.0, 2.0, s, 0.0);
      CheckSpecial("a = -0 (mass at 0)", q, -0.0, 2.0, s, 0.0);
      CheckSpecial("b = inf (mass at 0)", q, 2.0, kInf, s, 0.0);
      CheckSpecial("b = 0 (mass at 1)", q, 2.0, 0.0, s, 1.0);
      CheckSpecial("a = inf (mass at 1)", q, kInf, 2.0, s, 1.0);
      // Two degeneracies, and negatives.
      CheckSpecial("a = 0, b = 0", q, 0.0, 0.0, s, NAN);
      CheckSpecial("a = inf, b = inf", q, kInf, kInf, s, NAN);
      CheckSpecial("a = 0, b = inf", q, 0.0, kInf, s, NAN);
      CheckSpecial("a < 0", q, -1.0, 2.0, s, NAN);
      CheckSpecial("b < 0", q, 2.0, -1e-320, s, NAN);
      CheckSpecial("a = -inf", q, -kInf, 2.0, s, NAN);
    }
  }

  // NaN propagates WITH its payload, in any argument.
  uint64_t bits = 0x7FF8000000ABCDEFULL;
  double payload;
  std::memcpy(&payload, &bits, sizeof(payload));
  for (bool q : {false, true}) {
    const double g1 = q ? OneQ(payload, 2.0, 0.25) : OneP(payload, 2.0, 0.25);
    const double g2 = q ? OneQ(2.0, payload, 0.25) : OneP(2.0, payload, 0.25);
    const double g3 = q ? OneQ(2.0, 3.0, payload) : OneP(2.0, 3.0, payload);
    if (!SameBits(g1, payload) || !SameBits(g2, payload) ||
        !SameBits(g3, payload)) {
      std::fprintf(stderr,
                   "FAIL: beta_%c_inv did not propagate the NaN payload\n",
                   q ? 'q' : 'p');
      g_fail = 1;
    }
  }
}

// The tolerance is dominated by the REFERENCE, not by the kernel. Every
// closed form here is exp(L/e) with L a log, and L carries |L|*2^-53
// absolutely, so the reference's own relative error is |L|*2^-53/e -- which at
// s = 1e-300 and e = 2 is 173 ulp before corvus computes anything. (The
// kernel's answer at a = 2, p = 1e-300 is the exactly-representable 1e-150;
// libm's round trip through log and exp is the one that misses.)
double ClosedTolUlp(double lg, double e) {
  return 16.0 + std::fabs(lg) / e;
}

void CheckClosed(const char* what, double got, double want, double a, double b,
                 double s, double tol_ulp) {
  if (static_cast<double>(UlpDiff(got, want)) > tol_ulp) {
    std::fprintf(stderr,
                 "FAIL: %s at a=%.17g b=%.17g s=%.17g: got %.17g want ~%.17g "
                 "(%llu ulp)\n",
                 what, a, b, s, got, want,
                 static_cast<unsigned long long>(UlpDiff(got, want)));
    g_fail = 1;
  }
}

// b = 1: I_x(a,1) = x^a exactly, so beta_p_inv(a,1,p) = p^(1/a) and
// beta_q_inv(a,1,q) = (1-q)^(1/a).
// a = 1: I_x(1,b) = 1 - (1-x)^b, so beta_p_inv(1,b,p) = 1 - (1-p)^(1/b) and
// beta_q_inv(1,b,q) = 1 - q^(1/b).
// The near-one answers are skipped: there the libm reference is a subtraction
// of two nearly equal numbers and would be measuring itself.
void ClosedForm() {
  const double kExps[] = {0.25, 0.5, 1.0, 2.0, 7.0, 40.0};
  const double kSs[] = {1e-300, 1e-100, 1e-20, 1e-8, 0.001, 0.1,
                        0.25,   0.5,    0.75,  0.9,  0.99};
  for (double e : kExps) {
    for (double s : kSs) {
      {
        const double lg = std::log(s);
        const double want = std::exp(lg / e);
        if (want < 0.9 && want > 0.0) {
          CheckClosed("beta_p_inv(a,1,p) = p^(1/a)", OneP(e, 1.0, s), want, e,
                      1.0, s, ClosedTolUlp(lg, e));
        }
      }
      {
        const double lg = std::log1p(-s);
        const double want = std::exp(lg / e);
        if (want < 0.9 && want > 0.0) {
          CheckClosed("beta_q_inv(a,1,q) = (1-q)^(1/a)", OneQ(e, 1.0, s), want,
                      e, 1.0, s, ClosedTolUlp(lg, e));
        }
      }
      {
        const double lg = std::log1p(-s);
        const double want = -std::expm1(lg / e);
        if (want < 0.9 && want > 0.0) {
          CheckClosed("beta_p_inv(1,b,p) = 1-(1-p)^(1/b)", OneP(1.0, e, s),
                      want, 1.0, e, s, ClosedTolUlp(lg, e));
        }
      }
      {
        const double lg = std::log(s);
        const double want = -std::expm1(lg / e);
        if (want < 0.9 && want > 0.0) {
          CheckClosed("beta_q_inv(1,b,q) = 1-q^(1/b)", OneQ(1.0, e, s), want,
                      1.0, e, s, ClosedTolUlp(lg, e));
        }
      }
    }
  }
  // The uniform case, where the quantile is the probability itself and the
  // answer must be EXACT: both flips and both orientations collapse to
  // identity, so any stray rounding in the frame shows up here.
  for (double s : kSs) {
    const double gp = OneP(1.0, 1.0, s);
    const double gq = OneQ(1.0, 1.0, s);
    if (UlpDiff(gp, s) > 1 || UlpDiff(gq, 1.0 - s) > 1) {
      std::fprintf(stderr,
                   "FAIL: uniform a=b=1 at s=%.17g: p_inv %.17g (want %.17g), "
                   "q_inv %.17g (want %.17g)\n",
                   s, gp, s, gq, 1.0 - s);
      g_fail = 1;
    }
  }
}

// THE SWAP IDENTITY, which the public header documents as the lossless route
// for a near-one answer: I_x(a,b) = 1 - I_{1-x}(b,a), so
//     1 - beta_p_inv(a, b, p) == beta_p_inv(b, a, 1-p)
// up to the single rounding of the subtraction. Checked in the direction where
// BOTH sides are well conditioned (the answer away from either endpoint), so a
// disagreement means the frame is wrong rather than that a double ran out of
// bits.
void SwapIdentity() {
  const double kAs[] = {0.03, 0.5, 1.5, 4.0, 30.0, 1e4};
  const double kBs[] = {0.07, 0.9, 2.0, 11.0, 500.0, 1e6};
  const double kPs[] = {0.05, 0.2, 0.5, 0.8, 0.95};
  for (double a : kAs) {
    for (double b : kBs) {
      for (double p : kPs) {
        const double x = OneP(a, b, p);
        const double z = OneP(b, a, 1.0 - p);
        if (!(x > 1e-8 && x < 1.0 - 1e-8)) continue;
        // The tolerance carries the SUBTRACTION's own cost, which is the
        // point of the identity rather than an exception to it: 1 - x formed
        // from a double x is exact to half an ulp of 1, so it agrees with the
        // directly-solved z only to 1.1e-16/z relative. That gap IS the
        // lossless-near-1 story -- z is the number that keeps its digits.
        const double tol =
            1e-13 + 4.0 * std::numeric_limits<double>::epsilon() / z;
        const double rel = std::fabs((1.0 - x) - z) / z;
        if (!(rel < tol)) {
          std::fprintf(stderr,
                       "FAIL: swap identity a=%.17g b=%.17g p=%.17g: "
                       "1-x=%.17g vs %.17g (rel %.3e)\n",
                       a, b, p, 1.0 - x, z, rel);
          g_fail = 1;
        }
      }
    }
  }
}

// Round trip through the forward pair, as a wiring check -- the ULP gate is
// the accuracy claim. The tolerance is a measurement of the FORWARD's own
// condition number: one ulp of x maps to kappa = x g(x)/I ulps of p, which
// grows without bound near either endpoint and with the parameters, so a fixed
// bound would report on the conditioning rather than on either kernel.
void RoundTrip() {
  std::vector<double> as, bs, ss;
  for (double a : {1e-3, 0.1, 0.5, 1.0, 2.5, 19.0, 100.0, 1e4}) {
    for (double b : {1e-3, 0.1, 0.7, 1.0, 3.0, 21.0, 400.0, 1e6}) {
      for (double s : {1e-30, 1e-8, 1e-3, 0.1, 0.3, 0.5, 0.7, 0.9}) {
        as.push_back(a);
        bs.push_back(b);
        ss.push_back(s);
      }
    }
  }
  const size_t n = as.size();
  std::vector<double> x(n), back(n);

  for (bool q : {false, true}) {
    if (q) {
      corvus::beta_q_inv(as, bs, ss, x);
      corvus::beta_q(as, bs, x, back);
    } else {
      corvus::beta_p_inv(as, bs, ss, x);
      corvus::beta_p(as, bs, x, back);
    }
    for (size_t i = 0; i < n; ++i) {
      if (!(x[i] >= 0.0 && x[i] <= 1.0)) {
        std::fprintf(stderr,
                     "FAIL: beta_%c_inv(%.17g, %.17g, %.17g) = %.17g not in "
                     "[0, 1]\n",
                     q ? 'q' : 'p', as[i], bs[i], ss[i], x[i]);
        g_fail = 1;
        continue;
      }
      // No round trip to check where the answer saturates the representable
      // range: whole slabs of this domain have their true root far below the
      // smallest subnormal or within one ulp of 1.
      if (x[i] == 0.0 || x[i] == 1.0) continue;
      // NOR WHERE x IS NEAR 1, and that is the API's own statement rather than
      // a tolerance dodge. A double cannot hold 1 - 1e-12 to better than
      // 1e-16 absolute, so ANY x within ~1e-8 of the top carries a relative
      // error of at least 1e-8 in 1 - x -- and the forward's sensitivity there
      // is d I/d x, which for a shape parameter of 1e4 is enormous. The round
      // trip would be measuring the representation, not either kernel. This is
      // exactly the case the public header sends to the swap identity:
      // 1 - x at full precision is beta_p_inv(b, a, q).
      if (1.0 - x[i] < 1e-8) continue;
      const double kappa = std::sqrt(std::fmax(as[i], bs[i]));
      // The second term is the NEAR-ONE representation cost, not slack: an x
      // reported as 1 - d carries an absolute error of half an ulp of 1, so
      // d itself is only good to 1.1e-16/d relative, and the forward's
      // dependence on d goes like b. Nothing either kernel does can improve
      // it -- that is the whole reason the public header documents the swap
      // identity for callers who need 1 - x.
      const double tol = 1e-11 +
                         4096.0 * kappa *
                             std::numeric_limits<double>::epsilon() +
                         8.0 * std::fmax(bs[i], 1.0) *
                             std::numeric_limits<double>::epsilon() /
                             (1.0 - x[i]);
      const double rel = std::fabs(back[i] - ss[i]) / ss[i];
      if (!(rel < tol)) {
        std::fprintf(stderr,
                     "FAIL: %c round trip a=%.17g b=%.17g s=%.17g -> x=%.17g "
                     "-> %.17g (rel %.3e, tol %.3e)\n",
                     q ? 'q' : 'p', as[i], bs[i], ss[i], x[i], back[i], rel,
                     tol);
        g_fail = 1;
      }
    }
  }
}

// Every probe alone, then at each lane offset inside a vector of unrelated
// points, then interleaved with points from other regimes. All bit-identical.
void LaneMix() {
  struct Pt {
    double a, b, s;
  };
  const Pt pts[] = {
      {1e-4, 1e-4, 0.3},    {1e-4, 5.0, 1e-40},  {0.1, 0.1, 0.5},
      {0.5, 3.0, 1e-8},     {1.0, 1.0, 0.25},    {3.0, 7.0, 0.75},
      {19.0, 21.0, 1e-12},  {40.0, 40.0, 0.5},   {2.0, 1e6, 0.999},
      {1e4, 1e4, 0.5},      {1e4, 2.0, 1e-100},  {1e8, 3.0, 0.3},
      {1e20, 1e20, 0.5},    {1e100, 2.0, 0.5},   {1e-300, 0.7, 0.7},
      {2.0, 2.0, 1e-300},   {1e-20, 1e-20, 0.5}, {5.0, 1e300, 0.4},
  };
  constexpr size_t kNPts = sizeof(pts) / sizeof(pts[0]);
  const Pt filler[] = {{7.0, 2.0, 0.4},
                       {1e-2, 1e-2, 1e-9},
                       {1e6, 1e6, 0.6},
                       {30.0, 0.5, 1e-30}};

  double alone_p[kNPts], alone_q[kNPts];
  for (size_t i = 0; i < kNPts; ++i) {
    alone_p[i] = OneP(pts[i].a, pts[i].b, pts[i].s);
    alone_q[i] = OneQ(pts[i].a, pts[i].b, pts[i].s);
  }

  constexpr size_t kN = 13;  // not a multiple of any lane count in the fleet
  for (size_t i = 0; i < kNPts; ++i) {
    for (size_t off = 0; off < kN; ++off) {
      std::vector<double> a(kN), b(kN), s(kN), gp(kN), gq(kN);
      for (size_t j = 0; j < kN; ++j) {
        const Pt& f = filler[(j + i) % 4];
        a[j] = f.a;
        b[j] = f.b;
        s[j] = f.s;
      }
      a[off] = pts[i].a;
      b[off] = pts[i].b;
      s[off] = pts[i].s;
      corvus::beta_p_inv(a, b, s, gp);
      corvus::beta_q_inv(a, b, s, gq);
      if (!SameBits(gp[off], alone_p[i]) || !SameBits(gq[off], alone_q[i])) {
        std::fprintf(stderr,
                     "FAIL: lane mix changed the answer at a=%.17g b=%.17g "
                     "s=%.17g (offset %zu): p %.17g vs %.17g, q %.17g vs "
                     "%.17g\n",
                     pts[i].a, pts[i].b, pts[i].s, off, gp[off], alone_p[i],
                     gq[off], alone_q[i]);
        g_fail = 1;
      }
    }
  }

  // One vector holding all the probes at once (a non-lane-multiple length).
  std::vector<double> a(kNPts), b(kNPts), s(kNPts), gp(kNPts), gq(kNPts);
  for (size_t i = 0; i < kNPts; ++i) {
    a[i] = pts[i].a;
    b[i] = pts[i].b;
    s[i] = pts[i].s;
  }
  corvus::beta_p_inv(a, b, s, gp);
  corvus::beta_q_inv(a, b, s, gq);
  for (size_t i = 0; i < kNPts; ++i) {
    if (!SameBits(gp[i], alone_p[i]) || !SameBits(gq[i], alone_q[i])) {
      std::fprintf(stderr,
                   "FAIL: all-probe vector changed the answer at a=%.17g "
                   "b=%.17g s=%.17g\n",
                   pts[i].a, pts[i].b, pts[i].s);
      g_fail = 1;
    }
  }
}

}  // namespace

int main() {
  if (!corvus_test::ReportAndCheckTarget()) return 2;
  Specials();
  ClosedForm();
  SwapIdentity();
  RoundTrip();
  LaneMix();
  if (g_fail == 0) std::printf("PASS: betainv smoke\n");
  return g_fail;
}
