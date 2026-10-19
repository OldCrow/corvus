// Guard harness (#53): cross-function checks that need no reference file, so
// they cover the bands the reference sets never sampled. The 2026-10-10
// defect hunt found every one of #47-#61 in such a band; the ULP gates pass
// there because they have no rows to fail on.
//
// FIVE CHECKS, every one over the public API on the dispatched tier:
//   1. harness self-test -- UlpDiff across opposite signs (ulp_utils.h UB).
//   2. scalar = span     -- the eight scalar entry points are bit-identical
//      to a length-N span call over an IEEE-specials product grid (the span
//      call mixes every special into every lane position).
//   3. sign/range        -- every forward and inverse output of the
//      incomplete gamma/beta families lies in its range with a clear sign
//      bit, over subnormal and huge shapes; NaN only for a NaN contract.
//   4. beta symmetry     -- beta_p(a, b, x) against beta_q(b, a, 1 - x) for
//      x in [1/2, 1), where 1 - x is exact (Sterbenz), at x = 1/2 exactly,
//      at the R4 window boundary fl((a+1)/(a+b+2)) and its neighbours, and
//      for tiny shape pairs.
//   5. inverse bracket   -- x = inv(p) must satisfy F(x - k ulp) <= p <=
//      F(x + k ulp) (reversed for the Q side), with F's own error allowed
//      for. Coarse by design: k and the F slack sit well above every gate,
//      so a failure here is a wrong root, not a last-bit miss. The ULP
//      gates own the last bits.
//
// No check here has an oracle beyond corvus's own forward functions, so a
// pass proves consistency, not accuracy. What they buy is coverage: seeded
// grids reach the subnormal/huge corners at negligible cost on every tier.
//
// Open defects the grids reach are strict known-defect bands (below): the
// test fails outside them, and fails when one stops failing. `--all` prints
// every failure, banded ones included.
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
const double kNaN = std::numeric_limits<double>::quiet_NaN();
const double kDmin = std::numeric_limits<double>::denorm_min();
const double kMax = std::numeric_limits<double>::max();

int g_fail = 0;
size_t g_print_limit = 8;  // --all prints every failure

// ---------------------------------------------------------------------------
// KNOWN-DEFECT BANDS [DECIDED 2026-10-17, user]: STRICT, ONE PER OPEN ISSUE.
//
// The guards landed before the fixes they guard (v1.1.0 order row 2; the
// fixes are rows 3-10). Each open defect the grids reach is a band below,
// tagged with its issue and pinned to the extent measured on AVX3_ZEN4 and
// SSE4 (clang-cl, 2026-10-17; identical failing points on both). The test
// fails on
//   * any failure outside every band -- a new defect or a regression, and
//   * any band with NO failure -- the fix has landed, so the band must go.
// A fix commit therefore deletes its band: the band's presence is the
// fail-first record, its deletion is the guard landing with the fix. Never
// widen a band to absorb a new failure; file it.
//
// Arguments are the failing call's: `x` is the probability for an inverse.

struct Band {
  int issue;
  const char* check;
  bool (*match)(const char* fn, double a, double b, double x);
  size_t hits = 0;
};

bool Is(const char* fn, const char* name) { return std::strcmp(fn, name) == 0; }

Band g_bands[] = {
    // #48: gammainv diverges for small a in the upper tail.
    {48, "inverse bracket",
     [](const char* fn, double a, double, double x) {
       return a < 0.5 && ((Is(fn, "gamma_p_inv") && x > 0.5) ||
                          (Is(fn, "gamma_q_inv") && x < 0.5));
     }},
    // #50: betainv off by up to ~2e9 ulp for a moderate shape and a tiny one,
    // p below the tiny one (beta_q_inv: the shapes swapped). Measured to big
    // in [20, 1.8e4], small <= 5.4e-6 -- wider than the issue's b <~ 1e-8.
    {50, "inverse bracket",
     [](const char* fn, double a, double b, double x) {
       const bool p = Is(fn, "beta_p_inv");
       if (!p && !Is(fn, "beta_q_inv")) return false;
       const double big = p ? a : b, small = p ? b : a;
       return big >= 10.0 && big <= 3e4 && small <= 1e-5 && x < small;
     }},
    // #55: betainv returns ~0.5 for a >~ 1e17 with the root near 1. Capped
    // at the grid's 1e20: c85e4de's witness points at a = 1e100 must stay
    // outside every band, or a revert of its clamp would pass as #55.
    {55, "inverse bracket",
     [](const char* fn, double a, double, double) {
       return (Is(fn, "beta_p_inv") || Is(fn, "beta_q_inv")) && a >= 1e16 &&
              a <= 1e20;
     }},
    // #55 (found under #47): both shapes tiny -- the backward-error regime
    // -- with a >> b, where the orientation probe falls back to 0.5
    // (beta_q_inv(1e-20, 1e-70, q) = 0.5; the root is 1). Measured on
    // AVX3_ZEN4/SSE4/SSE2, identical: a <= 2.8e-17, b/a <= 1.6e-33, never
    // a < b. Row 9's orientation pass deletes it.
    {55, "inverse bracket",
     [](const char* fn, double a, double b, double) {
       return (Is(fn, "beta_p_inv") || Is(fn, "beta_q_inv")) && a <= 1e-16 &&
              b < a * 1e-32;
     }},
};

// ---------------------------------------------------------------------------
// Failure accounting. Every check reports through Fail(); a failure inside a
// known-defect band counts against the band, any other against the check.
// The first few of each print in full, the rest only count.

struct Check {
  const char* name;
  size_t n = 0;
  size_t fails = 0;
  size_t known = 0;
};

void Fail(Check* c, const char* fn, double a, double b, double x, double got,
          const char* why) {
  for (Band& band : g_bands) {
    if (std::strcmp(band.check, c->name) != 0 || !band.match(fn, a, b, x))
      continue;
    ++band.hits;
    ++c->known;
    if (g_print_limit == SIZE_MAX) {
      std::printf("  known #%d %-12s a=%a b=%a x=%a -> %a\n", band.issue, fn,
                  a, b, x, got);
    }
    return;
  }
  ++c->fails;
  if (c->fails <= g_print_limit) {
    std::printf("  FAIL %-12s %-12s a=%a b=%a x=%a -> %a  (%s)\n", c->name, fn,
                a, b, x, got, why);
  }
}

void Report(const Check& c) {
  std::printf("%-16s %8zu points, %6zu failures, %6zu in known bands\n", c.name,
              c.n, c.fails, c.known);
  if (c.fails > 0) ++g_fail;
}

void ReportBands() {
  for (const Band& band : g_bands) {
    std::printf("known-defect band #%d (%s): %zu failures\n", band.issue,
                band.check, band.hits);
    if (band.hits == 0) {
      std::printf(
          "  FAIL band #%d no longer fails here: if #%d is fixed, delete its "
          "band in the same commit\n",
          band.issue, band.issue);
      ++g_fail;
    }
  }
}

// Deterministic generator (splitmix64): a fixed seed keeps every grid
// identical across tiers and runs, so a failure always names the same point.
struct Rng {
  uint64_t s;
  uint64_t Next() {
    uint64_t z = (s += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
  }
  double Unit() { return static_cast<double>(Next() >> 11) * 0x1p-53; }
  // log-uniform in [lo, hi], lo > 0
  double LogU(double lo, double hi) {
    return std::exp(std::log(lo) + Unit() * (std::log(hi) - std::log(lo)));
  }
};

// x moved k ulps in direction dir (+1/-1), clamped to [lo, hi].
double Step(double x, int k, int dir, double lo, double hi) {
  for (int i = 0; i < k; ++i) x = std::nextafter(x, dir > 0 ? kInf : -kInf);
  return x < lo ? lo : (x > hi ? hi : x);
}

// a <= b within `slack` ulps (a may exceed b by up to slack).
bool LeSlack(double a, double b, uint64_t slack) {
  return a <= b || UlpDiff(a, b) <= slack;
}

// ---------------------------------------------------------------------------
// 1. Harness self-test.

void CheckHarness() {
  Check c{"harness"};
  struct {
    double a, b;
    uint64_t want;
  } rows[] = {
      {-kDmin, kDmin, 2},
      {-0.0, 0.0, 0},
      {-1.0, 1.0, 2 * 0x3FF0000000000000ULL},
      {-2.0, 2.0, 2 * 0x4000000000000000ULL},
      {-3.0, 3.0, 2 * 0x4008000000000000ULL},
      {-kMax, kMax, 2 * 0x7FEFFFFFFFFFFFFFULL},
      {1.0, std::nextafter(1.0, 2.0), 1},
  };
  for (const auto& r : rows) {
    ++c.n;
    const uint64_t got = UlpDiff(r.a, r.b);
    if (got != r.want || UlpDiff(r.b, r.a) != r.want) {
      ++c.fails;
      std::printf("  FAIL harness UlpDiff(%a, %a) = %llu, want %llu\n", r.a,
                  r.b, static_cast<unsigned long long>(got),
                  static_cast<unsigned long long>(r.want));
    }
  }
  Report(c);
}

// ---------------------------------------------------------------------------
// 2. Scalar entry points = span path, bit for bit.

std::vector<double> Specials() {
  double payload;
  const uint64_t bits = 0x7FF8000000000ABCULL;  // quiet NaN with a payload
  std::memcpy(&payload, &bits, sizeof(payload));
  return {0.0,     -0.0,    kDmin,  -kDmin,   0x1p-1022, 1e-300, 1e-17,
          1e-8,    0.25,    0.5,    std::nextafter(1.0, 0.0), 1.0,
          std::nextafter(1.0, 2.0), 2.0,    15.2,   1e5,      0x1p59,
          1e300,   kMax,    kInf,   -1.0,    -kInf,  kNaN,     payload};
}

template <class Span2, class Scalar2>
void ScalarSpan2(Check* c, const char* fn, Span2 span_fn, Scalar2 scalar_fn) {
  const std::vector<double> s = Specials();
  std::vector<double> a, x;
  for (double ai : s)
    for (double xi : s) {
      a.push_back(ai);
      x.push_back(xi);
    }
  std::vector<double> out(a.size());
  span_fn(a, x, out);
  for (size_t i = 0; i < a.size(); ++i) {
    ++c->n;
    const double one = scalar_fn(a[i], x[i]);
    if (!SameBits(one, out[i])) Fail(c, fn, a[i], 0.0, x[i], one, "scalar != span");
  }
}

template <class Span3, class Scalar3>
void ScalarSpan3(Check* c, const char* fn, Span3 span_fn, Scalar3 scalar_fn) {
  const std::vector<double> s = Specials();
  std::vector<double> a, b, x;
  for (double ai : s)
    for (double bi : s)
      for (double xi : s) {
        a.push_back(ai);
        b.push_back(bi);
        x.push_back(xi);
      }
  std::vector<double> out(a.size());
  span_fn(a, b, x, out);
  for (size_t i = 0; i < a.size(); ++i) {
    ++c->n;
    const double one = scalar_fn(a[i], b[i], x[i]);
    if (!SameBits(one, out[i])) Fail(c, fn, a[i], b[i], x[i], one, "scalar != span");
  }
}

void CheckScalarSpan() {
  Check c{"scalar=span"};
  using S = std::span<const double>;
  using O = std::span<double>;
  ScalarSpan2(&c, "gamma_p", [](S a, S x, O o) { corvus::gamma_p(a, x, o); },
              [](double a, double x) { return corvus::gamma_p(a, x); });
  ScalarSpan2(&c, "gamma_q", [](S a, S x, O o) { corvus::gamma_q(a, x, o); },
              [](double a, double x) { return corvus::gamma_q(a, x); });
  ScalarSpan2(&c, "gamma_p_inv",
              [](S a, S x, O o) { corvus::gamma_p_inv(a, x, o); },
              [](double a, double x) { return corvus::gamma_p_inv(a, x); });
  ScalarSpan2(&c, "gamma_q_inv",
              [](S a, S x, O o) { corvus::gamma_q_inv(a, x, o); },
              [](double a, double x) { return corvus::gamma_q_inv(a, x); });
  ScalarSpan3(&c, "beta_p",
              [](S a, S b, S x, O o) { corvus::beta_p(a, b, x, o); },
              [](double a, double b, double x) { return corvus::beta_p(a, b, x); });
  ScalarSpan3(&c, "beta_q",
              [](S a, S b, S x, O o) { corvus::beta_q(a, b, x, o); },
              [](double a, double b, double x) { return corvus::beta_q(a, b, x); });
  ScalarSpan3(&c, "beta_p_inv",
              [](S a, S b, S x, O o) { corvus::beta_p_inv(a, b, x, o); },
              [](double a, double b, double x) {
                return corvus::beta_p_inv(a, b, x);
              });
  ScalarSpan3(&c, "beta_q_inv",
              [](S a, S b, S x, O o) { corvus::beta_q_inv(a, b, x, o); },
              [](double a, double b, double x) {
                return corvus::beta_q_inv(a, b, x);
              });
  Report(c);
}

// ---------------------------------------------------------------------------
// 3. Sign/range over subnormal and huge shapes.

std::vector<double> Shapes() {
  std::vector<double> v;
  for (int k = 1; k <= 9; ++k) v.push_back(k * kDmin);  // odd and even k (#51)
  for (double s : {0x1p-1060, 1e-310, 0x1p-1022, 1e-300, 1e-200, 1e-140,
                   1e-100, 1e-40, 1e-20, 1e-12, 1e-8, 1e-4, 0.01, 0.1, 0.47,
                   0.5, 1.0, 1.5, 2.0, 7.5, 15.2, 20.0, 26.0, 100.0, 1e4, 1e8,
                   1e17, 0x1p59, 1e30, 1e100, 1e300, kMax})
    v.push_back(s);
  return v;
}

std::vector<double> GammaXs() {
  std::vector<double> v = {0.0, kDmin, 1e-300, 1e-20, 1e-4, 0.5, 1.0};
  for (double x = 1.6; x <= 2.5; x += 0.05) v.push_back(x);  // #51 window
  for (double s : {3.0, 15.0, 30.0, 100.0, 1e4, 1e17, 1e300, kMax, kInf})
    v.push_back(s);
  return v;
}

std::vector<double> UnitXs() {
  std::vector<double> v = {0.0, kDmin, 1e-300, 1e-100, 1e-20, 1e-8, 0.01,
                           0.25, std::nextafter(0.5, 0.0), 0.5,
                           std::nextafter(0.5, 1.0), 0.75, 0.99,
                           std::nextafter(1.0, 0.0), 1.0};
  return v;
}

bool InUnit(double r) { return !std::isnan(r) && !std::signbit(r) && r <= 1.0; }

void CheckRange() {
  Check c{"sign/range"};
  const std::vector<double> sh = Shapes();
  // gamma_p / gamma_q: in [0, 1] for every finite positive a, x >= 0.
  {
    std::vector<double> a, x;
    for (double ai : sh)
      for (double xi : GammaXs()) {
        a.push_back(ai);
        x.push_back(xi);
      }
    std::vector<double> p(a.size()), q(a.size());
    corvus::gamma_p(a, x, p);
    corvus::gamma_q(a, x, q);
    for (size_t i = 0; i < a.size(); ++i) {
      c.n += 2;
      if (!InUnit(p[i])) Fail(&c, "gamma_p", a[i], 0.0, x[i], p[i], "range");
      if (!InUnit(q[i])) Fail(&c, "gamma_q", a[i], 0.0, x[i], q[i], "range");
    }
  }
  // beta_p / beta_q: in [0, 1] for finite positive a, b and x in [0, 1].
  {
    std::vector<double> a, b, x;
    for (double ai : sh)
      for (double bi : sh)
        for (double xi : UnitXs()) {
          a.push_back(ai);
          b.push_back(bi);
          x.push_back(xi);
        }
    std::vector<double> p(a.size()), q(a.size());
    corvus::beta_p(a, b, x, p);
    corvus::beta_q(a, b, x, q);
    for (size_t i = 0; i < a.size(); ++i) {
      c.n += 2;
      if (!InUnit(p[i])) Fail(&c, "beta_p", a[i], b[i], x[i], p[i], "range");
      if (!InUnit(q[i])) Fail(&c, "beta_q", a[i], b[i], x[i], q[i], "range");
    }
  }
  // Inverses: gamma in [0, +inf], beta in [0, 1], for p in [0, 1].
  {
    std::vector<double> a, s;
    for (double ai : sh)
      for (double si : UnitXs()) {
        a.push_back(ai);
        s.push_back(si);
      }
    std::vector<double> xp(a.size()), xq(a.size());
    corvus::gamma_p_inv(a, s, xp);
    corvus::gamma_q_inv(a, s, xq);
    for (size_t i = 0; i < a.size(); ++i) {
      c.n += 2;
      if (std::isnan(xp[i]) || std::signbit(xp[i]))
        Fail(&c, "gamma_p_inv", a[i], 0.0, s[i], xp[i], "range");
      if (std::isnan(xq[i]) || std::signbit(xq[i]))
        Fail(&c, "gamma_q_inv", a[i], 0.0, s[i], xq[i], "range");
    }
  }
  {
    std::vector<double> a, b, s;
    for (double ai : sh)
      for (double bi : sh)
        for (double si : UnitXs()) {
          a.push_back(ai);
          b.push_back(bi);
          s.push_back(si);
        }
    std::vector<double> xp(a.size()), xq(a.size());
    corvus::beta_p_inv(a, b, s, xp);
    corvus::beta_q_inv(a, b, s, xq);
    for (size_t i = 0; i < a.size(); ++i) {
      c.n += 2;
      if (!InUnit(xp[i])) Fail(&c, "beta_p_inv", a[i], b[i], s[i], xp[i], "range");
      if (!InUnit(xq[i])) Fail(&c, "beta_q_inv", a[i], b[i], s[i], xq[i], "range");
    }
  }
  Report(c);
}

// ---------------------------------------------------------------------------
// 4. Beta reflection symmetry I_x(a, b) = 1 - I_{1-x}(b, a).
//
// Both sides are corvus values of the same number, each within its region
// gate (<= 3 ulp, test_beta_ulp.cpp), so they may differ by the sum.

constexpr uint64_t kSymUlp = 6;

void CheckSymmetry() {
  Check c{"beta symmetry"};
  const double tiny[] = {1e-300, 1e-200, 1e-140, 1e-100, 1e-40, 1e-20,
                         1e-12,  1e-8,   1e-4,   1e-3,   0.01,  0.1,
                         0.5,    1.0,    3.0,    40.0};
  std::vector<double> a, b, x;
  for (double ai : tiny)
    for (double bi : tiny) {
      std::vector<double> xs = {0.5, 0.75, 0.9, 0.999};
      const double w = (ai + 1.0) / (ai + bi + 2.0);  // R4 window boundary
      for (int k = -3; k <= 3; ++k) xs.push_back(Step(w, k < 0 ? -k : k, k < 0 ? -1 : 1, 0.0, 1.0));
      for (int k = 1; k <= 3; ++k) xs.push_back(Step(0.5, k, 1, 0.0, 1.0));
      for (double xi : xs) {
        if (!(xi >= 0.5 && xi < 1.0)) continue;  // 1 - x exact only here
        a.push_back(ai);
        b.push_back(bi);
        x.push_back(xi);
      }
    }
  std::vector<double> y(x.size());
  for (size_t i = 0; i < x.size(); ++i) y[i] = 1.0 - x[i];
  std::vector<double> p(x.size()), q(x.size());
  corvus::beta_p(a, b, x, p);
  corvus::beta_q(b, a, y, q);
  for (size_t i = 0; i < x.size(); ++i) {
    ++c.n;
    if (std::isnan(p[i]) || UlpDiff(p[i], q[i]) > kSymUlp)
      Fail(&c, "beta_p/q", a[i], b[i], x[i], p[i] - q[i], "P(a,b,x) != Q(b,a,1-x)");
  }
  Report(c);
}

// ---------------------------------------------------------------------------
// 5. Inverse bracket.
//
// If x is within k ulp of the true root, the true root lies in
// [x - k ulp, x + k ulp], and monotone F puts p between F at the two ends
// -- exactly, before F's own rounding. kBracketUlp covers the inverses'
// 1-2 ulp gates with room; kForwardSlack covers F's own error (gamma <= 4,
// beta <= 3 ulp per region gate) with room. A wrong root is off by millions
// of ulp (#48-#50, #55), so the margins cost no sensitivity.

constexpr int kBracketUlp = 4;
constexpr uint64_t kForwardSlack = 8;

void BracketOne(Check* c, const char* fn, bool upper, double a, double b,
                double s, double x, double f_lo, double f_hi) {
  ++c->n;
  if (std::isnan(x)) {
    Fail(c, fn, a, b, s, x, "NaN root");
    return;
  }
  // P side: F(x-) <= p <= F(x+). Q side: Q(x-) >= q >= Q(x+).
  const bool ok = upper ? (LeSlack(s, f_lo, kForwardSlack) &&
                           LeSlack(f_hi, s, kForwardSlack))
                        : (LeSlack(f_lo, s, kForwardSlack) &&
                           LeSlack(s, f_hi, kForwardSlack));
  if (!ok) Fail(c, fn, a, b, s, x, "p outside [F(x-k), F(x+k)]");
}

void BracketGamma(Check* c, const std::vector<double>& a,
                  const std::vector<double>& s) {
  const size_t n = a.size();
  std::vector<double> xp(n), xq(n);
  corvus::gamma_p_inv(a, s, xp);
  corvus::gamma_q_inv(a, s, xq);
  std::vector<double> lo(n), hi(n), flo(n), fhi(n);
  for (int side = 0; side < 2; ++side) {
    const std::vector<double>& x = side == 0 ? xp : xq;
    for (size_t i = 0; i < n; ++i) {
      lo[i] = Step(x[i], kBracketUlp, -1, 0.0, kInf);
      hi[i] = Step(x[i], kBracketUlp, 1, 0.0, kInf);
    }
    if (side == 0) {
      corvus::gamma_p(a, lo, flo);
      corvus::gamma_p(a, hi, fhi);
    } else {
      corvus::gamma_q(a, lo, flo);
      corvus::gamma_q(a, hi, fhi);
    }
    for (size_t i = 0; i < n; ++i)
      BracketOne(c, side == 0 ? "gamma_p_inv" : "gamma_q_inv", side == 1, a[i],
                 0.0, s[i], x[i], flo[i], fhi[i]);
  }
}

void BracketBeta(Check* c, const std::vector<double>& a,
                 const std::vector<double>& b, const std::vector<double>& s) {
  const size_t n = a.size();
  std::vector<double> xp(n), xq(n);
  corvus::beta_p_inv(a, b, s, xp);
  corvus::beta_q_inv(a, b, s, xq);
  std::vector<double> lo(n), hi(n), flo(n), fhi(n);
  for (int side = 0; side < 2; ++side) {
    const std::vector<double>& x = side == 0 ? xp : xq;
    for (size_t i = 0; i < n; ++i) {
      lo[i] = Step(x[i], kBracketUlp, -1, 0.0, 1.0);
      hi[i] = Step(x[i], kBracketUlp, 1, 0.0, 1.0);
    }
    if (side == 0) {
      corvus::beta_p(a, b, lo, flo);
      corvus::beta_p(a, b, hi, fhi);
    } else {
      corvus::beta_q(a, b, lo, flo);
      corvus::beta_q(a, b, hi, fhi);
    }
    for (size_t i = 0; i < n; ++i)
      BracketOne(c, side == 0 ? "beta_p_inv" : "beta_q_inv", side == 1, a[i],
                 b[i], s[i], x[i], flo[i], fhi[i]);
  }
}

void CheckBracket() {
  Check c{"inverse bracket"};
  Rng rng{53};
  {
    std::vector<double> a, s;
    // Seeded log grid over the whole shape range, both tails.
    for (int i = 0; i < 4000; ++i) {
      a.push_back(rng.LogU(1e-20, 1e20));
      const double t = rng.LogU(1e-300, 0.5);
      s.push_back(rng.Unit() < 0.5 ? t : 1.0 - t);
    }
    // Headline points of the gammainv issues (#48, #49).
    for (double ai : {0.1, 1e-5, 0.0928, 0.3, 0.47})
      for (double si : {0.5, 0.9, 0.99, 0.0081835764484083223, 5.8252e-7}) {
        a.push_back(ai);
        s.push_back(si);
      }
    for (double ai : {20.0, 22.5, 26.0})
      for (double si : {1e-300, 1e-250, 1e-220}) {
        a.push_back(ai);
        s.push_back(si);
      }
    // #49's band, seeded: the big-a forward's E with u -> -1 (a just above
    // a_T = 20, tiny p, down to denorm_min).
    for (int i = 0; i < 600; ++i) {
      a.push_back(19.0 + 9.0 * rng.Unit());
      s.push_back(rng.LogU(kDmin, 1e-150));
    }
    BracketGamma(&c, a, s);
  }
  {
    std::vector<double> a, b, s;
    for (int i = 0; i < 6000; ++i) {
      a.push_back(rng.LogU(1e-20, 1e20));
      b.push_back(rng.LogU(1e-20, 1e20));
      const double t = rng.LogU(1e-300, 0.5);
      s.push_back(rng.Unit() < 0.5 ? t : 1.0 - t);
    }
    // A moderate shape against a small one (#50; the small shape reaches
    // 0.1 so the sweep also shows where the band ends) and the huge-a band
    // (#55), seeded.
    for (int i = 0; i < 4000; ++i) {
      a.push_back(rng.LogU(1.0, 1e5));
      b.push_back(rng.LogU(1e-200, 1e-1));
      s.push_back(rng.LogU(1e-300, 1e-8));
    }
    for (int i = 0; i < 500; ++i) {
      a.push_back(rng.LogU(1e16, 1e20));
      b.push_back(rng.LogU(1e-200, 1e-6));
      s.push_back(rng.Unit() < 0.5 ? rng.LogU(1e-300, 0.5) : 0.9);
    }
    // c85e4de's two betainv fixes (huge-parameter prescale, exp-floor
    // clamp): the hunt's mutation test showed no gate failed with either
    // reverted. These are its witness points.
    const double pts[][3] = {{0.001, 1e307, 0.9999999999},
                             {0.001, 1e307, 1e-300},
                             {1e100, 0.001, 1e-20},
                             {1e100, 0.001, 1e-300},
                             {1e100, 0.001, 1e-100}};
    for (const auto& p : pts) {
      a.push_back(p[0]);
      b.push_back(p[1]);
      s.push_back(p[2]);
    }
    // Both shapes tiny -- betainv's backward-error regime and its
    // orientation probe (#55, found under #47). Drawn LAST so the grids
    // above keep their points.
    for (int i = 0; i < 1500; ++i) {
      a.push_back(rng.LogU(1e-300, 1e-16));
      b.push_back(rng.LogU(1e-300, 1e-16));
      const double t = rng.LogU(1e-300, 0.5);
      s.push_back(rng.Unit() < 0.5 ? t : 1.0 - t);
    }
    BracketBeta(&c, a, b, s);
  }
  Report(c);
}

}  // namespace

int main(int argc, char** argv) {
  if (argc > 1 && std::strcmp(argv[1], "--all") == 0) g_print_limit = SIZE_MAX;
  if (!corvus_test::ReportAndCheckTarget()) return 1;
  CheckHarness();
  CheckScalarSpan();
  CheckRange();
  CheckSymmetry();
  CheckBracket();
  ReportBands();
  if (g_fail) {
    std::printf("FAIL: %d check(s) failed\n", g_fail);
    return 1;
  }
  std::printf("PASS\n");
  return 0;
}
