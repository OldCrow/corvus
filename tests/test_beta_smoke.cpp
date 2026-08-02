// Smoke test for corvus::beta_p / corvus::beta_q: the full specials table,
// P + Q = 1, the exact symmetric value I_{1/2}(a,a) = 1/2, and lane-mix
// determinism.
//
// The last one is the interesting one, and it is why this file exists
// separately from the ULP gate. Both summed regions (R1's power series and
// R4's alpha-scaled series) stop per lane on a freeze mask and let the whole
// vector break out once every lane has stopped, so a lane's iteration count is
// decided by its neighbours; the region cores additionally run on ALL lanes
// with their inactive ones scrubbed. If a frozen lane's accumulator were
// merely having zeros added to it rather than being left alone, the extra dd
// renormalizations would be visible in the last bits -- and only ever in mixed
// vectors, which is exactly the case a single-point test never builds. So
// every probe point here is evaluated alone, at several lane offsets, and
// interleaved with points from the other regions, and all of it must be
// bit-identical.
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <vector>

#include "corvus/corvus.h"
#include "expect_target.h"

namespace {

const double kInf = std::numeric_limits<double>::infinity();
const double kNan = std::numeric_limits<double>::quiet_NaN();

int g_fail = 0;

int64_t OrderedBits(double x) {
  int64_t b;
  std::memcpy(&b, &x, sizeof(b));
  return b < 0 ? (INT64_MIN - b) : b;
}
uint64_t UlpDiff(double a, double b) {
  return static_cast<uint64_t>(std::llabs(OrderedBits(a) - OrderedBits(b)));
}
bool SameBits(double a, double b) {
  uint64_t ba, bb;
  std::memcpy(&ba, &a, sizeof(ba));
  std::memcpy(&bb, &b, sizeof(bb));
  return ba == bb;
}

using Fn = void (*)(std::span<const double>, std::span<const double>,
                    std::span<const double>, std::span<double>);

double One(Fn fn, double a, double b, double x) {
  double out = 0.0;
  fn(std::span<const double>(&a, 1), std::span<const double>(&b, 1),
     std::span<const double>(&x, 1), std::span<double>(&out, 1));
  return out;
}

// A NaN expectation means "any NaN"; a zero expectation additionally requires
// the sign to be positive (the design pins +0, matching gamma and C99).
void CheckSpecial(const char* fname, Fn fn, double a, double b, double x,
                  double want) {
  const double got = One(fn, a, b, x);
  bool ok;
  if (std::isnan(want)) {
    ok = std::isnan(got);
  } else if (want == 0.0) {
    ok = got == 0.0 && !std::signbit(got);
  } else {
    ok = got == want;
  }
  if (!ok) {
    std::fprintf(stderr, "FAIL: %s(%g, %g, %g) = %.17g, want %.17g\n", fname, a,
                 b, x, got, want);
    g_fail = 1;
  }
}

void Specials() {
  struct Case {
    double a, b, x, p, q;
    const char* why;
  };
  // PLAN.md, "Specials [pinned now; gamma-consistent doctrine: one degenerate
  // parameter gets its limit, two degeneracies (or a degenerate parameter
  // meeting the x-boundary its mass sits on) -> NaN]".
  const Case cases[] = {
      // x at the ends, ordinary parameters.
      {1.0, 1.0, 0.0, 0.0, 1.0, "x = 0"},
      {0.5, 3.0, 0.0, 0.0, 1.0, "x = 0"},
      {1.0, 1.0, 1.0, 1.0, 0.0, "x = 1"},
      {0.5, 3.0, 1.0, 1.0, 0.0, "x = 1"},
      {1e300, 1e-300, 0.0, 0.0, 1.0, "x = 0, extreme params"},
      {1e-300, 1e300, 1.0, 1.0, 0.0, "x = 1, extreme params"},

      // a = 0: all the mass sits at 0.
      {0.0, 1.0, 0.5, 1.0, 0.0, "a = 0"},
      {0.0, 1.0, 1e-300, 1.0, 0.0, "a = 0"},
      {0.0, 1.0, 1.0, 1.0, 0.0, "a = 0 at x = 1"},
      {0.0, 1.0, 0.0, kNan, kNan, "a = 0 meets its own mass point"},

      // b = +inf: also all the mass at 0.
      {1.0, kInf, 0.5, 1.0, 0.0, "b = +inf"},
      {1.0, kInf, 1.0, 1.0, 0.0, "b = +inf at x = 1"},
      {1.0, kInf, 0.0, kNan, kNan, "b = +inf meets its own mass point"},

      // b = 0: all the mass sits at 1.
      {1.0, 0.0, 0.5, 0.0, 1.0, "b = 0"},
      {1.0, 0.0, 1.0 - 1e-16, 0.0, 1.0, "b = 0"},
      {1.0, 0.0, 0.0, 0.0, 1.0, "b = 0 at x = 0"},
      {1.0, 0.0, 1.0, kNan, kNan, "b = 0 meets its own mass point"},

      // a = +inf: also all the mass at 1.
      {kInf, 1.0, 0.5, 0.0, 1.0, "a = +inf"},
      {kInf, 1.0, 0.0, 0.0, 1.0, "a = +inf at x = 0"},
      {kInf, 1.0, 1.0, kNan, kNan, "a = +inf meets its own mass point"},

      // Two degeneracies.
      {0.0, 0.0, 0.5, kNan, kNan, "a = b = 0"},
      {kInf, kInf, 0.5, kNan, kNan, "a = b = +inf"},
      {0.0, kInf, 0.5, kNan, kNan, "a = 0, b = +inf"},
      {kInf, 0.0, 0.5, kNan, kNan, "a = +inf, b = 0"},

      // Out of domain.
      {-1.0, 1.0, 0.5, kNan, kNan, "a < 0"},
      {1.0, -1.0, 0.5, kNan, kNan, "b < 0"},
      {-1.0, -1.0, 0.5, kNan, kNan, "a, b < 0"},
      {1.0, 1.0, -0.5, kNan, kNan, "x < 0"},
      {1.0, 1.0, 1.5, kNan, kNan, "x > 1"},
      {1.0, 1.0, kInf, kNan, kNan, "x = +inf"},
      {1.0, 1.0, -kInf, kNan, kNan, "x = -inf"},
      {-kInf, 1.0, 0.5, kNan, kNan, "a = -inf"},
      {1.0, -kInf, 0.5, kNan, kNan, "b = -inf"},

      // NaN anywhere.
      {kNan, 1.0, 0.5, kNan, kNan, "a NaN"},
      {1.0, kNan, 0.5, kNan, kNan, "b NaN"},
      {1.0, 1.0, kNan, kNan, kNan, "x NaN"},
      {kNan, kNan, kNan, kNan, kNan, "all NaN"},
      {kNan, 0.0, 0.0, kNan, kNan, "NaN beats every other rule"},
  };
  for (const Case& c : cases) {
    CheckSpecial("beta_p", corvus::beta_p, c.a, c.b, c.x, c.p);
    CheckSpecial("beta_q", corvus::beta_q, c.a, c.b, c.x, c.q);
  }

  // NaN payload propagation, the erf/erfc/lgamma/gamma convention.
  const double payload = std::nan("42");
  const double probes[3][3] = {
      {payload, 1.0, 0.5}, {1.0, payload, 0.5}, {1.0, 1.0, payload}};
  for (const auto& pr : probes) {
    const double got = One(corvus::beta_p, pr[0], pr[1], pr[2]);
    if (!SameBits(got, payload)) {
      std::fprintf(stderr,
                   "FAIL: beta_p did not propagate the NaN payload for "
                   "(%g, %g, %g)\n",
                   pr[0], pr[1], pr[2]);
      g_fail = 1;
    }
  }
}

// A spread of in-domain points covering all four regions, both orientations,
// and both sides of every boundary the router tests.
std::vector<std::array<double, 3>> Probes() {
  return {
      // R4, tiny-first native (a <= b) and swapped (b < a).
      {1e-300, 1.0, 0.3},   {1e-20, 1.0, 0.4},    {0.015625, 4.0, 0.4},
      {0.001, 2.0, 0.1},    {1.0, 1e-300, 0.7},   {5.0, 0.0078125, 0.9},
      {8.0, 0.015625, 1.0 - 9.5e-7},
      // R1, native and swapped.
      {0.5, 3.0, 0.2},      {2.0, 5.0, 0.3},      {1.0, 1.0, 0.25},
      {8.0, 2.0, 0.4},      {3.0, 100.0, 0.05},   {100.0, 3.0, 0.95},
      {2.0, 20.0, 0.35},    {1e-8, 3.0, 0.44},
      // R2, both orientations.
      {1.0, 1.0, 0.6},      {3.0, 30.0, 0.5},     {0.5, 40.0, 0.9},
      {10.0, 10.0, 0.9},    {2.0, 2.0, 0.8},      {30.0, 3.0, 0.5},
      {1e6, 1.0, 1.0 - 1e-5}, {1.0, 1e10, 8e-8},
      // R3 ridge, on it and both sides.
      {100.0, 100.0, 0.5},  {64.0, 64.0, 0.5},    {1e4, 1e4, 0.5},
      {200.0, 100.0, 0.6},  {1e6, 1e6, 0.5000001},{50.0, 150.0, 0.25},
      {1e10, 1e10, 0.5},    {33.0, 40.0, 0.45},   {1000.0, 250.0, 0.75},
      // Boundary brackets: xi1, B1, T_ridge, eps_R4.
      {2.0, 4.0, 0.45},     {2.0, 4.0, 0.4500001},{4.0, 2.0, 0.55},
      {16.0, 16.0, 0.5},    {16.5, 16.5, 0.5},    {0.015625, 8.0, 0.4},
  };
}

void PlusQ() {
  const auto probes = Probes();
  std::vector<double> a, b, x;
  for (const auto& pr : probes) {
    a.push_back(pr[0]);
    b.push_back(pr[1]);
    x.push_back(pr[2]);
  }
  std::vector<double> p(a.size()), q(a.size());
  corvus::beta_p(a, b, x, p);
  corvus::beta_q(a, b, x, q);
  for (size_t i = 0; i < a.size(); ++i) {
    const uint64_t u = UlpDiff(p[i] + q[i], 1.0);
    if (u > 1) {
      std::fprintf(stderr,
                   "FAIL: P + Q off by %llu ULP at a=%.17g b=%.17g x=%.17g "
                   "(P=%.17g Q=%.17g)\n",
                   static_cast<unsigned long long>(u), a[i], b[i], x[i], p[i],
                   q[i]);
      g_fail = 1;
    }
  }
}

// I_{1/2}(a, a) = 1/2 EXACTLY, for every a: the integrand is symmetric about
// 1/2. This is the kernel's most brutal single check -- it says the routed
// evaluation is correctly rounded at a point where the true value is a
// representable double, in every region the diagonal crosses.
void SymmetricHalf() {
  const double as[] = {1e-300, 1e-8, 0.0078125, 0.015625, 0.5,  1.0,
                       2.0,    5.0,  20.0,      63.0,     64.0, 100.0,
                       1000.0, 1e6,  1e10,      1e100};
  std::vector<double> a, b, x;
  for (double v : as) {
    a.push_back(v);
    b.push_back(v);
    x.push_back(0.5);
  }
  std::vector<double> p(a.size()), q(a.size());
  corvus::beta_p(a, b, x, p);
  corvus::beta_q(a, b, x, q);
  for (size_t i = 0; i < a.size(); ++i) {
    if (p[i] != 0.5 || q[i] != 0.5) {
      std::fprintf(stderr,
                   "FAIL: I_1/2(%.17g, %.17g) = P %.17g / Q %.17g, want 0.5 "
                   "exactly (%llu / %llu ULP off)\n",
                   a[i], a[i], p[i], q[i],
                   static_cast<unsigned long long>(UlpDiff(p[i], 0.5)),
                   static_cast<unsigned long long>(UlpDiff(q[i], 0.5)));
      g_fail = 1;
    }
  }
}

// Evaluate every probe alone, then in batches that place it at a different
// lane offset each time and surround it with points from other regions.
void LaneMix() {
  const auto probes = Probes();
  const size_t n = probes.size();

  for (int fn_idx = 0; fn_idx < 2; ++fn_idx) {
    Fn fn = fn_idx == 0 ? corvus::beta_p : corvus::beta_q;
    const char* name = fn_idx == 0 ? "beta_p" : "beta_q";

    std::vector<double> alone(n);
    for (size_t i = 0; i < n; ++i) {
      alone[i] = One(fn, probes[i][0], probes[i][1], probes[i][2]);
    }

    // Offsets 0..8 cover every lane position on every tier up to AVX-512, and
    // the leading filler is drawn from the far end of the probe list so each
    // probe sits next to a different region's point every time.
    for (size_t off = 0; off < 9; ++off) {
      std::vector<double> a, b, x;
      for (size_t k = 0; k < off; ++k) {
        a.push_back(probes[(n - 1 - k) % n][0]);
        b.push_back(probes[(n - 1 - k) % n][1]);
        x.push_back(probes[(n - 1 - k) % n][2]);
      }
      for (const auto& pr : probes) {
        a.push_back(pr[0]);
        b.push_back(pr[1]);
        x.push_back(pr[2]);
      }
      std::vector<double> out(a.size());
      fn(a, b, x, out);
      for (size_t i = 0; i < n; ++i) {
        if (!SameBits(out[off + i], alone[i])) {
          std::fprintf(stderr,
                       "FAIL: %s not lane-position invariant at a=%.17g "
                       "b=%.17g x=%.17g (alone %.17g, offset %zu %.17g)\n",
                       name, probes[i][0], probes[i][1], probes[i][2], alone[i],
                       off, out[off + i]);
          g_fail = 1;
        }
      }
    }

    // Reversed order: same points, every neighbour different.
    std::vector<double> ra, rb, rx;
    for (size_t i = n; i-- > 0;) {
      ra.push_back(probes[i][0]);
      rb.push_back(probes[i][1]);
      rx.push_back(probes[i][2]);
    }
    std::vector<double> rout(n);
    fn(ra, rb, rx, rout);
    for (size_t i = 0; i < n; ++i) {
      if (!SameBits(rout[n - 1 - i], alone[i])) {
        std::fprintf(stderr,
                     "FAIL: %s not order invariant at a=%.17g b=%.17g x=%.17g "
                     "(alone %.17g, reversed %.17g)\n",
                     name, probes[i][0], probes[i][1], probes[i][2], alone[i],
                     rout[n - 1 - i]);
        g_fail = 1;
      }
    }
  }
}

// Aliasing is part of the public contract, and the length is deliberately not
// a multiple of any lane count so the masked tail runs.
void Aliasing() {
  const auto probes = Probes();
  std::vector<double> a, b, x;
  for (const auto& pr : probes) {
    a.push_back(pr[0]);
    b.push_back(pr[1]);
    x.push_back(pr[2]);
  }
  std::vector<double> want(a.size());
  corvus::beta_p(a, b, x, want);
  std::vector<double> buf = x;
  corvus::beta_p(a, b, buf, buf);
  for (size_t i = 0; i < a.size(); ++i) {
    if (!SameBits(buf[i], want[i])) {
      std::fprintf(stderr,
                   "FAIL: beta_p in-place differs at a=%.17g b=%.17g x=%.17g\n",
                   a[i], b[i], x[i]);
      g_fail = 1;
    }
  }
}

// The three analytic lines the design names, as a coarse correctness net that
// does not depend on the reference file: I_x(a,1) = x^a, I_x(1,b) = 1-(1-x)^b,
// I_x(1/2,1/2) = (2/pi) asin(sqrt(x)). Tolerance is loose (4 ULP) because the
// libm right-hand sides carry their own error.
void AnalyticLines() {
  struct Case {
    double a, b, x, want;
  };
  std::vector<Case> cs;
  for (double xx : {0.05, 0.25, 0.5, 0.75, 0.99}) {
    for (double aa : {0.5, 1.0, 2.5, 10.0}) {
      cs.push_back({aa, 1.0, xx, std::pow(xx, aa)});
      cs.push_back({1.0, aa, xx, -std::expm1(aa * std::log1p(-xx))});
    }
    cs.push_back({0.5, 0.5, xx,
                  2.0 / 3.14159265358979323846 * std::asin(std::sqrt(xx))});
  }
  std::vector<double> a, b, x;
  for (const Case& c : cs) {
    a.push_back(c.a);
    b.push_back(c.b);
    x.push_back(c.x);
  }
  std::vector<double> got(a.size());
  corvus::beta_p(a, b, x, got);
  for (size_t i = 0; i < cs.size(); ++i) {
    const uint64_t u = UlpDiff(got[i], cs[i].want);
    if (u > 4) {
      std::fprintf(stderr,
                   "FAIL: analytic line off by %llu ULP at a=%.17g b=%.17g "
                   "x=%.17g (got %.17g want %.17g)\n",
                   static_cast<unsigned long long>(u), cs[i].a, cs[i].b,
                   cs[i].x, got[i], cs[i].want);
      g_fail = 1;
    }
  }
}

}  // namespace

int main() {
  if (!corvus_test::ReportAndCheckTarget()) return 2;
  Specials();
  PlusQ();
  SymmetricHalf();
  AnalyticLines();
  LaneMix();
  Aliasing();
  if (g_fail == 0) std::printf("PASS: beta_p/beta_q smoke tests\n");
  return g_fail;
}
