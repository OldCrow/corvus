// Smoke test for corvus::gamma_p / corvus::gamma_q: the specials table, the
// P + Q = 1 identity, and lane-mix determinism.
//
// The third one is the interesting one. Both summed regions (R1's series and
// R4's alternating series) stop per lane on a freeze mask and let the whole
// vector break out once every lane has stopped, so a lane's iteration count
// is decided by its neighbours. If a frozen lane's accumulator were merely
// having zeros added to it rather than being left alone, the extra dd
// renormalizations would be visible in the last bits -- and only ever in
// mixed vectors, which is exactly the case a single-point test never builds.
// So every probe point here is evaluated alone, at several lane offsets, and
// interleaved with points from the other regions, and all of it must be
// bit-identical.
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

double One(void (*fn)(std::span<const double>, std::span<const double>,
                      std::span<double>),
           double a, double x) {
  double out = 0.0;
  fn(std::span<const double>(&a, 1), std::span<const double>(&x, 1),
     std::span<double>(&out, 1));
  return out;
}

// Expected value for a specials entry. A NaN expectation means "any NaN".
void CheckSpecial(const char* fname,
                  void (*fn)(std::span<const double>, std::span<const double>,
                             std::span<double>),
                  double a, double x, double want, bool want_pos_zero) {
  const double got = One(fn, a, x);
  bool ok;
  if (std::isnan(want)) {
    ok = std::isnan(got);
  } else if (want == 0.0 && want_pos_zero) {
    ok = got == 0.0 && !std::signbit(got);
  } else {
    ok = got == want;
  }
  if (!ok) {
    std::fprintf(stderr, "FAIL: %s(%g, %g) = %g, want %g\n", fname, a, x, got,
                 want);
    g_fail = 1;
  }
}

void Specials() {
  struct Case {
    double a, x, p, q;
  };
  // Specials table (SciPy limits).
  const Case cases[] = {
      {1.0, 0.0, 0.0, 1.0},     // no mass yet
      {0.5, 0.0, 0.0, 1.0},
      {0.0, 1.0, 1.0, 0.0},     // Gamma(0) is infinite: all the mass is below
      {0.0, 1e-30, 1.0, 0.0},
      {0.0, 0.0, kNan, kNan},
      {1.0, kInf, 1.0, 0.0},
      {1e250, kInf, 1.0, 0.0},
      {kInf, 1.0, 0.0, 1.0},
      {kInf, 1e250, 0.0, 1.0},
      {kInf, kInf, kNan, kNan},
      {-1.0, 1.0, kNan, kNan},
      {-0.5, 2.0, kNan, kNan},
      {1.0, -1.0, kNan, kNan},
      {-1.0, -1.0, kNan, kNan},
      {kNan, 1.0, kNan, kNan},
      {1.0, kNan, kNan, kNan},
      {kNan, kNan, kNan, kNan},
      {-kInf, 1.0, kNan, kNan},
      {1.0, -kInf, kNan, kNan},
  };
  for (const Case& c : cases) {
    CheckSpecial("gamma_p", corvus::gamma_p, c.a, c.x, c.p, true);
    CheckSpecial("gamma_q", corvus::gamma_q, c.a, c.x, c.q, true);
  }

  // NaN payload propagation, the erf/erfc/lgamma convention.
  const double payload = std::nan("42");
  for (int which = 0; which < 2; ++which) {
    const double got = which == 0 ? One(corvus::gamma_p, payload, 1.0)
                                  : One(corvus::gamma_p, 1.0, payload);
    if (!SameBits(got, payload)) {
      std::fprintf(stderr, "FAIL: gamma_p did not propagate the NaN payload\n");
      g_fail = 1;
    }
  }
}

// A spread of in-domain points covering all four regions and both sides of
// every boundary the router tests.
std::vector<std::pair<double, double>> Probes() {
  return {
      // R4 / R1 small-a corner (a <= 3/2, x <= 4), incl. the e^-gamma line.
      {1e-300, 1.0},   {1e-30, 0.5},    {0.25, 0.1},
      {0.5, 0.5},      {0.5, 4.0},      {1.0, 0.5},
      {1.0, 2.0},      {1.5, 4.0},      {1.5, 2.5},
      {1.0, 0.5614594835668851},  // x = e^-gamma, the R4 cancellation line
      // R1 small-a beyond the box.
      {3.0, 4.0},      {7.5, 8.5},      {19.0, 20.0},
      // R2 small-a.
      {1.0, 8.0},      {3.0, 30.0},     {19.0, 60.0},   {0.5, 40.0},
      // R1 / R2 large-a wings.
      {25.0, 10.0},    {100.0, 20.0},   {1000.0, 400.0},
      {25.0, 60.0},    {100.0, 250.0},  {1000.0, 2200.0},
      // R3 ridge, both sides and on it.
      {20.0, 20.0},    {20.0, 15.0},    {20.0, 30.0},
      {100.0, 99.0},   {100.0, 101.0},  {1e6, 1e6},
      {1e6, 1000000.5},{1e10, 1e10},    {1e10, 1.0000001e10},
      // R3 boundary brackets.
      {20.0, 10.0},    {20.0, 40.0},    {1000.0, 1900.0},
  };
}

void PlusQ() {
  const auto probes = Probes();
  std::vector<double> a, x;
  for (const auto& pr : probes) {
    a.push_back(pr.first);
    x.push_back(pr.second);
  }
  std::vector<double> p(a.size()), q(a.size());
  corvus::gamma_p(a, x, p);
  corvus::gamma_q(a, x, q);
  for (size_t i = 0; i < a.size(); ++i) {
    const uint64_t u = UlpDiff(p[i] + q[i], 1.0);
    if (u > 1) {
      std::fprintf(stderr,
                   "FAIL: P + Q off by %llu ULP at a=%.17g x=%.17g "
                   "(P=%.17g Q=%.17g)\n",
                   static_cast<unsigned long long>(u), a[i], x[i], p[i], q[i]);
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
    auto fn = fn_idx == 0 ? corvus::gamma_p : corvus::gamma_q;
    const char* name = fn_idx == 0 ? "gamma_p" : "gamma_q";

    std::vector<double> alone(n);
    for (size_t i = 0; i < n; ++i) {
      alone[i] = One(fn, probes[i].first, probes[i].second);
    }

    // Offsets 0..8 cover every lane position on every tier up to AVX-512,
    // and the leading filler is drawn from the far end of the probe list so
    // each probe sits next to a different region's point every time.
    for (size_t off = 0; off < 9; ++off) {
      std::vector<double> a, x;
      for (size_t k = 0; k < off; ++k) {
        a.push_back(probes[(n - 1 - k) % n].first);
        x.push_back(probes[(n - 1 - k) % n].second);
      }
      for (const auto& pr : probes) {
        a.push_back(pr.first);
        x.push_back(pr.second);
      }
      std::vector<double> out(a.size());
      fn(a, x, out);
      for (size_t i = 0; i < n; ++i) {
        if (!SameBits(out[off + i], alone[i])) {
          std::fprintf(stderr,
                       "FAIL: %s not lane-position invariant at a=%.17g "
                       "x=%.17g (alone %.17g, offset %zu %.17g)\n",
                       name, probes[i].first, probes[i].second, alone[i], off,
                       out[off + i]);
          g_fail = 1;
        }
      }
    }

    // Reversed order: same points, every neighbour different.
    std::vector<double> ra, rx;
    for (size_t i = n; i-- > 0;) {
      ra.push_back(probes[i].first);
      rx.push_back(probes[i].second);
    }
    std::vector<double> rout(n);
    fn(ra, rx, rout);
    for (size_t i = 0; i < n; ++i) {
      if (!SameBits(rout[n - 1 - i], alone[i])) {
        std::fprintf(stderr,
                     "FAIL: %s not order invariant at a=%.17g x=%.17g "
                     "(alone %.17g, reversed %.17g)\n",
                     name, probes[i].first, probes[i].second, alone[i],
                     rout[n - 1 - i]);
        g_fail = 1;
      }
    }
  }
}

// The huge-a diagonal (#12): P(a, a) = Q(a, a) ~ 1/2 all the way to
// DBL_MAX. Above 2^996 the non-FMA tiers' Dekker split of a overflows in
// GammaTemme's DdRecip/DdMulD unless the operands are prescaled, and the
// laundered NaN came back as an exact, plausible-looking (1, 0) pair. Only
// the exact diagonal is reachable at these magnitudes (x = a +/- 1 ulp is
// genuinely saturated), so the interval assert is deliberately loose: the
// failure mode is a wrong exact 1/0, not a ULP slip. 2^996 sits on the
// last correct-before-the-fix boundary and serves as the control point.
void HugeDiagonal() {
  const double diag[] = {0x1.0p+996, 0x1.0p+997, 0x1.0p+1000,
                         std::numeric_limits<double>::max()};
  for (double a : diag) {
    const double p = One(corvus::gamma_p, a, a);
    const double q = One(corvus::gamma_q, a, a);
    if (!(p > 0.4 && p < 0.6) || !(q > 0.4 && q < 0.6)) {
      std::fprintf(stderr,
                   "FAIL: gamma diagonal a=%.17g: P=%.17g Q=%.17g "
                   "(want both in (0.4, 0.6))\n",
                   a, p, q);
      g_fail = 1;
    }
  }
}

// Aliasing is part of the public contract, and the length is deliberately
// not a multiple of any lane count so the masked tail runs.
void Aliasing() {
  const auto probes = Probes();
  std::vector<double> a, x;
  for (const auto& pr : probes) {
    a.push_back(pr.first);
    x.push_back(pr.second);
  }
  std::vector<double> want(a.size());
  corvus::gamma_p(a, x, want);
  std::vector<double> buf = x;
  corvus::gamma_p(a, buf, buf);
  for (size_t i = 0; i < a.size(); ++i) {
    if (!SameBits(buf[i], want[i])) {
      std::fprintf(stderr, "FAIL: gamma_p in-place differs at a=%.17g x=%.17g\n",
                   a[i], x[i]);
      g_fail = 1;
    }
  }
}

}  // namespace

int main() {
  if (!corvus_test::ReportAndCheckTarget()) return 2;
  Specials();
  PlusQ();
  LaneMix();
  HugeDiagonal();
  Aliasing();
  if (g_fail == 0) std::printf("PASS: gamma_p/gamma_q smoke tests\n");
  return g_fail;
}
