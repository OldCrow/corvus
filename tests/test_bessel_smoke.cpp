// Smoke test for corvus::i0/i1/i0e/i1e: the full specials table (including
// the sign of i1/i1e at -0 and at +-inf), overflow-boundary exactness (last
// finite double -> finite, one ULP above -> inf, signed for i1), a handful
// of known values spanning both regions and both signs, the even/odd
// identities on a small grid, and lane-mix determinism. The ULP gate lives
// in test_bessel_ulp.
//
// EVERY VECTOR BUILT HERE HAS A LENGTH THAT IS NOT A MULTIPLE OF 2, 4 OR 8
// (the masked-tail rule, AGENTS.md): the widest validated tier is 8 lanes
// (AVX-512), so a length divisible by 2/4/8 could complete without ever
// exercising the masked tail path on some tier.
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <span>
#include <vector>

#include "corvus/corvus.h"
#include "expect_target.h"

namespace {

const double kInf = std::numeric_limits<double>::infinity();

// Mirrors src/bessel_data.h's pinned overflow boundaries (last FINITE
// double). Hardcoded here rather than including the internal data header --
// the house smoke-test convention (test_trigamma_smoke.cpp,
// test_betainv_smoke.cpp) keeps smoke tests to the PUBLIC surface; the ULP
// gate is what re-derives routing from the internal header.
constexpr double kI0OverflowX = 0x1.64fe5304e83e4p+9;
constexpr double kI1OverflowX = 0x1.64fe69ff9fec7p+9;

int g_fail = 0;

int64_t OrderedBits(double x) {
  int64_t b;
  std::memcpy(&b, &x, sizeof(b));
  return b < 0 ? (INT64_MIN - b) : b;
}
uint64_t UlpDiff(double a, double b) {
  if (std::isnan(a) || std::isnan(b)) return UINT64_MAX;
  return static_cast<uint64_t>(std::llabs(OrderedBits(a) - OrderedBits(b)));
}
bool SameBits(double a, double b) {
  uint64_t ba, bb;
  std::memcpy(&ba, &a, sizeof(ba));
  std::memcpy(&bb, &b, sizeof(bb));
  return ba == bb;
}

double I0(double x) {
  double out = 0.0;
  corvus::i0(std::span<const double>(&x, 1), std::span<double>(&out, 1));
  return out;
}
double I1(double x) {
  double out = 0.0;
  corvus::i1(std::span<const double>(&x, 1), std::span<double>(&out, 1));
  return out;
}
double I0e(double x) {
  double out = 0.0;
  corvus::i0e(std::span<const double>(&x, 1), std::span<double>(&out, 1));
  return out;
}
double I1e(double x) {
  double out = 0.0;
  corvus::i1e(std::span<const double>(&x, 1), std::span<double>(&out, 1));
  return out;
}

// A NaN expectation means "any NaN"; everything else is matched WITH its
// sign, which is the whole content of i1/i1e's +-0 and +-inf specials.
void CheckSpecial(const char* what, double got, double want) {
  bool ok;
  if (std::isnan(want)) {
    ok = std::isnan(got);
  } else {
    ok = got == want && std::signbit(got) == std::signbit(want);
  }
  if (!ok) {
    std::fprintf(stderr, "FAIL: %s: got %.17g, want %.17g\n", what, got,
                want);
    g_fail = 1;
  }
}

void Specials() {
  // x = +-0.
  CheckSpecial("i0(+0)", I0(0.0), 1.0);
  CheckSpecial("i0(-0)", I0(-0.0), 1.0);
  CheckSpecial("i0e(+0)", I0e(0.0), 1.0);
  CheckSpecial("i0e(-0)", I0e(-0.0), 1.0);
  CheckSpecial("i1(+0)", I1(0.0), 0.0);
  CheckSpecial("i1(-0)", I1(-0.0), -0.0);
  CheckSpecial("i1e(+0)", I1e(0.0), 0.0);
  CheckSpecial("i1e(-0)", I1e(-0.0), -0.0);

  // x = +-inf.
  CheckSpecial("i0(+inf)", I0(kInf), kInf);
  CheckSpecial("i0(-inf)", I0(-kInf), kInf);
  CheckSpecial("i0e(+inf)", I0e(kInf), 0.0);
  CheckSpecial("i0e(-inf)", I0e(-kInf), 0.0);
  CheckSpecial("i1(+inf)", I1(kInf), kInf);
  CheckSpecial("i1(-inf)", I1(-kInf), -kInf);
  CheckSpecial("i1e(+inf)", I1e(kInf), 0.0);
  CheckSpecial("i1e(-inf)", I1e(-kInf), -0.0);

  // NaN propagates WITH its payload, matching every other corvus kernel.
  uint64_t bits = 0x7FF8000000ABCDEFULL;
  double payload;
  std::memcpy(&payload, &bits, sizeof(payload));
  if (!SameBits(I0(payload), payload)) {
    std::fprintf(stderr, "FAIL: i0 did not propagate the NaN payload\n");
    g_fail = 1;
  }
  if (!SameBits(I1(payload), payload)) {
    std::fprintf(stderr, "FAIL: i1 did not propagate the NaN payload\n");
    g_fail = 1;
  }
  if (!SameBits(I0e(payload), payload)) {
    std::fprintf(stderr, "FAIL: i0e did not propagate the NaN payload\n");
    g_fail = 1;
  }
  if (!SameBits(I1e(payload), payload)) {
    std::fprintf(stderr, "FAIL: i1e did not propagate the NaN payload\n");
    g_fail = 1;
  }
}

// Last finite double -> finite; one ULP above -> inf (signed for i1/i1e's
// odd pair, which saturates by sign; i0/i0e are even and saturate to +inf
// regardless of the sign of x).
void OverflowBoundary() {
  const double i0_last = I0(kI0OverflowX);
  const double i0_next = I0(std::nextafter(kI0OverflowX, kInf));
  if (!std::isfinite(i0_last) || i0_next != kInf) {
    std::fprintf(stderr,
                 "FAIL: i0 overflow boundary: i0(last)=%.17g (want finite), "
                 "i0(next)=%.17g (want inf)\n",
                 i0_last, i0_next);
    g_fail = 1;
  }
  // Same magnitude boundary on the negative side; i0 is even, so it also
  // saturates to +inf (never -inf) there.
  const double i0_last_neg = I0(-kI0OverflowX);
  const double i0_next_neg = I0(std::nextafter(-kI0OverflowX, -kInf));
  if (!std::isfinite(i0_last_neg) || i0_next_neg != kInf) {
    std::fprintf(stderr,
                 "FAIL: i0 negative overflow boundary: i0(last)=%.17g, "
                 "i0(next)=%.17g (want +inf)\n",
                 i0_last_neg, i0_next_neg);
    g_fail = 1;
  }

  const double i1_last = I1(kI1OverflowX);
  const double i1_next = I1(std::nextafter(kI1OverflowX, kInf));
  if (!std::isfinite(i1_last) || i1_next != kInf) {
    std::fprintf(stderr,
                 "FAIL: i1 overflow boundary: i1(last)=%.17g (want finite), "
                 "i1(next)=%.17g (want inf)\n",
                 i1_last, i1_next);
    g_fail = 1;
  }
  const double i1_last_neg = I1(-kI1OverflowX);
  const double i1_next_neg = I1(std::nextafter(-kI1OverflowX, -kInf));
  if (!std::isfinite(i1_last_neg) || i1_next_neg != -kInf) {
    std::fprintf(stderr,
                 "FAIL: i1 negative overflow boundary: i1(last)=%.17g, "
                 "i1(next)=%.17g (want -inf)\n",
                 i1_last_neg, i1_next_neg);
    g_fail = 1;
  }
}

// Known values pulled directly from the trusted mpmath-generated reference
// rows (tests/data/i{0,1}{,e}_reference.txt), spanning series/tail and both
// signs. Not an accuracy gate -- test_bessel_ulp owns that -- just a check
// that each region is wired to the right formula and sign.
void KnownValues() {
  struct Point {
    double x;
    double want;
    const char* what;
  };
  const Point kPoints[] = {
      // Series region, including the exact split seam x = 8.
      {0x1.0000000000000p+3, 0x1.ab9069e3504fap+8, "i0 seam x=8"},
      {-0x1.0000000000000p+3, 0x1.ab9069e3504fap+8, "i0 seam x=-8 (even)"},
      {0x1.0000000000000p+3, 0x1.8fdf85e46607cp+8, "i1 seam x=8"},
      {0x1.01371c9a75f8ap+2, 0x1.a6d87ab067e66p-3, "i0e series"},
      // Tiny x: i1 -> round(x/2), including subnormals.
      {-0x1.8a46e4a15fb8ap-794, -0x1.8a46e4a15fb8ap-795, "i1 tiny -> x/2"},
      // Tail region, both signs.
      {0x1.1184775f691a1p+3, 0x1.65240223d4dfap+9, "i0 tail"},
      {-0x1.01f00f80245b2p+9, 0x1.568b79ea0efbfp+738, "i0 deep tail (neg)"},
      {0x1.61ca0ba64e8dep+3, 0x1.ca3409c34325bp+12, "i1 tail"},
      {0x1.3be49a0933c7dp+3, 0x1.f3a35f6a25e23p-4, "i1e tail"},
      {-0x1.2ff6ca48103d8p+992, 0x1.76e72ad308e60p-498,
       "i0e deep tail (neg, near DBL_MAX)"},
      {-0x1.4e718d7d76254p+664, -0x1.65695b56fc951p-334,
       "i1e deep tail (neg)"},
  };
  for (const Point& p : kPoints) {
    // Dispatch on the leading character of `what` ("i0 "/"i1 "/"i0e"/"i1e")
    // is fragile; instead each row is checked against ALL FOUR functions
    // is wasteful -- so the table above is grouped by which function it
    // names, matched by substring.
    double got;
    if (std::strstr(p.what, "i0e") == p.what) {
      got = I0e(p.x);
    } else if (std::strstr(p.what, "i1e") == p.what) {
      got = I1e(p.x);
    } else if (std::strstr(p.what, "i0") == p.what) {
      got = I0(p.x);
    } else {
      got = I1(p.x);
    }
    const uint64_t u = UlpDiff(got, p.want);
    if (u > 2) {
      std::fprintf(stderr, "FAIL: %s: got %.17g, want %.17g (%llu ULP)\n",
                   p.what, got, p.want, static_cast<unsigned long long>(u));
      g_fail = 1;
    }
  }
}

// Even/odd identities on a small grid spanning both regions, and the finite
// von Mises mean-resultant ratio i1e/i0e. Length 13: not a multiple of
// 2/4/8.
void Identities() {
  const double kXs[] = {0.001, 0.25, 1.0,  2.5, 5.0,   7.999, 8.0,
                        8.001, 15.0, 50.0, 200, 500.0, 700.0};
  static_assert(sizeof(kXs) / sizeof(kXs[0]) == 13, "length check");
  for (double x : kXs) {
    const double i0p = I0(x), i0n = I0(-x);
    const double i1p = I1(x), i1n = I1(-x);
    const double i0ep = I0e(x), i0en = I0e(-x);
    const double i1ep = I1e(x), i1en = I1e(-x);
    if (!SameBits(i0p, i0n)) {
      std::fprintf(stderr, "FAIL: i0(-%.4g) != i0(%.4g): %.17g vs %.17g\n", x,
                   x, i0n, i0p);
      g_fail = 1;
    }
    if (!SameBits(i0ep, i0en)) {
      std::fprintf(stderr, "FAIL: i0e(-%.4g) != i0e(%.4g): %.17g vs %.17g\n",
                   x, x, i0en, i0ep);
      g_fail = 1;
    }
    if (!SameBits(i1n, -i1p)) {
      std::fprintf(stderr, "FAIL: i1(-%.4g) != -i1(%.4g): %.17g vs %.17g\n",
                   x, x, i1n, -i1p);
      g_fail = 1;
    }
    if (!SameBits(i1en, -i1ep)) {
      std::fprintf(stderr, "FAIL: i1e(-%.4g) != -i1e(%.4g): %.17g vs %.17g\n",
                   x, x, i1en, -i1ep);
      g_fail = 1;
    }
    const double ratio = i1ep / i0ep;  // A(kappa): finite, in [0, 1)
    if (!(std::isfinite(ratio) && ratio >= 0.0 && ratio < 1.0)) {
      std::fprintf(stderr,
                   "FAIL: i1e/i0e ratio at x=%.4g not in [0,1): %.17g\n", x,
                   ratio);
      g_fail = 1;
    }
  }
}

// Every probe alone, then at each lane offset inside a vector of unrelated
// points, then interleaved with points from other regimes. All
// bit-identical. Lengths 11 and 9: not multiples of 2/4/8.
void LaneMix() {
  const double pts[] = {
      -700.5,   -8.0000001, -8.0, -7.9999999, -0.5,      -1e-300,
      0.0,      1e-300,     0.5,  7.9999999,  8.0,
  };
  const double filler[] = {3.0, -3.0, 600.0, -600.0, 0.0, 1e-10};
  constexpr size_t kNPts = sizeof(pts) / sizeof(pts[0]);
  static_assert(kNPts == 11, "length check");

  double alone_i0[kNPts], alone_i1[kNPts], alone_i0e[kNPts], alone_i1e[kNPts];
  for (size_t i = 0; i < kNPts; ++i) {
    alone_i0[i] = I0(pts[i]);
    alone_i1[i] = I1(pts[i]);
    alone_i0e[i] = I0e(pts[i]);
    alone_i1e[i] = I1e(pts[i]);
  }

  constexpr size_t kN = 9;  // not a multiple of any lane count in the fleet
  for (size_t i = 0; i < kNPts; ++i) {
    for (size_t off = 0; off < kN; ++off) {
      std::vector<double> in(kN), g0(kN), g1(kN), g0e(kN), g1e(kN);
      for (size_t j = 0; j < kN; ++j) in[j] = filler[(j + i) % 6];
      in[off] = pts[i];
      corvus::i0(in, g0);
      corvus::i1(in, g1);
      corvus::i0e(in, g0e);
      corvus::i1e(in, g1e);
      if (!SameBits(g0[off], alone_i0[i]) || !SameBits(g1[off], alone_i1[i]) ||
          !SameBits(g0e[off], alone_i0e[i]) ||
          !SameBits(g1e[off], alone_i1e[i])) {
        std::fprintf(stderr,
                     "FAIL: lane mix changed the answer at x=%.17g "
                     "(offset %zu)\n",
                     pts[i], off);
        g_fail = 1;
      }
    }
  }
}

}  // namespace

int main() {
  if (!corvus_test::ReportAndCheckTarget()) return 2;
  Specials();
  OverflowBoundary();
  KnownValues();
  Identities();
  LaneMix();
  if (g_fail == 0) std::printf("PASS: bessel smoke\n");
  return g_fail;
}
