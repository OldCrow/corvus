// Smoke: lbeta domain policy, specials, symmetry, exact values, the -inf
// saturation edge, and the masked-tail path (odd vector lengths). The
// dense accuracy work lives in test_lbeta_ulp.cpp.
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <vector>

#include "corvus/corvus.h"
#include "expect_target.h"

namespace {

int g_fail = 0;

void Check(bool ok, const char* what) {
  if (!ok) {
    std::fprintf(stderr, "FAIL: %s\n", what);
    ++g_fail;
  }
}

bool SameBits(double x, double y) {
  return std::memcmp(&x, &y, sizeof(x)) == 0;
}

double One(double a, double b) {
  // Length-3 vectors (not a lane-count multiple) so the masked tail runs.
  std::vector<double> va{a, 1.0, a}, vb{b, 1.0, b}, out(3);
  corvus::lbeta(va, vb, out);
  // Bitwise: NaN lanes must agree too (== would reject NaN == NaN).
  Check(SameBits(out[0], out[2]) || (std::isnan(out[0]) && std::isnan(out[2])),
        "masked-tail lane disagrees with full lane");
  return out[0];
}

}  // namespace

int main() {
  if (!corvus_test::ReportAndCheckTarget()) return 2;

  const double inf = std::numeric_limits<double>::infinity();
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double dmax = (std::numeric_limits<double>::max)();

  // B(1,1) = 1: ln B = 0 is ON the zero curve, where the contract is the
  // 2^-53 absolute band, not exactness (the dd assembly leaves ~1e-32).
  Check(std::fabs(One(1.0, 1.0)) <= 0x1.0p-53, "lbeta(1,1) within 2^-53 of 0");
  Check(One(1.0, 2.0) == -0x1.62e42fefa39efp-1, "lbeta(1,2) == ln(1/2)");

  // Symmetry on a small grid spanning the bands.
  const double grid[] = {1e-300, 0.5, 1.0, 7.5, 256.5, 1e5, 1e100, 0x1.1p+990};
  for (double a : grid) {
    for (double b : grid) {
      if (One(a, b) != One(b, a)) {
        Check(false, "symmetry lbeta(a,b) == lbeta(b,a)");
      }
    }
  }

  // Domain policy: non-positive, inf, NaN -> NaN.
  Check(std::isnan(One(0.0, 1.0)), "lbeta(0,b) is NaN");
  Check(std::isnan(One(-1.0, 1.0)), "lbeta(-1,b) is NaN");
  Check(std::isnan(One(1.0, -0.5)), "lbeta(a,<0) is NaN");
  Check(std::isnan(One(inf, 1.0)), "lbeta(inf,b) is NaN");
  Check(std::isnan(One(1.0, inf)), "lbeta(a,inf) is NaN");
  Check(std::isnan(One(nan, 1.0)), "lbeta(NaN,b) is NaN");
  Check(std::isnan(One(1.0, nan)), "lbeta(a,NaN) is NaN");

  // Tiny parameter: lbeta(m, b) ~ -ln m, positive and finite.
  Check(One(1e-300, 2.0) > 690.0 && One(1e-300, 2.0) < 692.0,
        "lbeta(1e-300, 2) ~ 690.8");

  // Huge parameters: finite deep in the big band, -inf at the far edge
  // (true ln B(DBL_MAX, DBL_MAX) ~ -2.5e308 < -DBL_MAX).
  const double big = One(0x1.1p+990, 0x1.1p+990);
  Check(std::isfinite(big) && big < -0x1.0p+989, "big band finite and large");
  Check(One(dmax, dmax) == -inf, "lbeta(DBL_MAX, DBL_MAX) == -inf");

  // Monotonic sanity: increasing b decreases ln B (fixed a).
  Check(One(2.5, 3.0) > One(2.5, 4.0), "monotone decreasing in b");

  if (g_fail == 0) {
    std::printf("PASS: lbeta smoke\n");
    return 0;
  }
  return 1;
}
