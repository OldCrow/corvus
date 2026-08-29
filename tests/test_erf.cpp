#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>

#include "corvus/corvus.h"
#include "expect_target.h"
#include "ulp_utils.h"

namespace {

using corvus_test::SameBits;

int failures = 0;

void Check(bool ok, const char* what) {
  if (!ok) {
    std::fprintf(stderr, "FAIL: %s\n", what);
    ++failures;
  }
}

}  // namespace

int main() {
  if (!corvus_test::ReportAndCheckTarget()) return 2;

  // Table kernel is 1-ULP class; std::erf itself may be ~1 ULP off, so allow
  // a few ULP at |erf| <= 1. The strict gate is test_erf_ulp vs mpmath.
  constexpr double kTol = 5e-16;

  // Dense sweep across the interesting range, deliberately not a multiple of
  // any lane count so the masked tail path is exercised.
  constexpr size_t kN = 12007;
  std::vector<double> in(kN), out(kN);
  for (size_t i = 0; i < kN; ++i) {
    in[i] = -6.0 + 12.0 * static_cast<double>(i) / (kN - 1);
  }
  corvus::erf(in, out);

  double max_err = 0.0;
  for (size_t i = 0; i < kN; ++i) {
    max_err = std::max(max_err, std::abs(out[i] - std::erf(in[i])));
  }
  std::printf("max abs error vs std::erf over [-6,6]: %.3e\n", max_err);
  Check(max_err < kTol, "sweep accuracy within provisional tolerance");

  // Specials.
  const double inf = std::numeric_limits<double>::infinity();
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double dbl_max = std::numeric_limits<double>::max();
  // #14 N12: erf(+/-DBL_MAX), folded into the same batch call (length 7,
  // still not a lane multiple) so the masked-lane path carries more than one
  // special at once. corvus.h documents "Max 1 ULP over the full domain" and
  // erf(+/-inf) = +/-1; true erf(x) is within 2^-1074 of 1 (resp. -1) once
  // |1 - erf(x)| underflows the smallest denormal, which happens long before
  // x = DBL_MAX (erf(6) is already 1 - 2.15e-17), so the correctly-rounded
  // double result at DBL_MAX is exactly +/-1, not merely close to it.
  std::vector<double> sp_in = {0.0, -0.0, inf, -inf, nan, dbl_max, -dbl_max};
  std::vector<double> sp_out(sp_in.size());
  corvus::erf(sp_in, sp_out);
  // #14 N6: exact-bits + signbit, not just `== 0.0` -- IEEE == treats +0 and
  // -0 as equal, so a plain `== 0.0` cannot tell +0 from -0 and would pass
  // even if the sign of erf(-0) regressed to +0.
  Check(sp_out[0] == 0.0 && !std::signbit(sp_out[0]), "erf(+0) == +0");
  Check(sp_out[1] == 0.0 && std::signbit(sp_out[1]), "erf(-0) == -0");
  Check(sp_out[2] == 1.0, "erf(inf) == 1");
  Check(sp_out[3] == -1.0, "erf(-inf) == -1");
  Check(std::isnan(sp_out[4]), "erf(nan) is nan");
  Check(sp_out[5] == 1.0, "erf(DBL_MAX) == 1");
  Check(sp_out[6] == -1.0, "erf(-DBL_MAX) == -1");

  // #14 N7: NaN propagates WITH its payload (corvus.h: "NaN propagates
  // (payload preserved)"), not just as some quiet NaN.
  {
    const double payload = std::nan("42");
    std::vector<double> pin = {payload};
    std::vector<double> pout(1);
    corvus::erf(pin, pout);
    Check(SameBits(pout[0], payload), "erf(nan) preserves its NaN payload");
  }

  // Exact aliasing (in-place).
  std::vector<double> buf(in);
  corvus::erf(buf, buf);
  for (size_t i = 0; i < kN; i += 1000) {
    Check(buf[i] == out[i], "in-place matches out-of-place");
  }

  if (failures == 0) {
    std::printf("all checks passed\n");
    return 0;
  }
  std::fprintf(stderr, "%d check(s) failed\n", failures);
  return 1;
}
