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

  // Smoke tolerance vs std::erfc (itself only ~1-2 ULP): relative, loose.
  // The strict gate is test_erfc_ulp vs mpmath.
  constexpr double kRelTol = 5e-15;

  // Dense sweep across core + tail; length deliberately not a multiple of
  // any lane count so the masked tail path is exercised.
  constexpr size_t kN = 24007;
  std::vector<double> in(kN), out(kN);
  for (size_t i = 0; i < kN; ++i) {
    in[i] = -6.5 + 33.0 * static_cast<double>(i) / (kN - 1);  // [-6.5, 26.5]
  }
  corvus::erfc(in, out);

  double max_rel = 0.0;
  for (size_t i = 0; i < kN; ++i) {
    const double want = std::erfc(in[i]);
    if (want > 0.0) {
      max_rel = std::max(max_rel, std::abs(out[i] - want) / want);
    }
  }
  std::printf("max rel error vs std::erfc over [-6.5,26.5]: %.3e\n", max_rel);
  Check(max_rel < kRelTol, "sweep relative accuracy vs libm smoke bound");

  // Specials.
  const double inf = std::numeric_limits<double>::infinity();
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double dbl_max = std::numeric_limits<double>::max();
  // #14 N12: erfc(+/-DBL_MAX), folded into the same batch call (length 9,
  // still not a lane multiple). corvus.h documents "erfc(-inf) = 2,
  // erfc(+inf) = 0, results underflow gradually past x ~ 26.5"; DBL_MAX is
  // far past that underflow point (erfc(30) is already the hard-zero row
  // below), so the correctly-rounded result there is exactly +0, and by
  // erfc(-x) = 2 - erfc(x), exactly 2 at -DBL_MAX.
  std::vector<double> sp_in = {0.0, -0.0, inf, -inf, nan, 30.0, -30.0,
                               dbl_max, -dbl_max};
  std::vector<double> sp_out(sp_in.size());
  corvus::erfc(sp_in, sp_out);
  Check(sp_out[0] == 1.0, "erfc(0) == 1");
  Check(sp_out[1] == 1.0, "erfc(-0) == 1");
  Check(sp_out[2] == 0.0, "erfc(inf) == 0");
  Check(sp_out[3] == 2.0, "erfc(-inf) == 2");
  Check(std::isnan(sp_out[4]), "erfc(nan) is nan");
  Check(sp_out[5] == 0.0, "erfc(30) == 0 (underflow)");
  Check(sp_out[6] == 2.0, "erfc(-30) == 2");
  Check(sp_out[7] == 0.0, "erfc(DBL_MAX) == 0 (underflow)");
  Check(sp_out[8] == 2.0, "erfc(-DBL_MAX) == 2");

  // #14 N7: NaN propagates WITH its payload (corvus.h: "NaN propagates
  // (payload preserved)"), not just as some quiet NaN.
  {
    const double payload = std::nan("42");
    std::vector<double> pin = {payload};
    std::vector<double> pout(1);
    corvus::erfc(pin, pout);
    Check(SameBits(pout[0], payload), "erfc(nan) preserves its NaN payload");
  }

  // Region boundaries: each side of the core/tail split (6.0) and the
  // polynomial interval splits (10, 17) stays close to libm. (Do not test
  // step-to-step jumps: erfc's own relative slope is ~2a per ULP of x, so
  // the true function moves ~1e-13 across two ULPs at x = 17.)
  for (double b : {6.0, 10.0, 17.0}) {
    std::vector<double> bx = {std::nextafter(b, 0.0), b, std::nextafter(b, 99.0)};
    std::vector<double> by(bx.size());
    corvus::erfc(bx, by);
    for (size_t k = 0; k < bx.size(); ++k) {
      const double want = std::erfc(bx[k]);
      Check(std::abs(by[k] - want) / want < 2e-15, "boundary point near libm");
    }
  }

  // Exact aliasing (in-place).
  std::vector<double> buf(in);
  corvus::erfc(buf, buf);
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
