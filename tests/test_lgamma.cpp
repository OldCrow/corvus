#include <algorithm>
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

  // Smoke tolerance vs std::lgamma (itself only a few ULP, and on some libms
  // wrong outright for subnormal arguments). The strict gate is
  // test_lgamma_ulp vs mpmath.
  constexpr double kRelTol = 5e-14;

  // Dense sweep over every positive region; length deliberately not a
  // multiple of any lane count so the masked tail path is exercised.
  constexpr size_t kN = 24007;
  std::vector<double> in(kN), out(kN);
  for (size_t i = 0; i < kN; ++i) {
    in[i] = 0.001 + 40.0 * static_cast<double>(i) / (kN - 1);  // (0, 40]
  }
  corvus::lgamma(in, out);

  // Scaled by max(|want|, 1): relative where lgamma is large, absolute where
  // it is small. Not a softening -- libm is simply not a usable relative
  // reference near x = 1 and x = 2, where its own error is ~1e-16 absolute
  // against a result of ~1e-3, i.e. hundreds of ULP. Relative accuracy at the
  // zeros is what test_lgamma_ulp checks, against mpmath.
  double max_err = 0.0;
  double worst = 0.0;
  for (size_t i = 0; i < kN; ++i) {
    const double want = std::lgamma(in[i]);
    const double err =
        std::abs(out[i] - want) / std::max(std::abs(want), 1.0);
    if (err > max_err) {
      max_err = err;
      worst = in[i];
    }
  }
  std::printf("max scaled error vs std::lgamma over (0,40]: %.3e at x=%.17g\n",
              max_err, worst);
  Check(max_err < kRelTol, "sweep accuracy vs libm smoke bound");

  // The zeros are exact, and positively signed. C99 requires +0, and the
  // t*B(t) form naturally produces the sign of B -- negative at x = 1 -- so
  // this is a real check, not a tautology.
  const double inf = std::numeric_limits<double>::infinity();
  const double nan = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> sp_in = {1.0, 2.0,  0.0, -0.0, -1.0, -2.0,
                               -40.0, inf, -inf, nan,  1e308, 3e305};
  std::vector<double> sp_out(sp_in.size());
  corvus::lgamma(sp_in, sp_out);
  Check(sp_out[0] == 0.0 && !std::signbit(sp_out[0]), "lgamma(1) == +0");
  Check(sp_out[1] == 0.0 && !std::signbit(sp_out[1]), "lgamma(2) == +0");
  Check(sp_out[2] == inf, "lgamma(+0) == +inf");
  Check(sp_out[3] == inf, "lgamma(-0) == +inf");
  Check(sp_out[4] == inf, "lgamma(-1) == +inf (pole)");
  Check(sp_out[5] == inf, "lgamma(-2) == +inf (pole)");
  Check(sp_out[6] == inf, "lgamma(-40) == +inf (pole)");
  Check(sp_out[7] == inf, "lgamma(+inf) == +inf");
  Check(sp_out[8] == inf, "lgamma(-inf) == +inf");
  Check(std::isnan(sp_out[9]), "lgamma(nan) is nan");
  Check(sp_out[10] == inf, "lgamma(1e308) == +inf (overflow)");
  Check(sp_out[11] == inf, "lgamma(3e305) == +inf (overflow)");

  // #14 N7: NaN propagates WITH its payload (corvus.h: "NaN propagates
  // (payload preserved)"), not just as some quiet NaN.
  {
    const double payload = std::nan("42");
    std::vector<double> pin = {payload};
    std::vector<double> pout(1);
    corvus::lgamma(pin, pout);
    Check(SameBits(pout[0], payload), "lgamma(nan) preserves its NaN payload");
  }

  // The last finite result, and the first infinite one. This straddles the
  // point where the Stirling product itself would overflow if it were
  // grouped the textbook way; see src/lgamma-inl.h.
  {
    constexpr double kMaxArg = 0x1.754d9278b51a7p+1014;
    std::vector<double> bx = {kMaxArg, std::nextafter(kMaxArg, inf)};
    std::vector<double> by(bx.size());
    corvus::lgamma(bx, by);
    Check(std::isfinite(by[0]), "lgamma is finite at the overflow threshold");
    Check(by[1] == inf, "lgamma overflows one ulp above it");
  }

  // Subnormal arguments: lgamma(x) = -log x to well within a ULP, and the
  // kernel has to prescale to take that log at all.
  {
    std::vector<double> tx = {5e-324, 1e-320, 1e-310, 2.2250738585072014e-308};
    std::vector<double> ty(tx.size());
    corvus::lgamma(tx, ty);
    for (size_t k = 0; k < tx.size(); ++k) {
      const double want = -std::log(tx[k]);
      Check(std::isfinite(ty[k]) && std::abs(ty[k] - want) / want < 1e-14,
            "subnormal argument gives -log x");
    }
  }

  // Region boundaries: both sides of every split stay close to libm.
  for (double b : {0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.0}) {
    std::vector<double> bx = {std::nextafter(b, 0.0), b, std::nextafter(b, 99.0)};
    std::vector<double> by(bx.size());
    corvus::lgamma(bx, by);
    for (size_t k = 0; k < bx.size(); ++k) {
      const double want = std::lgamma(bx[k]);
      Check(std::abs(by[k] - want) / std::abs(want) < 1e-14,
            "boundary point near libm");
    }
  }

  // Exact aliasing (in-place).
  std::vector<double> buf(in);
  corvus::lgamma(buf, buf);
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
