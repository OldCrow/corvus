// Smoke tests for corvus::erfinv / corvus::erfcinv: round-trip against
// corvus's own erf/erfc (there is no std::erfinv to compare against), every
// documented special value, in-place aliasing, and a non-lane-multiple
// length so the masked tail path is exercised. The strict accuracy claim
// against the mpmath oracle is test_erfinv_ulp, not this file.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>

#include "corvus/corvus.h"
#include "expect_target.h"

namespace {

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

  // --- erfinv round-trip: erfinv(erf(x)) ~= x -------------------------
  // erf itself is max 1 ULP; inverting it amplifies that error by the local
  // derivative d(erfinv)/dy = 1/erf'(x) = (sqrt(pi)/2)*exp(x^2), which is
  // where the tolerance's steep x-dependence comes from -- at x = 5.8 the
  // amplification is already ~7 orders of magnitude past the centre.
  {
    constexpr size_t kN = 10007;  // not a multiple of any lane count
    std::vector<double> xs(kN), y(kN), rt(kN);
    for (size_t i = 0; i < kN; ++i) {
      xs[i] = -5.8 + 11.6 * static_cast<double>(i) / (kN - 1);
    }
    corvus::erf(xs, y);
    corvus::erfinv(y, rt);
    double max_err = 0.0, worst_x = 0.0;
    for (size_t i = 0; i < kN; ++i) {
      const double deriv = 0.8862269254527579 * std::exp(xs[i] * xs[i]);
      const double tol = 50.0 * 0x1p-52 * deriv;  // 50x safety margin
      const double err = std::abs(rt[i] - xs[i]);
      if (err > tol) {
        Check(false, "erfinv(erf(x)) round-trip within tolerance");
      }
      if (err > max_err) {
        max_err = err;
        worst_x = xs[i];
      }
    }
    std::printf("erfinv(erf(x)) max abs round-trip error over [-5.8,5.8]: "
                "%.3e at x=%.17g\n",
                max_err, worst_x);
  }

  // --- erfcinv round-trip: erfcinv(erfc(x)) ~= x, x >= 0 ---------------
  // Unlike erfinv's round trip, this one does NOT blow up with x: erfc's
  // error is RELATIVE (bounded ULP of erfc(x) itself, however tiny), and
  // erfcinv's condition number in that relative error is ~2x^2 (PLAN.md,
  // "Phase C" condition analysis) -- the whole reason erfcinv was scheduled
  // right after erfc's tail was flagged as a possible open item. Capped at
  // 26 rather than erfc's full domain: past that erfc(x) itself is
  // subnormal and loses significant bits to gradual underflow, so a
  // round-trip through it is testing subnormal precision loss rather than
  // erfcinv -- test_erfinv_ulp's T-far region (root-found reference, not a
  // round trip) is what actually covers that band, down to 2^-1074.
  {
    constexpr size_t kN = 8009;
    std::vector<double> xs(kN), y(kN), rt(kN);
    for (size_t i = 0; i < kN; ++i) {
      xs[i] = 26.0 * static_cast<double>(i) / (kN - 1);
    }
    corvus::erfc(xs, y);
    corvus::erfcinv(y, rt);
    double max_err = 0.0, worst_x = 0.0;
    for (size_t i = 0; i < kN; ++i) {
      if (y[i] == 0.0) continue;  // erfc underflowed; erfcinv(0) == +inf
      const double tol =
          50.0 * 0x1p-52 * std::max(1.0, 2.0 * xs[i] * xs[i]);
      const double err = std::abs(rt[i] - xs[i]);
      if (err > tol) {
        Check(false, "erfcinv(erfc(x)) round-trip within tolerance");
      }
      if (err > max_err) {
        max_err = err;
        worst_x = xs[i];
      }
    }
    std::printf("erfcinv(erfc(x)) max abs round-trip error over [0,27.2]: "
                "%.3e at x=%.17g\n",
                max_err, worst_x);
  }

  const double inf = std::numeric_limits<double>::infinity();
  const double nan = std::numeric_limits<double>::quiet_NaN();

  // --- erfinv specials --------------------------------------------------
  {
    std::vector<double> in = {0.0, -0.0, 1.0, -1.0,
                              std::nextafter(1.0, 2.0),
                              std::nextafter(-1.0, -2.0),
                              2.0, -2.0, nan};
    std::vector<double> out(in.size());
    corvus::erfinv(in, out);
    Check(out[0] == 0.0 && !std::signbit(out[0]), "erfinv(+0) == +0");
    Check(out[1] == 0.0 && std::signbit(out[1]), "erfinv(-0) == -0");
    Check(out[2] == inf, "erfinv(1) == +inf");
    Check(out[3] == -inf, "erfinv(-1) == -inf");
    Check(std::isnan(out[4]), "erfinv(1+ulp) is NaN");
    Check(std::isnan(out[5]), "erfinv(-1-ulp) is NaN");
    Check(std::isnan(out[6]), "erfinv(2) is NaN");
    Check(std::isnan(out[7]), "erfinv(-2) is NaN");
    Check(std::isnan(out[8]), "erfinv(nan) is NaN");
  }

  // --- erfcinv specials ---------------------------------------------------
  {
    std::vector<double> in = {0.0, 2.0, 1.0, -0.0, std::nextafter(0.0, -1.0),
                              std::nextafter(2.0, 3.0), -0.5, 2.5, nan};
    std::vector<double> out(in.size());
    corvus::erfcinv(in, out);
    Check(out[0] == inf, "erfcinv(0) == +inf");
    Check(out[1] == -inf, "erfcinv(2) == -inf");
    Check(out[2] == 0.0 && !std::signbit(out[2]), "erfcinv(1) == +0");
    Check(out[3] == inf, "erfcinv(-0) == +inf");
    Check(std::isnan(out[4]), "erfcinv(0-ulp) is NaN");
    Check(std::isnan(out[5]), "erfcinv(2+ulp) is NaN");
    Check(std::isnan(out[6]), "erfcinv(-0.5) is NaN");
    Check(std::isnan(out[7]), "erfcinv(2.5) is NaN");
    Check(std::isnan(out[8]), "erfcinv(nan) is NaN");
  }

  // --- routing-boundary continuity (both sides of every split) -----------
  for (double b : {0.5, 1.5}) {
    std::vector<double> bx = {std::nextafter(b, 0.0), b, std::nextafter(b, 3.0)};
    std::vector<double> by(bx.size());
    corvus::erfcinv(bx, by);
    for (size_t k = 0; k + 1 < bx.size(); ++k) {
      Check(std::abs(by[k] - by[k + 1]) < 1e-9,
            "erfcinv continuous across a routing boundary");
    }
  }
  for (double b : {0.5, -0.5}) {
    std::vector<double> bx = {std::nextafter(b, 0.0), b, std::nextafter(b, 2.0)};
    std::vector<double> by(bx.size());
    corvus::erfinv(bx, by);
    for (size_t k = 0; k + 1 < bx.size(); ++k) {
      Check(std::abs(by[k] - by[k + 1]) < 1e-9,
            "erfinv continuous across a routing boundary");
    }
  }

  // --- in-place aliasing --------------------------------------------------
  {
    constexpr size_t kN = 4001;
    std::vector<double> in(kN), out(kN);
    for (size_t i = 0; i < kN; ++i) {
      in[i] = -0.999 + 1.998 * static_cast<double>(i) / (kN - 1);
    }
    corvus::erfinv(in, out);
    std::vector<double> buf(in);
    corvus::erfinv(buf, buf);
    for (size_t i = 0; i < kN; i += 137) {
      Check(buf[i] == out[i], "erfinv in-place matches out-of-place");
    }
  }
  {
    constexpr size_t kN = 4003;
    std::vector<double> in(kN), out(kN);
    for (size_t i = 0; i < kN; ++i) {
      in[i] = 0.001 + 1.998 * static_cast<double>(i) / (kN - 1);
    }
    corvus::erfcinv(in, out);
    std::vector<double> buf(in);
    corvus::erfcinv(buf, buf);
    for (size_t i = 0; i < kN; i += 137) {
      Check(buf[i] == out[i], "erfcinv in-place matches out-of-place");
    }
  }

  if (failures == 0) {
    std::printf("all checks passed\n");
    return 0;
  }
  std::fprintf(stderr, "%d check(s) failed\n", failures);
  return 1;
}
