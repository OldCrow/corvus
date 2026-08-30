#include <cmath>
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

int failures = 0;

void Check(bool ok, const char* what) {
  if (!ok) {
    std::fprintf(stderr, "FAIL: %s\n", what);
    ++failures;
  }
}

// A quiet NaN with a specific payload, so payload propagation is observable
// bit-for-bit rather than "some NaN came out".
double PayloadNan() {
  const uint64_t bits = 0x7FF8DEADBEEF1234ull;
  double v;
  std::memcpy(&v, &bits, sizeof(v));
  return v;
}

}  // namespace

int main() {
  if (!corvus_test::ReportAndCheckTarget()) return 2;

  const double inf = std::numeric_limits<double>::infinity();
  const double pnan = PayloadNan();

  // --- 1. Specials in one batch call. Length 19 is deliberately a
  // non-multiple of every lane count (2/4/8) so the masked-tail path runs
  // in the same call that carries the specials (#14).
  {
    std::vector<double> in = {0.0,   -0.0, inf,  -inf, pnan,
                              0.5,   1.0,  -2.5, 1e9,  1e300,
                              2.5e-310,  // subnormal
                              -0.75, 3.0,  -1e9, 42.0, -300.25,
                              1e-8,  2.0,  -6.5};
    std::vector<double> cout_(in.size()), sout(in.size());
    corvus::cos(in, cout_);
    corvus::sin(in, sout);

    Check(cout_[0] == 1.0, "cos(+0) == 1 exactly");
    Check(cout_[1] == 1.0, "cos(-0) == 1 exactly");
    // Signed zero cannot be told apart by IEEE ==; assert the exact bits
    // (#14 N6).
    Check(SameBits(sout[0], 0.0), "sin(+0) is +0 with signbit clear");
    Check(SameBits(sout[1], -0.0), "sin(-0) is -0 with signbit set");
    Check(std::isnan(cout_[2]), "cos(+inf) is nan");
    Check(std::isnan(cout_[3]), "cos(-inf) is nan");
    Check(std::isnan(sout[2]), "sin(+inf) is nan");
    Check(std::isnan(sout[3]), "sin(-inf) is nan");
    // NaN propagates WITH its payload, not just as some quiet NaN.
    Check(SameBits(cout_[4], pnan), "cos(nan) preserves its NaN payload");
    Check(SameBits(sout[4], pnan), "sin(nan) preserves its NaN payload");
    for (size_t i = 5; i < in.size(); ++i) {
      Check(std::isfinite(cout_[i]) && cout_[i] >= -1.0 && cout_[i] <= 1.0,
            "cos of finite input is finite and in [-1,1]");
      Check(std::isfinite(sout[i]) && sout[i] >= -1.0 && sout[i] <= 1.0,
            "sin of finite input is finite and in [-1,1]");
    }
  }

  // --- 2. Payload NaN in every lane position, so no lane of any tier's
  // vector (or its masked tail) can drop or requantize a NaN. Length 17 is
  // again a non-multiple of 2/4/8.
  for (size_t p = 0; p < 16; ++p) {
    std::vector<double> in(17, 0.5);
    in[p] = pnan;
    std::vector<double> cout_(in.size()), sout(in.size());
    corvus::cos(in, cout_);
    corvus::sin(in, sout);
    Check(SameBits(cout_[p], pnan), "cos: NaN payload preserved at position p");
    Check(SameBits(sout[p], pnan), "sin: NaN payload preserved at position p");
    for (size_t i = 0; i < in.size(); ++i) {
      if (i == p) continue;
      Check(std::isfinite(cout_[i]), "cos: non-NaN lanes finite around NaN");
      Check(std::isfinite(sout[i]), "sin: non-NaN lanes finite around NaN");
    }
  }

  // --- 3. Exact aliasing (in.data() == out.data()) must be bit-identical to
  // the out-of-place call. Length 21: non-multiple of 2/4/8 again.
  {
    std::vector<double> in = {0.5,    -0.5,  1.0,   -2.5,   3.75,  1e300,
                              -1e300, 1e9,   -42.0, 0.125,  6.0,   -6.0,
                              100.5,  -0.01, 2.5,   -3.125, 7.0,   -8.5,
                              0.25,   1e-3,  -1e6};
    std::vector<double> expect_cos(in.size()), expect_sin(in.size());
    corvus::cos(in, expect_cos);
    corvus::sin(in, expect_sin);

    std::vector<double> buf(in);
    corvus::cos(buf, buf);
    for (size_t i = 0; i < in.size(); ++i) {
      Check(SameBits(buf[i], expect_cos[i]),
            "cos in-place matches out-of-place bit-for-bit");
    }
    buf = in;
    corvus::sin(buf, buf);
    for (size_t i = 0; i < in.size(); ++i) {
      Check(SameBits(buf[i], expect_sin[i]),
            "sin in-place matches out-of-place bit-for-bit");
    }
  }

  // --- 4. Parity consistency between the two entry points: sin is odd and
  // cos is even, and both must hold BIT-exactly (same reduction path for x
  // and -x), including at huge arguments and subnormals.
  {
    std::vector<double> xs = {0.5,     1.0,     2.5,      3.14159, 1e9,
                              1e300,   8388608.0 /* 2^23 */, 1e-300,
                              5e-324 /* min subnormal */,   2.5e-310,
                              0.001,   6.25,    77.0,     1e6,     123456.789,
                              0.75,    1.5,     4.0,      9.0,     16.0,
                              1e-8,    1e-16,   3e5,      2.2e10,  7.7e100,
                              1e200,   0.1,     0.2,      0.3,     0.4,
                              10.0,    20.0,    30.0,     100.0,   1000.0,
                              1e4,     3.3e7,   4.4e8,    5.5e11,  6.6e13,
                              1.25e-5, 2.5e-7,  9.9e-100, 1.1e50,  2.2e150,
                              0.9,     1.1,     2.9,      3.1,     6.28};
    std::vector<double> neg(xs.size());
    for (size_t i = 0; i < xs.size(); ++i) neg[i] = -xs[i];

    std::vector<double> sp(xs.size()), sn(xs.size());
    std::vector<double> cp(xs.size()), cn(xs.size());
    corvus::sin(xs, sp);
    corvus::sin(neg, sn);
    corvus::cos(xs, cp);
    corvus::cos(neg, cn);
    for (size_t i = 0; i < xs.size(); ++i) {
      Check(SameBits(sn[i], -sp[i]), "sin(-x) is the exact negation of sin(x)");
      Check(SameBits(cn[i], cp[i]), "cos(-x) is bit-identical to cos(x)");
    }
  }

  if (failures == 0) {
    std::printf("all checks passed\n");
    return 0;
  }
  std::fprintf(stderr, "%d check(s) failed\n", failures);
  return 1;
}
