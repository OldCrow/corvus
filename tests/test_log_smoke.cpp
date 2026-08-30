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

  // --- 1a. log specials in one batch call. Length 19 is deliberately a
  // non-multiple of every lane count (2/4/8) so the masked-tail path runs
  // in the same call that carries the specials (#14).
  {
    std::vector<double> in = {0.0,   -0.0,  -1.0,  -1e300, pnan,
                              1.0,   inf,
                              5e-324,     // min subnormal
                              2.5e-310,   // subnormal
                              1e300, 0.5,  2.0,    1e-8,
                              3.0,   1e6,  0.001,  42.0,   7.5,
                              123456.789};
    std::vector<double> out(in.size());
    corvus::log(in, out);

    Check(out[0] == -inf, "log(+0) is -inf");
    Check(out[1] == -inf, "log(-0) is -inf");
    // Negative inputs: any NaN is acceptable, payload unspecified.
    Check(std::isnan(out[2]), "log(-1) is nan");
    Check(std::isnan(out[3]), "log(-1e300) is nan");
    // NaN propagates WITH its payload, not just as some quiet NaN.
    Check(SameBits(out[4], pnan), "log(nan) preserves its NaN payload");
    // Signed zero cannot be told apart by IEEE ==; assert the exact bits
    // (#14 N6).
    Check(SameBits(out[5], 0.0), "log(1) is +0 with signbit clear");
    Check(out[6] == inf, "log(+inf) is +inf");
    for (size_t i = 7; i < in.size(); ++i) {
      Check(std::isfinite(out[i]), "log of ordinary positive input is finite");
    }
  }

  // --- 1b. log1p specials in one batch call, same length-19 shape (#14).
  {
    std::vector<double> in = {-1.0,  -1.5,  -1e300, 0.0,   -0.0,
                              pnan,  inf,
                              -0.9999999999999999,  // finite, near the pole
                              1e300, -0.5,  1e-300, -1e-300,
                              5e-324,     // min subnormal
                              -5e-324,
                              0.5,   2.0,   42.0,   1e-8,  -0.25};
    std::vector<double> out(in.size());
    corvus::log1p(in, out);

    Check(out[0] == -inf, "log1p(-1) is -inf");
    // Inputs below -1: any NaN is acceptable, payload unspecified.
    Check(std::isnan(out[1]), "log1p(-1.5) is nan");
    Check(std::isnan(out[2]), "log1p(-1e300) is nan");
    // Signed zero cannot be told apart by IEEE ==; assert the exact bits
    // (#14 N6).
    Check(SameBits(out[3], 0.0), "log1p(+0) is +0 with signbit clear");
    Check(SameBits(out[4], -0.0), "log1p(-0) is -0 with signbit set");
    // NaN propagates WITH its payload, not just as some quiet NaN.
    Check(SameBits(out[5], pnan), "log1p(nan) preserves its NaN payload");
    Check(out[6] == inf, "log1p(+inf) is +inf");
    for (size_t i = 7; i < in.size(); ++i) {
      Check(std::isfinite(out[i]), "log1p of ordinary in-domain input is finite");
    }
  }

  // --- 2. Payload NaN in every lane position, so no lane of any tier's
  // vector (or its masked tail) can drop or requantize a NaN. Length 17 is
  // again a non-multiple of 2/4/8. 0.5 is in-domain for both functions.
  for (size_t p = 0; p < 16; ++p) {
    std::vector<double> in(17, 0.5);
    in[p] = pnan;
    std::vector<double> lout(in.size()), l1out(in.size());
    corvus::log(in, lout);
    corvus::log1p(in, l1out);
    Check(SameBits(lout[p], pnan), "log: NaN payload preserved at position p");
    Check(SameBits(l1out[p], pnan),
          "log1p: NaN payload preserved at position p");
    for (size_t i = 0; i < in.size(); ++i) {
      if (i == p) continue;
      Check(std::isfinite(lout[i]), "log: non-NaN lanes finite around NaN");
      Check(std::isfinite(l1out[i]), "log1p: non-NaN lanes finite around NaN");
    }
  }

  // --- 3. Exact aliasing (in.data() == out.data()) must be bit-identical to
  // the out-of-place call. Length 21: non-multiple of 2/4/8 again.
  {
    // In-domain for log: strictly positive.
    std::vector<double> lin = {0.5,   1.0,    2.0,      3.75,   1e300, 1e-300,
                               1e9,   42.0,   0.125,    6.0,    100.5, 0.01,
                               2.5,   3.125,  7.0,      8.5,    0.25,  1e-3,
                               1e6,   2.5e-310, 5e-324};
    std::vector<double> expect_log(lin.size());
    corvus::log(lin, expect_log);

    std::vector<double> buf(lin);
    corvus::log(buf, buf);
    for (size_t i = 0; i < lin.size(); ++i) {
      Check(SameBits(buf[i], expect_log[i]),
            "log in-place matches out-of-place bit-for-bit");
    }

    // In-domain for log1p: strictly greater than -1.
    std::vector<double> l1in = {-0.999, -0.5,  -0.25,  -0.01, -1e-300, 0.5,
                                1.0,    2.0,   3.75,   1e300, 1e-300,  1e9,
                                42.0,   0.125, 6.0,    100.5, 0.01,    2.5,
                                10.0,   1e-8,  5e-324};
    std::vector<double> expect_log1p(l1in.size());
    corvus::log1p(l1in, expect_log1p);

    buf = l1in;
    corvus::log1p(buf, buf);
    for (size_t i = 0; i < l1in.size(); ++i) {
      Check(SameBits(buf[i], expect_log1p[i]),
            "log1p in-place matches out-of-place bit-for-bit");
    }
  }

  if (failures == 0) {
    std::printf("all checks passed\n");
    return 0;
  }
  std::fprintf(stderr, "%d check(s) failed\n", failures);
  return 1;
}
