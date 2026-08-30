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
    std::vector<double> in = {0.0,    -0.0,   inf,    -inf,  pnan,
                              1.0,    -1.0,   709.0,  -745.0,
                              -708.5,     // subnormal result
                              710.0,      // overflow -> +inf
                              -746.0,     // underflow -> +0
                              1e-300, -1e-300, 0.5,   -0.5,
                              2.5e-310,   // subnormal input
                              100.0,  -100.0};
    std::vector<double> out(in.size());
    corvus::exp(in, out);

    Check(out[0] == 1.0, "exp(+0) == 1 exactly");
    Check(out[1] == 1.0, "exp(-0) == 1 exactly");
    Check(out[2] == inf, "exp(+inf) is +inf");
    // Signed zero cannot be told apart by IEEE ==; assert the exact bits
    // (#14 N6).
    Check(SameBits(out[3], 0.0), "exp(-inf) is +0 with signbit clear");
    // NaN propagates WITH its payload, not just as some quiet NaN.
    Check(SameBits(out[4], pnan), "exp(nan) preserves its NaN payload");
    // Overflow and underflow lanes have exact expected results.
    Check(out[10] == inf, "exp(710) overflows to +inf");
    Check(SameBits(out[11], 0.0), "exp(-746) underflows to +0 with signbit clear");
    for (size_t i = 5; i < in.size(); ++i) {
      if (i == 10 || i == 11) continue;
      Check(std::isfinite(out[i]) && out[i] > 0.0,
            "exp of ordinary finite input is finite and > 0");
    }
  }

  // --- 2. Payload NaN in every lane position, so no lane of any tier's
  // vector (or its masked tail) can drop or requantize a NaN. Length 17 is
  // again a non-multiple of 2/4/8.
  for (size_t p = 0; p < 16; ++p) {
    std::vector<double> in(17, 0.5);
    in[p] = pnan;
    std::vector<double> out(in.size());
    corvus::exp(in, out);
    Check(SameBits(out[p], pnan), "exp: NaN payload preserved at position p");
    for (size_t i = 0; i < in.size(); ++i) {
      if (i == p) continue;
      Check(std::isfinite(out[i]) && out[i] > 0.0,
            "exp: non-NaN lanes finite and > 0 around NaN");
    }
  }

  // --- 3. Exact aliasing (in.data() == out.data()) must be bit-identical to
  // the out-of-place call. Length 21: non-multiple of 2/4/8 again.
  {
    std::vector<double> in = {0.5,    -0.5,   1.0,     -2.5,   3.75,   -745.1,
                              709.5,  100.0,  -100.0,  0.125,  6.0,    -6.0,
                              -300.25, -0.01, 2.5,     -3.125, 7.0,    -8.5,
                              0.25,   1e-3,   -650.0};
    std::vector<double> expect(in.size());
    corvus::exp(in, expect);

    std::vector<double> buf(in);
    corvus::exp(buf, buf);
    for (size_t i = 0; i < in.size(); ++i) {
      Check(SameBits(buf[i], expect[i]),
            "exp in-place matches out-of-place bit-for-bit");
    }
  }

  if (failures == 0) {
    std::printf("all checks passed\n");
    return 0;
  }
  std::fprintf(stderr, "%d check(s) failed\n", failures);
  return 1;
}
