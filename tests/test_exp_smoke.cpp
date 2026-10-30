#include <array>
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

  // --- 1b. expm1 (#45) specials, same length-19 shape (#14).
  {
    std::vector<double> in = {0.0,    -0.0,   inf,    -inf,  pnan,
                              1e-300, -1e-300,
                              5e-324,     // min subnormal: expm1(x) = x
                              -5e-324,
                              710.0,      // overflow -> +inf
                              -38.0,      // exactly -1
                              -1000.0,    // exactly -1
                              0x1p-6, -0x1p-6, 0.5,   -0.5,
                              709.0,  -37.0,  2.5e-310};
    std::vector<double> out(in.size());
    corvus::expm1(in, out);

    // Signed zero cannot be told apart by IEEE ==; assert the exact bits
    // (#14 N6).
    Check(SameBits(out[0], 0.0), "expm1(+0) is +0");
    Check(SameBits(out[1], -0.0), "expm1(-0) is -0");
    Check(out[2] == inf, "expm1(+inf) is +inf");
    Check(out[3] == -1.0, "expm1(-inf) is -1");
    Check(SameBits(out[4], pnan), "expm1(nan) preserves its NaN payload");
    Check(out[5] == 1e-300 && out[6] == -1e-300, "expm1(+-1e-300) = x");
    Check(SameBits(out[7], 5e-324) && SameBits(out[8], -5e-324),
          "expm1 of a subnormal is the subnormal");
    Check(out[9] == inf, "expm1(710) overflows to +inf");
    Check(out[10] == -1.0 && out[11] == -1.0, "expm1(x <= -38) is -1");
    Check(out[17] > -1.0 && out[17] < -0.999999999999999,
          "expm1(-37) sits just above -1");
    Check(SameBits(out[18], 2.5e-310), "expm1(2.5e-310) = x");
    for (size_t i = 12; i < 17; ++i) {
      Check(std::isfinite(out[i]), "expm1 of ordinary input is finite");
    }
  }

  // --- 2. Payload NaN in every lane position, so no lane of any tier's
  // vector (or its masked tail) can drop or requantize a NaN. Length 17 is
  // again a non-multiple of 2/4/8.
  for (size_t p = 0; p < 16; ++p) {
    std::vector<double> in(17, 0.5);
    in[p] = pnan;
    std::vector<double> out(in.size()), m1out(in.size());
    corvus::exp(in, out);
    corvus::expm1(in, m1out);
    Check(SameBits(out[p], pnan), "exp: NaN payload preserved at position p");
    Check(SameBits(m1out[p], pnan),
          "expm1: NaN payload preserved at position p");
    for (size_t i = 0; i < in.size(); ++i) {
      if (i == p) continue;
      Check(std::isfinite(out[i]) && out[i] > 0.0,
            "exp: non-NaN lanes finite and > 0 around NaN");
      Check(std::isfinite(m1out[i]) && m1out[i] > 0.0,
            "expm1: non-NaN lanes finite and > 0 around NaN");
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

    std::vector<double> expect_m1(in.size());
    corvus::expm1(in, expect_m1);
    buf = in;
    corvus::expm1(buf, buf);
    for (size_t i = 0; i < in.size(); ++i) {
      Check(SameBits(buf[i], expect_m1[i]),
            "expm1 in-place matches out-of-place bit-for-bit");
    }
  }

  // --- 4. expm1 lane mix (#45): the kernel skips the series or the
  // exponential branch only when no lane of the vector needs it, so a
  // point's bits must not depend on its neighbours. Probes on both sides of
  // the 2^-6 cut and the -38 shortcut, specials included, against every
  // lane position.
  {
    const std::vector<std::array<double, 1>> probes = {
        {1e-300}, {-0x1.8p-40}, {0x1p-6}, {-0x1p-6}, {0x1.0000000000001p-6},
        {-0x1.0000000000001p-6}, {0.3}, {-0.7}, {5.5}, {-20.0}, {-37.9},
        {-38.0}, {-500.0}, {700.0}, {709.9}, {inf}, {-inf}, {0.0}, {-0.0},
        {pnan}, {2.5e-310}};
    failures += corvus_test::LaneMix(
        "expm1", probes,
        [](const auto& c, auto& o) { corvus::expm1(c[0], o); });
  }

  if (failures == 0) {
    std::printf("all checks passed\n");
    return 0;
  }
  std::fprintf(stderr, "%d check(s) failed\n", failures);
  return 1;
}
