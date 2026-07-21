#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>

#include "corvus/corvus.h"

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
  std::printf("corvus active SIMD target: %s\n", corvus::active_target());

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
  std::vector<double> sp_in = {0.0, -0.0, inf, -inf, nan};
  std::vector<double> sp_out(sp_in.size());
  corvus::erf(sp_in, sp_out);
  Check(sp_out[0] == 0.0, "erf(0) == 0");
  Check(sp_out[1] == 0.0 && std::signbit(sp_out[1]), "erf(-0) == -0");
  Check(sp_out[2] == 1.0, "erf(inf) == 1");
  Check(sp_out[3] == -1.0, "erf(-inf) == -1");
  Check(std::isnan(sp_out[4]), "erf(nan) is nan");

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
