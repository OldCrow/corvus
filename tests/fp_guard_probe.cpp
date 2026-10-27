// #78: one TU, compiled two ways by tests/CMakeLists.txt, both with a
// fast-math flag first on the line, as a parent project's CMAKE_CXX_FLAGS
// would put it.
// - fp_guard_override adds CORVUS_FP_FLAGS after it. The guard must stay
//   silent, and at run time NaN and signed zero must behave as IEEE says.
// - fp_guard_tripwire does not. src/fp_guard.h must stop the build.
#include <cmath>
#include <cstdio>

#include "src/fp_guard.h"

int main() {
  // volatile: the values must come from run-time arithmetic, not folding.
  volatile double zero = 0.0;
  const double nan = zero / zero;
  const double neg_zero = -zero;
  int failures = 0;
  if (nan == nan || !std::isnan(nan)) {
    std::puts("FAIL: NaN compares equal or is not NaN (finite-math folding)");
    ++failures;
  }
  if (!std::signbit(neg_zero) || std::signbit(zero + neg_zero)) {
    std::puts("FAIL: signed zero not preserved (-0, or 0 + -0 = +0)");
    ++failures;
  }
  if (failures == 0) std::puts("fp guard override: IEEE NaN and signed zero");
  return failures;
}
