// Reports the runtime-dispatched SIMD target and, when CORVUS_EXPECT_TARGET is
// set in the environment, fails if dispatch did not land on it.
//
// Why this exists: tier sweeps cap targets with CORVUS_DISABLED_TARGETS and
// then *assume* the cap bit. A cap that silently fails to take effect leaves
// the suite green while measuring a tier nobody asked for — not hypothetical,
// Highway's HWY_BROKEN_* blocklists do exactly that (every AVX3* target
// vanishes under MSVC with no diagnostic), and an uncapped run on unknown CI
// hardware has the same failure mode. Accuracy bounds and benchmark numbers
// are only meaningful next to the tier that produced them, so the tier is
// asserted rather than assumed.
//
// Unset or empty CORVUS_EXPECT_TARGET means "report only", so local runs and
// native-dispatch runs need no configuration.
#ifndef CORVUS_TESTS_EXPECT_TARGET_H_
#define CORVUS_TESTS_EXPECT_TARGET_H_

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "corvus/corvus.h"

namespace corvus_test {

// Always prints the active target. Returns false only on an explicit mismatch.
// Safe to call before any kernel call: active_target() goes through the same
// dynamic-dispatch table, so it resolves the target the kernels will use.
inline bool ReportAndCheckTarget() {
  const char* active = corvus::active_target();
  std::printf("corvus active SIMD target: %s\n", active);

  // MSVC deprecates std::getenv in favour of _dupenv_s (C4996), which /WX
  // promotes to an error. The concern does not apply here: one read at
  // startup, no writes to the environment, and the pointer is not retained
  // past the comparisons below. Suppressed locally rather than with
  // _CRT_SECURE_NO_WARNINGS so the deprecation keeps biting everywhere else.
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 4996)
#endif
  const char* want = std::getenv("CORVUS_EXPECT_TARGET");
#if defined(_MSC_VER)
#pragma warning(pop)
#endif
  if (want == nullptr || want[0] == '\0') return true;
  if (std::strcmp(active, want) == 0) return true;

  std::fprintf(stderr,
               "FAIL: expected SIMD target '%s' but dispatch selected '%s'.\n"
               "A CORVUS_DISABLED_TARGETS cap did not take effect, or this "
               "compiler's Highway blocklist removed the intended tier.\n",
               want, active);
  return false;
}

}  // namespace corvus_test

#endif  // CORVUS_TESTS_EXPECT_TARGET_H_
