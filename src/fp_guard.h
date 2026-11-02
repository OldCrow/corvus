// Refuses a value-unsafe floating-point mode at compile time (#78).
//
// The double-double layer (src/dd-inl.h) rests on TwoSum, Fast2Sum and
// Dekker, which hold only when every IEEE operation rounds as written:
// reassociation or a reciprocal in place of a division breaks them. The NaN
// scrubs ahead of every value-derived gather, and the NaN/inf specials, need
// IsNaN and IsInf to mean what they say: finite-math folds them to false.
//
// CORVUS_FP_FLAGS (CMakeLists.txt) turns these modes off again after any
// flags a parent project sets -- -fno-fast-math on GCC, Clang and clang-cl,
// /fp:precise on MSVC -- so a fast-math consumer still builds a correct
// corvus. This guard fires only when something re-enables them after those
// flags, or the compiler is not one CORVUS_FP_FLAGS knows.
//
// What it can see, measured 2026-10-26 on GCC 16.2 (mingw), clang-cl 22.1
// and cl 19.51:
// - GCC: __FAST_MATH__, __FINITE_MATH_ONLY__ = 1, __ASSOCIATIVE_MATH__,
//   __RECIPROCAL_MATH__ and __NO_SIGNED_ZEROS__, for -ffast-math and for each
//   flag alone.
// - Clang and clang-cl: __FAST_MATH__ and __FINITE_MATH_ONLY__ = 1 (plus
//   _M_FP_FAST on clang-cl) for -ffast-math or /fp:fast and
//   -ffinite-math-only. Clang defines no macro for -fassociative-math or
//   -freciprocal-math given alone, so only the -fno-fast-math override
//   covers those.
// - MSVC: _M_FP_FAST for /fp:fast.
// __FINITE_MATH_ONLY__ is defined as 0 in the default mode, so it is tested
// by value.
#ifndef CORVUS_FP_GUARD_H_
#define CORVUS_FP_GUARD_H_

#if defined(__FAST_MATH__) || \
    (defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__) || \
    defined(__ASSOCIATIVE_MATH__) || defined(__RECIPROCAL_MATH__) || \
    defined(__NO_SIGNED_ZEROS__) || defined(_M_FP_FAST)
#error "corvus: value-unsafe floating-point mode (fast-math, finite-math, associative/reciprocal math, no-signed-zeros or /fp:fast) -- see src/fp_guard.h"
#endif

#endif  // CORVUS_FP_GUARD_H_
