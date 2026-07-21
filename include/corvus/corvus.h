#ifndef CORVUS_CORVUS_H_
#define CORVUS_CORVUS_H_

/// \file corvus.h
/// \brief Public API: SIMD-vectorized statistical special functions.
///
/// This is the only installed header. The SIMD backend is an implementation
/// detail: no Highway types appear here, and consumers need only link
/// libcorvus. All functions dispatch at runtime to the best SIMD tier the
/// CPU supports (SSE2..AVX-512, NEON).
///
/// Common contract for all batch functions:
///  - `in` and `out` must have the same length.
///  - Exact aliasing is allowed (`in.data() == out.data()`); partial overlap
///    is undefined behavior.
///  - No allocation, no exceptions, thread-safe (stateless; the dispatch
///    pointer is resolved on first call).
///  - Accuracy bounds are measured against a correctly-rounded mpmath
///    oracle and enforced per-tier by the test suite; see docs/ACCURACY.md.

#include <cstddef>
#include <span>

namespace corvus {

inline constexpr int kVersionMajor = 0;
inline constexpr int kVersionMinor = 1;
inline constexpr int kVersionPatch = 0;

/// \brief Name of the SIMD target selected by runtime dispatch.
/// \return A static string such as "AVX2", "SSE4", or "NEON". Note Highway
///   names AVX-512 tiers "AVX3" (see README naming note).
const char* active_target();

/// \brief out[i] = erf(in[i]).
///
/// Max 1 ULP over the full domain on validated tiers. Specials: erf(+/-0)
/// = +/-0, erf(+/-inf) = +/-1, NaN propagates (payload preserved).
void erf(std::span<const double> in, std::span<double> out);

/// \brief out[i] = erfc(in[i]).
///
/// Accuracy on validated tiers: max 1 ULP for |x| <= 6 and for subnormal
/// results; max 5 ULP for normal-result x > 6 (bounded by the backend
/// vector exp; see docs/ACCURACY.md). Specials: erfc(-inf) = 2,
/// erfc(+inf) = 0, results underflow gradually past x ~ 26.5, NaN
/// propagates (payload preserved).
void erfc(std::span<const double> in, std::span<double> out);

}  // namespace corvus

#endif  // CORVUS_CORVUS_H_
