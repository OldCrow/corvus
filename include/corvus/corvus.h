#ifndef CORVUS_CORVUS_H_
#define CORVUS_CORVUS_H_

#include <cstddef>
#include <span>

namespace corvus {

inline constexpr int kVersionMajor = 0;
inline constexpr int kVersionMinor = 1;
inline constexpr int kVersionPatch = 0;

// Name of the SIMD target selected by runtime dispatch (e.g. "AVX2", "NEON").
const char* ActiveTarget();

// out[i] = erf(in[i]). Spans must be the same length; may alias exactly
// (in.data() == out.data()) but must not partially overlap.
void Erf(std::span<const double> in, std::span<double> out);

// out[i] = erfc(in[i]). Same aliasing rules as Erf.
void Erfc(std::span<const double> in, std::span<double> out);

}  // namespace corvus

#endif  // CORVUS_CORVUS_H_
