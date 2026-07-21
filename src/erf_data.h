#ifndef CORVUS_ERF_DATA_H_
#define CORVUS_ERF_DATA_H_

#include <cstdint>

namespace corvus::detail {

// Flat table, 3 doubles per grid point {E_hi, S, E_lo} on r_j = j/256,
// j = 0..kErfTableLastIndex. Defined once in erf_data.cpp so every SIMD
// target compiled by foreach_target shares a single copy.
extern const double kErfTable[];
extern const double kErfAMax;
extern const int64_t kErfTableLastIndex;

}  // namespace corvus::detail

#endif  // CORVUS_ERF_DATA_H_
