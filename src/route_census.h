// Route census (#82 T7): which kernel routes a lane takes.
//
// Kernels mark each route -- a router region, a sub-map branch, a rare-path
// block behind an AllFalse skip, a saturation or special-value fixup -- with
// CORVUS_ROUTE(d, "family/route", mask) (src/ops-inl.h). A region block that
// runs when ANY lane is in the region opens CORVUS_ROUTE_SCOPE(d, region):
// routes inside it, the cores it calls included, count only while lane 0 is
// in every enclosing scope. In a normal build both macros are empty and cost
// nothing. With -DCORVUS_ROUTE_CENSUS=ON, tests/hunt_probe evaluates one row
// per call in lane 0 and prints the routes it took; tools/hunt/
// route_census.py takes the route list from the source and reports every
// route no reference row reaches and every row whose route disagrees with
// the hunt's router replica.
//
// Not part of the public surface; the census build is a test tool only.
#ifndef CORVUS_SRC_ROUTE_CENSUS_H_
#define CORVUS_SRC_ROUTE_CENSUS_H_

#if defined(CORVUS_ROUTE_CENSUS)
#include <vector>

namespace corvus::census {

// Records `name` (a string literal) when `lane0` is set and lane 0 is inside
// every open scope; consecutive repeats (iteration loops) are kept once.
void Hit(const char* name, bool lane0) noexcept;

// The routes recorded since the last call, in order; clears the record.
std::vector<const char*> Take();

void Enter(bool lane0) noexcept;
void Leave(bool lane0) noexcept;

// A region block: lane 0 in the region or not, for the block's lifetime.
class Scope {
 public:
  explicit Scope(bool lane0) noexcept : lane0_(lane0) { Enter(lane0_); }
  ~Scope() { Leave(lane0_); }
  Scope(const Scope&) = delete;
  Scope& operator=(const Scope&) = delete;

 private:
  bool lane0_;
};

}  // namespace corvus::census
#endif  // CORVUS_ROUTE_CENSUS

#endif  // CORVUS_SRC_ROUTE_CENSUS_H_
