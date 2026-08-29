// Measures max ULP deviation of corvus::erfinv / corvus::erfcinv against the
// mpmath-generated correctly-rounded (erfinv) / root-found (erfcinv)
// reference, with a per-region breakdown matching the kernel's own routing
// (see src/erfinv-inl.h): C (the central polynomial), T-mid (seed + Halley
// in residual space, |result| < 6) and T-far (seed + Halley in log space,
// |result| >= 6 -- reachable only through erfcinv; erfinv's own range never
// exceeds ~5.86).
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "corvus/corvus.h"
#include "expect_target.h"
#include "ulp_utils.h"

namespace {

using corvus_test::UlpDiff;

// Gates, set from measured values with no margin (house rule): regressions
// should trip them.
constexpr uint64_t kMaxUlpErfinvC = 1;
constexpr uint64_t kMaxUlpErfinvT = 1;
constexpr uint64_t kMaxUlpErfcinvC = 1;
constexpr uint64_t kMaxUlpErfcinvTMid = 1;
constexpr uint64_t kMaxUlpErfcinvTFar = 1;

struct Region {
  const char* name;
  uint64_t bound;
  uint64_t max_ulp = 0;
  size_t n = 0;
  size_t miss = 0;
  double worst_x = 0.0;
};

bool LoadReference(const char* path, std::vector<double>* in,
                   std::vector<double>* want) {
  const auto rows = corvus_test::LoadRef(path, 2, 5000);
  for (const auto& row : rows) {
    in->push_back(corvus_test::ParseDouble(row.tok[0], path, row.line));
    want->push_back(corvus_test::ParseDouble(row.tok[1], path, row.line));
  }
  return true;
}

int ReportRegions(const char* label, const Region* regions, int n_regions) {
  int rc = 0;
  for (int i = 0; i < n_regions; ++i) {
    const Region& r = regions[i];
    std::printf(
        "%-8s %-13s n=%6zu  max ULP=%3llu (gate %llu)  not-CR: %zu (%.2f%%)  "
        "worst x=%.17g\n",
        label, r.name, r.n, static_cast<unsigned long long>(r.max_ulp),
        static_cast<unsigned long long>(r.bound), r.miss,
        r.n ? 100.0 * static_cast<double>(r.miss) / static_cast<double>(r.n)
            : 0.0,
        r.worst_x);
    if (r.n > 0 && r.max_ulp > r.bound) {
      std::fprintf(stderr, "FAIL: %s %s exceeds gate\n", label, r.name);
      rc = 1;
    }
  }
  return rc;
}

}  // namespace

int main(int argc, char** argv) {
  if (!corvus_test::ReportAndCheckTarget()) return 2;

  const char* erfinv_path =
      argc > 1 ? argv[1] : "tests/data/erfinv_reference.txt";
  const char* erfcinv_path =
      argc > 2 ? argv[2] : "tests/data/erfcinv_reference.txt";

  int rc = 0;

  // --- erfinv: C (|y| <= 1/2) vs T (1/2 < |y| < 1); T never reaches "far" ---
  {
    std::vector<double> in, want;
    if (!LoadReference(erfinv_path, &in, &want)) return 2;
    std::vector<double> got(in.size());
    corvus::erfinv(in, got);

    Region regions[] = {{"C", kMaxUlpErfinvC}, {"T", kMaxUlpErfinvT}};
    for (size_t i = 0; i < in.size(); ++i) {
      Region& r = std::abs(in[i]) <= 0.5 ? regions[0] : regions[1];
      const uint64_t u = UlpDiff(got[i], want[i]);
      ++r.n;
      if (u > 0) ++r.miss;
      if (u > r.max_ulp) {
        r.max_ulp = u;
        r.worst_x = in[i];
      }
    }
    rc |= ReportRegions("erfinv", regions, 2);
  }

  // --- erfcinv: C ([1/2,3/2]) vs T-mid (|x|<6) vs T-far (|x|>=6) ---
  {
    std::vector<double> in, want;
    if (!LoadReference(erfcinv_path, &in, &want)) return 2;
    std::vector<double> got(in.size());
    corvus::erfcinv(in, got);

    Region regions[] = {{"C", kMaxUlpErfcinvC},
                        {"T-mid", kMaxUlpErfcinvTMid},
                        {"T-far", kMaxUlpErfcinvTFar}};
    for (size_t i = 0; i < in.size(); ++i) {
      const double z = in[i];
      Region& r = (z >= 0.5 && z <= 1.5)   ? regions[0]
                  : (std::abs(want[i]) < 6.0) ? regions[1]
                                              : regions[2];
      const uint64_t u = UlpDiff(got[i], want[i]);
      ++r.n;
      if (u > 0) ++r.miss;
      if (u > r.max_ulp) {
        r.max_ulp = u;
        r.worst_x = z;
      }
    }
    rc |= ReportRegions("erfcinv", regions, 3);
  }

  if (rc == 0) std::printf("PASS: all regions within gates\n");
  return rc;
}
