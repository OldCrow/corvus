// Measures max ULP deviation of corvus::erfc against the mpmath-generated
// correctly-rounded reference, with a per-region breakdown (the tail is
// gated by the backend Exp accuracy, the core by the table method).
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "corvus/corvus.h"
#include "expect_target.h"
#include "ulp_utils.h"

namespace {

using corvus_test::UlpDiff;

// Gates, set from measured values with no margin:
// regressions should trip them. The tail-normal 2-ULP bound is fit-limited:
// what remains is the tail polynomial G, not the exponential (a dd Horner
// would buy 2 -> 1 at a poor speed trade -- documented, accepted).
constexpr uint64_t kMaxUlpCore = 1;          // |x| <= 6
constexpr uint64_t kMaxUlpTailNormal = 2;    // |x| > 6, normal results
constexpr uint64_t kMaxUlpTailSubnormal = 1; // |x| > 6, subnormal results

struct Region {
  const char* name;
  uint64_t bound;
  uint64_t max_ulp = 0;
  size_t n = 0;
  size_t miss = 0;  // not correctly rounded
  double worst_x = 0.0;
};

}  // namespace

int main(int argc, char** argv) {
  // Before doing any work: a wrong tier makes every number below meaningless.
  if (!corvus_test::ReportAndCheckTarget()) return 2;

  const char* path = argc > 1 ? argv[1] : "tests/data/erfc_reference.txt";
  const auto rows = corvus_test::LoadRef(path, 2, 10000);

  std::vector<double> in, want;
  for (const auto& row : rows) {
    in.push_back(corvus_test::ParseDouble(row.tok[0], path, row.line));
    want.push_back(corvus_test::ParseDouble(row.tok[1], path, row.line));
  }

  std::vector<double> got(in.size());
  corvus::erfc(in, got);

  Region regions[] = {
      {"core |x|<=6", kMaxUlpCore},
      {"tail normal", kMaxUlpTailNormal},
      {"tail subnormal", kMaxUlpTailSubnormal},
  };
  constexpr double kMinNormal = 2.2250738585072014e-308;

  for (size_t i = 0; i < in.size(); ++i) {
    Region& r = std::abs(in[i]) <= 6.0 ? regions[0]
                : (want[i] >= kMinNormal ? regions[1] : regions[2]);
    const uint64_t u = UlpDiff(got[i], want[i]);
    ++r.n;
    if (u > 0) {
      ++r.miss;
    }
    if (u > r.max_ulp) {
      r.max_ulp = u;
      r.worst_x = in[i];
    }
  }

  int rc = 0;
  for (const Region& r : regions) {
    std::printf("%-15s n=%6zu  max ULP=%3llu (gate %llu)  not-CR: %zu (%.2f%%)  worst x=%.17g\n",
                r.name, r.n, static_cast<unsigned long long>(r.max_ulp),
                static_cast<unsigned long long>(r.bound), r.miss,
                r.n ? 100.0 * static_cast<double>(r.miss) / static_cast<double>(r.n) : 0.0,
                r.worst_x);
    if (r.max_ulp > r.bound) {
      std::fprintf(stderr, "FAIL: %s exceeds gate\n", r.name);
      rc = 1;
    }
  }
  if (rc == 0) {
    std::printf("PASS: all regions within gates\n");
  }
  return rc;
}
