// Measures max ULP deviation of corvus::trigamma against the mpmath-generated
// correctly-rounded reference, split by the kernel's regions.
//
// ONE METRIC, EVERYWHERE. This test is structurally simpler than
// test_digamma_ulp, and deliberately so: psi_1(x) = sum over n >= 0 of
// 1/(x + n)^2 is a sum of squares, hence strictly positive wherever finite, so
// it has NO zeros -- not on the positive axis, and not on the negative axis
// either, where its per-interval minima rise monotonically from 8.933 to
// pi^2. digamma and lgamma both need an absolute band because they pass
// through zero at points with no closed form, and a fixed absolute error near
// such a point is an unbounded relative one. Nothing here does. A single
// relative (ULP) metric therefore covers the whole real line, and adding an
// absolute band would report a number about a region that does not exist.
//
// The buckets follow the kernel's own region split, so a regression names the
// branch that moved. The asymptotic bucket spans [X0, inf) as one gate: the
// 1/x-only path above kTrigammaAsymCut is the same branch with its corrections
// selected away, not a separate approximation. Its numbers are also printed
// separately, report-only, because the two sub-paths fail for different
// reasons and the split costs nothing.
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "corvus/corvus.h"
#include "expect_target.h"
#include "src/trigamma_data.h"
#include "ulp_utils.h"

namespace {

using corvus_test::SameBits;
using corvus_test::UlpDiff;

// Gate PINNED to measured, no margin. Identical cells on every validated
// leg — AVX3_ZEN4 native, AVX2/SSE4/SSSE3/SSE2 capped (Ryzen), Linux CI
// sweep, NEON (CI, identical to native including not-CR counts), Windows
// MSVC: (0,1) correctly rounded, all other buckets 1 ULP max. The
// down-walk's predicted 12.36x amplification cost did not consume a bit on
// this reference set — it shows as the walk bucket's elevated not-CR rate
// (2.97% FMA / 4.44% non-FMA) instead.
constexpr uint64_t kMaxUlp = 1;

struct Region {
  const char* name;
  uint64_t max_ulp = 0;
  size_t n = 0;
  size_t miss = 0;  // not correctly rounded
  double worst_x = 0.0;
};

void Accumulate(Region& r, double x, uint64_t u) {
  ++r.n;
  if (u > 0) ++r.miss;
  if (u > r.max_ulp) {
    r.max_ulp = u;
    r.worst_x = x;
  }
}

void Report(const Region& r, bool gated) {
  char gate[24];
  if (gated) {
    std::snprintf(gate, sizeof(gate), "gate %llu",
                  static_cast<unsigned long long>(kMaxUlp));
  } else {
    std::snprintf(gate, sizeof(gate), "report only");
  }
  std::printf(
      "%-20s n=%6zu  max ULP=%3llu (%s)  not-CR: %zu (%.2f%%)  "
      "worst x=%.17g\n",
      r.name, r.n, static_cast<unsigned long long>(r.max_ulp), gate, r.miss,
      r.n ? 100.0 * static_cast<double>(r.miss) / static_cast<double>(r.n)
          : 0.0,
      r.worst_x);
}

}  // namespace

int main(int argc, char** argv) {
  // Before doing any work: a wrong tier makes every number below meaningless.
  if (!corvus_test::ReportAndCheckTarget()) return 2;

  const char* path = argc > 1 ? argv[1] : "tests/data/trigamma_reference.txt";
  const auto rows = corvus_test::LoadRef(path, 2, 10000);

  std::vector<double> in, want;
  for (const auto& row : rows) {
    in.push_back(corvus_test::ParseDouble(row.tok[0], path, row.line));
    want.push_back(corvus_test::ParseDouble(row.tok[1], path, row.line));
  }

  std::vector<double> got(in.size());
  {
    // #14 N4: split the whole-set call in two so the masked LoadN/StoreN
    // tail path runs even when the reference set's length happens to be a
    // multiple of every SIMD tier's lane count. The final call's length (3)
    // is below every lane count, and N-3 is odd whenever N is even, so
    // neither call can land on a lane boundary either.
    const std::span<const double> in_s(in);
    const std::span<double> got_s(got);
    const size_t n = in_s.size();
    corvus::trigamma(in_s.first(n - 3), got_s.first(n - 3));
    corvus::trigamma(in_s.last(3), got_s.last(3));
  }

  // Routing mirrors src/trigamma_data.h's own constants, so the test and the
  // kernel cannot drift apart.
  Region regions[] = {
      {"pos (0,1)"}, {"pos zone [1,2)"}, {"pos walk [2,8)"},
      {"pos asym [8,inf)"}, {"negative"},
  };
  // Report-only split of the asymptotic bucket at the 1/x-only cut.
  Region asym_sub[] = {{"  ... [8,2^89)"}, {"  ... [2^89,inf)"}};

  for (size_t i = 0; i < in.size(); ++i) {
    const double x = in[i];
    // #14 N6: a reference value of exactly 0.0 must be checked bit-exact --
    // UlpDiff maps +0 and -0 to the same point, so it cannot see a sign
    // regression there. psi_1 is a sum of squares and has no zeros (see the
    // file header), so this reference set has no such row today; the check
    // exists for the next regen.
    const bool zero_ref = want[i] == 0.0;
    const uint64_t u = zero_ref ? (SameBits(got[i], want[i]) ? 0 : UINT64_MAX)
                                 : UlpDiff(got[i], want[i]);
    if (zero_ref && u != 0) {
      std::fprintf(stderr,
                   "FAIL: signed-zero mismatch at x=%.17g: got=%.17g want=%.17g\n",
                   x, got[i], want[i]);
    }
    if (x < 0.0) {
      Accumulate(regions[4], x, u);
    } else if (x < corvus::detail::kTrigammaZoneLo) {
      Accumulate(regions[0], x, u);
    } else if (x < corvus::detail::kTrigammaZoneHi) {
      Accumulate(regions[1], x, u);
    } else if (x < corvus::detail::kTrigammaX0) {
      Accumulate(regions[2], x, u);
    } else {
      Accumulate(regions[3], x, u);
      Accumulate(asym_sub[x < corvus::detail::kTrigammaAsymCut ? 0 : 1], x, u);
    }
  }

  // #14 N10: a routing-threshold edit that empties a bucket (including a
  // report-only one, e.g. one side of the asym_sub 1/x-only split) must fail
  // the gate, not silently disable it by reporting on zero rows.
  for (const Region& r : regions) {
    if (r.n == 0) {
      std::fprintf(stderr, "FAIL: bucket '%s' is empty (0 rows)\n", r.name);
      return 1;
    }
  }
  for (const Region& r : asym_sub) {
    if (r.n == 0) {
      std::fprintf(stderr, "FAIL: bucket '%s' is empty (0 rows)\n", r.name);
      return 1;
    }
  }

  int rc = 0;
  for (const Region& r : regions) {
    Report(r, true);
    if (r.max_ulp > kMaxUlp) {
      std::fprintf(stderr, "FAIL: %s exceeds gate\n", r.name);
      rc = 1;
    }
    if (&r == &regions[3]) {
      for (const Region& s : asym_sub) Report(s, false);
    }
  }

  if (rc == 0) std::printf("PASS: all regions within gates\n");
  return rc;
}
