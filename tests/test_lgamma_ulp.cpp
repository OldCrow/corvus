// Measures max ULP deviation of corvus::lgamma against the mpmath-generated
// correctly-rounded reference, split by the kernel's regions.
//
// The negative axis is reported in two bands, and that split is the point
// rather than a convenience. lgamma has zeros there wherever |Gamma(x)| = 1 --
// infinitely many, none with a closed form, so none reproducible by an
// exact-argument form the way x = 1 and x = 2 are. Near them a fixed absolute
// error is an unbounded relative one, so the |result| < 1 band is gated on
// absolute error and the rest on ULP. Reporting one blended number would hide
// which of the two is actually moving.
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "corvus/corvus.h"
#include "expect_target.h"
#include "ulp_utils.h"

namespace {

using corvus_test::SameBits;
using corvus_test::UlpDiff;

// Gates, set from measured values (see docs/ACCURACY.md) with no
// margin: regressions should trip them. Measured identical on AVX3_ZEN4,
// AVX2, SSE4, SSSE3 and SSE2.
constexpr uint64_t kMaxUlpSmall = 1;   // 0 < x < 1/2   (shift by one log)
constexpr uint64_t kMaxUlpZone = 1;    // 1/2 <= x <= 5/2 (the zeros live here)
constexpr uint64_t kMaxUlpMid = 1;     // 5/2 < x < 8   (product recurrence)
// Stirling is correctly rounded on every reference point, so its gate is 0 --
// the tightest in the project, and deliberately so: this region has no fitted
// approximation carrying error, only log_dd and dd assembly, and a single ULP
// appearing here would mean one of those moved.
constexpr uint64_t kMaxUlpBig = 0;     // x >= 8        (Stirling)
constexpr uint64_t kMaxUlpNeg = 1;     // x < 0, |lgamma| >= 1
constexpr double kMaxAbsNegSmall = 0x1p-53;  // x < 0, |lgamma| < 1

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

  const char* path = argc > 1 ? argv[1] : "tests/data/lgamma_reference.txt";
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
    corvus::lgamma(in_s.first(n - 3), got_s.first(n - 3));
    corvus::lgamma(in_s.last(3), got_s.last(3));
  }

  Region regions[] = {
      {"pos small", kMaxUlpSmall}, {"pos zone", kMaxUlpZone},
      {"pos mid", kMaxUlpMid},     {"pos big", kMaxUlpBig},
      {"neg |lg|>=1", kMaxUlpNeg},
  };
  double neg_small_abs = 0.0;
  double neg_small_worst_x = 0.0;
  size_t neg_small_n = 0, neg_small_miss = 0, neg_small_zero_bad = 0;

  for (size_t i = 0; i < in.size(); ++i) {
    const double x = in[i];
    if (x < 0.0 && std::abs(want[i]) < 1.0) {
      ++neg_small_n;
      if (want[i] == 0.0) {
        // #14 N6: a plain subtraction maps +0/-0 to the same magnitude just
        // like UlpDiff does (see ulp_utils.h's policy comment), so a
        // signed-zero regression on this branch would be just as invisible
        // to `e` below. Check bit-exact instead whenever the reference is
        // exactly zero.
        if (!SameBits(got[i], want[i])) {
          ++neg_small_miss;
          ++neg_small_zero_bad;
          std::fprintf(stderr,
                       "FAIL: neg |lg|<1 signed-zero mismatch at x=%.17g: "
                       "got=%.17g want=%.17g\n",
                       x, got[i], want[i]);
        }
        continue;
      }
      const double e = std::abs(got[i] - want[i]);
      if (got[i] != want[i]) ++neg_small_miss;
      if (e > neg_small_abs) {
        neg_small_abs = e;
        neg_small_worst_x = x;
      }
      continue;
    }
    Region& r = x < 0.0        ? regions[4]
                : x < 0.5      ? regions[0]
                : x <= 2.5     ? regions[1]
                : x < 8.0      ? regions[2]
                                : regions[3];
    // #14 N6: a reference value of exactly 0.0 (the zeros at x = 1 and
    // x = 2) must be checked bit-exact -- UlpDiff maps +0 and -0 to the
    // same point, so it cannot see a sign regression there.
    const bool zero_ref = want[i] == 0.0;
    const uint64_t u = zero_ref ? (SameBits(got[i], want[i]) ? 0 : UINT64_MAX)
                                 : UlpDiff(got[i], want[i]);
    if (zero_ref && u != 0) {
      std::fprintf(stderr,
                   "FAIL: signed-zero mismatch at x=%.17g: got=%.17g want=%.17g\n",
                   x, got[i], want[i]);
    }
    ++r.n;
    if (u > 0) ++r.miss;
    if (u > r.max_ulp) {
      r.max_ulp = u;
      r.worst_x = x;
    }
  }

  // #14 N10: a routing-threshold edit that empties a bucket must fail the
  // gate, not silently disable it by reporting on zero rows.
  for (const Region& r : regions) {
    if (r.n == 0) {
      std::fprintf(stderr, "FAIL: bucket '%s' is empty (0 rows)\n", r.name);
      return 1;
    }
  }
  if (neg_small_n == 0) {
    std::fprintf(stderr, "FAIL: bucket 'neg |lg|<1' is empty (0 rows)\n");
    return 1;
  }

  int rc = neg_small_zero_bad > 0 ? 1 : 0;
  for (const Region& r : regions) {
    std::printf(
        "%-13s n=%6zu  max ULP=%3llu (gate %llu)  not-CR: %zu (%.2f%%)  "
        "worst x=%.17g\n",
        r.name, r.n, static_cast<unsigned long long>(r.max_ulp),
        static_cast<unsigned long long>(r.bound), r.miss,
        r.n ? 100.0 * static_cast<double>(r.miss) / static_cast<double>(r.n)
            : 0.0,
        r.worst_x);
    if (r.max_ulp > r.bound) {
      std::fprintf(stderr, "FAIL: %s exceeds gate\n", r.name);
      rc = 1;
    }
  }
  std::printf(
      "%-13s n=%6zu  max abs=%.3e (gate %.1e)  not-CR: %zu (%.2f%%)  "
      "worst x=%.17g\n",
      "neg |lg|<1", neg_small_n, neg_small_abs, kMaxAbsNegSmall, neg_small_miss,
      neg_small_n ? 100.0 * static_cast<double>(neg_small_miss) /
                        static_cast<double>(neg_small_n)
                  : 0.0,
      neg_small_worst_x);
  if (neg_small_abs > kMaxAbsNegSmall) {
    std::fprintf(stderr, "FAIL: neg |lg|<1 exceeds absolute gate\n");
    rc = 1;
  }

  if (rc == 0) std::printf("PASS: all regions within gates\n");
  return rc;
}
