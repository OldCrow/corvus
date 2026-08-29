// Measures max ULP deviation of corvus::digamma against the mpmath-generated
// correctly-rounded reference, split by the kernel's regions.
//
// The negative axis is reported in two bands, and that split is the point
// rather than a convenience. psi has a zero between every consecutive pair of
// negative poles -- infinitely many, none with a closed form, so none
// reproducible by an exact-argument form the way the positive root at
// x0 ~ 1.4616 is. Near them a fixed absolute error is an unbounded relative
// one, so the |psi| < 1 band is gated on ABSOLUTE error (reported in units of
// 2^-53) and the rest on ULP. That is lgamma's convention verbatim, and it is
// the family's dual metric: relative where |psi| >= 1, 2^-53-class
// absolute inside the zero bands. Reporting one blended number would hide
// which of the two is actually moving.
//
// The positive-axis buckets follow the kernel's own region split, so a
// regression names the branch that moved.
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

// Gates PINNED to measured, no margin. Identical cells on
// every validated leg: AVX3_ZEN4 native, AVX2/SSE4/SSSE3/SSE2 capped
// (Ryzen), Linux CI sweep, NEON (CI), Windows MSVC — 1 ULP max in all five
// relative buckets, 1.00 x 2^-53 absolute in the negative zero bands (NEON
// matches native down to the not-CR counts).
constexpr uint64_t kMaxUlpRel = 1;
constexpr double kMaxAbsUnits = 1.0;  // in units of 2^-53
constexpr double kAbsUnit = 0x1p-53;

struct Region {
  const char* name;
  uint64_t max_ulp = 0;
  size_t n = 0;
  size_t miss = 0;  // not correctly rounded
  double worst_x = 0.0;
};

}  // namespace

int main(int argc, char** argv) {
  // Before doing any work: a wrong tier makes every number below meaningless.
  if (!corvus_test::ReportAndCheckTarget()) return 2;

  const char* path = argc > 1 ? argv[1] : "tests/data/digamma_reference.txt";
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
    corvus::digamma(in_s.first(n - 3), got_s.first(n - 3));
    corvus::digamma(in_s.last(3), got_s.last(3));
  }

  Region regions[] = {
      {"pos (0,1)"},   {"pos zone"},     {"pos walk"},
      {"pos asym"},    {"neg |psi|>=1"},
  };
  double neg_small_abs = 0.0;
  double neg_small_worst_x = 0.0;
  size_t neg_small_n = 0, neg_small_miss = 0, neg_small_zero_bad = 0;

  for (size_t i = 0; i < in.size(); ++i) {
    const double x = in[i];
    // The lgamma convention: on the negative axis, the metric is chosen by
    // |psi| alone. |psi| < 1 happens only inside the zero bands (|psi| grows
    // past 1 within ~1/8 of the way to either neighbouring pole), so this one
    // test reproduces the per-zero band membership the design describes
    // without carrying a table of zeros the kernel does not have either.
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
                       "FAIL: neg |psi|<1 signed-zero mismatch at x=%.17g: "
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
    Region& r = x < 0.0   ? regions[4]
                : x < 1.0 ? regions[0]
                : x < 2.0 ? regions[1]
                : x < 8.0 ? regions[2]
                          : regions[3];
    // #14 N6: a reference value of exactly 0.0 must be checked bit-exact --
    // UlpDiff maps +0 and -0 to the same point, so it cannot see a sign
    // regression there. This reference set has no such row today (digamma's
    // only positive zero, x ~ 1.4616, is irrational, so no sampled reference
    // row lands exactly on 0.0); the check exists for the next regen.
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
    std::fprintf(stderr, "FAIL: bucket 'neg |psi|<1' is empty (0 rows)\n");
    return 1;
  }

  int rc = neg_small_zero_bad > 0 ? 1 : 0;
  for (const Region& r : regions) {
    std::printf(
        "%-13s n=%6zu  max ULP=%3llu (gate %llu)  not-CR: %zu (%.2f%%)  "
        "worst x=%.17g\n",
        r.name, r.n, static_cast<unsigned long long>(r.max_ulp),
        static_cast<unsigned long long>(kMaxUlpRel), r.miss,
        r.n ? 100.0 * static_cast<double>(r.miss) / static_cast<double>(r.n)
            : 0.0,
        r.worst_x);
    if (r.max_ulp > kMaxUlpRel) {
      std::fprintf(stderr, "FAIL: %s exceeds gate\n", r.name);
      rc = 1;
    }
  }
  std::printf(
      "%-13s n=%6zu  max abs=%.3e = %.2f x 2^-53 (gate %.1f)  not-CR: %zu "
      "(%.2f%%)  worst x=%.17g\n",
      "neg |psi|<1", neg_small_n, neg_small_abs, neg_small_abs / kAbsUnit,
      kMaxAbsUnits, neg_small_miss,
      neg_small_n ? 100.0 * static_cast<double>(neg_small_miss) /
                        static_cast<double>(neg_small_n)
                  : 0.0,
      neg_small_worst_x);
  if (neg_small_abs > kMaxAbsUnits * kAbsUnit) {
    std::fprintf(stderr, "FAIL: neg |psi|<1 exceeds absolute gate\n");
    rc = 1;
  }

  if (rc == 0) std::printf("PASS: all regions within gates\n");
  return rc;
}
