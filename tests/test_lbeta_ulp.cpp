// Measures corvus::lbeta against the mpmath-generated correctly-rounded
// reference (tests/data/lbeta_reference.txt, rows "a b lbeta(a,b)" in hex),
// bucketed BY FORMULA:
//   * relative rows, |ln B| >= 1: ULP metric, gate pinned to measured;
//   * the zero-manifold band, |ln B| < 1: ln B has a zero curve through
//     (1,1) where the result is inherently relative-ill-conditioned (the
//     lgamma-negative-axis precedent), so the contract there is ABSOLUTE,
//     in units of 2^-53, gate pinned to measured;
//   * the big band, min(a,b) > kLbetaBigMin = 2^990, reported as its own
//     bucket so the Stirling-direct path's quality names itself;
//   * -inf rows (true ln B below the double range): exact bit match, no
//     ULP distance (same rationale as the bessel boundary rows).
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#include "corvus/corvus.h"
#include "expect_target.h"
#include "ulp_utils.h"

namespace {

using corvus_test::OrderedBits;
using corvus_test::UlpDiff;

// Gates PINNED to measured, no margin: lbeta
// measured CORRECTLY ROUNDED on every row of every band -- relative and
// big-band max 0 ULP (4569 + 165 rows), absolute band max 0.500 x 2^-53
// (the final rounding itself), -inf boundary 44/44 exact (see
// docs/ACCURACY.md).
constexpr uint64_t kMaxUlp = 0;         // relative rows and big-band rows
constexpr double kMaxAbs53 = 0.5;       // |lnB| < 1 band, units of 2^-53

constexpr double kBigMin = 0x1.0p+990;  // mirrors kLbetaBigMin (beta-inl.h)

}  // namespace

int main(int argc, char** argv) {
  if (!corvus_test::ReportAndCheckTarget()) return 2;

  const char* path = argc > 1 ? argv[1] : "tests/data/lbeta_reference.txt";
  const auto rows = corvus_test::LoadRef(path, 3, 2000);

  std::vector<double> a, b, want;
  a.reserve(rows.size());
  b.reserve(rows.size());
  want.reserve(rows.size());
  for (const auto& row : rows) {
    a.push_back(corvus_test::ParseDouble(row.tok[0], path, row.line));
    b.push_back(corvus_test::ParseDouble(row.tok[1], path, row.line));
    want.push_back(corvus_test::ParseDouble(row.tok[2], path, row.line));
  }

  std::vector<double> got(a.size());
  corvus::lbeta(a, b, got);

  uint64_t rel_max = 0, big_max = 0;
  size_t rel_n = 0, rel_miss = 0, big_n = 0, big_miss = 0;
  double rel_worst_a = 0, rel_worst_b = 0, big_worst_a = 0, big_worst_b = 0;
  double abs_max = 0;
  size_t abs_n = 0;
  double abs_worst_a = 0, abs_worst_b = 0;
  size_t inf_n = 0, inf_fail = 0;

  for (size_t i = 0; i < a.size(); ++i) {
    if (std::isinf(want[i])) {
      ++inf_n;
      if (got[i] != want[i]) {
        ++inf_fail;
        std::fprintf(stderr, "FAIL: lbeta(%.17g, %.17g) = %.17g, want %.17g\n",
                     a[i], b[i], got[i], want[i]);
      }
      continue;
    }
    if (std::min(a[i], b[i]) > kBigMin) {
      const uint64_t u = UlpDiff(got[i], want[i]);
      ++big_n;
      if (u > 0) ++big_miss;
      if (u > big_max) {
        big_max = u;
        big_worst_a = a[i];
        big_worst_b = b[i];
      }
    } else if (std::fabs(want[i]) >= 1.0) {
      const uint64_t u = UlpDiff(got[i], want[i]);
      ++rel_n;
      if (u > 0) ++rel_miss;
      if (u > rel_max) {
        rel_max = u;
        rel_worst_a = a[i];
        rel_worst_b = b[i];
      }
    } else {
      const double e = std::fabs(got[i] - want[i]) * 0x1.0p+53;
      ++abs_n;
      if (e > abs_max) {
        abs_max = e;
        abs_worst_a = a[i];
        abs_worst_b = b[i];
      }
    }
  }

  std::printf("lbeta (%s):\n", path);
  std::printf("  relative (|lnB|>=1) n=%6zu  max ULP=%3llu (gate %llu)  "
              "not-CR: %zu (%.2f%%)  worst a=%.17g b=%.17g\n",
              rel_n, static_cast<unsigned long long>(rel_max),
              static_cast<unsigned long long>(kMaxUlp), rel_miss,
              rel_n ? 100.0 * static_cast<double>(rel_miss) /
                          static_cast<double>(rel_n)
                    : 0.0,
              rel_worst_a, rel_worst_b);
  std::printf("  big band (m>2^990)  n=%6zu  max ULP=%3llu (gate %llu)  "
              "not-CR: %zu  worst a=%.17g b=%.17g\n",
              big_n, static_cast<unsigned long long>(big_max),
              static_cast<unsigned long long>(kMaxUlp), big_miss, big_worst_a,
              big_worst_b);
  std::printf("  abs band (|lnB|<1)  n=%6zu  max abs=%.3f x 2^-53 (gate "
              "%.1f)  worst a=%.17g b=%.17g\n",
              abs_n, abs_max, kMaxAbs53, abs_worst_a, abs_worst_b);
  std::printf("  -inf boundary       n=%6zu  exact-match fails=%zu\n", inf_n,
              inf_fail);

  int rc = 0;
  if (rel_max > kMaxUlp) {
    std::fprintf(stderr, "FAIL: relative band exceeds gate\n");
    rc = 1;
  }
  if (big_max > kMaxUlp) {
    std::fprintf(stderr, "FAIL: big band exceeds gate\n");
    rc = 1;
  }
  if (abs_max > kMaxAbs53) {
    std::fprintf(stderr, "FAIL: absolute band exceeds gate\n");
    rc = 1;
  }
  if (inf_fail != 0) rc = 1;

  if (rc == 0) std::printf("PASS: all regions within gates\n");
  return rc;
}
