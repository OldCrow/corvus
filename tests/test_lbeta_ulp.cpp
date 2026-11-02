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

using corvus_test::SameBits;
using corvus_test::UlpDiff;

// Gates: lbeta was documented and gated as CORRECTLY ROUNDED (relative and
// big-band 0 ULP, absolute band 0.5 x 2^-53). #76 (order row 12j, review
// M-7) found near-tie misses: lbeta(0x1.5c788bfb641a0p+1,
// 0x1.5f0d82e6d1d90p+0) is 1 ULP and lbeta(0x1.62f92a1b0bb25p+0,
// 0x1.9546088f2e3ccp+0) 0.5012 x 2^-53 absolute, both 0.5012 ulp real. The
// contract is now 1 ULP / 1 x 2^-53, measured CR on the reference set
// except such ties (see docs/ACCURACY.md); the -inf boundary stays exact.
constexpr uint64_t kMaxUlp = 1;         // relative rows and big-band rows
constexpr double kMaxAbs53 = 1.0;       // |lnB| < 1 band, units of 2^-53
// a = 1 or b = 1: -ln x from one log (#80), judged relatively at every
// magnitude. The absolute band would hide the lgamma pair's ~2^-74 floor here.
constexpr uint64_t kMaxUlpUnit = 0;

constexpr double kBigMin = 0x1.0p+990;  // mirrors kLbetaBigMin (beta-inl.h)

// #14 N4: masked-tail split. Splitting into subspans [0, n-3) and [n-3, n)
// makes neither half a lane multiple regardless of what n is, so the
// masked-tail path always runs; each call is independently masked and
// stateless, so this is row-for-row identical to one whole-set call.
void SplitCall(const std::vector<double>& a, const std::vector<double>& b,
               std::vector<double>& out) {
  const size_t n = out.size();
  const size_t split = n - 3;
  corvus::lbeta(std::span<const double>(a).subspan(0, split),
                std::span<const double>(b).subspan(0, split),
                std::span<double>(out).subspan(0, split));
  corvus::lbeta(std::span<const double>(a).subspan(split),
                std::span<const double>(b).subspan(split),
                std::span<double>(out).subspan(split));
}

}  // namespace

int main(int argc, char* const* argv) {
  if (!corvus_test::ReportAndCheckTarget()) return 2;

  const char* path = argc > 1 ? argv[1] : "tests/data/lbeta_reference.txt";
  const auto rows = corvus_test::LoadRef(path, 3, 2000);

  std::vector<double> a, b, want, res;
  a.reserve(rows.size());
  b.reserve(rows.size());
  want.reserve(rows.size());
  for (const auto& row : rows) {
    a.push_back(corvus_test::ParseDouble(row.tok[0], path, row.line));
    b.push_back(corvus_test::ParseDouble(row.tok[1], path, row.line));
    want.push_back(corvus_test::ParseDouble(row.tok[2], path, row.line));
    res.push_back(corvus_test::Residual(row, 0));
  }

  std::vector<double> got(a.size());
  SplitCall(a, b, got);

  uint64_t rel_max = 0, big_max = 0;
  size_t rel_n = 0, rel_miss = 0, big_n = 0, big_miss = 0;
  double rel_worst_a = 0, rel_worst_b = 0, big_worst_a = 0, big_worst_b = 0;
  double abs_max = 0;
  size_t abs_n = 0;
  double abs_worst_a = 0, abs_worst_b = 0;
  size_t inf_n = 0, inf_fail = 0;
  corvus_test::RealMax rel_real, big_real, unit_real;
  uint64_t unit_max = 0;
  size_t unit_n = 0, unit_miss = 0;
  double unit_worst_a = 0, unit_worst_b = 0;

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
    if (a[i] == 1.0 || b[i] == 1.0) {
      // ln B(1, 1) = +0 is judged bit-exact: UlpDiff cannot see the sign.
      const bool zero_ref = want[i] == 0.0;
      const uint64_t u = zero_ref ? (SameBits(got[i], want[i]) ? 0 : UINT64_MAX)
                                  : UlpDiff(got[i], want[i]);
      ++unit_n;
      if (u > 0) ++unit_miss;
      if (!zero_ref) unit_real.Add(corvus_test::RealErr(got[i], want[i], res[i]), i);
      if (u > unit_max) {
        unit_max = u;
        unit_worst_a = a[i];
        unit_worst_b = b[i];
      }
    } else if (std::min(a[i], b[i]) > kBigMin) {
      const uint64_t u = UlpDiff(got[i], want[i]);
      ++big_n;
      if (u > 0) ++big_miss;
      big_real.Add(corvus_test::RealErr(got[i], want[i], res[i]), i);
      if (u > big_max) {
        big_max = u;
        big_worst_a = a[i];
        big_worst_b = b[i];
      }
    } else if (std::fabs(want[i]) >= 1.0) {
      const uint64_t u = UlpDiff(got[i], want[i]);
      ++rel_n;
      if (u > 0) ++rel_miss;
      rel_real.Add(corvus_test::RealErr(got[i], want[i], res[i]), i);
      if (u > rel_max) {
        rel_max = u;
        rel_worst_a = a[i];
        rel_worst_b = b[i];
      }
    } else {
      // The absolute band is judged in absolute units, so it carries no
      // real gate (nor do the exact-match -inf rows above).
      const double e = std::fabs(got[i] - want[i]) * 0x1.0p+53;
      ++abs_n;
      if (e > abs_max) {
        abs_max = e;
        abs_worst_a = a[i];
        abs_worst_b = b[i];
      }
    }
  }

  int rc = 0;
  std::printf("lbeta (%s):\n", path);
  std::printf("  relative (|lnB|>=1) n=%6zu  max ULP=%3llu (gate %llu)  "
              "not-CR: %zu (%.2f%%)  worst a=%.17g b=%.17g\n",
              rel_n, static_cast<unsigned long long>(rel_max),
              static_cast<unsigned long long>(kMaxUlp), rel_miss,
              rel_n ? 100.0 * static_cast<double>(rel_miss) /
                          static_cast<double>(rel_n)
                    : 0.0,
              rel_worst_a, rel_worst_b);
  // #82 T1: real error per ULP-judged band.
  if (rel_n != 0) {
    rc |= corvus_test::CheckReal(path, "lbeta", "relative", kMaxUlp, rel_real);
  }
  std::printf("  unit line (a|b=1)   n=%6zu  max ULP=%3llu (gate %llu)  "
              "not-CR: %zu  worst a=%.17g b=%.17g\n",
              unit_n, static_cast<unsigned long long>(unit_max),
              static_cast<unsigned long long>(kMaxUlpUnit), unit_miss,
              unit_worst_a, unit_worst_b);
  if (unit_n != 0) {
    rc |= corvus_test::CheckReal(path, "lbeta", "unit line", kMaxUlpUnit,
                                 unit_real);
  }
  std::printf("  big band (m>2^990)  n=%6zu  max ULP=%3llu (gate %llu)  "
              "not-CR: %zu  worst a=%.17g b=%.17g\n",
              big_n, static_cast<unsigned long long>(big_max),
              static_cast<unsigned long long>(kMaxUlp), big_miss, big_worst_a,
              big_worst_b);
  if (big_n != 0) {
    rc |= corvus_test::CheckReal(path, "lbeta", "big band", kMaxUlp, big_real);
  }
  std::printf("  abs band (|lnB|<1)  n=%6zu  max abs=%.3f x 2^-53 (gate "
              "%.1f)  worst a=%.17g b=%.17g\n",
              abs_n, abs_max, kMaxAbs53, abs_worst_a, abs_worst_b);
  std::printf("  -inf boundary       n=%6zu  exact-match fails=%zu\n", inf_n,
              inf_fail);

  // #14 N10: a band that never accumulated a row is a gate that never ran --
  // every one of these five bands is populated by design in the certified
  // reference set (docs/ACCURACY.md), so n=0 here means a coverage gap.
  if (rel_n == 0) {
    std::fprintf(stderr, "FAIL: relative band is empty (n=0) -- vacuous gate\n");
    rc = 1;
  }
  if (unit_n == 0) {
    std::fprintf(stderr, "FAIL: unit line is empty (n=0) -- vacuous gate\n");
    rc = 1;
  }
  if (big_n == 0) {
    std::fprintf(stderr, "FAIL: big band is empty (n=0) -- vacuous gate\n");
    rc = 1;
  }
  if (abs_n == 0) {
    std::fprintf(stderr,
                 "FAIL: absolute band is empty (n=0) -- vacuous gate\n");
    rc = 1;
  }
  if (inf_n == 0) {
    std::fprintf(stderr,
                 "FAIL: -inf boundary is empty (n=0) -- vacuous gate\n");
    rc = 1;
  }
  if (rel_max > kMaxUlp) {
    std::fprintf(stderr, "FAIL: relative band exceeds gate\n");
    rc = 1;
  }
  if (unit_max > kMaxUlpUnit) {
    std::fprintf(stderr, "FAIL: unit line exceeds gate\n");
    rc = 1;
  }
  if (big_max > kMaxUlp) {
    std::fprintf(stderr, "FAIL: big band exceeds gate\n");
    rc = 1;
  }
  rc |= corvus_test::CheckBound(path, "lbeta", "absolute", "abs53", kMaxAbs53);
  if (abs_max > kMaxAbs53) {
    std::fprintf(stderr, "FAIL: absolute band exceeds gate\n");
    rc = 1;
  }
  if (inf_fail != 0) rc = 1;

  if (rc == 0) std::printf("PASS: all regions within gates\n");
  return rc;
}
