// Measures max ULP deviation of corvus::gamma_p / corvus::gamma_q against
// the mpmath / exact-Temme reference, broken down the way the kernel is
// actually built: by region AND by whether the value under test is the side
// the kernel computed DIRECTLY or the 1 (-) direct complement.
//
// The split matters more here than in the other ULP tests. Every region
// computes the smaller of the pair directly and gets the other by
// subtraction from one; the direct side is the one carrying a relative
// accuracy claim down to the subnormal band, while the complement is a
// different statement: usually large, but as small as ~0.13 in R1 (#56),
// where it inherits the direct side's error amplified ~5x. Reporting them
// together would average two different claims.
//
// The routing below is re-derived from the same src/gamma_data.h constants
// the kernel reads, so the two cannot drift apart silently.
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <span>
#include <string>
#include <vector>

#include "corvus/corvus.h"
#include "expect_target.h"
#include "src/gamma_data.h"
#include "src/lgamma_data.h"  // R4's a-bound is lgamma's centre-2 zone, shifted
#include "ulp_utils.h"

namespace {

using corvus_test::LoadRef;
using corvus_test::ParseDouble;
using corvus_test::SameBits;
using corvus_test::UlpDiff;

// Gates pinned to the values measured, no margin, per house rule. Indexed
// [region][0 = direct, 1 = complement]; cells that hold no points for a
// function carry 0 so any routing change that populates them trips the gate
// and forces a re-measure. R1 and R2 re-pinned 2026-10-19 (#56: dd series
// terms and dd CF tail; the reference set gained 1,024 ridge rows for
// a < 20 and 302 Temme-tail rows for #57), measured on Zen 4 clang-cl at
// AVX3_ZEN4, SSE4 and SSE2: R1 P 2 -> 1, R1 Q complement 4 -> 1 (was 6 on
// the new rows), R2 Q 2 -> 1, R2 P complement 1 -> 0. R3 Q stays 2: the
// new tail rows read 3 before #57's z.lo fix and 2 after (#64 tracks the
// remaining real-error margin). Earlier pins: Kaby Lake, AppleClang, AVX2
// plus the SSE4/SSSE3/SSE2 caps. Row 12a (#67, 2026-10-20): 1,700 rows in
// R4's small-a band; R4 Q 4 -> 1 with lgamma(1+a) from the centre-2 zone,
// and R4 P complement 0 -> 1 (2 rows, 0.03%; the published bound is 1 ULP
// since row 12's erratum). R2 P complement still measures 0 on these rows,
// but its published bound is also 1 ULP (dense seam rows miss by 1).
constexpr uint64_t kGateP[4][2] = {
    {1, 0},  // R1: direct P max 1
    {0, 0},  // R2: complement P max 0 on these rows (bound 1 ULP)
    {2, 1},  // R3
    {0, 1},  // R4: complement P max 1
};
constexpr uint64_t kGateQ[4][2] = {
    {0, 1},  // R1: complement Q max 1
    {1, 0},  // R2: direct Q max 1
    {2, 1},  // R3
    {1, 0},  // R4: direct Q max 1
};

struct Region {
  const char* name;
  uint64_t bound;
  uint64_t max_ulp = 0;
  size_t n = 0;
  size_t miss = 0;
  double worst_a = 0.0;
  double worst_x = 0.0;
};

// Region codes, in the kernel's own order.
enum : int { kR1 = 0, kR2 = 1, kR3 = 2, kR4 = 3 };

// Re-derivation of GammaVec's router. `want_p` selects which function's
// overlap rule applies; `direct_is_p` reports which side that region
// computes directly at this (a, x).
int Route(double a, double x, bool want_p, bool* direct_is_p) {
  const bool small = a < corvus::detail::kGammaAT;
  const bool xle = x <= a + 1.0;
  const bool lo = 2.0 * x <= a;
  const bool hi = x >= 2.0 * a;
  const bool box = a <= corvus::detail::kLgammaZoneHi - 1.0 && x <= 4.0;

  const bool r1 = small ? xle : lo;
  const bool r3 = !small && !lo && !hi;

  const bool r4 = want_p ? (box && !xle) : box;
  if (r4) {
    *direct_is_p = false;
    return kR4;
  }
  if (r1) {
    *direct_is_p = true;
    return kR1;
  }
  if (r3) {
    *direct_is_p = x < a;
    return kR3;
  }
  *direct_is_p = false;
  return kR2;
}

// LoadRef exits on a missing, malformed or short file, so this cannot fail.
void LoadReference(const char* path, std::vector<double>* a,
                   std::vector<double>* x, std::vector<double>* p,
                   std::vector<double>* q) {
  const auto rows = LoadRef(path, 4, 10000);
  for (const auto& row : rows) {
    a->push_back(ParseDouble(row.tok[0], path, row.line));
    x->push_back(ParseDouble(row.tok[1], path, row.line));
    p->push_back(ParseDouble(row.tok[2], path, row.line));
    q->push_back(ParseDouble(row.tok[3], path, row.line));
  }
}

// `live[i]` says whether region-cell i CAN receive any row for the function
// under measurement. R1 is always P-direct and R2/R4 are always Q-direct
// (Route()'s fixed orientation, independent of (a, x)), so for any single
// call exactly one of each region's {dir, cmp} cells is structurally unfired
// -- that is the "cells that hold no points for a function carry 0" case the
// gate table above already documents, not a regression. R3's direct side
// varies with (a, x), so both of its cells are always live.
int ReportRegions(const char* label, const Region* r, const bool* live,
                  int n) {
  int rc = 0;
  for (int i = 0; i < n; ++i) {
    std::printf(
        "%-8s %-14s n=%6zu  max ULP=%3llu (gate %llu)  not-CR: %zu (%.2f%%)  "
        "worst a=%.17g x=%.17g\n",
        label, r[i].name, r[i].n,
        static_cast<unsigned long long>(r[i].max_ulp),
        static_cast<unsigned long long>(r[i].bound), r[i].miss,
        r[i].n ? 100.0 * static_cast<double>(r[i].miss) /
                     static_cast<double>(r[i].n)
               : 0.0,
        r[i].worst_a, r[i].worst_x);
    // N10: a LIVE bucket that never received a row passes every gate above
    // vacuously -- a routing or reference-set regression, not a clean run.
    if (live[i] && r[i].n == 0) {
      std::fprintf(stderr, "FAIL: %s %s bucket is empty (n=0)\n", label,
                   r[i].name);
      std::exit(1);
    }
    if (r[i].n > 0 && r[i].max_ulp > r[i].bound) {
      std::fprintf(stderr, "FAIL: %s %s exceeds gate\n", label, r[i].name);
      rc = 1;
    }
  }
  return rc;
}

// N4: split every whole-reference-set call into [0, n-3) and [n-3, n) so the
// trailing 3-row group always runs a masked tail, independent of whether n
// itself happens to be a lane multiple on the tier under test. when #14 added it,
// gamma's 16734 rows were even (no tail on any 2-lane tier) and gammainv-q's
// 6520 rows were a multiple of 8 (no tail on any tier at all) -- exactly the
// coverage gap this closes. Downstream per-row loops are untouched: they
// only see the fully populated `got` vector.
void CallSplit(void (*fn)(std::span<const double>, std::span<const double>,
                          std::span<double>),
               const std::vector<double>& a, const std::vector<double>& x,
               std::vector<double>& out) {
  const size_t n = a.size();
  const size_t split = n - 3;
  fn(std::span<const double>(a).first(split),
     std::span<const double>(x).first(split),
     std::span<double>(out).first(split));
  fn(std::span<const double>(a).subspan(split),
     std::span<const double>(x).subspan(split),
     std::span<double>(out).subspan(split));
}

int Measure(const char* label, bool want_p, const std::vector<double>& a,
            const std::vector<double>& x, const std::vector<double>& got,
            const std::vector<double>& want) {
  // Index [region][0 = direct side, 1 = complement].
  const auto& g = want_p ? kGateP : kGateQ;
  Region reg[8] = {
      {"R1 series dir", g[0][0]},   {"R1 series cmp", g[0][1]},
      {"R2 cf     dir", g[1][0]},   {"R2 cf     cmp", g[1][1]},
      {"R3 temme  dir", g[2][0]},   {"R3 temme  cmp", g[2][1]},
      {"R4 smalla dir", g[3][0]},   {"R4 smalla cmp", g[3][1]},
  };
  int signed_zero_fail = 0;
  size_t zero_ref_rows = 0;
  for (size_t i = 0; i < a.size(); ++i) {
    bool direct_is_p = false;
    const int code = Route(a[i], x[i], want_p, &direct_is_p);
    const bool is_direct = (direct_is_p == want_p);
    Region& r = reg[2 * code + (is_direct ? 0 : 1)];
    uint64_t u;
    // N6: a reference value of exactly +/-0 has no ULP neighbourhood --
    // UlpDiff maps both zeros to the same point (ulp_utils.h's policy
    // comment), so it cannot see a sign regression. Check bit-exactness
    // instead and fail with a dedicated message on any mismatch.
    if (want[i] == 0.0) {
      ++zero_ref_rows;
      if (!SameBits(got[i], want[i])) {
        std::fprintf(stderr,
                     "FAIL: %s signed-zero mismatch at a=%.17g x=%.17g "
                     "got=%a want=%a\n",
                     label, a[i], x[i], got[i], want[i]);
        signed_zero_fail = 1;
      }
      u = 0;
    } else {
      u = UlpDiff(got[i], want[i]);
    }
    ++r.n;
    if (u > 0) ++r.miss;
    if (u > r.max_ulp) {
      r.max_ulp = u;
      r.worst_a = a[i];
      r.worst_x = x[i];
    }
  }
  std::printf("%-8s %zu exact-zero reference rows (checked via SameBits)\n",
              label, zero_ref_rows);
  // R1's dir/cmp split and R2/R4's are each fixed by Route() regardless of
  // (a, x) -- see the comment on ReportRegions -- so exactly one cell of
  // each is dead for this want_p. R3 varies with (a, x); both its cells
  // (indices 4, 5) stay live.
  bool live[8] = {true, true, true, true, true, true, true, true};
  live[want_p ? 1 : 0] = false;  // R1: the side Route() never fires here
  live[want_p ? 2 : 3] = false;  // R2: ditto
  live[want_p ? 6 : 7] = false;  // R4: ditto
  return ReportRegions(label, reg, live, 8) | signed_zero_fail;
}

}  // namespace

int main(int argc, char* const* argv) {
  if (!corvus_test::ReportAndCheckTarget()) return 2;

  const char* p_path = argc > 1 ? argv[1] : "tests/data/gamma_p_reference.txt";
  const char* q_path = argc > 2 ? argv[2] : "tests/data/gamma_q_reference.txt";

  int rc = 0;
  {
    std::vector<double> a, x, p, q;
    LoadReference(p_path, &a, &x, &p, &q);
    std::vector<double> got(a.size());
    CallSplit(corvus::gamma_p, a, x, got);
    rc |= Measure("gamma_p", true, a, x, got, p);
  }
  {
    std::vector<double> a, x, p, q;
    LoadReference(q_path, &a, &x, &p, &q);
    std::vector<double> got(a.size());
    CallSplit(corvus::gamma_q, a, x, got);
    rc |= Measure("gamma_q", false, a, x, got, q);
  }

  if (rc == 0) std::printf("PASS: all regions within gates\n");
  return rc;
}
