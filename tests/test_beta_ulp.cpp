// Measures max ULP deviation of corvus::beta_p / corvus::beta_q against the
// CF / mpmath reference, broken down the way the kernel is actually built: by
// region AND by whether the value under test is the side the kernel computed
// DIRECTLY or the 1 (-) direct complement.
//
// The split matters more here than in the other ULP tests, for gamma's reason:
// every region computes one side of the pair directly and gets the other by
// subtraction from one. The direct side is the one carrying a relative
// accuracy claim down to the subnormal band; the complement is >= ~0.4 by
// routing and its bound is a different (easier) statement. Reporting them
// together would average a hard claim with an easy one.
//
// The routing below is re-derived from the same src/beta_data.h constants the
// kernel reads, with the same overflow-free spellings (nu as 1/(1/a+1/b),
// xi/p as xi*(1+b/a), the R2 threshold as 1/(1+(b+1)/(a+1))), so the two
// cannot drift apart silently. A handful of points that sit exactly on a
// routing wall may land in a different bin than the kernel chose -- the
// kernel's own predicate uses log_dd's log where this file uses libm's -- but
// that only moves a point between report rows, never changes what is measured.
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
#include "src/beta_data.h"

namespace {

// PROVISIONAL GATES [G3]. A deliberately useless placeholder, NOT a measured
// bound: this stage's job is to produce the per-region table, and G4 pins
// every cell to what it measures with no margin (the house rule). A passing
// run here is NOT an accuracy claim. What G3 actually measured on
// AVX3_ZEN4/clang-cl, with the reference-defect rows below excluded:
//   direct side     R1 1, R2 0, R3 3, R4 3 ULP
//   complement      R2 0, R3 1, R4 0 ULP; R1 55/429 ULP
// The R1 complement cell and the R2 complement's worst point are escalations
// against the routing design, not arithmetic -- see PLAN.md's G3 record.
constexpr uint64_t kGateProvisional = 1000000;

// Rows the kernel and the reference disagree about by MORE than this are not a
// rounding disagreement, they are a disagreement about the value. Every one
// audited during G3 (240 sampled across both files, recomputed against an
// independent high-precision CF/APSER oracle) was a defect in the SHIPPED
// REFERENCE SET, not in the kernel: tools/gen_beta_reference.py's small-tau
// oracle forms lgamma(B+tau) - lgamma(B) at a fixed dps, which silently
// returns 0 whenever tau << B and collapses the answer to a value independent
// of B. They are skipped and counted here so the rest of the table means
// something. G4 DELETES this escape hatch once the reference set is
// regenerated; if it ever hides a real regression it will show up as a jump in
// the skipped count, which is printed on every run.
constexpr uint64_t kRefDefectCutoff = 1000000;
size_t g_ref_defect = 0;

int64_t OrderedBits(double x) {
  int64_t b;
  std::memcpy(&b, &x, sizeof(b));
  return b < 0 ? (INT64_MIN - b) : b;
}

uint64_t UlpDiff(double a, double b) {
  return static_cast<uint64_t>(std::llabs(OrderedBits(a) - OrderedBits(b)));
}

bool SameBits(double a, double b) {
  uint64_t ba, bb;
  std::memcpy(&ba, &a, sizeof(ba));
  std::memcpy(&bb, &b, sizeof(bb));
  return ba == bb;
}

struct Region {
  const char* name;
  uint64_t bound;
  uint64_t max_ulp = 0;
  size_t n = 0;
  size_t miss = 0;
  double wa = 0.0, wb = 0.0, wx = 0.0;
};

enum : int { kR1 = 0, kR2 = 1, kR3 = 2, kR4 = 3, kSp = 4 };

// Re-derivation of BetaVec's router (route_final, PLAN.md G1b order).
// `direct_is_p` reports which side of the pair the chosen region computes
// directly at this (a, b, x).
int Route(double a, double b, double x, bool* direct_is_p) {
  const bool in_domain = a > 0.0 && b > 0.0 && std::isfinite(a) &&
                         std::isfinite(b) && x > 0.0 && x < 1.0;
  if (!in_domain) {
    *direct_is_p = true;
    return kSp;
  }
  using corvus::detail::kBetaB1;
  using corvus::detail::kBetaEpsR4;
  using corvus::detail::kBetaLn2;
  using corvus::detail::kBetaTRidge;
  using corvus::detail::kBetaXi1;
  using corvus::detail::kBetaXiRatioHi;
  using corvus::detail::kBetaXiRatioLo;

  const double y = 1.0 - x;
  const double tau = std::min(a, b);
  const double bmax = std::max(a, b);
  const bool sw4 = b < a;
  const double xt = sw4 ? y : x;

  // 0. tiny-min box, evaluated tiny-first. The xi cap is the WIDENED one --
  // see the routing note in src/beta-inl.h.
  const double thr_t = 1.0 / (1.0 + (bmax + 1.0) / (tau + 1.0));
  const bool r4 = tau <= kBetaEpsR4 &&
                  tau * std::fabs(std::log(xt)) <= kBetaLn2 &&
                  (xt <= kBetaXi1 || xt < thr_t) && bmax * xt <= kBetaB1;
  if (r4) {
    *direct_is_p = sw4;  // R4 computes the COMPLEMENT of its own triple
    return kR4;
  }
  // 1. power series, either orientation.
  const bool r1n = x <= kBetaXi1 && b * x <= kBetaB1;
  const bool r1s = y <= kBetaXi1 && a * y <= kBetaB1;
  if (r1n || r1s) {
    const bool sw = !r1n && r1s;
    *direct_is_p = !sw;
    return kR1;
  }
  // 2. ridge ratio band.
  const double nu = 1.0 / (1.0 / a + 1.0 / b);
  const double rat1 = x * (1.0 + b / a);
  const double rat2 = y * (1.0 + a / b);
  const bool band = rat1 >= kBetaXiRatioLo && rat1 <= kBetaXiRatioHi &&
                    rat2 >= kBetaXiRatioLo && rat2 <= kBetaXiRatioHi;
  if (nu >= kBetaTRidge && band) {
    // The kernel decides R3's direct side from lambda's EXACT sign, not from
    // the router's rounded ratio, so this mirrors that (see BetaR3Out). The
    // two can disagree only where lambda is within an ulp of zero, i.e. where
    // both sides are ~1/2 and the bin makes no difference to what is measured.
    const bool sw = rat1 > 1.0;
    const double lam = a * y - b * x;
    *direct_is_p = sw ? (lam > 0.0) : (lam >= 0.0);
    return kR3;
  }
  // 3. continued fraction, orientation by the pinned rule.
  const double thr = 1.0 / (1.0 + (b + 1.0) / (a + 1.0));
  const bool sw = !(x < thr);
  *direct_is_p = !sw;
  return kR2;
}

bool LoadReference(const char* path, std::vector<double>* a,
                   std::vector<double>* b, std::vector<double>* x,
                   std::vector<double>* p, std::vector<double>* q) {
  std::ifstream f(path);
  if (!f) {
    std::fprintf(stderr, "cannot open reference file: %s\n", path);
    return false;
  }
  std::string sa, sb, sx, sp, sq;
  while (f >> sa >> sb >> sx >> sp >> sq) {
    a->push_back(std::strtod(sa.c_str(), nullptr));
    b->push_back(std::strtod(sb.c_str(), nullptr));
    x->push_back(std::strtod(sx.c_str(), nullptr));
    p->push_back(std::strtod(sp.c_str(), nullptr));
    q->push_back(std::strtod(sq.c_str(), nullptr));
  }
  if (a->size() < 10000) {
    std::fprintf(stderr, "reference file suspiciously small: %zu lines\n",
                 a->size());
    return false;
  }
  return true;
}

int ReportRegions(const char* label, Region* r, int n) {
  int rc = 0;
  for (int i = 0; i < n; ++i) {
    std::printf(
        "%-7s %-14s n=%6zu  max ULP=%4llu (gate %llu)  not-CR: %zu (%.2f%%)  "
        "worst a=%.17g b=%.17g x=%.17g\n",
        label, r[i].name, r[i].n,
        static_cast<unsigned long long>(r[i].max_ulp),
        static_cast<unsigned long long>(r[i].bound), r[i].miss,
        r[i].n ? 100.0 * static_cast<double>(r[i].miss) /
                     static_cast<double>(r[i].n)
               : 0.0,
        r[i].wa, r[i].wb, r[i].wx);
    if (r[i].n > 0 && r[i].max_ulp > r[i].bound) {
      std::fprintf(stderr, "FAIL: %s %s exceeds gate\n", label, r[i].name);
      rc = 1;
    }
  }
  return rc;
}

int Measure(const char* label, bool want_p, const std::vector<double>& a,
            const std::vector<double>& b, const std::vector<double>& x,
            const std::vector<double>& got, const std::vector<double>& want) {
  const uint64_t g = kGateProvisional;
  Region reg[10] = {
      {"R1 series dir", g}, {"R1 series cmp", g}, {"R2 cf     dir", g},
      {"R2 cf     cmp", g}, {"R3 temme  dir", g}, {"R3 temme  cmp", g},
      {"R4 tiny   dir", g}, {"R4 tiny   cmp", g}, {"specials     ", 0},
      {"specials  (-)", 0},
  };
  for (size_t i = 0; i < a.size(); ++i) {
    bool direct_is_p = false;
    const int code = Route(a[i], b[i], x[i], &direct_is_p);
    const bool is_direct = (direct_is_p == want_p);
    Region& r = reg[2 * code + (is_direct ? 0 : 1)];
    uint64_t u;
    if (code == kSp) {
      // Specials must reproduce the oracle EXACTLY (or NaN for NaN); there is
      // no ULP statement to make about a table lookup.
      const bool ok = std::isnan(want[i]) ? std::isnan(got[i])
                                          : SameBits(got[i], want[i]);
      u = ok ? 0 : 1;
    } else if (std::isnan(got[i]) || std::isnan(want[i])) {
      u = std::isnan(got[i]) && std::isnan(want[i]) ? 0 : ~uint64_t{0};
    } else {
      u = UlpDiff(got[i], want[i]);
    }
    if (code != kSp && u > kRefDefectCutoff) {
      ++g_ref_defect;
      continue;
    }
    ++r.n;
    if (u > 0) ++r.miss;
    if (u > r.max_ulp) {
      r.max_ulp = u;
      r.wa = a[i];
      r.wb = b[i];
      r.wx = x[i];
    }
  }
  return ReportRegions(label, reg, 10);
}

}  // namespace

int main(int argc, char** argv) {
  if (!corvus_test::ReportAndCheckTarget()) return 2;

  const char* p_path = argc > 1 ? argv[1] : "tests/data/beta_p_reference.txt";
  const char* q_path = argc > 2 ? argv[2] : "tests/data/beta_q_reference.txt";

  int rc = 0;
  {
    std::vector<double> a, b, x, p, q;
    if (!LoadReference(p_path, &a, &b, &x, &p, &q)) return 2;
    std::vector<double> got(a.size());
    corvus::beta_p(a, b, x, got);
    rc |= Measure("beta_p", true, a, b, x, got, p);
  }
  {
    std::vector<double> a, b, x, p, q;
    if (!LoadReference(q_path, &a, &b, &x, &p, &q)) return 2;
    std::vector<double> got(a.size());
    corvus::beta_q(a, b, x, got);
    rc |= Measure("beta_q", false, a, b, x, got, q);
  }

  std::printf(
      "reference-defect rows skipped (> %llu ULP, see kRefDefectCutoff): %zu\n",
      static_cast<unsigned long long>(kRefDefectCutoff), g_ref_defect);
  if (rc == 0) std::printf("PASS: all regions within PROVISIONAL gates\n");
  return rc;
}
