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
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "corvus/corvus.h"
#include "expect_target.h"
#include "src/beta_data.h"
#include "ulp_utils.h"

namespace {

using corvus_test::OrderedBits;
using corvus_test::SameBits;
using corvus_test::UlpDiff;

// PINNED GATES: every cell is the value measured on AVX3_ZEN4/clang-cl
// against the certified reference set (max over the beta_p and beta_q
// runs), pinned with NO margin (the house rule). The R1-cmp and
// R4-postroute cells are wider because their double-class lgamma
// components are amplified by cancellation (src/beta-inl.h derives both
// sites). A tier that
// measures above any cell is an escalation, not a gate bump.
//
// The pb-corner reference rows stress the gammalim slice's E_g at
// b ~ 1e307 with t ~ w*m in the deep tail; 6 of 430 rows measure exactly
// 1 ULP there (worst at (40, 1e307, 1.6e-306), verified against mpmath
// betainc at dps 380 -- the kernel is one ulp off there). Every other cell
// is unaffected.
//
// R2-dir 0 -> 1 at the #13 regeneration: the single-rounding reference fix
// moved ONE subnormal-Q row (1.97e-300, 24.68, 0.474) one ulp closer to
// truth (re-verified against mpmath betainc at dps 400; the true value
// sits 0.60 ulp above the old double-rounded row). The kernel agreed with
// the OLD row, so its true error there is 1 ulp -- the old 0 was an
// artifact of the reference's own double rounding, not a kernel change.
constexpr uint64_t kGateR1Dir = 1, kGateR1Cmp = 1;
constexpr uint64_t kGateR2Dir = 1, kGateR2Cmp = 0;
constexpr uint64_t kGateR3Dir = 3, kGateR3Cmp = 1;
constexpr uint64_t kGateR4Dir = 2, kGateR4Cmp = 0;
constexpr uint64_t kGatePrDir = 1, kGatePrCmp = 0;
constexpr uint64_t kGateGlDir = 1, kGateGlCmp = 1;

struct Region {
  const char* name;
  uint64_t bound;
  uint64_t max_ulp = 0;
  size_t n = 0;
  size_t miss = 0;
  double wa = 0.0, wb = 0.0, wx = 0.0;
};

enum : int { kR1 = 0, kR2 = 1, kR3 = 2, kR4 = 3, kSp = 4, kPr = 5, kGl = 6 };

// Re-derivation of BetaVec's router (route_final, including the (C)
// gamma-limit slice). `direct_is_p` reports which
// side of the pair the chosen region computes directly at this (a, b, x).
// pref/qref are the REFERENCE values -- the post-route decision in the
// kernel is on its own dd R1 value, and the reference is the same number
// to far more bits than the 2^-11 bar needs.
int Route(double a, double b, double x, double pref, double qref,
          bool* direct_is_p) {
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
  // 1. power series, either orientation -- with a near-one post-route: a
  // fired R1 orientation whose evaluated value
  // exceeds kBetaNearOne folds into R4's analytic assembly (same
  // orientation) instead, provided the fired first parameter is at or
  // below the kBetaPrTauMax zone ceiling.
  const bool r1n = x <= kBetaXi1 && b * x <= kBetaB1;
  const bool r1s = y <= kBetaXi1 && a * y <= kBetaB1;
  if (r1n || r1s) {
    const bool sw = !r1n && r1s;
    const double eval_v = sw ? qref : pref;  // I of the fired triple
    const double fired_alpha = sw ? b : a;
    if (eval_v > corvus::detail::kBetaNearOne &&
        fired_alpha <= corvus::detail::kBetaPrTauMax) {
      *direct_is_p = sw;  // val is the fired orientation's COMPLEMENT
      return kPr;
    }
    *direct_is_p = !sw;
    return kR1;
  }
  // 2. ridge ratio band, floor lowered to kBetaGlRidgeMin above the
  // gamma-limit threshold [(C) slice, ridge part].
  const double nu = 1.0 / (1.0 / a + 1.0 / b);
  const double rat1 = x * (1.0 + b / a);
  const double rat2 = y * (1.0 + a / b);
  const bool band = rat1 >= kBetaXiRatioLo && rat1 <= kBetaXiRatioHi &&
                    rat2 >= kBetaXiRatioLo && rat2 <= kBetaXiRatioHi;
  const bool gl_hi = bmax >= corvus::detail::kBetaGammaLim;
  if (band && (nu >= kBetaTRidge ||
               (gl_hi && nu >= corvus::detail::kBetaGlRidgeMin))) {
    // The kernel decides R3's direct side from lambda's EXACT sign, not from
    // the router's rounded ratio, so this mirrors that (see BetaR3Out). The
    // two can disagree only where lambda is within an ulp of zero, i.e. where
    // both sides are ~1/2 and the bin makes no difference to what is measured.
    const bool sw = rat1 > 1.0;
    const double lam = a * y - b * x;
    *direct_is_p = sw ? (lam > 0.0) : (lam >= 0.0);
    return kR3;
  }
  // 3. continued fraction, orientation by the pinned rule -- except above
  // kBetaGammaLim, where the lane takes the (C) gamma-limit slice and val
  // holds the NATURALLY COMPUTED gamma side (series -> P_gamma, CF ->
  // Q_gamma; see the slice block in src/beta-inl.h for the mapping).
  const double thr = 1.0 / (1.0 + (b + 1.0) / (a + 1.0));
  const bool sw = !(x < thr);
  // Both-huge lanes are EXCLUDED from the slice (kernel mirror): the
  // small/huge mapping needs exactly one parameter above the limit;
  // both-huge off-band lanes stay plain R2 and saturate via the PB
  // E-clamp.
  if (gl_hi && std::min(a, b) < corvus::detail::kBetaGammaLim) {
    const double ra = sw ? b : a;   // routed alpha
    const double rxi = sw ? y : x;  // routed xi
    const bool hf = ra >= corvus::detail::kBetaGammaLim;
    const double s = hf ? (sw ? a : b) : ra;  // small routed parameter
    const double huge = hf ? ra : (sw ? a : b);
    const double t = -huge * std::log(hf ? rxi : 1.0 - rxi);
    const bool ser = (s < 20.0 && t <= s + 1.0) || (s >= 20.0 && s >= 2.0 * t);
    const bool agree = (ser == !hf);  // val == routed value?
    *direct_is_p = agree ? !sw : sw;
    return kGl;
  }
  *direct_is_p = !sw;
  return kR2;
}

// #14 N4: masked-tail split. beta_p_reference.txt/beta_q_reference.txt are
// both 38100 rows -- a multiple of 4, so a single whole-set call never
// exercises the masked-tail path below AVX-512. Splitting into subspans
// [0, n-3) and [n-3, n) makes neither half a lane multiple regardless of
// what n is, so the tail path always runs; each call is independently
// masked and stateless, so this is row-for-row identical to one call.
template <typename Fn, typename... In>
void SplitCall(Fn fn, std::vector<double>& out, const In&... in) {
  const size_t n = out.size();
  const size_t split = n - 3;
  fn(std::span<const double>(in).subspan(0, split)...,
     std::span<double>(out).subspan(0, split));
  fn(std::span<const double>(in).subspan(split)...,
     std::span<double>(out).subspan(split));
}

bool LoadReference(const char* path, std::vector<double>* a,
                   std::vector<double>* b, std::vector<double>* x,
                   std::vector<double>* p, std::vector<double>* q) {
  const auto rows = corvus_test::LoadRef(path, 5, 10000);
  a->reserve(rows.size());
  b->reserve(rows.size());
  x->reserve(rows.size());
  p->reserve(rows.size());
  q->reserve(rows.size());
  for (const auto& row : rows) {
    a->push_back(corvus_test::ParseDouble(row.tok[0], path, row.line));
    b->push_back(corvus_test::ParseDouble(row.tok[1], path, row.line));
    x->push_back(corvus_test::ParseDouble(row.tok[2], path, row.line));
    p->push_back(corvus_test::ParseDouble(row.tok[3], path, row.line));
    q->push_back(corvus_test::ParseDouble(row.tok[4], path, row.line));
  }
  return true;
}

int ReportRegions(const char* label, const Region* r, int n, bool want_p) {
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
    // #14 N10: a bucket that never accumulated a row is a gate that never
    // actually ran. The one structural exception is the specials
    // complement cell: Route() hard-codes direct_is_p = true for kSp (see
    // its body), so exactly one of "specials"/"specials (-)" (indices 8
    // and 9 below) is unreachable in any single beta_p or beta_q run --
    // that is routing, not a coverage gap, and is the only cell skipped.
    const bool structurally_empty = (i == 8 && !want_p) || (i == 9 && want_p);
    if (r[i].n == 0 && !structurally_empty) {
      std::fprintf(stderr, "FAIL: %s %s bucket is empty (n=0) -- vacuous gate\n",
                   label, r[i].name);
      rc = 1;
    }
    if (r[i].n > 0 && r[i].max_ulp > r[i].bound) {
      std::fprintf(stderr, "FAIL: %s %s exceeds gate\n", label, r[i].name);
      rc = 1;
    }
  }
  return rc;
}

int Measure(const char* label, bool want_p, const std::vector<double>& a,
            const std::vector<double>& b, const std::vector<double>& x,
            const std::vector<double>& pref, const std::vector<double>& qref,
            const std::vector<double>& got, const std::vector<double>& want) {
  Region reg[14] = {
      {"R1 series dir", kGateR1Dir}, {"R1 series cmp", kGateR1Cmp},
      {"R2 cf     dir", kGateR2Dir}, {"R2 cf     cmp", kGateR2Cmp},
      {"R3 temme  dir", kGateR3Dir}, {"R3 temme  cmp", kGateR3Cmp},
      {"R4 tiny   dir", kGateR4Dir}, {"R4 tiny   cmp", kGateR4Cmp},
      {"specials     ", 0},          {"specials  (-)", 0},
      {"R4 postrt dir", kGatePrDir}, {"R4 postrt cmp", kGatePrCmp},
      {"R2 gammalim d", kGateGlDir}, {"R2 gammalim c", kGateGlCmp},
  };
  for (size_t i = 0; i < a.size(); ++i) {
    bool direct_is_p = false;
    const int code = Route(a[i], b[i], x[i], pref[i], qref[i], &direct_is_p);
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
    } else if (want[i] == 0.0) {
      // #14 N6: a reference of exactly zero (saturated rows carry an exact
      // 0/1 pair by construction) has no ULP neighbourhood -- UlpDiff maps
      // +0 and -0 to the same point (ulp_utils.h's policy), so a
      // signed-zero regression needs an exact-bits check instead.
      const bool ok = SameBits(got[i], want[i]);
      if (!ok) {
        std::fprintf(stderr,
                     "FAIL: %s %s signed-zero mismatch at a=%.17g b=%.17g "
                     "x=%.17g got=%.17g want=%.17g\n",
                     label, r.name, a[i], b[i], x[i], got[i], want[i]);
      }
      u = ok ? 0 : ~uint64_t{0};
    } else {
      u = UlpDiff(got[i], want[i]);
    }
    ++r.n;
    if (u > 0) ++r.miss;
    // Row-level dump for gate pinning: BETA_ULP_DUMP=<min ULP> prints
    // every not-correctly-rounded row at or above the threshold.
    // getenv deprecation suppressed locally, expect_target.h's pattern.
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 4996)
#endif
    static const char* dump_env = std::getenv("BETA_ULP_DUMP");
#if defined(_MSC_VER)
#pragma warning(pop)
#endif
    if (u > 0 && dump_env && u >= std::strtoull(dump_env, nullptr, 10)) {
      std::fprintf(stderr,
                   "DUMP %s %s ulp=%llu a=%.17g b=%.17g x=%.17g got=%a "
                   "want=%a\n",
                   label, r.name, static_cast<unsigned long long>(u), a[i],
                   b[i], x[i], got[i], want[i]);
    }
    if (u > r.max_ulp) {
      r.max_ulp = u;
      r.wa = a[i];
      r.wb = b[i];
      r.wx = x[i];
    }
  }
  return ReportRegions(label, reg, 14, want_p);
}

// ---- Post-pass (i): monotonicity in x over the reference set -------------
// P(a, b, x) is strictly increasing in x. The REFERENCE values must be
// non-decreasing with NO slack (the independent harness certified exact
// monotonicity over every (a, b) group; any violation here is a regen
// regression). The KERNEL values are each within their bucket's ULP bound
// of a monotone function, so adjacent values may legally dip by the sum of
// two bounds; a dip beyond kMonoSlackUlp is a seam discontinuity the
// pointwise gates cannot see (both sides individually in-budget).
constexpr uint64_t kMonoSlackUlp = 4;

int MonoPostPass(const std::vector<double>& a, const std::vector<double>& b,
                 const std::vector<double>& x, const std::vector<double>& pref,
                 const std::vector<double>& got) {
  std::map<std::pair<double, double>, std::vector<size_t>> groups;
  for (size_t i = 0; i < a.size(); ++i) {
    if (a[i] > 0.0 && b[i] > 0.0 && std::isfinite(a[i]) &&
        std::isfinite(b[i]) && x[i] > 0.0 && x[i] < 1.0) {
      groups[{a[i], b[i]}].push_back(i);
    }
  }
  size_t n_groups = 0, ref_bad = 0, ker_bad = 0;
  uint64_t worst = 0;
  double wa = 0.0, wb = 0.0, wx = 0.0;
  for (auto& kv : groups) {
    auto& idx = kv.second;
    if (idx.size() < 3) continue;
    ++n_groups;
    std::sort(idx.begin(), idx.end(),
              [&](size_t i, size_t j) { return x[i] < x[j]; });
    for (size_t k = 1; k < idx.size(); ++k) {
      const size_t i0 = idx[k - 1], i1 = idx[k];
      if (x[i1] == x[i0]) continue;
      if (pref[i1] < pref[i0]) {
        ++ref_bad;
        std::fprintf(stderr,
                     "MONO ref: a=%.17g b=%.17g x=%.17g->%.17g P %.17g->%.17g\n",
                     a[i0], b[i0], x[i0], x[i1], pref[i0], pref[i1]);
      }
      if (got[i1] < got[i0]) {
        const uint64_t u = UlpDiff(got[i1], got[i0]);
        if (u > worst) {
          worst = u;
          wa = a[i0];
          wb = b[i0];
          wx = x[i1];
        }
        if (u > kMonoSlackUlp) {
          ++ker_bad;
          std::fprintf(
              stderr,
              "MONO kernel: a=%.17g b=%.17g x=%.17g->%.17g P dips %llu ulp\n",
              a[i0], b[i0], x[i0], x[i1], static_cast<unsigned long long>(u));
        }
      }
    }
  }
  std::printf(
      "beta_p  monotonicity-in-x: %zu groups; ref violations=%zu; kernel "
      "dips > %llu ulp: %zu (worst dip %llu ulp at a=%.6g b=%.6g x=%.6g)\n",
      n_groups, ref_bad, static_cast<unsigned long long>(kMonoSlackUlp),
      ker_bad, static_cast<unsigned long long>(worst), wa, wb, wx);
  return (ref_bad || ker_bad) ? 1 : 0;
}

// ---- Post-pass (ii): dense seam-crossing sweeps ---------------------------
// One dense line per routing boundary, crossing the seam at fixed
// representative parameters. Monotonicity of the KERNEL output along the
// line (increasing in x and b, decreasing in a) is the continuity gate: a
// method-value mismatch at a seam shows as a wrong-direction step larger
// than the pointwise slack. Each line also reports the region sequence it
// actually crossed (via the router replica) so a constants drift that
// makes a sweep miss its seam is visible rather than silently green.
struct Seam {
  const char* name;
  int vary;  // 0 = x, 1 = a, 2 = b
  double fa, fb, fx;
  double lo, hi;
  bool geometric;
};

int SeamSweeps() {
  static const Seam kSeams[] = {
      {"R4->R2  bmax*xt=B1   [x]", 0, 0.01, 100.0, 0.0, 0.04, 0.12, false},
      {"R4->R1  tau=epsR4    [a]", 1, 0.0, 100.0, 0.05, 0.008, 0.03, true},
      {"R1->R2  b*x=B1       [x]", 0, 5.0, 30.0, 0.0, 0.20, 0.35, false},
      {"R2->R3->R2 band edges[x]", 0, 100.0, 100.0, 0.0, 0.15, 0.85, false},
      {"R2->R3  nu=TRidge    [a]", 1, 0.0, 100.0, 0.4, 40.0, 55.0, false},
      {"R1->Pr->R2 near-one  [x]", 0, 0.5, 100.0, 0.0, 0.03, 0.10, false},
      {"Pr->R1  bar in a     [a]", 1, 0.0, 20.0, 0.4, 0.2, 3.0, false},
      {"R2->Gl  bmax=2^59    [b]", 2, 0.5, 0.0, 1e-16, 0x1p58, 0x1p60, true},
      {"R1->Gl ser->Gl cf    [x]", 0, 10.0, 0x1p60, 0.0, 7.0 / 0x1p60,
       14.0 / 0x1p60, false},
      {"Gl->R3  nu=GlMin     [a]", 1, 0.0, 0x1p60, 1.735e-17, 15.0, 26.0,
       false},
  };
  constexpr int kN = 4001;
  int rc = 0;
  for (const Seam& s : kSeams) {
    std::vector<double> a(kN), b(kN), x(kN), p(kN);
    for (int i = 0; i < kN; ++i) {
      const double t = static_cast<double>(i) / (kN - 1);
      const double v = s.geometric
                           ? s.lo * std::exp(t * std::log(s.hi / s.lo))
                           : s.lo + t * (s.hi - s.lo);
      a[i] = s.vary == 1 ? v : s.fa;
      b[i] = s.vary == 2 ? v : s.fb;
      x[i] = s.vary == 0 ? v : s.fx;
    }
    corvus::beta_p(a, b, x, p);
    // Direction: P increases in x and b, decreases in a.
    const bool increasing = s.vary != 1;
    uint64_t worst = 0;
    size_t bad = 0;
    double wv = 0.0;
    for (int i = 1; i < kN; ++i) {
      const bool dip = increasing ? (p[i] < p[i - 1]) : (p[i] > p[i - 1]);
      if (!dip) continue;
      const uint64_t u = UlpDiff(p[i], p[i - 1]);
      if (u > worst) {
        worst = u;
        wv = s.vary == 0 ? x[i] : (s.vary == 1 ? a[i] : b[i]);
      }
      if (u > kMonoSlackUlp) ++bad;
    }
    // Region sequence actually crossed (router replica; complement from the
    // kernel value is more than accurate enough for the near-one bar).
    char seq[64];
    int sn = 0, last = -1;
    seq[0] = '\0';
    for (int i = 0; i < kN && sn < 56; ++i) {
      bool dp;
      const int code = Route(a[i], b[i], x[i], p[i], 1.0 - p[i], &dp);
      if (code != last) {
        static const char* names[] = {"R1", "R2", "R3", "R4",
                                      "Sp", "Pr", "Gl"};
        sn += std::snprintf(seq + sn, sizeof(seq) - sn, "%s%s",
                            last < 0 ? "" : ">", names[code]);
        last = code;
      }
    }
    std::printf(
        "beta_p  seam %-24s crossed %-14s wrong-dir > %llu ulp: %zu (worst "
        "%llu at %.9g)\n",
        s.name, seq, static_cast<unsigned long long>(kMonoSlackUlp), bad,
        static_cast<unsigned long long>(worst), wv);
    if (bad) rc = 1;
    if (last == -1 || std::strchr(seq, '>') == nullptr) {
      std::fprintf(stderr, "FAIL: seam sweep '%s' crossed no boundary (%s)\n",
                   s.name, seq);
      rc = 1;
    }
  }
  return rc;
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
    SplitCall(corvus::beta_p, got, a, b, x);
    rc |= Measure("beta_p", true, a, b, x, p, q, got, p);
    rc |= MonoPostPass(a, b, x, p, got);
    rc |= SeamSweeps();
  }
  {
    std::vector<double> a, b, x, p, q;
    if (!LoadReference(q_path, &a, &b, &x, &p, &q)) return 2;
    std::vector<double> got(a.size());
    SplitCall(corvus::beta_q, got, a, b, x);
    rc |= Measure("beta_q", false, a, b, x, p, q, got, q);
  }

  if (rc == 0) std::printf("PASS: all regions within G4 pinned gates\n");
  return rc;
}
