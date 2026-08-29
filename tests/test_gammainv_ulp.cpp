// Measures max ULP deviation of corvus::gamma_p_inv / corvus::gamma_q_inv
// against the bracket-certified reference set (tools/gen_gammainv_reference.py:
// every row's x is the double whose two half-ulp midpoints straddle the true
// root, certified at layered dps 60/100).
//
// THE BUCKETS ARE THE KERNEL'S OWN REGIMES, not the forward's regions. What
// decides the accuracy of an inverse is which SEED answered and how well
// conditioned the inversion is there, so the split is: the deep-small closed
// form, the two small-a seed territories either side of the shallow
// threshold, S1's Temme-quantile territory, and -- separately, which is the
// whole point -- the beyond-resolution rows.
//
// BEYOND-RESOLUTION MUST NOT BE POOLED WITH ANYTHING. For a
// above ~3e34 one ulp of x already moves a*phi past 800, so the entire
// P = 0 -> P = 1 transition happens inside a single ulp and the certified
// answer is x = a for every non-degenerate target. Those rows are correct and
// trivially so; they are about half the random grid, and pooling them would
// dilute every real-region statistic by a factor of two.
//
// The two cross-cuts at the end are report-only and deliberately overlap the
// partition: large-side inputs (s > 1/2) exercise the exact input flip and are
// spread across every bucket, and the subnormal/zero answers are the ones the
// deep-small branch's single-rounding scaling exists for.
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <span>
#include <string>
#include <vector>

#include "corvus/corvus.h"
#include "expect_target.h"
#include "src/gamma_data.h"      // kGammaAT: the S1 / small-a seed split
#include "src/gammainv_data.h"   // the pinned deep-small cut and thresholds
#include "ulp_utils.h"

namespace {

using corvus_test::LoadRef;
using corvus_test::ParseDouble;
using corvus_test::SameBits;
using corvus_test::UlpDiff;

// Gate PINNED to measured, no margin. Identical cells on
// every validated leg -- clang-cl AVX3_ZEN4 native, g++ SSE2-capped, MSVC
// AVX2, and the capped clang-cl sweep -- including not-CR counts and
// worst-case points: 1 ULP max in every bucket, with the deep-small
// closed-form, subnormal-x and x=0 bands correctly rounded (max 0).
constexpr uint64_t kMaxUlp = 1;

// Beyond-resolution floor. Not a kernel constant -- the kernel has no such
// branch, by design -- but a property of the domain: a*phi(1 +- ulp/a) > 800
// once a exceeds ~3e34, which is where the certified x collapses onto a.
constexpr double kBeyondResolutionA = 3e34;

struct Bucket {
  const char* name;
  uint64_t max_ulp = 0;
  size_t n = 0;
  size_t miss = 0;  // not correctly rounded
  double worst_a = 0.0;
  double worst_s = 0.0;
};

void Accumulate(Bucket& b, double a, double s, uint64_t u) {
  ++b.n;
  if (u > 0) ++b.miss;
  if (u > b.max_ulp) {
    b.max_ulp = u;
    b.worst_a = a;
    b.worst_s = s;
  }
}

void Report(const Bucket& b, bool gated) {
  char gate[24];
  if (gated) {
    std::snprintf(gate, sizeof(gate), "gate %llu",
                  static_cast<unsigned long long>(kMaxUlp));
  } else {
    std::snprintf(gate, sizeof(gate), "report only");
  }
  std::printf(
      "%-26s n=%6zu  max ULP=%3llu (%s)  not-CR: %zu (%.2f%%)  "
      "worst a=%.17g s=%.17g\n",
      b.name, b.n, static_cast<unsigned long long>(b.max_ulp), gate, b.miss,
      b.n ? 100.0 * static_cast<double>(b.miss) / static_cast<double>(b.n)
          : 0.0,
      b.worst_a, b.worst_s);
}

// Rows are three hex doubles: a, s, x.
bool LoadReference(const char* path, std::vector<double>* a,
                   std::vector<double>* s, std::vector<double>* x) {
  const auto rows = LoadRef(path, 3, 5000);
  for (const auto& row : rows) {
    a->push_back(ParseDouble(row.tok[0], path, row.line));
    s->push_back(ParseDouble(row.tok[1], path, row.line));
    x->push_back(ParseDouble(row.tok[2], path, row.line));
  }
  return true;
}

enum : int { kDeep = 0, kBeyond = 1, kS1 = 2, kSmallTail = 3, kSmallMid = 4 };

// The routing mirrors the kernel's own thresholds, read from the same
// generated header, so the two cannot drift apart silently.
int Bucketize(double a, double s, double x) {
  // The kernel cuts on x0*(1+a) against the pinned constant; x0 and the
  // certified x agree to well within the factor that decides this.
  if (x * (1.0 + a) < corvus::detail::kGammaInvDeepSmallCut) return kDeep;
  if (a >= kBeyondResolutionA) return kBeyond;
  if (a >= corvus::detail::kGammaAT) return kS1;
  const double t = s > 0.5 ? 1.0 - s : s;  // the side actually solved
  return t < corvus::detail::kGammaInvShallowThreshold ? kSmallTail : kSmallMid;
}

// N4: split every whole-reference-set call into [0, n-3) and [n-3, n) so the
// trailing 3-row group always runs a masked tail, independent of whether n
// itself happens to be a lane multiple on the tier under test. gamma's
// 16734 rows are even (no tail on any 2-lane tier today) and gammainv-q's
// 6520 rows are a multiple of 8 (no tail on any tier at all) -- exactly the
// coverage gap this closes. Downstream per-row loops are untouched: they
// only see the fully populated `got` vector.
void CallSplit(void (*fn)(std::span<const double>, std::span<const double>,
                          std::span<double>),
               const std::vector<double>& a, const std::vector<double>& s,
               std::vector<double>& out) {
  const size_t n = a.size();
  const size_t split = n - 3;
  fn(std::span<const double>(a).first(split),
     std::span<const double>(s).first(split),
     std::span<double>(out).first(split));
  fn(std::span<const double>(a).subspan(split),
     std::span<const double>(s).subspan(split),
     std::span<double>(out).subspan(split));
}

int Measure(const char* label, bool want_q, const std::vector<double>& a,
            const std::vector<double>& s, const std::vector<double>& want,
            const std::vector<double>& got) {
  Bucket b[5] = {
      {"deep-small closed form"}, {"beyond-resolution"}, {"S1 (a_T..3e34)"},
      {"small-a tail seed"},      {"small-a mid seed"},
  };
  Bucket cross[3] = {
      {"  ... s > 1/2 (flip)"}, {"  ... subnormal x"}, {"  ... x = 0"}};

  int signed_zero_fail = 0;
  for (size_t i = 0; i < a.size(); ++i) {
    uint64_t u;
    // N6: a reference value of exactly +/-0 has no ULP neighbourhood --
    // UlpDiff maps both zeros to the same point (ulp_utils.h's policy
    // comment), so it cannot see a sign regression. Check bit-exactness
    // instead and fail with a dedicated message on any mismatch; these are
    // exactly the rows the "x = 0" cross-cut below also counts.
    if (want[i] == 0.0) {
      if (!SameBits(got[i], want[i])) {
        std::fprintf(stderr,
                     "FAIL: %s signed-zero mismatch at a=%.17g s=%.17g "
                     "got=%a want=%a\n",
                     label, a[i], s[i], got[i], want[i]);
        signed_zero_fail = 1;
      }
      u = 0;
    } else {
      u = UlpDiff(got[i], want[i]);
    }
    Accumulate(b[Bucketize(a[i], s[i], want[i])], a[i], s[i], u);
    if (s[i] > 0.5) Accumulate(cross[0], a[i], s[i], u);
    if (want[i] > 0.0 && want[i] < (std::numeric_limits<double>::min)()) {
      Accumulate(cross[1], a[i], s[i], u);
    }
    if (want[i] == 0.0) Accumulate(cross[2], a[i], s[i], u);
  }

  int rc = 0;
  std::printf("--- %s (%s side) ---\n", label, want_q ? "q" : "p");
  std::printf("%-8s %zu exact-zero reference rows (checked via SameBits)\n",
              label, cross[2].n);
  for (const Bucket& r : b) {
    Report(r, true);
    // N10: a named bucket that never received a row passes its gate
    // vacuously -- a routing or reference-set regression, not a clean run.
    if (r.n == 0) {
      std::fprintf(stderr, "FAIL: %s %s bucket is empty (n=0)\n", label,
                   r.name);
      std::exit(1);
    }
    if (r.n > 0 && r.max_ulp > kMaxUlp) {
      std::fprintf(stderr, "FAIL: %s %s exceeds gate\n", label, r.name);
      rc = 1;
    }
  }
  for (const Bucket& r : cross) {
    Report(r, false);
    if (r.n == 0) {
      std::fprintf(stderr, "FAIL: %s %s bucket is empty (n=0)\n", label,
                   r.name);
      std::exit(1);
    }
  }
  return rc | signed_zero_fail;
}

}  // namespace

int main(int argc, char** argv) {
  // Before doing any work: a wrong tier makes every number below meaningless.
  if (!corvus_test::ReportAndCheckTarget()) return 2;

  const char* p_path =
      argc > 1 ? argv[1] : "tests/data/gammainv_p_reference.txt";
  const char* q_path =
      argc > 2 ? argv[2] : "tests/data/gammainv_q_reference.txt";

  int rc = 0;
  {
    std::vector<double> a, s, x;
    if (!LoadReference(p_path, &a, &s, &x)) return 2;
    std::vector<double> got(a.size());
    CallSplit(corvus::gamma_p_inv, a, s, got);
    rc |= Measure("gamma_p_inv", false, a, s, x, got);
  }
  {
    std::vector<double> a, s, x;
    if (!LoadReference(q_path, &a, &s, &x)) return 2;
    std::vector<double> got(a.size());
    CallSplit(corvus::gamma_q_inv, a, s, got);
    rc |= Measure("gamma_q_inv", true, a, s, x, got);
  }

  if (rc == 0) std::printf("PASS: all buckets within gates\n");
  return rc;
}
