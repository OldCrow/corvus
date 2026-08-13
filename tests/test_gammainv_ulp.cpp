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
#include <cstring>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

#include "corvus/corvus.h"
#include "expect_target.h"
#include "src/gamma_data.h"      // kGammaAT: the S1 / small-a seed split
#include "src/gammainv_data.h"   // the pinned deep-small cut and thresholds

namespace {

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

int64_t OrderedBits(double x) {
  int64_t b;
  std::memcpy(&b, &x, sizeof(b));
  return b < 0 ? (INT64_MIN - b) : b;
}

uint64_t UlpDiff(double a, double b) {
  if (std::isnan(a) || std::isnan(b)) return UINT64_MAX;
  if (std::isinf(a) || std::isinf(b)) return a == b ? 0 : UINT64_MAX;
  return static_cast<uint64_t>(std::llabs(OrderedBits(a) - OrderedBits(b)));
}

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
  std::ifstream f(path);
  if (!f) {
    std::fprintf(stderr, "cannot open reference file: %s\n", path);
    return false;
  }
  std::string sa, ss, sx;
  while (f >> sa >> ss >> sx) {
    a->push_back(std::strtod(sa.c_str(), nullptr));
    s->push_back(std::strtod(ss.c_str(), nullptr));
    x->push_back(std::strtod(sx.c_str(), nullptr));
  }
  if (a->size() < 5000) {
    std::fprintf(stderr, "reference file suspiciously small: %zu lines\n",
                 a->size());
    return false;
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

int Measure(const char* label, bool want_q, const std::vector<double>& a,
            const std::vector<double>& s, const std::vector<double>& want,
            const std::vector<double>& got) {
  Bucket b[5] = {
      {"deep-small closed form"}, {"beyond-resolution"}, {"S1 (a_T..3e34)"},
      {"small-a tail seed"},      {"small-a mid seed"},
  };
  Bucket cross[3] = {
      {"  ... s > 1/2 (flip)"}, {"  ... subnormal x"}, {"  ... x = 0"}};

  for (size_t i = 0; i < a.size(); ++i) {
    const uint64_t u = UlpDiff(got[i], want[i]);
    Accumulate(b[Bucketize(a[i], s[i], want[i])], a[i], s[i], u);
    if (s[i] > 0.5) Accumulate(cross[0], a[i], s[i], u);
    if (want[i] > 0.0 && want[i] < (std::numeric_limits<double>::min)()) {
      Accumulate(cross[1], a[i], s[i], u);
    }
    if (want[i] == 0.0) Accumulate(cross[2], a[i], s[i], u);
  }

  int rc = 0;
  std::printf("--- %s (%s side) ---\n", label, want_q ? "q" : "p");
  for (const Bucket& r : b) {
    Report(r, true);
    if (r.n > 0 && r.max_ulp > kMaxUlp) {
      std::fprintf(stderr, "FAIL: %s %s exceeds gate\n", label, r.name);
      rc = 1;
    }
  }
  for (const Bucket& r : cross) Report(r, false);
  return rc;
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
    corvus::gamma_p_inv(a, s, got);
    rc |= Measure("gamma_p_inv", false, a, s, x, got);
  }
  {
    std::vector<double> a, s, x;
    if (!LoadReference(q_path, &a, &s, &x)) return 2;
    std::vector<double> got(a.size());
    corvus::gamma_q_inv(a, s, got);
    rc |= Measure("gamma_q_inv", true, a, s, x, got);
  }

  if (rc == 0) std::printf("PASS: all buckets within gates\n");
  return rc;
}
