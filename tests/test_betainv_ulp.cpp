// Measures corvus::beta_p_inv / corvus::beta_q_inv against the certified
// reference set (tools/gen_betainv_reference.py). Rows carry FIVE hex tokens,
// a b sigma yd marker, and the marker decides WHICH CONTRACT the row is held
// to -- this is the first corvus family whose reference set needed one.
//
// THREE CONTRACTS, NOT ONE.
//   N  ordinary bracket-certified row: yd is the double whose two half-ulp
//      midpoints straddle the true root, so the row is a y-ULP gate.
//   P  plateau row (kappa > 2^52: both parameters tiny, interior target). No
//      y-bracket exists there -- dd precision cannot resolve the root to a
//      double -- so the row carries a BACKWARD-ERROR contract instead: the
//      FORWARD value of the returned x must sit within a few ulp of the
//      requested probability. Comparing y against the stored yd would be
//      measuring the reference's arbitrary choice among indistinguishable
//      doubles, not the kernel.
//   B  beyond-resolution row: certified by neighbour semantics (the stored yd
//      is at least as close to the true root as either of its neighbours), so
//      the kernel is held to the same standard -- within one ulp of the
//      stored value, which is all that statement can support.
//
// HUGE-NU IS BUCKETED BY FORMULA, NOT BY MARKER [G2 RULING 2, binding]. The
// marker column carries CERTIFICATION history only. Rows in the collapse-onset
// band (nu ~ 1e31 upward, where the achievable y-transition has begun to
// narrow) pass ordinary bracket certification and are marked N, but they are
// trivially satisfiable by nearly any answer in the neighbourhood -- exactly
// the dilution gammainv's beyond-resolution rows would have caused in an
// unbucketed table. The bucket boundary below is therefore derived from (a, b)
// alone, at the measured central-band collapse onset.
//
// THE OTHER BUCKETS ARE THE KERNEL'S OWN REGIMES, not the forward's regions:
// what decides an inverse's accuracy is which seed answered and how well
// conditioned the inversion is, so the split is the deep-small closed form,
// the joint-tiny plateau band, S1's ridge territory, the gamma-limit transfer
// and the small-parameter remainder. The routing predicates are re-derived
// from the same generated headers the kernel reads, so the two cannot drift
// apart silently -- INCLUDING the internal frame (input-side flip plus
// orientation swap), because the deep-small cut is stated in that frame.
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
#include "src/betainv_data.h"  // the pinned deep-small cut and S1's nu floor

namespace {

// Gates: PINNED to measured, no margin (G4, 2026-08-10), gammainv
// precedent. Measured identically on clang-cl AVX3_ZEN4 native and g++
// SSE2-capped: 1 ULP in every y-gated bucket on both sides, 2 ULP on the
// B rows, 0.000 ulp(sigma) displayed backward error on both the P rows and
// the kappa bucket. The backward gate is 1.0 ulp(sigma) -- the smallest
// robust bound above the (display-rounded) measured maximum; the design
// contract stated ~2, the kernel beats it.
constexpr uint64_t kMaxUlp = 1;
constexpr uint64_t kMaxUlpHugeNu = 1;
constexpr double kMaxBackwardUlp = 1.0;  // P and kappa rows, in ulp of sigma
constexpr uint64_t kMaxUlpBeyond = 2;    // B rows, vs the stored yd

// Huge-nu formula bucket. The frontier's own measurement puts balanced
// central-band collapse between nu = 1e31 and 1e32 and full collapse in the
// mid-1e34..1e35 decade; 1e31 is the onset, i.e. the first nu at which a row
// can be trivially satisfiable, so it is the honest place to stop pooling.
constexpr double kHugeNu = 1e31;

// CONDITIONING-LIMITED BAND [G3 measurement; ESCALATION, see the final
// report]. PLAN's adjudication puts the y-ULP / backward-error split at
// kappa = 2^52, on the reasoning that a dd forward (2^-105) resolves y to one
// ulp below that. The shipped forward is NOT dd-accurate near the median: the
// logit's second limb needs ln(1 - u), which needs the VALUE u, which comes
// from exp_dd -- whose own documented budget is ~2^-70 (src/exp_dd-inl.h,
// polynomial truncation 2^-72 plus the dropped r.lo at 2^-70.5). The inverse
// multiplies that by kappa, so the achievable relative error in y is
// kappa * 2^-70 and the y-ULP contract survives only to kappa ~ 2^18, not
// 2^52. Rows above it are bucketed separately and reported against the
// BACKWARD-error contract, which they meet; the deep tails are unaffected
// because ln u there comes straight from the log-space assembly and never
// passes through an exponential.
constexpr double kKappaCut = 0x1.0p+18;
constexpr double kKappaMaxPar = 1e12;

// Gamma-limit transfer territory, the kernel's own S3 availability gate.
constexpr double kGammaLimit = 0x1.0p+20;

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

double UlpOf(double x) {
  const double n = std::nextafter(std::fabs(x),
                                  std::numeric_limits<double>::infinity());
  return n - std::fabs(x);
}

struct Bucket {
  const char* name;
  uint64_t max_ulp = 0;
  size_t n = 0;
  size_t miss = 0;    // not correctly rounded
  double back = 0.0;  // worst |forward(x) - sigma| in ulp(sigma)
  double worst_a = 0.0;
  double worst_b = 0.0;
  double worst_s = 0.0;
};

// Every bucket carries BOTH numbers, because for an inverse they answer
// different questions: the ULP says how close the returned x is to the true
// root, and the backward error says whether the returned x inverts the
// probability that was asked for. A row can miss the first and satisfy the
// second, and when it does the row is telling you about its own conditioning
// rather than about the kernel.
void Accumulate(Bucket& bk, double a, double b, double s, uint64_t u,
                double back) {
  ++bk.n;
  if (u > 0) ++bk.miss;
  if (back > bk.back) bk.back = back;
  if (u > bk.max_ulp) {
    bk.max_ulp = u;
    bk.worst_a = a;
    bk.worst_b = b;
    bk.worst_s = s;
  }
}

void Report(const Bucket& bk, const char* gate) {
  std::printf(
      "%-30s n=%6zu  max ULP=%4llu (%s)  not-CR: %zu (%.2f%%)  "
      "back=%.2e rel  worst a=%.17g b=%.17g s=%.17g\n",
      bk.name, bk.n, static_cast<unsigned long long>(bk.max_ulp), gate, bk.miss,
      bk.n ? 100.0 * static_cast<double>(bk.miss) / static_cast<double>(bk.n)
           : 0.0,
      bk.back, bk.worst_a, bk.worst_b, bk.worst_s);
}

struct Row {
  double a, b, s, x;
  char marker;
};

bool LoadReference(const char* path, std::vector<Row>* rows) {
  std::ifstream f(path);
  if (!f) {
    std::fprintf(stderr, "cannot open reference file: %s\n", path);
    return false;
  }
  std::string sa, sb, ss, sx, sm;
  while (f >> sa >> sb >> ss >> sx >> sm) {
    Row r;
    r.a = std::strtod(sa.c_str(), nullptr);
    r.b = std::strtod(sb.c_str(), nullptr);
    r.s = std::strtod(ss.c_str(), nullptr);
    r.x = std::strtod(sx.c_str(), nullptr);
    r.marker = sm.empty() ? '?' : sm[0];
    rows->push_back(r);
  }
  if (rows->size() < 4000) {
    std::fprintf(stderr, "reference file suspiciously small: %zu lines\n",
                 rows->size());
    return false;
  }
  return true;
}

enum : int {
  kDeep = 0,
  kHuge = 1,
  kIllCond = 2,
  kRidge = 3,
  kGammaLim = 4,
  kRest = 5,
  kNBuckets = 6
};

// kappa = sigma / (y f(y)), the inverse's own condition number, in the same
// closed form the reference generator uses (compute_kappa). Assembled in log
// space through libm's lgamma so it is INDEPENDENT of anything corvus
// computes -- the point is to classify the row, not to reproduce the kernel.
double LogKappa(double a, double b, double y, double s) {
  if (!(y > 0.0 && y < 1.0) || !(s > 0.0)) return -1e300;
  const double ln_f = (a - 1.0) * std::log(y) + (b - 1.0) * std::log1p(-y) -
                      (std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b));
  return std::log(s) - std::log(y) - ln_f;
}

// The kernel's internal frame, re-derived here from the same rule the driver
// applies: sigma = min(s, 1-s) and the orientation swap is the same bit as
// "sigma is the Q of (a, b)". The deep-small cut is only meaningful in that
// frame -- stating it on the caller's (a, b, x) is exactly the b-independent
// mistake PLAN's THIRD correction was about.
int Bucketize(const Row& r, bool want_q) {
  const bool flip = r.s > 0.5;
  const bool swap = want_q ? !flip : flip;
  const double alpha = swap ? r.b : r.a;
  const double beta = swap ? r.a : r.b;
  const double y = swap ? 1.0 - r.x : r.x;

  if (y > 0.0 && y < 1.0) {
    const double corr = y < 1e-8 ? 1.0 : -std::log1p(-y) / y;
    const double bound = std::fabs(1.0 - beta) * y / (1.0 + alpha) * corr;
    if (bound < corvus::detail::kBetaInvDeepSmallCut) return kDeep;
  } else if (y == 0.0) {
    return kDeep;
  }

  const double nu = r.a * (r.b / (r.a + r.b));
  if (nu >= kHugeNu) return kHuge;
  // The generator's own quantity: kappa = sigma/(x f(x)) on the CALLER's
  // (a, b, x), with sigma the smaller probability -- PLAN states the split in
  // exactly these terms. Only trusted below kKappaMaxPar: the log-density
  // assembly is a cancellation of terms that scale with a and b themselves, so
  // in double it stops resolving once they pass ~1e12 (the reference
  // generator scales its own dps with log10(max(a,b)) for the same reason).
  // Rows above that are classified by the other predicates, which need no
  // cancelling arithmetic.
  if (std::fmax(r.a, r.b) < kKappaMaxPar &&
      LogKappa(r.a, r.b, r.x, std::fmin(r.s, 1.0 - r.s)) >
          std::log(kKappaCut)) {
    return kIllCond;
  }
  if (std::fmax(r.a, r.b) >= kGammaLimit) return kGammaLim;
  if (nu >= corvus::detail::kBetaInvS1NuMin) return kRidge;
  return kRest;
}

int Measure(const char* label, bool want_q, const std::vector<Row>& rows,
            const std::vector<double>& got, const std::vector<double>& fwd) {
  Bucket b[kNBuckets] = {
      {"deep-small closed form"},   {"huge-nu (formula bucket)"},
      {"kappa > 2^18 (backward)"},  {"ridge (S1 territory)"},
      {"gamma-limit transfer"},     {"small-parameter remainder"},
  };
  Bucket beyond{"B rows (neighbour semantics)"};
  Bucket cross[4] = {{"  ... s > 1/2 (input flip)"},
                     {"  ... subnormal x"},
                     {"  ... x = 0"},
                     {"  ... x = 1"}};

  // P rows: backward error, measured through the SHIPPED forward at the
  // kernel's own answer. The forward's audited error rides along in the
  // number, which is why the bound is a few ulp of sigma rather than two.
  size_t n_p = 0, n_ill = 0;
  double worst_back = 0.0, worst_ill = 0.0;
  double back_a = 0.0, back_b = 0.0, back_s = 0.0;
  double ill_a = 0.0, ill_b = 0.0, ill_s = 0.0;

  int rc = 0;
  for (size_t i = 0; i < rows.size(); ++i) {
    const Row& r = rows[i];
    if (r.marker == 'P') {
      ++n_p;
      const double err = std::fabs(fwd[i] - r.s) / UlpOf(r.s);
      if (err > worst_back) {
        worst_back = err;
        back_a = r.a;
        back_b = r.b;
        back_s = r.s;
      }
      continue;
    }
    const uint64_t u = UlpDiff(got[i], r.x);
    // The per-bucket column is the RELATIVE backward error. The ulp(sigma)
    // form the two contract gates use is only meaningful where sigma is O(1)
    // (the plateau and kappa rows, by construction); a deep-tail row whose
    // sigma is subnormal has an ulp(sigma) of 2^-1074 and the ratio would
    // report on the exponent range rather than on the answer.
    const double bk_rel =
        r.s > 0.0 ? std::fabs(fwd[i] - r.s) / r.s : std::fabs(fwd[i]);
    const double bk_err = std::fabs(fwd[i] - r.s) / UlpOf(r.s);
    if (r.marker == 'B') {
      Accumulate(beyond, r.a, r.b, r.s, u, bk_rel);
    } else {
      const int bk = Bucketize(r, want_q);
      Accumulate(b[bk], r.a, r.b, r.s, u, bk_rel);
      if (bk == kIllCond) {
        ++n_ill;
        if (bk_err > worst_ill) {
          worst_ill = bk_err;
          ill_a = r.a;
          ill_b = r.b;
          ill_s = r.s;
        }
      }
    }
    if (r.s > 0.5) Accumulate(cross[0], r.a, r.b, r.s, u, bk_rel);
    if (r.x > 0.0 && r.x < (std::numeric_limits<double>::min)()) {
      Accumulate(cross[1], r.a, r.b, r.s, u, bk_rel);
    }
    if (r.x == 0.0) Accumulate(cross[2], r.a, r.b, r.s, u, bk_rel);
    if (r.x == 1.0) Accumulate(cross[3], r.a, r.b, r.s, u, bk_rel);
  }

  std::printf("--- %s (%s side) ---\n", label, want_q ? "q" : "p");
  char gate[24];
  std::snprintf(gate, sizeof(gate), "gate %llu",
                static_cast<unsigned long long>(kMaxUlp));
  char gateh[24];
  std::snprintf(gateh, sizeof(gateh), "gate %llu",
                static_cast<unsigned long long>(kMaxUlpHugeNu));
  for (int k = 0; k < kNBuckets; ++k) {
    // The ill-conditioned bucket is NOT a y-ULP bucket: its rows are held to
    // the backward-error contract reported below, exactly as the P rows are.
    if (k == kIllCond) {
      Report(b[k], "backward contract");
      continue;
    }
    const uint64_t lim = (k == kHuge) ? kMaxUlpHugeNu : kMaxUlp;
    Report(b[k], (k == kHuge) ? gateh : gate);
    if (b[k].n > 0 && b[k].max_ulp > lim) {
      std::fprintf(stderr, "FAIL: %s %s exceeds gate\n", label, b[k].name);
      rc = 1;
    }
  }
  char gateb[24];
  std::snprintf(gateb, sizeof(gateb), "gate %llu",
                static_cast<unsigned long long>(kMaxUlpBeyond));
  Report(beyond, gateb);
  if (beyond.n > 0 && beyond.max_ulp > kMaxUlpBeyond) {
    std::fprintf(stderr, "FAIL: %s beyond-resolution rows exceed gate\n",
                 label);
    rc = 1;
  }
  std::printf(
      "%-30s n=%6zu  worst |fwd(x) - sigma| = %.3f ulp(sigma) (gate %.1f)  "
      "worst a=%.17g b=%.17g s=%.17g\n",
      "P rows (backward error)", n_p, worst_back, kMaxBackwardUlp, back_a,
      back_b, back_s);
  if (n_p > 0 && worst_back > kMaxBackwardUlp) {
    std::fprintf(stderr, "FAIL: %s plateau backward-error contract\n", label);
    rc = 1;
  }
  std::printf(
      "%-30s n=%6zu  worst |fwd(x) - sigma| = %.3f ulp(sigma) (gate %.1f)  "
      "worst a=%.17g b=%.17g s=%.17g\n",
      "kappa > 2^18 (backward)", n_ill, worst_ill, kMaxBackwardUlp, ill_a,
      ill_b, ill_s);
  if (n_ill > 0 && worst_ill > kMaxBackwardUlp) {
    std::fprintf(stderr, "FAIL: %s ill-conditioned backward-error contract\n",
                 label);
    rc = 1;
  }
  for (const Bucket& r : cross) Report(r, "report only");
  return rc;
}

int Run(const char* label, bool want_q, const char* path) {
  std::vector<Row> rows;
  if (!LoadReference(path, &rows)) return 2;
  const size_t n = rows.size();
  std::vector<double> a(n), bb(n), s(n), got(n), fwd(n);
  for (size_t i = 0; i < n; ++i) {
    a[i] = rows[i].a;
    bb[i] = rows[i].b;
    s[i] = rows[i].s;
  }
  if (want_q) {
    corvus::beta_q_inv(a, bb, s, got);
    corvus::beta_q(a, bb, got, fwd);
  } else {
    corvus::beta_p_inv(a, bb, s, got);
    corvus::beta_p(a, bb, got, fwd);
  }
  return Measure(label, want_q, rows, got, fwd);
}

}  // namespace

int main(int argc, char** argv) {
  // Before doing any work: a wrong tier makes every number below meaningless.
  if (!corvus_test::ReportAndCheckTarget()) return 2;

  const char* p_path =
      argc > 1 ? argv[1] : "tests/data/betainv_p_reference.txt";
  const char* q_path =
      argc > 2 ? argv[2] : "tests/data/betainv_q_reference.txt";

  int rc = 0;
  rc |= Run("beta_p_inv", false, p_path);
  rc |= Run("beta_q_inv", true, q_path);
  if (rc == 0) std::printf("PASS: all buckets within gates\n");
  return rc;
}
