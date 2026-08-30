// Measures max ULP deviation of corvus::cos/sin against the mpmath-generated
// correctly-rounded references (tests/data/{cos,sin}_reference.txt), split by
// the kernel's own regions: small (|x| <= kTrigDMax = 2^23, the exact-split
// reduction) and huge (|x| > 2^23, the Payne-Hanek reduction -- its rows
// include every exponent's CF worst reduction cancellation, so this gate is
// the end-to-end certification of the PH design), each split by sign (sin is
// odd by construction -- the |x|-symmetric kernel makes sin(-x) the exact
// negation of sin(x) -- so a sign-carrying regression names itself).
//
// The two files share one input set row-for-row (generator guarantee), so
// each is gated independently with the same region structure.
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "corvus/corvus.h"
#include "expect_target.h"
#include "src/trig_data.h"
#include "ulp_utils.h"

namespace {

using corvus_test::SameBits;
using corvus_test::UlpDiff;

using Fn = void (*)(std::span<const double>, std::span<double>);

// #14 N4: masked-tail split -- neither subspan is a lane multiple for any
// tier, so the masked LoadN/StoreN path always runs.
void SplitCall(Fn fn, const std::vector<double>& in, std::vector<double>& out) {
  const size_t n = in.size();
  const size_t split = n - 3;
  fn(std::span<const double>(in).subspan(0, split),
     std::span<double>(out).subspan(0, split));
  fn(std::span<const double>(in).subspan(split),
     std::span<double>(out).subspan(split));
}

// Gate PINNED to measured, no margin (docs/ACCURACY.md): 1 ULP in every
// region on every tier class. Measured 2026-08-30 on AVX3_ZEN4 (FMA) AND
// on the SSE2 cap (no FMA) -- the small region's reduction exactness is
// structural (30-bit split parts), not FMA-dependent, so unlike the donor
// (whose SSE2 tier gates at 2) the non-FMA tiers hold the same 1-ULP
// ceiling here. Huge region: the PH reduction carries r to ~2^-93 relative
// even at the deepest cancellation, so the cores' error dominates and the
// region pins at the same ceiling.
constexpr uint64_t kMaxUlp = 1;

struct Region {
  const char* name;
  uint64_t gate;
  uint64_t max_ulp = 0;
  size_t n = 0;
  size_t miss = 0;  // not correctly rounded
  double worst_x = 0.0;
};

void Accumulate(Region& r, double x, uint64_t u) {
  ++r.n;
  if (u > 0) ++r.miss;
  if (u > r.max_ulp) {
    r.max_ulp = u;
    r.worst_x = x;
  }
}

void Report(const Region& r) {
  std::printf(
      "  %-10s n=%6zu  max ULP=%3llu (gate %llu)  not-CR: %zu (%.2f%%)  "
      "worst x=%.17g\n",
      r.name, r.n, static_cast<unsigned long long>(r.max_ulp),
      static_cast<unsigned long long>(r.gate), r.miss,
      r.n ? 100.0 * static_cast<double>(r.miss) / static_cast<double>(r.n)
          : 0.0,
      r.worst_x);
}

int RunOne(const char* label, const char* path, Fn fn) {
  const auto rows = corvus_test::LoadRef(path, 2, 10000);

  std::vector<double> in;
  std::vector<double> want;
  in.reserve(rows.size());
  want.reserve(rows.size());
  for (const auto& row : rows) {
    in.push_back(corvus_test::ParseDouble(row.tok[0], path, row.line));
    want.push_back(corvus_test::ParseDouble(row.tok[1], path, row.line));
  }

  std::vector<double> got(in.size());
  SplitCall(fn, in, got);

  Region pos_small{"pos small", kMaxUlp};
  Region neg_small{"neg small", kMaxUlp};
  Region pos_huge{"pos huge", kMaxUlp};
  Region neg_huge{"neg huge", kMaxUlp};

  int rc = 0;
  for (size_t i = 0; i < in.size(); ++i) {
    const double x = in[i];
    uint64_t u;
    if (want[i] == 0.0) {
      // #14 N6: UlpDiff maps +0 and -0 to one point, so an exact-zero
      // reference row is held to bit identity. No row rounds to zero today
      // (specials live in the smoke test); this guards the next regen.
      const bool ok = SameBits(got[i], want[i]);
      if (!ok) {
        std::fprintf(stderr,
                     "FAIL: %s(%.17g) signed-zero mismatch: got=%.17g "
                     "want=%.17g\n",
                     label, x, got[i], want[i]);
        rc = 1;
      }
      u = ok ? 0 : ~uint64_t{0};
    } else {
      u = UlpDiff(got[i], want[i]);
    }
    const bool small = std::fabs(x) <= corvus::detail::kTrigDMax;
    if (x >= 0.0) {
      Accumulate(small ? pos_small : pos_huge, x, u);
    } else {
      Accumulate(small ? neg_small : neg_huge, x, u);
    }
  }

  std::printf("%s (%s):\n", label, path);
  for (const Region* r : {&pos_small, &neg_small, &pos_huge, &neg_huge}) {
    Report(*r);
    // #14 N10: an empty bucket is a gate that never ran; all four are
    // populated by generator design.
    if (r->n == 0) {
      std::fprintf(stderr, "FAIL: %s %s bucket is empty (n=0) -- vacuous gate\n",
                   label, r->name);
      rc = 1;
    }
    if (r->max_ulp > r->gate) {
      std::fprintf(stderr, "FAIL: %s %s max ULP %llu exceeds gate %llu\n",
                   label, r->name,
                   static_cast<unsigned long long>(r->max_ulp),
                   static_cast<unsigned long long>(r->gate));
      rc = 1;
    }
  }
  return rc;
}

}  // namespace

int main(int argc, char** argv) {
  if (!corvus_test::ReportAndCheckTarget()) return 2;

  const char* cos_path = argc > 1 ? argv[1] : "tests/data/cos_reference.txt";
  const char* sin_path = argc > 2 ? argv[2] : "tests/data/sin_reference.txt";

  int rc = 0;
  rc |= RunOne("cos", cos_path, corvus::cos);
  rc |= RunOne("sin", sin_path, corvus::sin);
  if (rc == 0) std::printf("PASS\n");
  return rc;
}
