// Measures max ULP deviation of corvus::i0/i1/i0e/i1e against the
// mpmath-generated correctly-rounded references (tests/data/i{0,1}{,e}_
// reference.txt), split by the kernel's own regions: series (|x| <=
// kBesselSplit), tail (|x| > kBesselSplit, finite expected result), and --
// for i0/i1 only -- the overflow BOUNDARY bracket, where the expected value
// is +-inf. Each of series/tail is further split by sign, so a sign-carrying
// regression in the odd i1/i1e pair or an even-function asymmetry in i0/i0e
// names itself rather than hiding inside a combined bucket.
//
// INF ROWS ARE NOT A ULP METRIC. ULP distance is undefined against an
// infinite reference (OrderedBits(inf) is a specific finite int64, so a
// naive UlpDiff would silently produce a huge-but-bounded number instead of
// flagging the comparison as ill-posed); those rows are instead held to an
// EXACT bit match (value and sign), separate from every other bucket.
//
// SUSPICIOUSLY-SMALL THRESHOLD. House ULP tests hard-code `size < 10000` as
// a sanity check against a truncated or empty reference file; the Bessel
// references are deliberately 2351-2515 rows, so this file uses 2000
// instead.
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#include "corvus/corvus.h"
#include "expect_target.h"
#include "src/bessel_data.h"
#include "ulp_utils.h"

namespace {

using corvus_test::OrderedBits;
using corvus_test::SameBits;
using corvus_test::UlpDiff;

using Fn = void (*)(std::span<const double>, std::span<double>);

// #14 N4: masked-tail split. Splitting into subspans [0, n-3) and [n-3, n)
// makes neither half a lane multiple regardless of what n is, so the
// masked-tail path always runs; each call is independently masked and
// stateless, so this is row-for-row identical to one whole-set call.
void SplitCall(Fn fn, const std::vector<double>& in, std::vector<double>& out) {
  const size_t n = in.size();
  const size_t split = n - 3;
  fn(std::span<const double>(in).subspan(0, split),
     std::span<double>(out).subspan(0, split));
  fn(std::span<const double>(in).subspan(split),
     std::span<double>(out).subspan(split));
}

// Gate PINNED to measured, no margin. Both regimes land at the design's
// expected 1 ULP ceiling on every validated tier.
constexpr uint64_t kMaxUlp = 1;

struct Region {
  const char* name;
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
      "  %-16s n=%6zu  max ULP=%3llu (gate %llu)  not-CR: %zu (%.2f%%)  "
      "worst x=%.17g\n",
      r.name, r.n, static_cast<unsigned long long>(r.max_ulp),
      static_cast<unsigned long long>(kMaxUlp), r.miss,
      r.n ? 100.0 * static_cast<double>(r.miss) / static_cast<double>(r.n)
          : 0.0,
      r.worst_x);
}

// Runs one reference file against `fn`. Returns 0 on success, 1 on any gate
// failure. `has_inf` selects whether an inf/-inf expected token is a
// meaningful (i0/i1, past the overflow boundary) or unexpected (i0e/i1e,
// which never overflow) row -- either way it is handled by the exact-match
// branch, but only i0/i1 are EXPECTED to exercise it.
int RunOne(const char* label, const char* path, Fn fn, bool has_inf) {
  const auto rows = corvus_test::LoadRef(path, 2, 2000);

  std::vector<double> in;
  std::vector<double> want;
  std::vector<bool> want_inf;
  in.reserve(rows.size());
  want.reserve(rows.size());
  want_inf.reserve(rows.size());
  for (const auto& row : rows) {
    in.push_back(corvus_test::ParseDouble(row.tok[0], path, row.line));
    const double w = corvus_test::ParseDouble(row.tok[1], path, row.line);
    want.push_back(w);
    want_inf.push_back(std::isinf(w));
  }

  std::vector<double> got(in.size());
  SplitCall(fn, in, got);

  Region pos_series{"pos series"};
  Region neg_series{"neg series"};
  Region pos_tail{"pos tail"};
  Region neg_tail{"neg tail"};
  size_t inf_n = 0, inf_fail = 0;

  for (size_t i = 0; i < in.size(); ++i) {
    const double x = in[i];
    if (want_inf[i]) {
      ++inf_n;
      // Exact match: value AND sign (i0's saturation is unsigned +inf; i1's
      // carries the sign of x). No ULP distance is computed here.
      const bool ok =
          got[i] == want[i] && std::signbit(got[i]) == std::signbit(want[i]);
      if (!ok) {
        ++inf_fail;
        std::fprintf(stderr,
                     "FAIL: %s(%.17g) = %.17g, want %.17g (inf boundary)\n",
                     label, x, got[i], want[i]);
      }
      continue;
    }
    uint64_t u;
    if (want[i] == 0.0) {
      // #14 N6: i1/i1e are odd, so their reference carries a real +-0 at
      // x = +-0. UlpDiff maps +0 and -0 to the same point (ulp_utils.h's
      // policy), so a signed-zero regression needs an exact-bits check.
      const bool ok = SameBits(got[i], want[i]);
      if (!ok) {
        std::fprintf(stderr,
                     "FAIL: %s(%.17g) signed-zero mismatch: got=%.17g "
                     "want=%.17g\n",
                     label, x, got[i], want[i]);
      }
      u = ok ? 0 : ~uint64_t{0};
    } else {
      u = UlpDiff(got[i], want[i]);
    }
    const bool series = std::fabs(x) <= corvus::detail::kBesselSplit;
    if (x >= 0.0) {
      Accumulate(series ? pos_series : pos_tail, x, u);
    } else {
      Accumulate(series ? neg_series : neg_tail, x, u);
    }
  }

  std::printf("%s (%s):\n", label, path);
  Report(pos_series);
  Report(neg_series);
  Report(pos_tail);
  Report(neg_tail);
  if (has_inf) {
    std::printf("  %-16s n=%6zu  exact-match fails=%zu\n", "boundary (inf)",
               inf_n, inf_fail);
  } else if (inf_n != 0) {
    // i0e/i1e never overflow -- an inf row here would mean the reference
    // itself is wrong for this file, not a kernel defect, but it is still
    // worth surfacing loudly rather than silently comparing it.
    std::fprintf(stderr,
                 "%s: %zu unexpected inf row(s) in a non-overflowing "
                 "function's reference\n",
                 label, inf_n);
  }

  int rc = 0;
  for (const Region* r : {&pos_series, &neg_series, &pos_tail, &neg_tail}) {
    // #14 N10: a bucket that never accumulated a row is a gate that never
    // ran -- the four series/tail x sign combinations are all populated by
    // design (the file banner), so n=0 here means a coverage gap.
    if (r->n == 0) {
      std::fprintf(stderr, "FAIL: %s %s bucket is empty (n=0) -- vacuous gate\n",
                   label, r->name);
      rc = 1;
    }
    if (r->max_ulp > kMaxUlp) {
      std::fprintf(stderr, "FAIL: %s %s exceeds gate\n", label, r->name);
      rc = 1;
    }
  }
  // has_inf functions (i0/i1) are documented to carry an overflow-boundary
  // bracket in their reference set; has_inf=false functions (i0e/i1e) are
  // already asserted not to, above.
  if (has_inf && inf_n == 0) {
    std::fprintf(stderr,
                 "FAIL: %s boundary (inf) bucket is empty (n=0) -- vacuous "
                 "gate\n",
                 label);
    rc = 1;
  }
  if (inf_fail != 0) rc = 1;
  return rc;
}

}  // namespace

int main(int argc, char** argv) {
  // Before doing any work: a wrong tier makes every number below meaningless.
  if (!corvus_test::ReportAndCheckTarget()) return 2;

  const char* path_i0 = argc > 1 ? argv[1] : "tests/data/i0_reference.txt";
  const char* path_i1 = argc > 2 ? argv[2] : "tests/data/i1_reference.txt";
  const char* path_i0e = argc > 3 ? argv[3] : "tests/data/i0e_reference.txt";
  const char* path_i1e = argc > 4 ? argv[4] : "tests/data/i1e_reference.txt";

  int rc = 0;
  rc |= RunOne("i0", path_i0, corvus::i0, /*has_inf=*/true);
  rc |= RunOne("i1", path_i1, corvus::i1, /*has_inf=*/true);
  rc |= RunOne("i0e", path_i0e, corvus::i0e, /*has_inf=*/false);
  rc |= RunOne("i1e", path_i1e, corvus::i1e, /*has_inf=*/false);

  if (rc == 0) std::printf("PASS: all regions within gates\n");
  return rc;
}
