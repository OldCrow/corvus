// Measures max ULP deviation of corvus::exp against the mpmath-generated
// correctly-rounded reference (tests/data/exp_reference.txt), split by
// RESULT class -- that is where exp's contract lives: normal results (the
// dd core + one rounding), subnormal results (the ScaleTwo one-effective-
// rounding claim), exact +0 rows (underflow past the last subnormal) and
// exact +inf rows (overflow) -- the reference brackets both crossover
// boundaries with +-160-ulp bit ladders, so these buckets certify that the
// kernel's saturation happens at the correctly rounded thresholds without
// any threshold blend in the kernel.
//
// INF AND ZERO ROWS ARE NOT A ULP METRIC (ulp_utils.h policy): both are
// held to exact bit matches, separate from the finite buckets.
//
// corvus::expm1 (#45) shares the TU and the exp_dd core, so it is gated
// here too, against tests/data/expm1_reference.txt: the series band
// (|x| <= 2^-6), the exponential branch (everything else with a finite
// result, including the exact -1 below x = -38), and exact +inf rows on
// the overflow ladder.
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "corvus/corvus.h"
#include "expect_target.h"
#include "ulp_utils.h"

namespace {

using corvus_test::SameBits;
using corvus_test::UlpDiff;

// #14 N4: masked-tail split -- neither subspan is a lane multiple for any
// tier, so the masked LoadN/StoreN path always runs.
using Fn = void (*)(std::span<const double>, std::span<double>);

void SplitCall(Fn fn, const std::vector<double>& in, std::vector<double>& out) {
  const size_t n = in.size();
  const size_t split = n - 3;
  fn(std::span<const double>(in).subspan(0, split),
     std::span<double>(out).subspan(0, split));
  fn(std::span<const double>(in).subspan(split),
     std::span<double>(out).subspan(split));
}

// Gates PINNED to measured, no margin (docs/ACCURACY.md), 2026-08-30 on
// AVX3_ZEN4 native AND the SSE2 cap (bucket-identical, same 26 rows):
// every NORMAL result in the set is correctly rounded (the core's ~2^-17
// ulp tie window never bites on these rows -- gate 0), and the subnormal
// band's fl(m.hi+m.lo)-then-ScaleTwo composition measures 1 (design bound
// ~0.75 ulp in the output's own ulp; 0.7495 measured in the top subnormal
// binade, see docs/ACCURACY.md).
constexpr uint64_t kGateNormal = 0;
constexpr uint64_t kGateSubnormal = 1;
// The contract's near-tie exception (corvus.h: "correctly rounded exp
// except within ~2^-15 ulp of a rounding tie"; #76, review L-6): a NORMAL
// result that close to a midpoint may land on either neighbour. Such a row
// -- stored residual |r| >= 1/2 - kTieWindow and the kernel 1 ULP off --
// counts as a tie, reported but not a gate miss. The witness
// exp(2^-53) = 1 + 2^-53 + 2^-107 (CR 1 + 2^-52, kernel 1.0) is in the set.
// Same window as the hunt family's TIE (tools/hunt/families/elementary.py).
constexpr double kTieWindow = 0x1p-15;
// expm1's exponential branch amplifies the core's ~2^-70 by up to ~65x just
// above the 2^-6 cut: ~2^-64 relative, so its documented near-tie window is
// ~2^-11 ulp (corvus.h). The series band keeps exp's window.
constexpr double kTieWindowExpm1 = 0x1p-11;

constexpr double kMinNormal = 0x1p-1022;

struct Region {
  const char* name;
  uint64_t gate;
  uint64_t max_ulp = 0;
  size_t n = 0;
  size_t miss = 0;
  double worst_x = 0.0;
  size_t ties = 0;  // near-tie rows under the contract's exception
  corvus_test::RealMax real{};
};

void Accumulate(Region& r, double x, uint64_t u, double real_err, size_t row,
                bool tie = false) {
  r.real.Add(real_err, row);
  ++r.n;
  if (u > 0) ++r.miss;
  if (tie && u == 1) {
    ++r.ties;
    return;
  }
  if (u > r.max_ulp) {
    r.max_ulp = u;
    r.worst_x = x;
  }
}

void Report(const Region& r) {
  std::printf(
      "  %-10s n=%6zu  max ULP=%3llu (gate %llu)  not-CR: %zu (%.2f%%), "
      "near-ties %zu  worst x=%.17g\n",
      r.name, r.n, static_cast<unsigned long long>(r.max_ulp),
      static_cast<unsigned long long>(r.gate), r.miss,
      r.n ? 100.0 * static_cast<double>(r.miss) / static_cast<double>(r.n)
          : 0.0,
      r.ties, r.worst_x);
}

int RunExp(const char* path) {
  const auto rows = corvus_test::LoadRef(path, 2, 10000);

  std::vector<double> in, want, res;
  in.reserve(rows.size());
  want.reserve(rows.size());
  res.reserve(rows.size());
  for (const auto& row : rows) {
    in.push_back(corvus_test::ParseDouble(row.tok[0], path, row.line));
    want.push_back(corvus_test::ParseDouble(row.tok[1], path, row.line));
    res.push_back(corvus_test::Residual(row, 0));
  }

  std::vector<double> got(in.size());
  SplitCall(corvus::exp, in, got);

  Region normal{"normal", kGateNormal};
  Region subnormal{"subnormal", kGateSubnormal};
  size_t inf_n = 0, inf_fail = 0;
  size_t zero_n = 0, zero_fail = 0;
  int rc = 0;

  for (size_t i = 0; i < in.size(); ++i) {
    const double x = in[i];
    const double w = want[i];
    if (std::isinf(w)) {
      // Overflow boundary: exact +inf, nothing else.
      ++inf_n;
      if (!SameBits(got[i], w)) {
        ++inf_fail;
        std::fprintf(stderr, "FAIL: exp(%.17g) = %.17g, want +inf\n", x,
                     got[i]);
      }
      continue;
    }
    if (w == 0.0) {
      // Underflow past the last subnormal: exact +0 (exp is never
      // negative, and -0 would be a sign regression UlpDiff cannot see --
      // #14 N6).
      ++zero_n;
      if (!SameBits(got[i], 0.0)) {
        ++zero_fail;
        std::fprintf(stderr, "FAIL: exp(%.17g) = %.17g, want exact +0\n", x,
                     got[i]);
      }
      continue;
    }
    const uint64_t u = UlpDiff(got[i], w);
    // The +inf and +0 rows above are exact-match checks and carry no real gate.
    const bool tie = w >= kMinNormal && std::fabs(res[i]) >= 0.5 - kTieWindow;
    Accumulate(w < kMinNormal ? subnormal : normal, x, u,
               corvus_test::RealErr(got[i], w, res[i]), i, tie);
  }

  std::printf("exp (%s):\n", path);
  Report(normal);
  Report(subnormal);
  std::printf("  %-10s n=%6zu  exact-match fails=%zu\n", "inf", inf_n,
              inf_fail);
  std::printf("  %-10s n=%6zu  exact-match fails=%zu\n", "zero", zero_n,
              zero_fail);

  for (const Region* r : {&normal, &subnormal}) {
    // #14 N10: an empty bucket is a gate that never ran; the reference's
    // subnormal-output band is a designed stratum.
    if (r->n == 0) {
      std::fprintf(stderr, "FAIL: exp %s bucket is empty (n=0) -- vacuous gate\n",
                   r->name);
      rc = 1;
    }
    if (r->max_ulp > r->gate) {
      std::fprintf(stderr, "FAIL: exp %s max ULP %llu exceeds gate %llu\n",
                   r->name, static_cast<unsigned long long>(r->max_ulp),
                   static_cast<unsigned long long>(r->gate));
      rc = 1;
    }
    // #82 T1: real error against the contract ceiling and the regression pin.
    if (corvus_test::CheckReal(path, "exp", r->name, r->gate, r->real) != 0) {
      rc = 1;
    }
  }
  if (inf_n == 0 || zero_n == 0) {
    std::fprintf(stderr,
                 "FAIL: exp inf/zero boundary bucket empty (inf n=%zu, zero "
                 "n=%zu) -- vacuous gate\n",
                 inf_n, zero_n);
    rc = 1;
  }
  if (inf_fail != 0 || zero_fail != 0) rc = 1;
  return rc;
}

int RunExpm1(const char* path) {
  const auto rows = corvus_test::LoadRef(path, 2, 10000);
  std::vector<double> in, want, res;
  for (const auto& row : rows) {
    in.push_back(corvus_test::ParseDouble(row.tok[0], path, row.line));
    want.push_back(corvus_test::ParseDouble(row.tok[1], path, row.line));
    res.push_back(corvus_test::Residual(row, 0));
  }
  std::vector<double> got(in.size());
  SplitCall(corvus::expm1, in, got);

  Region series{"series", 0};
  Region expo{"exp", 0};
  size_t inf_n = 0, inf_fail = 0;
  int rc = 0;
  for (size_t i = 0; i < in.size(); ++i) {
    const double x = in[i];
    const double w = want[i];
    if (std::isinf(w)) {
      // Overflow ladder: exact +inf, nothing else.
      ++inf_n;
      if (!SameBits(got[i], w)) {
        ++inf_fail;
        std::fprintf(stderr, "FAIL: expm1(%.17g) = %.17g, want +inf\n", x,
                     got[i]);
      }
      continue;
    }
    const uint64_t u = UlpDiff(got[i], w);
    const bool ser = std::fabs(x) <= 0x1p-6;
    const bool tie =
        std::fabs(res[i]) >= 0.5 - (ser ? kTieWindow : kTieWindowExpm1);
    Accumulate(ser ? series : expo, x, u,
               corvus_test::RealErr(got[i], w, res[i]), i, tie);
  }

  std::printf("expm1 (%s):\n", path);
  Report(series);
  Report(expo);
  std::printf("  %-10s n=%6zu  exact-match fails=%zu\n", "inf", inf_n,
              inf_fail);
  for (const Region* r : {&series, &expo}) {
    // #14 N10: an empty bucket is a gate that never ran.
    if (r->n == 0) {
      std::fprintf(stderr,
                   "FAIL: expm1 %s bucket is empty (n=0) -- vacuous gate\n",
                   r->name);
      rc = 1;
    }
    if (r->max_ulp > r->gate) {
      std::fprintf(stderr, "FAIL: expm1 %s max ULP %llu exceeds gate %llu\n",
                   r->name, static_cast<unsigned long long>(r->max_ulp),
                   static_cast<unsigned long long>(r->gate));
      rc = 1;
    }
    if (corvus_test::CheckReal(path, "expm1", r->name, r->gate, r->real) != 0) {
      rc = 1;
    }
  }
  if (inf_n == 0) {
    std::fprintf(stderr, "FAIL: expm1 inf bucket empty -- vacuous gate\n");
    rc = 1;
  }
  if (inf_fail != 0) rc = 1;
  return rc;
}

}  // namespace

int main(int argc, char* const* argv) {
  if (!corvus_test::ReportAndCheckTarget()) return 2;
  const char* exp_path = argc > 1 ? argv[1] : "tests/data/exp_reference.txt";
  const char* expm1_path =
      argc > 2 ? argv[2] : "tests/data/expm1_reference.txt";
  int rc = 0;
  rc |= RunExp(exp_path);
  rc |= RunExpm1(expm1_path);
  if (rc == 0) std::printf("PASS\n");
  return rc;
}
