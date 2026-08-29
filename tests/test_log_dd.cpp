// Accuracy gate for the internal double-double log kernel (src/log_dd-inl.h).
// Structure and rationale mirror test_exp_dd.cpp: the kernel has no public
// API, so this TU compiles the header through foreach_target and drives it
// directly, and the oracle carries double-double pairs so the gate can sit
// below the last bit of a double.
//
// log has no in-tree consumer yet (lgamma is Phase B), so this gate IS the
// acceptance test for the kernel.
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <vector>

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "tests/test_log_dd.cpp"
#include "hwy/foreach_target.h"  // IWYU pragma: keep

#include "src/log_dd-inl.h"
#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {

// Drives the dd-input overload, which subsumes the plain-double one (x_lo is
// 0 for most reference points). Masked tail, same as the shipped kernels.
void LogDdTestImpl(const double* xh, const double* xl, double* hi, double* lo,
                   size_t n) {
  // Named rather than decltype(d): d is const-qualified, and Dd<const Simd<>>
  // is a different type from the Dd<Simd<>> the kernel deduces.
  using Tag = op::ScalableTag<double>;
  const Tag d;
  const size_t N = op::Lanes(d);
  size_t i = 0;
  for (; i + N <= n; i += N) {
    const Dd<Tag> x{op::Load(d, xh + i), op::Load(d, xl + i)};
    const auto r = LogDd(d, x);
    op::Store(r.hi, d, hi + i);
    op::Store(r.lo, d, lo + i);
  }
  if (i < n) {
    const size_t m = n - i;
    const Dd<Tag> x{op::LoadN(d, xh + i, m), op::LoadN(d, xl + i, m)};
    const auto r = LogDd(d, x);
    op::StoreN(r.hi, d, hi + i, m);
    op::StoreN(r.lo, d, lo + i, m);
  }
}

const char* TestTargetNameImpl() { return hwy::TargetName(HWY_TARGET); }

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#if HWY_ONCE

#include "expect_target.h"
#include "ulp_utils.h"

namespace corvus {
HWY_EXPORT(LogDdTestImpl);
HWY_EXPORT(TestTargetNameImpl);

// Dispatch from inside namespace corvus: with a single compiled target (the
// SSE2 cap) Highway collapses HWY_DYNAMIC_DISPATCH to N_SSE2::FUNC, and a
// globally qualified call would name a namespace that does not exist.
void LogDdTest(const double* xh, const double* xl, double* hi, double* lo,
               size_t n) {
  HWY_DYNAMIC_DISPATCH(LogDdTestImpl)(xh, xl, hi, lo, n);
}
const char* TestTargetName() { return HWY_DYNAMIC_DISPATCH(TestTargetNameImpl)(); }
}  // namespace corvus

namespace {

// Gates set from measured values, no margin (see docs/ACCURACY.md and the
// note in test_exp_dd.cpp on why the dd gate is a floor, floored to 0.1 bits).
// 2^-67.88, identical on every tier and every compiler.
//
// It briefly read 2^-68.48 under GCC only: GCC's default -ffp-contract=fast
// was fusing a Mul into an adjacent Add inside the log1p path, which happened
// to be more accurate and was invisible to MSVC and to the no-FMA tiers. The
// build now sets -ffp-contract=off (see CORVUS_FP_FLAGS in the top-level
// CMakeLists and the rationale there); this gate is the honest, portable
// figure. If it ever passes by more than a rounding again, suspect the flag
// stopped being applied before believing the kernel improved.
constexpr double kMinRelExp = 67.8;
constexpr uint64_t kMaxUlp = 1;

using corvus_test::UlpDiff;

}  // namespace

int main(int argc, char** argv) {
  if (!corvus_test::ReportAndCheckTarget()) return 2;
  if (std::strcmp(corvus::TestTargetName(), corvus::active_target()) != 0) {
    std::fprintf(stderr, "FAIL: this TU and the library dispatched different targets\n");
    return 2;
  }

  const char* path = argc > 1 ? argv[1] : "tests/data/log_dd_reference.txt";
  const auto rows = corvus_test::LoadRef(path, 4, 5000);

  std::vector<double> xh, xl, want_hi, want_lo;
  for (const auto& row : rows) {
    xh.push_back(corvus_test::ParseDouble(row.tok[0], path, row.line));
    xl.push_back(corvus_test::ParseDouble(row.tok[1], path, row.line));
    want_hi.push_back(corvus_test::ParseDouble(row.tok[2], path, row.line));
    want_lo.push_back(corvus_test::ParseDouble(row.tok[3], path, row.line));
  }

  std::vector<double> hi(xh.size()), lo(xh.size());
  // Two calls, not one: the second covers only the last 3 rows, a length
  // below every lane count, so the masked LoadN/StoreN tail path inside
  // LogDdTestImpl runs on every tier regardless of whether the reference
  // set's own length happens to be a lane multiple (#14 N4).
  const size_t n = xh.size();
  const size_t split = n - 3;
  corvus::LogDdTest(xh.data(), xl.data(), hi.data(), lo.data(), split);
  corvus::LogDdTest(xh.data() + split, xl.data() + split, hi.data() + split,
                    lo.data() + split, n - split);

  double worst_rel = 0.0, worst_rel_x = 0.0;
  uint64_t worst_ulp = 0;
  double worst_ulp_x = 0.0;
  size_t miss = 0;
  // Reported separately: |x| within a factor 2 of 1, where log is small and
  // relative accuracy depends entirely on the centred mantissa and the exact
  // (R, L) = (1, 0) slots. A regression there would otherwise be diluted.
  double worst_near1 = 0.0, worst_near1_x = 0.0;
  size_t n_near1 = 0;

  for (size_t i = 0; i < xh.size(); ++i) {
    const double w = want_hi[i];
    const uint64_t u = UlpDiff(hi[i] + lo[i], w);
    if (u > 0) ++miss;
    if (u > worst_ulp) {
      worst_ulp = u;
      worst_ulp_x = xh[i];
    }
    const double rel = std::abs(((hi[i] - w) + (lo[i] - want_lo[i])) / w);
    if (rel > worst_rel) {
      worst_rel = rel;
      worst_rel_x = xh[i];
    }
    if (xh[i] > 0.5 && xh[i] < 2.0) {
      ++n_near1;
      if (rel > worst_near1) {
        worst_near1 = rel;
        worst_near1_x = xh[i];
      }
    }
  }

  int rc = 0;
  const double rel_bits = worst_rel > 0.0 ? -std::log2(worst_rel) : 1000.0;
  const double n1_bits = worst_near1 > 0.0 ? -std::log2(worst_near1) : 1000.0;
  std::printf("dd relative    n=%6zu  worst=2^-%.2f (gate 2^-%.2f)  worst x=%.17g\n",
              xh.size(), rel_bits, kMinRelExp, worst_rel_x);
  std::printf("  of which 0.5<x<2   n=%6zu  worst=2^-%.2f  worst x=%.17g\n",
              n_near1, n1_bits, worst_near1_x);
  std::printf("rounded        n=%6zu  max ULP=%3llu (gate %llu)  not-CR: %zu (%.2f%%)  worst x=%.17g\n",
              xh.size(), static_cast<unsigned long long>(worst_ulp),
              static_cast<unsigned long long>(kMaxUlp), miss,
              100.0 * static_cast<double>(miss) / static_cast<double>(xh.size()),
              worst_ulp_x);

  if (rel_bits < kMinRelExp) {
    std::fprintf(stderr, "FAIL: dd relative error exceeds gate\n");
    rc = 1;
  }
  if (worst_ulp > kMaxUlp) {
    std::fprintf(stderr, "FAIL: rounded ULP exceeds gate\n");
    rc = 1;
  }
  if (rc == 0) std::printf("PASS: all regions within gates\n");
  return rc;
}

#endif  // HWY_ONCE
