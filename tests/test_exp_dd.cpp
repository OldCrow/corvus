// Accuracy gate for the internal double-double exp kernel (src/exp_dd-inl.h).
//
// exp_dd has no public API -- it exists so that no accuracy-critical corvus
// kernel depends on the backend's exp -- so this test compiles the kernel
// itself through foreach_target rather than linking the library, and drives it
// directly. Two things are measured, because rounding to a double would hide
// most of what the kernel is for:
//   * the DD relative error against a dd oracle, which is the property the
//     erfc tail and (later) the incomplete gamma prefactor actually consume;
//   * the ULP error of the rounded result, including the subnormal band and
//     the flush-to-zero edge.
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <vector>

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "tests/test_exp_dd.cpp"
#include "hwy/foreach_target.h"  // IWYU pragma: keep

#include "src/exp_dd-inl.h"
#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {

// Same masked-tail discipline as the shipped kernels: one code path, no
// scalar fallback for the remainder.
void ExpDdTestImpl(const double* xh, const double* xl, double* hi, double* lo,
                   size_t n) {
  const op::ScalableTag<double> d;
  const size_t N = op::Lanes(d);
  size_t i = 0;
  for (; i + N <= n; i += N) {
    const auto r = ExpDd(d, op::Load(d, xh + i), op::Load(d, xl + i));
    op::Store(r.hi, d, hi + i);
    op::Store(r.lo, d, lo + i);
  }
  if (i < n) {
    const size_t m = n - i;
    const auto r = ExpDd(d, op::LoadN(d, xh + i, m), op::LoadN(d, xl + i, m));
    op::StoreN(r.hi, d, hi + i, m);
    op::StoreN(r.lo, d, lo + i, m);
  }
}

// This TU carries its own dispatch table. It is built with the same target
// defines as the library (CORVUS_HWY_TARGET_DEFS), and main() asserts the two
// tables agreed -- otherwise the reported tier would describe the library
// while the numbers came from here.
const char* TestTargetNameImpl() { return hwy::TargetName(HWY_TARGET); }

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#if HWY_ONCE

#include "expect_target.h"
#include "ulp_utils.h"

namespace corvus {
HWY_EXPORT(ExpDdTestImpl);
HWY_EXPORT(TestTargetNameImpl);

// The dispatch MUST happen inside namespace corvus, exactly as the shipped
// kernels do it. When a cap leaves a single target, Highway collapses
// HWY_DYNAMIC_DISPATCH to HWY_STATIC_DISPATCH, which prefixes the call with
// N_<target>:: -- so a qualified corvus::Foo call from global scope becomes
// N_SSE2::corvus::Foo and fails to compile. It builds fine at every other
// tier, which is what makes it worth a comment rather than a fix in silence.
void ExpDdTest(const double* xh, const double* xl, double* hi, double* lo,
               size_t n) {
  HWY_DYNAMIC_DISPATCH(ExpDdTestImpl)(xh, xl, hi, lo, n);
}
const char* TestTargetName() { return HWY_DYNAMIC_DISPATCH(TestTargetNameImpl)(); }
}  // namespace corvus

namespace {

// Gates set from measured values (see docs/ACCURACY.md), no margin.
// kMinRelExp is -log2 of the worst dd relative error: HIGHER is better, so
// this gate is a floor, unlike the ULP gates. It is the measurement (2^-68.45,
// identical on every validated tier) floored to 0.1 bits -- the same "round
// to the gate's own granularity" the integer ULP gates get for free. Do not
// round it to 68.5: that is above the measured value and fails immediately.
constexpr double kMinRelExp = 68.4;
constexpr uint64_t kMaxUlpNormal = 1;
constexpr uint64_t kMaxUlpSubnormal = 1;

using corvus_test::UlpDiff;

struct Region {
  const char* name;
  uint64_t bound;
  uint64_t max_ulp = 0;
  size_t n = 0;
  size_t miss = 0;
  double worst_x = 0.0;
};

}  // namespace

int main(int argc, char** argv) {
  if (!corvus_test::ReportAndCheckTarget()) return 2;

  const char* here = corvus::TestTargetName();
  if (std::strcmp(here, corvus::active_target()) != 0) {
    std::fprintf(stderr,
                 "FAIL: this TU dispatched '%s' but the library dispatched "
                 "'%s' — the two target sets are not the same build.\n",
                 here, corvus::active_target());
    return 2;
  }

  const char* path = argc > 1 ? argv[1] : "tests/data/exp_dd_reference.txt";
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
  // ExpDdTestImpl runs on every tier regardless of whether the reference
  // set's own length happens to be a lane multiple (#14 N4).
  const size_t n = xh.size();
  const size_t split = n - 3;
  corvus::ExpDdTest(xh.data(), xl.data(), hi.data(), lo.data(), split);
  corvus::ExpDdTest(xh.data() + split, xl.data() + split, hi.data() + split,
                    lo.data() + split, n - split);

  constexpr double kMinNormal = 2.2250738585072014e-308;
  Region regions[] = {
      {"normal", kMaxUlpNormal},
      {"subnormal", kMaxUlpSubnormal},
      {"zero/inf", 0},
  };

  double worst_rel = 0.0;      // dd relative error, normal results only
  double worst_rel_x = 0.0;
  size_t n_rel = 0;

  for (size_t i = 0; i < xh.size(); ++i) {
    const double w = want_hi[i];
    Region& r = (w >= kMinNormal && std::isfinite(w))
                    ? regions[0]
                    : (w > 0.0 ? regions[1] : regions[2]);
    ++r.n;

    // The kernel's answer as a double is the single rounding of hi + lo.
    const uint64_t u = UlpDiff(hi[i] + lo[i], w);
    if (u > 0) ++r.miss;
    if (u > r.max_ulp) {
      r.max_ulp = u;
      r.worst_x = xh[i];
    }

    // DD error. Both differences are exact (each pair agrees to within an ulp
    // of its own magnitude), so this measures the kernel below the last bit of
    // a double rather than at it.
    if (&r == &regions[0]) {
      const double err = (hi[i] - w) + (lo[i] - want_lo[i]);
      const double rel = std::abs(err / w);
      ++n_rel;
      if (rel > worst_rel) {
        worst_rel = rel;
        worst_rel_x = xh[i];
      }
    }
  }

  int rc = 0;
  const double rel_bits = worst_rel > 0.0 ? -std::log2(worst_rel) : 1000.0;
  std::printf("dd relative    n=%6zu  worst=2^-%.2f (gate 2^-%.2f)  worst x=%.17g\n",
              n_rel, rel_bits, kMinRelExp, worst_rel_x);
  if (rel_bits < kMinRelExp) {
    std::fprintf(stderr, "FAIL: dd relative error exceeds gate\n");
    rc = 1;
  }
  for (const Region& r : regions) {
    std::printf("%-14s n=%6zu  max ULP=%3llu (gate %llu)  not-CR: %zu (%.2f%%)  worst x=%.17g\n",
                r.name, r.n, static_cast<unsigned long long>(r.max_ulp),
                static_cast<unsigned long long>(r.bound), r.miss,
                r.n ? 100.0 * static_cast<double>(r.miss) / static_cast<double>(r.n) : 0.0,
                r.worst_x);
    if (r.max_ulp > r.bound) {
      std::fprintf(stderr, "FAIL: %s exceeds gate\n", r.name);
      rc = 1;
    }
  }
  if (rc == 0) std::printf("PASS: all regions within gates\n");
  return rc;
}

#endif  // HWY_ONCE
