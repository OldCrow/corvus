// Accuracy gate for Log1pmxDd (src/dd_special-inl.h), the shared primitive
// phi(u) = u - log1p(u) that the incomplete gamma's Temme region -- and
// later the incomplete beta -- rests on.
//
// phi is not a public function, so this test compiles the kernel header
// itself through foreach_target and drives it directly (the test_exp_dd /
// test_log_dd pattern). Rounding phi to a double would hide the property
// being measured: a*phi with a up to ~4e5 is the argument of an exponential,
// so what has to be small is phi's error BELOW the last bit of a double,
// which is why the reference carries a dd pair.
//
// The two gates split at |u| = 2^-40 for the reason the kernel's own branch
// exists. Above it the series' u.lo correction and, past the 1/16 cut, the
// amplification of log_dd's error by 2/u are what is being watched, and the
// meaningful measure is RELATIVE. Below it phi < 2^-81 and no amplification
// is possible, so an absolute bound is the honest statement.
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "tests/test_dd_special.cpp"
#include "hwy/foreach_target.h"  // IWYU pragma: keep

#include "src/dd_special-inl.h"
#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {

// Same masked-tail discipline as the shipped kernels: one code path. `ulo`
// carries the input's low word (all zero for dd_special_reference.txt).
static void Log1pmxTestImpl(const double* u, const double* ulo, double* hi,
                            double* lo, size_t n) {
  using DT = op::ScalableTag<double>;
  const DT d;
  const size_t N = op::Lanes(d);
  size_t i = 0;
  for (; i + N <= n; i += N) {
    const auto r = Log1pmxDd(
        d, Dd<DT>{op::Load(d, u + i), op::Load(d, ulo + i)});
    op::Store(r.hi, d, hi + i);
    op::Store(r.lo, d, lo + i);
  }
  if (i < n) {
    const size_t m = n - i;
    const auto r = Log1pmxDd(
        d, Dd<DT>{op::LoadN(d, u + i, m), op::LoadN(d, ulo + i, m)});
    op::StoreN(r.hi, d, hi + i, m);
    op::StoreN(r.lo, d, lo + i, m);
  }
}

// Expm1Dd's gate driver (#58), same shape.
static void Expm1TestImpl(const double* w, const double* wlo, double* hi,
                          double* lo, size_t n) {
  using DT = op::ScalableTag<double>;
  const DT d;
  const size_t N = op::Lanes(d);
  size_t i = 0;
  for (; i + N <= n; i += N) {
    const auto r =
        Expm1Dd(d, Dd<DT>{op::Load(d, w + i), op::Load(d, wlo + i)});
    op::Store(r.hi, d, hi + i);
    op::Store(r.lo, d, lo + i);
  }
  if (i < n) {
    const size_t m = n - i;
    const auto r = Expm1Dd(
        d, Dd<DT>{op::LoadN(d, w + i, m), op::LoadN(d, wlo + i, m)});
    op::StoreN(r.hi, d, hi + i, m);
    op::StoreN(r.lo, d, lo + i, m);
  }
}

static const char* UtilTargetNameImpl() { return hwy::TargetName(HWY_TARGET); }

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#if HWY_ONCE

#include "expect_target.h"
#include "ulp_utils.h"

namespace corvus {
HWY_EXPORT(Log1pmxTestImpl);
HWY_EXPORT(Expm1TestImpl);
HWY_EXPORT(UtilTargetNameImpl);

// Dispatch from INSIDE namespace corvus: with a single compiled target
// Highway collapses HWY_DYNAMIC_DISPATCH to N_<target>::FUNC, and a
// globally qualified call would then name a namespace that does not exist.
static void Log1pmxTest(const double* u, const double* ulo, double* hi,
                        double* lo, size_t n) {
  HWY_DYNAMIC_DISPATCH(Log1pmxTestImpl)(u, ulo, hi, lo, n);
}
static void Expm1Test(const double* w, const double* wlo, double* hi,
                      double* lo, size_t n) {
  HWY_DYNAMIC_DISPATCH(Expm1TestImpl)(w, wlo, hi, lo, n);
}
static const char* UtilTargetName() {
  return HWY_DYNAMIC_DISPATCH(UtilTargetNameImpl)();
}
}  // namespace corvus

namespace {

// Gates pinned to MEASURED (#17), floored to the print quantum, so a
// regression larger than the print quantum trips here rather than surfacing
// under a consumer's name. Absolute (|u| < 2^-40): 2^-187.00.
//
// Relative, ONE cell for every |u| >= 2^-40 row of both sets (u.lo = 0 here
// and #58's u.lo != 0 set) [row 12, 2026-10-20, user: option (i)]: 2^-67.78
// at u = 0.0625098, identical on AVX3_ZEN4, SSE4 and SSE2. The worst sits
// just above the 1/16 cut, where log_dd's error (~2^-72.6) is amplified by
// ~2/u. The earlier pins (71.64 here, 71.40 for u.lo != 0) held only
// because no row sampled (1/16, 0.064]; row 12 added 5,257 there. 100k-point
// probes reach 2^-67.63 (u ~ 0.0625809), so the primitive's true worst is
// ~2^-67.6, slightly below this pin's rows. A kernel fix near the cut is
// tracked for v1.2.0.
constexpr double kMinRelLog1pmx = 67.77;  // |u| >= 2^-40, relative, all u.lo
constexpr double kMaxAbs = 0x1.0p-186;  // |u| <  2^-40, absolute
constexpr double kSplit = 0x1.0p-40;

// #58's Expm1Dd set (tools/gen_dd_special_gates.py), relative, in bits. PIN:
// measured, floored to the print quantum, same rule as above; identical on
// AVX3_ZEN4, SSE4 and SSE2 (clang-cl, 2026-10-18): 2^-64.34 at w = -0.0187,
// the exponential branch just above kExpm1Cut = 2^-6 (before #58's fix, cut
// 2^-10: 2^-61.41 at w = -0.003).
constexpr double kMinRelExpm1 = 64.33;  // Expm1Dd, whole range

using DdFn = void (*)(const double*, const double*, double*, double*, size_t);

// Runs a four-column dd set (in.hi in.lo want.hi want.lo) through `fn` as two
// calls (n-3, 3: the masked tail runs on every tier) and returns the worst
// relative error in bits; `worst_in` gets its input's high word.
double WorstRelBits(const char* path, size_t min_rows, DdFn fn,
                    size_t* n_out, double* worst_in) {
  const auto rows = corvus_test::LoadRef(path, 4, min_rows);
  std::vector<double> in, inlo, want_hi, want_lo;
  for (const auto& row : rows) {
    in.push_back(corvus_test::ParseDouble(row.tok[0], path, row.line));
    inlo.push_back(corvus_test::ParseDouble(row.tok[1], path, row.line));
    want_hi.push_back(corvus_test::ParseDouble(row.tok[2], path, row.line));
    want_lo.push_back(corvus_test::ParseDouble(row.tok[3], path, row.line));
  }
  const size_t n = in.size();
  std::vector<double> hi(n), lo(n);
  fn(in.data(), inlo.data(), hi.data(), lo.data(), n - 3);
  fn(in.data() + n - 3, inlo.data() + n - 3, hi.data() + n - 3,
     lo.data() + n - 3, 3);
  double worst = 0.0;
  *worst_in = 0.0;
  for (size_t i = 0; i < n; ++i) {
    const double err = (hi[i] - want_hi[i]) + (lo[i] - want_lo[i]);
    const double rel = std::abs(err / want_hi[i]);
    if (rel > worst) {
      worst = rel;
      *worst_in = in[i];
    }
  }
  *n_out = n;
  return worst > 0.0 ? -std::log2(worst) : 1000.0;
}

}  // namespace

int main(int argc, char* const* argv) {
  if (!corvus_test::ReportAndCheckTarget()) return 2;

  const char* here = corvus::UtilTargetName();
  if (std::strcmp(here, corvus::active_target()) != 0) {
    std::fprintf(stderr,
                 "FAIL: this TU dispatched '%s' but the library dispatched "
                 "'%s' — the two target sets are not the same build.\n",
                 here, corvus::active_target());
    return 2;
  }

  const char* path =
      argc > 1 ? argv[1] : "tests/data/dd_special_reference.txt";
  const auto rows = corvus_test::LoadRef(path, 3, 1500);

  std::vector<double> u, want_hi, want_lo;
  for (const auto& row : rows) {
    u.push_back(corvus_test::ParseDouble(row.tok[0], path, row.line));
    want_hi.push_back(corvus_test::ParseDouble(row.tok[1], path, row.line));
    want_lo.push_back(corvus_test::ParseDouble(row.tok[2], path, row.line));
  }

  std::vector<double> hi(u.size()), lo(u.size()), ulo(u.size(), 0.0);
  // Two calls, not one: the second covers only the last 3 rows, a length
  // below every lane count, so the masked LoadN/StoreN tail path inside
  // Log1pmxTestImpl runs on every tier. Without this split the reference
  // set's length when #14 added it (3416, a multiple of 8) meant the tail
  // branch never ran on the AVX-512 tier at all (#14 N4).
  const size_t n = u.size();
  const size_t split = n - 3;
  corvus::Log1pmxTest(u.data(), ulo.data(), hi.data(), lo.data(), split);
  corvus::Log1pmxTest(u.data() + split, ulo.data() + split, hi.data() + split,
                      lo.data() + split, n - split);

  double worst_rel = 0.0, worst_rel_u = 0.0;
  double worst_abs = 0.0, worst_abs_u = 0.0;
  double tiny_rel = 0.0;
  size_t n_big = 0, n_small = 0;

  for (size_t i = 0; i < u.size(); ++i) {
    // Both differences are exact (each pair agrees with the other to within
    // an ulp of its own magnitude), so this measures the kernel below the
    // last bit of a double rather than at it.
    const double err = (hi[i] - want_hi[i]) + (lo[i] - want_lo[i]);
    const double phi = want_hi[i];
    const double rel = phi != 0.0 ? std::abs(err / phi) : 0.0;
    if (std::abs(u[i]) >= kSplit) {
      ++n_big;
      if (rel > worst_rel) {
        worst_rel = rel;
        worst_rel_u = u[i];
      }
    } else {
      ++n_small;
      if (std::abs(err) > worst_abs) {
        worst_abs = std::abs(err);
        worst_abs_u = u[i];
      }
      if (rel > tiny_rel) tiny_rel = rel;
    }
  }

  int rc = 0;
  const double rel_bits = worst_rel > 0.0 ? -std::log2(worst_rel) : 1000.0;
  std::printf(
      "Log1pmxDd |u|>=2^-40  n=%6zu  worst rel=2^-%.2f (gate 2^-%.2f)  "
      "worst u=%.17g\n",
      n_big, rel_bits, kMinRelLog1pmx, worst_rel_u);
  if (rel_bits < kMinRelLog1pmx) {
    std::fprintf(stderr, "FAIL: Log1pmxDd relative error exceeds gate\n");
    rc = 1;
  }
  const double abs_bits = worst_abs > 0.0 ? -std::log2(worst_abs) : 1000.0;
  const double tiny_bits = tiny_rel > 0.0 ? -std::log2(tiny_rel) : 1000.0;
  std::printf(
      "Log1pmxDd |u| <2^-40  n=%6zu  worst abs=2^-%.2f (gate 2^-%.2f)  "
      "worst u=%.17g  [rel there 2^-%.2f]\n",
      n_small, abs_bits, -std::log2(kMaxAbs), worst_abs_u, tiny_bits);
  if (worst_abs > kMaxAbs) {
    std::fprintf(stderr, "FAIL: Log1pmxDd absolute error exceeds gate\n");
    rc = 1;
  }

  // #58: nonzero u.lo (the correction term and the 2-arg path), and Expm1Dd.
  const char* ulo_path =
      argc > 2 ? argv[2] : "tests/data/log1pmx_ulo_reference.txt";
  const char* expm1_path =
      argc > 3 ? argv[3] : "tests/data/expm1_dd_reference.txt";
  size_t n_ulo = 0, n_e = 0;
  double w_ulo = 0.0, w_e = 0.0;
  const double ulo_bits =
      WorstRelBits(ulo_path, 1500, corvus::Log1pmxTest, &n_ulo, &w_ulo);
  std::printf(
      "Log1pmxDd u.lo != 0   n=%6zu  worst rel=2^-%.2f (gate 2^-%.2f)  "
      "worst u=%.17g\n",
      n_ulo, ulo_bits, kMinRelLog1pmx, w_ulo);
  if (ulo_bits < kMinRelLog1pmx) {
    std::fprintf(stderr, "FAIL: Log1pmxDd (u.lo != 0) exceeds gate\n");
    rc = 1;
  }
  const double e_bits =
      WorstRelBits(expm1_path, 2500, corvus::Expm1Test, &n_e, &w_e);
  std::printf(
      "Expm1Dd               n=%6zu  worst rel=2^-%.2f (gate 2^-%.2f)  "
      "worst w=%.17g\n",
      n_e, e_bits, kMinRelExpm1, w_e);
  if (e_bits < kMinRelExpm1) {
    std::fprintf(stderr, "FAIL: Expm1Dd relative error exceeds gate\n");
    rc = 1;
  }

  if (rc == 0) std::printf("PASS: all gates met\n");
  return rc;
}

#endif  // HWY_ONCE
