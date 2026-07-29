// Accuracy gate for Log1pmxDd (src/gamma-inl.h), the shared primitive
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
#include <fstream>
#include <string>
#include <vector>

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "tests/test_gamma_util.cpp"
#include "hwy/foreach_target.h"  // IWYU pragma: keep

#include "src/gamma-inl.h"
#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {

// Same masked-tail discipline as the shipped kernels: one code path.
void Log1pmxTestImpl(const double* u, double* hi, double* lo, size_t n) {
  using DT = op::ScalableTag<double>;
  const DT d;
  const size_t N = op::Lanes(d);
  size_t i = 0;
  for (; i + N <= n; i += N) {
    const auto r = Log1pmxDd(d, Dd<DT>{op::Load(d, u + i), op::Zero(d)});
    op::Store(r.hi, d, hi + i);
    op::Store(r.lo, d, lo + i);
  }
  if (i < n) {
    const size_t m = n - i;
    const auto r = Log1pmxDd(d, Dd<DT>{op::LoadN(d, u + i, m), op::Zero(d)});
    op::StoreN(r.hi, d, hi + i, m);
    op::StoreN(r.lo, d, lo + i, m);
  }
}

const char* UtilTargetNameImpl() { return hwy::TargetName(HWY_TARGET); }

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#if HWY_ONCE

#include "expect_target.h"

namespace corvus {
HWY_EXPORT(Log1pmxTestImpl);
HWY_EXPORT(UtilTargetNameImpl);

// Dispatch from INSIDE namespace corvus: with a single compiled target
// Highway collapses HWY_DYNAMIC_DISPATCH to N_<target>::FUNC, and a
// globally qualified call would then name a namespace that does not exist.
void Log1pmxTest(const double* u, double* hi, double* lo, size_t n) {
  HWY_DYNAMIC_DISPATCH(Log1pmxTestImpl)(u, hi, lo, n);
}
const char* UtilTargetName() {
  return HWY_DYNAMIC_DISPATCH(UtilTargetNameImpl)();
}
}  // namespace corvus

namespace {

// Gates from the design budget (PLAN.md, "Phase C part 2"), not from a
// measurement: the series truncates at 2^-75 and the log branch at the cut
// amplifies log_dd by 2/u. PROVISIONAL until the tier sweep pins them.
constexpr double kMinRelExp = 63.0;   // |u| >= 2^-40, relative
constexpr double kMaxAbs = 0x1.0p-100;  // |u| <  2^-40, absolute
constexpr double kSplit = 0x1.0p-40;

}  // namespace

int main(int argc, char** argv) {
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
      argc > 1 ? argv[1] : "tests/data/gamma_util_reference.txt";
  std::ifstream f(path);
  if (!f) {
    std::fprintf(stderr, "cannot open reference file: %s\n", path);
    return 2;
  }

  std::vector<double> u, want_hi, want_lo;
  std::string a, b, c;
  while (f >> a >> b >> c) {
    u.push_back(std::strtod(a.c_str(), nullptr));
    want_hi.push_back(std::strtod(b.c_str(), nullptr));
    want_lo.push_back(std::strtod(c.c_str(), nullptr));
  }
  if (u.size() < 1500) {
    std::fprintf(stderr, "reference file suspiciously small: %zu lines\n",
                 u.size());
    return 2;
  }

  std::vector<double> hi(u.size()), lo(u.size());
  corvus::Log1pmxTest(u.data(), hi.data(), lo.data(), u.size());

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
      n_big, rel_bits, kMinRelExp, worst_rel_u);
  if (rel_bits < kMinRelExp) {
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

  if (rc == 0) std::printf("PASS: all gates met\n");
  return rc;
}

#endif  // HWY_ONCE
