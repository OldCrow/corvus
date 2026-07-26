// Inverse error function and its complement. Per-target include guard
// (Highway -inl.h idiom).
//
// Both public functions route onto two shared cores, every routed argument
// exact by Sterbenz (PLAN.md, "Phase C part 1 -- erfinv/erfcinv design"):
//   erfinv:  |y| <= 1/2       -> C(y)
//            1/2 < |y| < 1    -> sign(y)*T(1 - |y|)
//   erfcinv: z in [1/2, 3/2]  -> C(1 - z)
//            z < 1/2          -> T(z)
//            z > 3/2          -> -T(2 - z)
// T(s) solves erfc(x) = s for s in (0, 1/2), returning the positive root
// (x in (~0.4769, ~27.217)); C(y) is the odd polynomial erfinv itself uses
// near 0. Both cores are shared verbatim between the two public functions.
//
// CORE C: x = y*Pc(y^2), Pc a single Chebyshev fit on v = y^2 in [0, 1/4],
// with 3 dd LEADING coefficients (the low-degree, dominant terms) and a
// plain-double tail (see src/erfinv_data.h / tools/gen_erfinv_data.py) --
// the same lead+tail split lgamma's zone polynomials use, for the same
// reason: the tail's rounding is attenuated by the v^3 it is multiplied by.
// NO Newton/Halley step here: the central condition number of x = erfinv(y)
// in y is ~1.0-1.2, so refining against corvus's own (already ~1-ULP) erf
// would pass that error straight through with nothing to gain, and using
// the DOUBLE-ROUNDED public erf/erfc would floor the step at ~1 ULP before
// it could even do that. Direct fit is both cheaper and more accurate.
// erfinv(+/-0) = +/-0 falls out of the odd form; subnormal y underflows y^2
// to 0, where Pc(0) = kErfinvCLeadHi[0] = erfinv'(0) = 2/sqrt(pi) is exact.
//
// CORE T: primary variable w = -log(s) (LogDdAny -- s can be subnormal,
// reachable only through erfcinv), t = sqrt(w.hi). A cheap seed
// x0 = t*Seed(t) (2^-19 relative; three intervals selected by t, coefficient
// select + one Horner pass, exactly the erfc tail's 3-interval pattern) is
// refined by ONE dd Halley step, split by the SAME threshold erfc's own
// core/tail split uses (t >= kErfinvTFar <=> x >= 6):
//
//   MID (t < kErfinvTFar): residual space, f(x) = erfc(x) - s. Uses
//   ErfcCoreDd (src/erfc_core-inl.h) -- the erfc kernel's own pre-rounding
//   compensated dd pair -- rather than the public erfc, so the ~2^-13
//   relative cancellation between erfc(x0) and s near the root is absorbed
//   by a dd subtraction (~2^-77 to spare) instead of being floored by a
//   rounding that already happened. f is then reduced to one double (a
//   single rounding is more than enough precision for a term that only
//   needs ~2^-40 relative accuracy to correct a 2^-19-accurate seed --
//   see the Halley derivation below). f' = erfc'(x0) = -2/sqrt(pi)*e^{-x0^2}
//   via the backend Exp: NOT accuracy-critical here, its few-ULP error is
//   attenuated by the small step to a few 2^-69ths.
//
//   FAR (t >= kErfinvTFar): log space, F(x) = log(erfc(x)) - log(s), using
//   corvus's own tail model erfc(x) = e^{-x^2}*G(1/x)/x (REUSING
//   erfc_tail_data.h via ErfcTailGFromU -- no new tail fit needed):
//     F = w - x0^2 - log(x0/G(1/x0))
//   x0^2 is split exactly via ops::ProdLow (x0 <= 28, nowhere near the
//   2^996 Dekker-overflow precondition); the log term only needs
//   double-accuracy ABSOLUTE (LogDd's hi word alone suffices -- see the
//   budget note in gen_erfinv_data.py). This form needs NO exponential at
//   all: F' = erfc'/erfc rewritten through the same G model is
//   F' = -2*x0/(sqrt(pi)*G), where the e^{-x^2} factor cancels exactly.
//   Log space is also the ONLY viable form for subnormal s: erfc(x0) - s in
//   residual space would underflow to zero long before x reaches its true
//   root, while w = -log(s) stays finite and well-scaled.
//
// HALLEY, both branches. The general formula for a root of f is
//   x1 = x0 - 2*f*f' / (2*f'^2 - f*f'')
// and both branches' f''/F'' reduce it to a form with one fewer op:
//   mid: f'' = -2*x0*f'                       => delta = -f / (f' + x0*f)
//   far: F'' = -2*x0*F' - F'^2                 => delta = -2*F / (2*F' + F*(2*x0 + F'))
// (Derivations: mid's f = erfc(x)-s has f'' = erfc''(x) = -2x*erfc'(x) =
// -2x*f' directly. Far's F = log(erfc(x))-log(s) has F' = erfc'/erfc and
// F'' = (erfc''*erfc - erfc'^2)/erfc^2 = erfc''/erfc - F'^2 = -2x*F' - F'^2,
// using erfc'' = -2x*erfc' again.) x1 = fl(x0 + delta) is the kernel's one
// unavoidable rounding; tools/gen_erfinv_data.py's self-check verifies the
// PRE-rounding value is within 2^-56 of the true root over dense samples
// including both region boundaries and the full subnormal-s range, so this
// rounding is the only source of the measured ULP bound.
#if defined(CORVUS_ERFINV_INL_H_) == defined(HWY_TARGET_TOGGLE)
#ifdef CORVUS_ERFINV_INL_H_
#undef CORVUS_ERFINV_INL_H_
#else
#define CORVUS_ERFINV_INL_H_
#endif

#include <limits>

#include "src/dd-inl.h"
#include "src/erfc_core-inl.h"
#include "src/erfinv_data.h"
#include "src/log_dd-inl.h"
#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {

namespace op = ops;

// 2/sqrt(pi) and sqrt(pi): the erfc derivative constant and its reciprocal
// partner used by the Halley derivatives below. Not stored in erfinv_data.h
// since they are closed-form, not fitted.
inline constexpr double kErfinvTwoOverSqrtPi = 0x1.20dd750429b6dp+0;  // 2/sqrt(pi)
inline constexpr double kErfinvSqrtPi = 0x1.c5bf891b4ef6bp+0;        // sqrt(pi)

// Pc(v) = L0 + v*(L1 + v*(L2 + v*S(v))), v = y^2, as a dd. See the file
// header for why the leading terms are dd and the tail is plain double.
template <class D>
HWY_INLINE Dd<D> ErfinvCentralDd(D d, op::V<D> v) {
  const auto* co = detail::kErfinvCCoef;
  auto s = op::Set(d, co[detail::kErfinvCNCoef - 1]);
  for (int k = detail::kErfinvCNCoef - 2; k >= 0; --k) {
    s = op::MulAdd(s, v, op::Set(d, co[k]));
  }

  const auto* lh = detail::kErfinvCLeadHi;
  const auto* ll = detail::kErfinvCLeadLo;
  constexpr int kLead = detail::kErfinvCLead;
  Dd<D> acc{op::Set(d, lh[kLead - 1]), op::Set(d, ll[kLead - 1])};
  acc = DdAddD(d, acc, op::Mul(s, v));
  for (int k = kLead - 2; k >= 0; --k) {
    acc = DdAdd(d, Dd<D>{op::Set(d, lh[k]), op::Set(d, ll[k])},
                DdMulD(d, acc, v));
  }
  return acc;
}

// x = y*Pc(y^2), rounded once. The explicit CopySign is load-bearing, not
// decorative: erfinv is odd on this whole domain so sign(result) == sign(y)
// mathematically, but at y = +/-0 the dd assembly's internal Fast2Sum adds
// a +0 and a -0 partial sum together -- and IEEE 754 defines (-0)+(+0) as
// +0 in round-to-nearest, which silently turns erfinv(-0) into +0 without
// this. CopySign is a no-op everywhere else since the sign already matches.
template <class D>
HWY_INLINE op::V<D> ErfinvCentral(D d, op::V<D> y) {
  const auto v = op::Mul(y, y);
  const auto pc = ErfinvCentralDd(d, v);
  return op::CopySign(DdToDouble(DdMulD(d, pc, y)), y);
}

// mid: delta = -f / (f' + x0*f), f = erfc(x0) - s, f' = erfc'(x0).
template <class D>
HWY_INLINE op::V<D> HalleyMid(D d, op::V<D> x0, op::V<D> s) {
  // Discarded lanes can carry a NaN x0 (erfinv(NaN)/erfcinv(z < 0) route
  // s <= 0 here, and sqrt(-log s) goes NaN). ErfcCoreDd's table index is
  // round(ac*256), NOT masked, so NaN must be scrubbed before the gather --
  // the same precondition erfc.cpp satisfies with its nan mask. Without
  // this, a NaN reaches the index or not purely by platform accident
  // (x86 minpd returns the non-NaN second operand; ARM fcvtzs(NaN) = 0):
  // neither is a guarantee, and Highway's debug-mode gather bounds assert
  // would trip. The scrubbed lanes' results are garbage and discarded.
  const auto x0s = op::IfThenElse(op::IsNaN(x0), op::Zero(d), x0);
  const auto e = ErfcCoreDd(d, x0s, x0s);  // x0 > 0 on live lanes
  const auto f = DdToDouble(DdAddD(d, e, op::Neg(s)));
  const auto fp = op::Mul(op::Set(d, -kErfinvTwoOverSqrtPi),
                          op::Exp(d, op::Neg(op::Mul(x0, x0))));
  return op::Neg(op::Div(f, op::MulAdd(x0, f, fp)));
}

// far: delta = -2F / (2F' + F*(2x0 + F')), F = w - x0^2 - log(x0/G(1/x0)).
template <class D>
HWY_INLINE op::V<D> HalleyFar(D d, op::V<D> x0, Dd<D> wneg) {
  const auto ssq = op::Mul(x0, x0);
  const auto slo = op::ProdLow(d, x0, x0, ssq);  // exact: x0 <= ~28 << 2^996

  const auto u = op::Div(op::Set(d, 1.0), x0);   // double accuracy suffices
  const auto g = ErfcTailGFromU(d, x0, u);
  const auto log_term = op::Sub(LogDd(d, x0).hi, LogDd(d, g).hi);

  auto fdd = DdAdd(d, wneg, Dd<D>{op::Neg(ssq), op::Neg(slo)});
  fdd = DdAddD(d, fdd, op::Neg(log_term));
  const auto f = DdToDouble(fdd);

  const auto fp = op::Neg(op::Div(op::Mul(op::Set(d, 2.0), x0),
                                  op::Mul(op::Set(d, kErfinvSqrtPi), g)));
  const auto denom = op::MulAdd(f, op::MulAdd(op::Set(d, 2.0), x0, fp),
                                op::Mul(op::Set(d, 2.0), fp));
  return op::Neg(op::Div(op::Mul(op::Set(d, 2.0), f), denom));
}

// T(s): the positive root of erfc(x) = s, for s in (0, 1/2). See the file
// header for the seed + Halley derivation. Not domain-checked: callers only
// ever route a valid s here on the SELECTED lane; other lanes may carry
// garbage (e.g. s outside (0, 1/2) from a discarded branch of an outer
// IfThenElse) that must not trap -- and IEEE arithmetic does not trap on
// log/sqrt of a non-positive value, it quietly produces NaN, which the
// caller's select discards.
template <class D>
HWY_INLINE op::V<D> ErfcInvCore(D d, op::V<D> s) {
  const auto w = LogDdAny(d, s);
  const Dd<D> wneg{op::Neg(w.hi), op::Neg(w.lo)};  // w = -log(s) > 0
  const auto t = op::Sqrt(wneg.hi);

  const auto m1 = op::Lt(t, op::Set(d, detail::kErfinvTSplit));  // t < 2
  const auto m2 = op::Lt(t, op::Set(d, detail::kErfinvTFar));    // mid overall
  const auto far = op::Ge(t, op::Set(d, detail::kErfinvTFar));

  const auto u = op::Div(op::Set(d, 1.0), t);
  const auto var = op::IfThenElse(far, u, t);

  const auto* c = detail::kErfinvSeedCoef;
  constexpr int kN = detail::kErfinvSeedNCoef;
  auto seed = Sel3(d, m1, m2, c[0][kN - 1], c[1][kN - 1], c[2][kN - 1]);
  for (int k = kN - 2; k >= 0; --k) {
    seed = op::MulAdd(seed, var, Sel3(d, m1, m2, c[0][k], c[1][k], c[2][k]));
  }
  const auto x0 = op::Mul(t, seed);

  op::V<D> x1;
  if (op::AllFalse(d, far)) {
    x1 = op::Add(x0, HalleyMid(d, x0, s));
  } else if (op::AllTrue(d, far)) {
    x1 = op::Add(x0, HalleyFar(d, x0, wneg));
  } else {
    x1 = op::Add(x0, op::IfThenElse(far, HalleyFar(d, x0, wneg),
                                    HalleyMid(d, x0, s)));
  }
  return x1;
}

template <class D>
HWY_INLINE op::V<D> ErfinvVec(D d, op::V<D> y) {
  const auto ay = op::Abs(y);
  const auto nan = op::IsNaN(y);
  const auto one = op::Set(d, 1.0);
  const auto inf = op::Set(d, std::numeric_limits<double>::infinity());
  const auto nan_out = op::Set(d, std::numeric_limits<double>::quiet_NaN());

  const auto c_m = op::Ge(op::Set(d, detail::kErfinvYSplit), ay);  // |y|<=1/2

  op::V<D> res;
  if (op::AllTrue(d, c_m)) {
    res = ErfinvCentral(d, y);
  } else if (op::AllFalse(d, c_m)) {
    const auto s = op::Sub(one, ay);  // exact by Sterbenz, ay in [1/2, inf)
    res = op::CopySign(ErfcInvCore(d, s), y);
  } else {
    const auto s = op::Sub(one, ay);
    res = op::IfThenElse(c_m, ErfinvCentral(d, y),
                         op::CopySign(ErfcInvCore(d, s), y));
  }

  res = op::IfThenElse(op::Eq(ay, one), op::CopySign(inf, y), res);
  res = op::IfThenElse(op::Gt(ay, one), nan_out, res);
  res = op::IfThenElse(nan, y, res);  // propagate NaN (payload preserved)
  return res;
}

// z < 1/2 -> T(z); z > 3/2 -> -T(2 - z), 2 - z exact by Sterbenz there.
// Computes ONE selected argument and calls ErfcInvCore once, rather than
// evaluating both branches -- unlike the C/T outer split, there is no
// "erfc tail is cheap to skip" asymmetry to trade against a second
// evaluation, and s is exact either way so the select costs nothing extra.
template <class D>
HWY_INLINE op::V<D> ErfcinvTail(D d, op::V<D> z) {
  const auto lo_m = op::Lt(z, op::Set(d, 0.5));
  const auto s = op::IfThenElse(lo_m, z, op::Sub(op::Set(d, 2.0), z));
  const auto x = ErfcInvCore(d, s);
  return op::IfThenElse(lo_m, x, op::Neg(x));
}

template <class D>
HWY_INLINE op::V<D> ErfcinvVec(D d, op::V<D> z) {
  const auto one = op::Set(d, 1.0);
  const auto two = op::Set(d, 2.0);
  const auto half = op::Set(d, 0.5);
  const auto inf = op::Set(d, std::numeric_limits<double>::infinity());
  const auto nan_out = op::Set(d, std::numeric_limits<double>::quiet_NaN());
  const auto nan_in = op::IsNaN(z);

  // |z - 1| <= 1/2 <=> z in [1/2, 3/2]: routes to C(1 - z), 1 - z exact by
  // Sterbenz there. Written this way (rather than two chained comparisons
  // ANDed together) because the ops facade has no mask-AND primitive. The
  // ARGUMENT to C is computed as the DIRECT Sub(one, z), not Neg(z - one):
  // at z = 1 those are NOT bit-identical -- Sub(one, z) is +0 (x - x is
  // always +0 in round-to-nearest), while Neg(z - one) is -0. The design's
  // "erfcinv(1) = +0 falls out via C(1-1)" claim depends on using the form
  // that actually produces +0.
  const auto zm1 = op::Sub(z, one);
  const auto c_m = op::Ge(half, op::Abs(zm1));
  const auto one_minus_z = op::Sub(one, z);

  op::V<D> res;
  if (op::AllTrue(d, c_m)) {
    res = ErfinvCentral(d, one_minus_z);
  } else if (op::AllFalse(d, c_m)) {
    res = ErfcinvTail(d, z);
  } else {
    res = op::IfThenElse(c_m, ErfinvCentral(d, one_minus_z), ErfcinvTail(d, z));
  }

  // z = 0 and z = 2 are exact. The out-of-domain check uses two DIRECT
  // comparisons against z itself, rather than the |z - 1| > 1 trick used
  // for c_m: for z a tiny negative subnormal, z - 1 rounds to exactly -1
  // (the deviation is many orders of magnitude below 1's ulp), which would
  // silently swallow the very case this check exists to catch.
  res = op::IfThenElse(op::Eq(z, op::Zero(d)), inf, res);
  res = op::IfThenElse(op::Eq(z, two), op::Neg(inf), res);
  res = op::IfThenElse(op::Lt(z, op::Zero(d)), nan_out, res);
  res = op::IfThenElse(op::Gt(z, two), nan_out, res);
  res = op::IfThenElse(nan_in, z, res);  // propagate NaN (payload preserved)
  return res;
}

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#endif  // include guard
