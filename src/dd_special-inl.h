// Shared dd special-function primitives: phi(u) = u - log1p(u) and
// e^w - 1, both to dd precision relative to the result. Per-target include
// guard (Highway -inl.h idiom).
//
// These sit one layer above the dd cores (they consume LogDdAny and ExpDd)
// and one layer below the function families: the incomplete gamma's Temme
// ridge rests on Log1pmxDd, its R4 corner on Expm1Dd, and the incomplete
// beta reuses both -- which is why they live here rather than in any one
// family's TU. Byte-identity of every consumer's ULP tables across a hoist
// into this file is the regression guard for such moves.
// Tables: src/dd_special_data.h
// (tools/gen_dd_special_data.py). Accuracy gate: tests/test_dd_special.cpp.
#if defined(CORVUS_DD_SPECIAL_INL_H_) == defined(HWY_TARGET_TOGGLE)
#ifdef CORVUS_DD_SPECIAL_INL_H_
#undef CORVUS_DD_SPECIAL_INL_H_
#else
#define CORVUS_DD_SPECIAL_INL_H_
#endif

#include "src/dd-inl.h"
#include "src/dd_special_data.h"
#include "src/exp_dd-inl.h"
#include "src/log_dd-inl.h"
#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {

namespace op = ops;

// Expm1Dd series/exponential crossover (#58: was 2^-10). Above it e^w (-) 1
// inherits ExpDd's ~2^-70 relative error on e^w -- not the dd's 2^-105:
// ExpDd's p2 = rh^2*q is a plain double -- amplified by 1/|w|, so at the
// old 2^-10 cut it reached 2^-61.4 (measured). Moving the cut to 2^-6 caps
// the amplification at 2^6 (~2^-64 just above the cut); below it the series
// runs to w^11/11!, whose truncation w^11/12! is <= 2^-94.8 relative at the
// cut, and whose plain-double tail from w^3 rounds at |w^2/6| * 2^-53 <=
// 2^-67.6 relative.
inline constexpr double kExpm1Cut = 0x1.0p-6;

// Log1pmxDd's series: how many leading terms of T get dd treatment. A
// rounding introduced at Horner step k reaches T attenuated by u^k, so with
// |u| <= 1/16 the plain-double tail from k = 6 contributes
// 2^-53 * u^6 / (T * 8) ~ 2^-79 relative -- comfortably under the 2^-74
// series truncation the generator's self-check (a) certifies.
inline constexpr int kPhiLead = 6;

// The leading six coefficients are carried as dd pairs
// (detail::kPhiCoefLo holds the low words), and this is not
// belt-and-braces: rounding 1/3 to a double is 2^-55.9 absolute, which is
// 2^-58.5 RELATIVE on phi at the cut -- fine for phi, not for a*phi, the
// argument of an exponential (a relative eps on phi is an a*phi*eps
// relative error on the result, and a*phi reaches ~740 with |u| still at
// the cut when a ~ 4e5, where an all-double table costs double-digit
// ULP). The dd leads put the coefficient
// contribution at 2^-79, leaving the 2^-75 truncation dominant -- the same
// lead/tail split as lgamma's zone polynomials, for the same reason. The
// generator self-checks each pair against the exact rational.

// ------------------------------------------------------------------------
// phi(u) = u - log1p(u), to dd precision RELATIVE to phi.
//
// This is the shared primitive the whole ridge rests on (the incomplete beta
// reuses it). Note what it is NOT: u (-) LogDd(fl(1+u)). For small u,
// phi ~ u^2/2 while log1p(u) ~ u, so that spelling amplifies log_dd's ~2^-68
// by 2/u -- at u = 2^-20 that is 2^-49 relative, and multiplied by an a of
// 1e6 it is a wrong answer, not a rounding.
//
// SMALL |u| (<= kPhiCut = 1/16): phi = u^2 * T(u) with
// T(u) = sum_k (-1)^k u^k / (k+2), the series in src/dd_special_data.h. No
// cancellation anywhere -- u^2 is formed exactly (TwoProd) and T is O(1).
// Truncation at 18 terms is u^18/10 <= 2^-75 relative at the cut.
//
// The u.lo handling is the subtle part. Evaluating the series at u.hi and
// stopping there costs (2/3)*|u.lo| <= 2^-58 relative at the cut -- fine for
// phi itself, fatal once multiplied by an a large enough to make a*phi ~ 700.
// The fix is one first-order correction, and it is written on phi rather
// than on T because phi' = u/(1+u) is closed form:
//     phi(u.hi + u.lo) = phi(u.hi) + u.lo * u.hi/(1 + u.hi) + O(u.lo^2).
// The neglected term is (u.lo^2/2)/(1+u)^2, which against phi ~ u^2/2 is
// u.lo^2/u^2 <= 2^-106 relative -- uniformly, at every u.
//
// LARGE |u|: no cancellation to protect against (phi >= 2^-9 while u <= 1),
// so the direct difference is used, with 1 + u handed to LogDdAny as a dd.
// At the cut this inherits log_dd's error amplified by 2/u; log_dd near
// argument 1 is better than its ~2^-68 global bound (the table's L_j is
// exact to 2^-107 and log1p(r) is relative to itself), but not enough to
// make the seam continuous in accuracy. Measured worst for the whole
// primitive: ~2^-67.6 relative, just above the cut (2^-67.78 at
// u = 0.0625098 on the reference rows, 2^-67.63 near u = 0.0625809 in
// 100k-point probes; tests/test_dd_special.cpp, row 12). The 2^-71.41
// recorded at #58 came from rows that never sampled (1/16, 0.064]. Kept as
// is (#68, closed): every consumer takes phi into an exponent capped near
// 800 or a term <= ~8 (LgammaDiffDd, m <= 256), so the seam is <= ~0.03 ULP
// downstream. Measured with a shadow build whose cut sat at 9/128: of
// 18,000 consumer outputs placed in the band (gamma, gammainv, beta,
// betainv), one changed, by a 1-ulp near-tie flip.
//
// TWO SPELLINGS OF 1 + u. The 2-arg overload takes w = 1 + u from the
// CALLER, for use where the caller owns an exact closed form (beta's
// BetaPsiCore: 1 + u = c*xi/alpha with no subtraction anywhere). The 1-arg
// form derives w from u itself -- TwoSum(1, u.hi) plus u.lo folded into the
// low word -- which is exact as an unevaluated sum, BUT NOT a normalized dd
// when u is near -1: for 1 + u ~ 2^-53 the folded u.lo is comparable to (or
// larger than) the high word, and LogDdAny's expansion in lo/hi drops the
// cubic (~1e-4-class error in a consumer's E at u.lo/w.hi ~ 1); at u.hi == -1
// exactly the high word is zero and LogDdAny has no valid path at all.
// HAZARD RULE: a caller whose u can approach -1 closer than ~2^-8 must use
// the 2-arg overload with an independently computed w. Every 1-arg call
// site today satisfies this (gamma's Temme band has u in [-1/2, 1], beta's
// R3 ratio band likewise, Log1pDdWide's arguments are >= -1 + 2^-12-class).
//
// FLOOR: |u| < 2^-484. The t*uh steps lose <= 2^-1073 against T >= 0.47.
// u^2 and u^2*T leave phi off by <= min(2^-1072, 1.1 u^2) ABSOLUTE (not
// relative to phi), which each caller absorbs:
//  * Log1pDdWide (w - phi): <= 2^-535 |w|.
//  * LgammaDiffDd (z*phi(m/z)): <= 1.1 m*w, under 2^-483 of m ln z (z >= Z0).
//  * BetaPsiCore (p*phi): |u| < 2^-484 needs y < 2^-326 (lambda is 0 or
//    >= 2^-157 alpha*y; x + y = 1), so alpha/beta > 2^326 and the other
//    term dominates cs: error <= 2^-324 cs (mirrored for v).
//  * gamma Temme, gammainv forward: u = x/a - 1 of doubles is 0 or >= 2^-54.
//  * BetaInvPhiOverSq: S1 seed only; BetaInvConsider scores it on the forward
//    and Newton polishes, so the error cannot reach a result.
//
// The block above documents Log1pmxDd as a whole. Its two branches are
// split into Log1pmxDdSeries (|u| <= kPhiCut) and Log1pmxDdDirect (#45), so
// a consumer that knows its lanes' band can skip the branch it does not
// need; Log1pmxDd computes both and selects, with the same bits as before.
template <class D>
CORVUS_OUTLINE Dd<D> Log1pmxDdSeries(D d, Dd<D> u) {
  const auto one = op::Set(d, 1.0);
  const auto uh = u.hi;
  constexpr int kN =
      static_cast<int>(sizeof(detail::kPhiCoef) / sizeof(double));
  const auto* c = detail::kPhiCoef;

  auto s = op::Set(d, c[kN - 1]);
  for (int k = kN - 2; k >= kPhiLead; --k) {
    s = op::MulAdd(s, uh, op::Set(d, c[k]));
  }
  Dd<D> t{op::Set(d, c[kPhiLead - 1]),
          op::Set(d, detail::kPhiCoefLo[kPhiLead - 1])};
  t = DdAddD(d, t, op::Mul(s, uh));
  for (int k = kPhiLead - 2; k >= 0; --k) {
    t = DdAdd(d, DdMulD(d, t, uh),
              Dd<D>{op::Set(d, c[k]), op::Set(d, detail::kPhiCoefLo[k])});
  }
  auto ser = DdMul(d, TwoProd(d, uh, uh), t);
  return DdAddD(d, ser, op::Mul(u.lo, op::Div(uh, op::Add(one, uh))));
}

// phi(u) = u - log1p(u) by the direct difference, for |u| > kPhiCut, with
// w = 1 + u supplied as a dd (see TWO SPELLINGS OF 1 + u below).
template <class D>
HWY_INLINE Dd<D> Log1pmxDdDirect(D d, Dd<D> u, Dd<D> w) {
  const auto lg = LogDdAny(d, w);
  return DdAdd(d, u, Dd<D>{op::Neg(lg.hi), op::Neg(lg.lo)});
}

template <class D>
HWY_INLINE Dd<D> Log1pmxDd(D d, Dd<D> u, Dd<D> w) {
  const auto ser = Log1pmxDdSeries(d, u);
  const auto big = Log1pmxDdDirect(d, u, w);

  const auto m = op::Ge(op::Set(d, detail::kPhiCut), op::Abs(u.hi));
  CORVUS_ROUTE(d, "dd_special/log1pmx_series", m);
  CORVUS_ROUTE(d, "dd_special/log1pmx_direct", op::MaskNot(m));
  return Dd<D>{op::IfThenElse(m, ser.hi, big.hi),
               op::IfThenElse(m, ser.lo, big.lo)};
}

template <class D>
HWY_INLINE Dd<D> Log1pmxDd(D d, Dd<D> u) {
  auto w = TwoSum(d, op::Set(d, 1.0), u.hi);  // exact
  w.lo = op::Add(w.lo, u.lo);
  return Log1pmxDd(d, u, w);
}

// e^w - 1 for a dd w, relative to the result: ~2^-64 worst (just above
// kExpm1Cut, where ExpDd's own error is amplified), ~2^-67.6 in the series
// -- NOT the dd's 2^-105 (#58; tests/test_dd_special.cpp gates it).
//
// Small |w| is the only interesting case: e^w (-) 1 there is a subtraction
// of two numbers within |w| of each other, which throws away log2(1/|w|)
// bits of the exponential's accuracy. The series keeps them. Structured
// like Log1pDd: w and w^2/2 in dd (w^2/2 is 2^-7 relative to w at the cut,
// so a rounded w^2 would show), everything from w^3 on in plain double,
// where |w^3/6| is 2^-14.6 relative at the cut and its own rounding lands
// at 2^-67.6.
template <class D>
CORVUS_OUTLINE Dd<D> Expm1DdSeries(D d, Dd<D> w) {
  const auto wh = w.hi;

  // q = 1/6 + w/24 + ... + w^8/11!, so w^3*q runs the series to w^11/11!.
  auto q = op::Set(d, 1.0 / 39916800.0);
  q = op::MulAdd(q, wh, op::Set(d, 1.0 / 3628800.0));
  q = op::MulAdd(q, wh, op::Set(d, 1.0 / 362880.0));
  q = op::MulAdd(q, wh, op::Set(d, 1.0 / 40320.0));
  q = op::MulAdd(q, wh, op::Set(d, 1.0 / 5040.0));
  q = op::MulAdd(q, wh, op::Set(d, 1.0 / 720.0));
  q = op::MulAdd(q, wh, op::Set(d, 1.0 / 120.0));
  q = op::MulAdd(q, wh, op::Set(d, 1.0 / 24.0));
  q = op::MulAdd(q, wh, op::Set(d, 1.0 / 6.0));
  // FLOOR: |wh| < 2^-484 (series lanes); sq.lo is off by <= min(2^-1073,
  // 1.01 wh^2) < 2^-535 |result|, inside the relative contract every caller
  // (gamma R4, beta, betainv) relies on.
  const auto sq = TwoProd(d, wh, wh);
  const auto half = op::Set(d, 0.5);
  // w^2/2 as a dd; halving is exact and the 2*w.hi*w.lo cross term is kept
  // for the same reason Log1pDd keeps it.
  const Dd<D> half_sq{op::Mul(sq.hi, half),
                      op::MulAdd(wh, w.lo, op::Mul(sq.lo, half))};
  auto ser = DdAdd(d, w, half_sq);
  return DdAddD(d, ser, op::Mul(op::Mul(sq.hi, wh), q));
}

// Expm1Dd = Expm1DdSeries below kExpm1Cut, e^w (-) 1 above (the block
// above documents both). The series is split out (#45) so the public
// expm1 can skip whichever branch no lane needs; the bits are unchanged.
template <class D>
CORVUS_OUTLINE Dd<D> Expm1Dd(D d, Dd<D> w) {
  const auto one = op::Set(d, 1.0);
  const auto wh = w.hi;
  const auto ser = Expm1DdSeries(d, w);

  const auto ex = ExpDd(d, wh, w.lo);
  const auto big = DdAddD(d, ex, op::Neg(one));

  const auto m = op::Ge(op::Set(d, kExpm1Cut), op::Abs(wh));
  CORVUS_ROUTE(d, "dd_special/expm1_series", m);
  CORVUS_ROUTE(d, "dd_special/expm1_exp", op::MaskNot(m));
  return Dd<D>{op::IfThenElse(m, ser.hi, big.hi),
               op::IfThenElse(m, ser.lo, big.lo)};
}

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#endif  // include guard
