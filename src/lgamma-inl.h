// log |Gamma(x)| over the whole real axis. Per-target include guard (Highway
// -inl.h idiom). The kernel is split out of lgamma.cpp because the incomplete
// gamma/beta prefactor exp(a*log x - x - lgamma(a)) will consume it directly,
// the same way erfc consumes erf_core-inl.h.
//
// EVERYTHING BELOW IS ACCUMULATED IN DOUBLE-DOUBLE AND ROUNDED EXACTLY ONCE.
// That is not belt-and-braces: three of the five regions end in a subtraction
// of two same-signed quantities, and carrying ~106 bits is what makes those
// cancellations cost bits instead of digits.
//
// REGIONS (x > 0; the negative axis reflects onto these -- see LgammaVec)
//   (0, 1/2)     lgamma(x) = lgamma(1 + x) - log x, and t = x is EXACT, so
//                this is the centre-1 polynomial at t = x with a log
//                subtracted. Note the shift is never formed: writing y = x + 1
//                and then t = y - 1 would round y and put ~1 ULP of psi(y)
//                into the result near x = 1/2.
//   [1/2, 3/2)   t = x - 1, exact by Sterbenz; centre-1 polynomial.
//   [3/2, 5/2]   t = x - 2, exact by Sterbenz; centre-2 polynomial.
//   (5/2, X0)    walk down by Gamma(z) = Gamma(z-1)*(z-1) until the argument
//                lands in (3/2, 5/2], then the centre-2 polynomial plus
//                log P, P the product of the factors walked past.
//   [X0, inf)    Stirling.
//
// WHY THE RECURRENCE GOES DOWN, NOT UP
//   x - k is EXACTLY representable for every integer k < x < 2^52: x is a
//   multiple of ulp(x) <= 1, and |x - k| < |x| so the result needs no finer
//   grid than x already sits on. x + k is a multiple of the same grid but
//   LARGER, so it can need a bit the format does not have -- x + 8 for
//   x in (5/2, 8) loses up to two. Walking up would therefore have to carry
//   every shifted argument as a dd pair; walking down keeps them all exact,
//   and lands in (3/2, 5/2] where one polynomial serves every lane.
//
// ZEROS
//   lgamma vanishes at x = 1 and x = 2 and nowhere else on the positive axis.
//   Both are reproduced by construction: the region containing them evaluates
//   t*B(t) with t exact, so a zero t gives a zero result no matter what B is,
//   and relative accuracy survives arbitrarily close to them. This is the
//   whole reason the zone polynomials are fitted to lgamma(c + t)/t rather
//   than to lgamma itself.
//
// THE NEGATIVE AXIS
//   The usual reflection log(pi/|sin pi x|) - lgamma(1 - x) is not used as
//   written, for two reasons.
//   * 1 - x is inexact, and lgamma(1 - x) would then need a dd argument
//     threaded through every region above. Gamma(1-x) = -x*Gamma(-x) trades
//     that for one extra log with an argument, -x, that is exact:
//         lgamma(x) = -log|u| - log(sin(pi u)/(pi u)) - log(-x) - lgamma(-x)
//     so the positive-axis pipeline runs once, on |x|, for lanes of either
//     sign, and the reflection is a post-correction.
//   * pi/|sin pi x| overflows as x approaches an integer. Factoring
//     sin(pi u) = pi*u*(sin(pi u)/(pi u)) cancels pi analytically and leaves
//     -log|u| -- which is exactly the term that should diverge at a pole --
//     plus a bounded correction fitted on u^2.
//   u = x - round(x) is exact for every x: for |x| < 1/2 it is x itself, and
//   otherwise |round(x) - x| <= 1/2 <= |x|/2 puts Sterbenz in range.
//   Relative accuracy is NOT claimed near the |Gamma| = 1 crossings, where
//   lgamma passes through zero at points with no closed form; see
//   docs/ACCURACY.md for the measured split.
#if defined(CORVUS_LGAMMA_INL_H_) == defined(HWY_TARGET_TOGGLE)
#ifdef CORVUS_LGAMMA_INL_H_
#undef CORVUS_LGAMMA_INL_H_
#else
#define CORVUS_LGAMMA_INL_H_
#endif

#include <limits>

#include "src/dd-inl.h"
#include "src/lgamma_data.h"
#include "src/log_dd-inl.h"
#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {

namespace op = ops;

// Pick between the two zone coefficient sets. One Horner pass over per-lane
// selected coefficients beats two Horners plus a select on the result: the
// centre-1 fit needs degree 34 and the centre-2 fit 21, so running both would
// cost the longer one twice over.
template <class D, class M>
HWY_INLINE op::V<D> Sel2(D d, M m, double v0, double v1) {
  return op::IfThenElse(m, op::Set(d, v0), op::Set(d, v1));
}

// B(t) = L0 + t*(L1 + t*(L2 + t*S(t))) for the selected centre, as a dd.
// Only S and the single product t*S are evaluated in double; because they
// enter multiplied by t^3 (|t| <= 1/2), their rounding reaches the result
// attenuated ~8x. The generator measures what survives and refuses to emit a
// table that misses budget -- see tools/gen_lgamma_data.py.
template <class D, class M>
HWY_INLINE Dd<D> ZoneBracket(D d, op::V<D> t, M c1) {
  const auto* co = detail::kLgammaZoneCoef;
  auto s = Sel2(d, c1, co[0][detail::kLgammaZoneNCoef - 1],
                co[1][detail::kLgammaZoneNCoef - 1]);
  for (int k = detail::kLgammaZoneNCoef - 2; k >= 0; --k) {
    s = op::MulAdd(s, t, Sel2(d, c1, co[0][k], co[1][k]));
  }

  const auto* lh = detail::kLgammaZoneLeadHi;
  const auto* ll = detail::kLgammaZoneLeadLo;
  Dd<D> acc{Sel2(d, c1, lh[0][2], lh[1][2]), Sel2(d, c1, ll[0][2], ll[1][2])};
  acc = DdAddD(d, acc, op::Mul(s, t));
  for (int k = 1; k >= 0; --k) {
    acc = DdAdd(d, Dd<D>{Sel2(d, c1, lh[0][k], lh[1][k]),
                         Sel2(d, c1, ll[0][k], ll[1][k])},
                DdMulD(d, acc, t));
  }
  return acc;
}

// Stirling for x >= X0. The remainder phi is the only part in plain double:
// it is ~1/(12x), so a full ULP of it is ~2^-60 of lgamma(X0), and X0 is the
// tightest point in the region (lgamma increases from there).
//
// Grouped as x*(log x - 1) - log(x)/2 + ..., NOT as the textbook
// (x - 1/2)*log x - x. The two are algebraically identical but the textbook
// form computes a product that exceeds the result by x, and near the overflow
// threshold (x ~ 2.556e305) that is enough to send the intermediate to
// infinity -- and its dd residual to NaN -- while the true lgamma is still
// finite, over a band ~0.1% wide in x. Grouped this way the intermediate
// exceeds the result by only ~log(x)/2, which is ~351 against an ulp of
// 2e292 at that magnitude: the product and the result overflow at the same
// double. log x - 1 needs no protection of its own, since x >= 8 keeps it
// above 1.
template <class D>
HWY_INLINE Dd<D> LgammaStirling(D d, op::V<D> x) {
  const auto l = LogDd(d, x);

  // The product is formed on x scaled down by 2^200 and scaled back after.
  // ops::ProdLow's non-FMA path is Dekker's split, whose intermediate
  // a*(2^27+1) overflows once |a| exceeds 2^996 ~ 6.7e299 -- and this is the
  // one place in corvus where an operand gets that large. Unguarded, SSE2
  // through SSE4 returned an infinite residual for x above that, while every
  // FMA target was correct; the tier sweep is what surfaced it.
  //
  // Scaling by a power of two is exact and every step below is linear in x,
  // so the scaled and unscaled computations agree bit for bit -- this is not
  // a separate code path for the non-FMA targets, and the results stay
  // identical across tiers. 2^200 sits between the two limits with room at
  // both ends: x >= 8 keeps the scaled product's residual normal (~2^-250 at
  // the low end), and x <= kLgammaMaxArg keeps the split under its ceiling.
  const auto down = op::Set(d, 0x1p-200);
  const auto up = op::Set(d, 0x1p200);
  auto st = DdMulD(d, DdAddD(d, l, op::Set(d, -1.0)), op::Mul(x, down));
  st = Dd<D>{op::Mul(st.hi, up), op::Mul(st.lo, up)};

  st = DdAdd(d, st, Dd<D>{op::Set(d, detail::kLgammaHalfLog2PiHi),
                          op::Set(d, detail::kLgammaHalfLog2PiLo)});
  // -log(x)/2; halving is exact, so the pair stays normalized.
  const auto half = op::Set(d, -0.5);
  st = DdAdd(d, st, Dd<D>{op::Mul(l.hi, half), op::Mul(l.lo, half)});

  const auto w = op::Div(op::Set(d, 1.0), op::Mul(x, x));  // 1/x^2
  auto p = op::Set(d, detail::kLgammaStirCoef[detail::kLgammaStirNCoef - 1]);
  for (int k = detail::kLgammaStirNCoef - 2; k >= 0; --k) {
    p = op::MulAdd(p, w, op::Set(d, detail::kLgammaStirCoef[k]));
  }
  return DdAddD(d, st, op::Div(p, x));
}

// lgamma on (0, X0): the two zone polynomials plus, for the outer bands, one
// logarithm.
//
// The recurrence and the (0, 1/2) shift are folded into a single log call.
// Region (0, 1/2) wants -log x and region (5/2, X0) wants +log P; every other
// lane wants neither, and P is exactly 1 there because no recurrence step
// fires. So one LogDdAny of "x or P" plus a sign select covers all three.
template <class D>
HWY_INLINE Dd<D> LgammaLow(D d, op::V<D> x) {
  const auto one = op::Set(d, 1.0);
  const auto small = op::Lt(x, op::Set(d, detail::kLgammaZoneLo));

  // All-zone fast path. When every lane sits in [kLgammaZoneLo, kLgammaZoneHi]
  // the recurrence multiplies P by exact ones and the log runs on P == 1,
  // whose table slots carry R = 1, L = 0 exactly -- both contribute exactly
  // zero, and DdAdd with an exact zero is a value-preserving renormalization.
  // Skipping them returns the identical rounded result while sparing the
  // vector kLgammaMidSteps dd multiplies and a full LogDdAny. This is the
  // same AllTrue/AllFalse region split LgammaPosDd already performs, one
  // level down; mixed vectors fall through to the full path unchanged.
  if (op::AllFalse(d, small) &&
      op::AllFalse(d, op::Gt(x, op::Set(d, detail::kLgammaZoneHi)))) {
    const auto c1 = op::Lt(x, op::Set(d, detail::kLgammaZoneMid));
    const auto t =
        op::IfThenElse(c1, op::Sub(x, one), op::Sub(x, op::Set(d, 2.0)));
    return DdMulD(d, ZoneBracket(d, t, c1), t);
  }

  // Walk down to (3/2, 5/2]. Step k fires when the argument has not yet
  // arrived, i.e. when x still exceeds kLgammaZoneHi + (k-1). Clamping the
  // driver to X0 keeps the loop finite for the Stirling lanes, whose y and P
  // are discarded -- an unclamped x = 1e308 would otherwise build an infinite
  // P and poison its dd residual with a NaN.
  const auto xr = op::Min(x, op::Set(d, detail::kLgammaX0));
  auto y = x;
  Dd<D> prod{one, op::Zero(d)};
  for (int k = 1; k <= detail::kLgammaMidSteps; ++k) {
    const auto fire =
        op::Gt(xr, op::Set(d, detail::kLgammaZoneHi + (k - 1)));
    // The fire masks shrink monotonically in k, so once no lane fires every
    // remaining step multiplies P by an exact one and leaves y unchanged --
    // breaking out returns the identical result.
    if (op::AllFalse(d, fire)) break;
    const auto step = op::Sub(xr, op::Set(d, static_cast<double>(k)));  // exact
    prod = DdMulD(d, prod, op::IfThenElse(fire, step, one));
    y = op::IfThenElse(fire, step, y);
  }

  // t = x - centre, exact in every band: x for (0, 1/2), x - 1 by Sterbenz on
  // [1/2, 3/2), and y - 2 for the rest (y is exact and y - 2 shrinks it).
  const auto c1 = op::Lt(x, op::Set(d, detail::kLgammaZoneMid));
  const auto t = op::IfThenElse(
      small, x,
      op::IfThenElse(c1, op::Sub(x, one), op::Sub(y, op::Set(d, 2.0))));
  const auto base = DdMulD(d, ZoneBracket(d, t, c1), t);

  // One log for both outer bands: x when shifting, P otherwise (P == 1, hence
  // log == 0, throughout the zone). LogDdAny rather than LogDd because x may
  // be subnormal here.
  const Dd<D> arg{op::IfThenElse(small, x, prod.hi),
                  op::IfThenElse(small, op::Zero(d), prod.lo)};
  const auto lg = LogDdAny(d, arg);
  const Dd<D> corr{op::IfThenElse(small, op::Neg(lg.hi), lg.hi),
                   op::IfThenElse(small, op::Neg(lg.lo), lg.lo)};
  return DdAdd(d, base, corr);
}

// lgamma for x > 0, as a dd. The caller rounds -- the reflection needs the
// full precision, since on the negative axis this value is subtracted from
// three others of comparable size.
template <class D>
HWY_INLINE Dd<D> LgammaPosDd(D d, op::V<D> x) {
  const auto big = op::Ge(x, op::Set(d, detail::kLgammaX0));

  // Both paths are domain-clamped so the unused one cannot fault or produce a
  // NaN that a select would then have to launder: Stirling needs x >= X0 for
  // LogDd's normal-argument contract, LgammaLow needs a bounded x for the
  // recurrence.
  if (op::AllTrue(d, big)) {
    return LgammaStirling(d, x);
  }
  if (op::AllFalse(d, big)) {
    return LgammaLow(d, x);
  }
  const auto st = LgammaStirling(d, op::Max(x, op::Set(d, detail::kLgammaX0)));
  const auto lo = LgammaLow(d, x);
  return Dd<D>{op::IfThenElse(big, st.hi, lo.hi),
               op::IfThenElse(big, st.lo, lo.lo)};
}

// log(sin(pi u)/(pi u)) for |u| <= 1/2, as a dd. Even in u, so it is fitted in
// v = u^2; v is formed exactly (TwoProd) because the leading term is
// -pi^2 v/6, large enough that v's own rounding would show.
template <class D>
HWY_INLINE Dd<D> LogSinc(D d, op::V<D> u) {
  const auto v = TwoProd(d, u, u);
  auto s = op::Set(d, detail::kLgammaSinCoef[detail::kLgammaSinNCoef - 1]);
  for (int k = detail::kLgammaSinNCoef - 2; k >= 0; --k) {
    s = op::MulAdd(s, v.hi, op::Set(d, detail::kLgammaSinCoef[k]));
  }
  Dd<D> acc{op::Set(d, detail::kLgammaSinLeadHi[1]),
            op::Set(d, detail::kLgammaSinLeadLo[1])};
  acc = DdAddD(d, acc, op::Mul(s, v.hi));
  acc = DdAdd(d, Dd<D>{op::Set(d, detail::kLgammaSinLeadHi[0]),
                       op::Set(d, detail::kLgammaSinLeadLo[0])},
              DdMul(d, acc, v));
  return DdMul(d, acc, v);
}

template <class D>
HWY_INLINE op::V<D> LgammaVec(D d, op::V<D> x) {
  const auto zero = op::Zero(d);
  const auto inf = op::Set(d, std::numeric_limits<double>::infinity());
  const auto neg = op::Lt(x, zero);

  // The positive pipeline runs on |x| for every lane: for x < 0 the value it
  // produces, lgamma(-x), is exactly the term the reflection needs. x = 0 is
  // NOT clamped away here -- it walks the ordinary path, where LogDdAny reads
  // a zero exponent field and returns a finite nonsense value that the pole
  // override below replaces. Clamping instead would have to pick a floor, and
  // any floor above the smallest subnormal silently truncates real arguments.
  const auto ax = op::Abs(x);
  const auto g = LgammaPosDd(d, ax);

  Dd<D> res = g;
  if (!op::AllFalse(d, neg)) {
    // u = x - round(x), exact; |u| <= 1/2. u == 0 marks a pole.
    const auto u = op::Sub(x, op::Round(x));
    const auto au = op::Abs(u);
    // lgamma(x) = -log|u| - logsinc(u) - log(-x) - lgamma(-x).
    auto r = DdAdd(d, LogDdAny(d, au), LogSinc(d, u));
    r = DdAdd(d, r, DdAdd(d, LogDdAny(d, ax), g));
    const Dd<D> refl{op::Neg(r.hi), op::Neg(r.lo)};
    const auto pole = op::Eq(u, zero);
    res = Dd<D>{op::IfThenElse(neg, op::IfThenElse(pole, inf, refl.hi), res.hi),
                op::IfThenElse(neg, op::IfThenElse(pole, zero, refl.lo), res.lo)};
  }

  // Adding +0 is the identity on every finite nonzero value and turns -0 into
  // +0. lgamma(1) and lgamma(2) are exactly zero and C99 requires the sign to
  // be positive; t*B(t) at t = 0 produces the sign of B, which is negative at
  // x = 1. One op, and it does not depend on how any target's Fast2Sum happens
  // to sign a zero.
  auto out = op::Add(DdToDouble(res), zero);

  // lgamma(+-0) = +inf, and so is every pole. The magnitude test covers three
  // cases at once and is written on |x| deliberately: +inf and overflow above
  // ~2.556e305, and the whole negative range below -2.556e305, where every
  // double is an integer and therefore a pole (so -inf lands right too). NaN
  // propagates with its payload, matching erf/erfc.
  out = op::IfThenElse(op::Eq(x, zero), inf, out);
  out = op::IfThenElse(op::Gt(ax, op::Set(d, detail::kLgammaMaxArg)), inf, out);
  return op::IfThenElse(op::IsNaN(x), x, out);
}

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#endif  // include guard
