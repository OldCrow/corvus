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
//   (5/2, X0)    Taylor table on a 1/32 grid (see LgammaTable).
//   [X0, inf)    Stirling.
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
//     sign, and the reflection is a post-correction. The two logs are one
//     (#31): -log|u| - log(-x) = -log(|u|*|x|), the product exact by
//     TwoProd. Tiny |x|, where u = x and x^2 could underflow, scales both
//     factors by 2^600 and takes 1200 ln 2 back out.
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

// ------------------------------------------------------------------------
// OUTLINED log: shared HWY_NOINLINE wrapper, see OutlinedLogDd in
// log_dd-inl.h. This file reaches it from three call sites -- the
// zone/Stirling seam in LgammaLow and both reflection-term logs in
// LgammaVec.

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
CORVUS_OUTLINE Dd<D> ZoneBracket(D d, op::V<D> t, M c1) {
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
                // FLOOR: |acc*t| < 2^-968 needs |t| < 2^-964 (|acc| >= 0.067);
                // the <= 2^-1072 loss lands on L1 or L0 (>= 0.32), 2^-1070
                // relative in B for every caller (LgammaLow, R4, gammainv).
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
// threshold (x ~ 2.5599833e305) that is enough to send the intermediate to
// infinity -- and its dd residual to NaN -- while the true lgamma is still
// finite, over a band ~0.1% wide in x. Grouped this way the intermediate
// exceeds the result by only ~log(x)/2, which is ~351 against an ulp of
// 2e292 at that magnitude: the product and the result overflow at the same
// double. log x - 1 needs no protection of its own, since x >= 8 keeps it
// above 1.
template <class D>
CORVUS_OUTLINE Dd<D> LgammaStirling(D d, op::V<D> x) {
  const auto l = LogDd(d, x);

  // The product is formed on x scaled down by 2^200 and scaled back after.
  // ops::ProdLow's non-FMA path is Dekker's split, whose intermediate
  // a*(2^27+1) overflows once |a| exceeds 2^996 ~ 6.7e299 -- and this is the
  // one place in corvus where an operand gets that large. Unguarded, SSE2
  // through SSE4 return an infinite residual for x above that threshold,
  // while every FMA target is correct.
  //
  // Scaling by a power of two is exact and every step below is linear in x,
  // so the scaled and unscaled computations agree bit for bit -- this is not
  // a separate code path for the non-FMA targets, and the results stay
  // identical across tiers. 2^200 sits between the two limits with room at
  // both ends: x >= 8 keeps the scaled product's residual normal (~2^-250 at
  // the low end), and x <= kLgammaMaxArg keeps the split under its ceiling.
  const auto down = op::Set(d, 0x1p-200);
  const auto up = op::Set(d, 0x1p200);
  // FLOOR: x in [8, kLgammaMaxArg ~ 2^1014] (Stirling lanes; the driver
  // passes Max(x, X0) on mixed vectors), so log x - 1 is in [1.07, 703] and
  // x*2^-200 in [2^-197, 2^814]: the product is >= 2^-197.
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

// lgamma on (5/2, X0) from the Taylor table (#31): with x0 = 5/2 + i/32 the
// nearest node and r = x - x0,
//   lgamma(x) = c0 + c1*r + r^2*(c2 + c3*r + ... + c8*r^6).
// x - 5/2 and r are exact: x, 5/2 and x0 all lie on x's grid (ulp(x) <=
// 2^-50 here) and neither difference leaves it, so |r| <= 1/64 with no
// rounding; scaling by 32 is exact, and Round's tie at a cell edge only
// picks either neighbour, with |r| = 1/64 on both. c0 and c1 are dd, the last
// places where a double's half ulp would show: the tail is <= ~2^-16, so its
// rounding lands near 2^-69 absolute. tools/gen_lgamma_data.py replays this
// arithmetic, fused and unfused, on bit-stepped cell edges: 2^-65.0 relative
// to lgamma (>= 0.28 here). Truncation (the r^9 term) is below that.
//
// A row holds c0.hi, c0.lo, c1.hi, c1.lo, then c8 .. c2, so each lane's eleven
// gathers touch one or two cache lines.
//
// Lanes outside the band must arrive clamped into it: the index is clamped
// too, but a far-off x would still make r huge.
template <class D>
CORVUS_OUTLINE Dd<D> LgammaTable(D d, op::V<D> x) {
  const op::SignedTag<D> di;
  const auto hi = op::Set(d, detail::kLgammaZoneHi);
  const auto invh = op::Set(d, detail::kLgammaTabInvH);
  const auto h = op::Set(d, 1.0 / detail::kLgammaTabInvH);
  auto kf = op::Round(op::Mul(op::Sub(x, hi), invh));
  kf = op::Min(op::Max(kf, op::Zero(d)),
               op::Set(d, static_cast<double>(detail::kLgammaTabN - 1)));
  const auto r = op::Sub(op::Sub(x, hi), op::Mul(kf, h));  // exact
  const auto i = op::ConvertToInt(
      di,
      op::Mul(kf, op::Set(d, static_cast<double>(detail::kLgammaTabStride))));
  const double* tab = detail::kLgammaTab;
  auto p = op::GatherIndex(d, tab + 4, i);
  for (int k = 5; k < detail::kLgammaTabStride; ++k) {
    p = op::MulAdd(p, r, op::GatherIndex(d, tab + k, i));
  }
  const Dd<D> c0{op::GatherIndex(d, tab + 0, i), op::GatherIndex(d, tab + 1, i)};
  const Dd<D> c1{op::GatherIndex(d, tab + 2, i), op::GatherIndex(d, tab + 3, i)};
  // FLOOR: r is 0 (exact) or a multiple of ulp(x) >= 2^-52, and |c1| >= 0.70
  // on the band.
  const auto lead = DdAdd(d, c0, DdMulD(d, c1, r));
  return DdAddD(d, lead, op::Mul(p, op::Mul(r, r)));
}

// lgamma on (0, X0). Three lane classes, each computed only when some lane
// needs it: the zone polynomials on [1/2, 5/2], the same polynomial at t = x
// less log x on (0, 1/2), and the table on (5/2, X0). Stirling lanes, which
// LgammaPosDd may also pass through here, take a discarded table value.
template <class D>
CORVUS_OUTLINE Dd<D> LgammaLow(D d, op::V<D> x) {
  const auto one = op::Set(d, 1.0);
  const auto zhi = op::Set(d, detail::kLgammaZoneHi);
  const auto small = op::Lt(x, op::Set(d, detail::kLgammaZoneLo));
  const auto over = op::Gt(x, zhi);
  {
    // Region routes; each mask is bounded on both sides so Stirling lanes
    // and NaN match none.
    const auto pos = op::Gt(x, op::Zero(d));
    const auto below_x0 = op::Lt(x, op::Set(d, detail::kLgammaX0));
    CORVUS_ROUTE(d, "lgamma/small", op::And(pos, small));
    CORVUS_ROUTE(d, "lgamma/small_subnormal",
                 op::And(pos, op::Lt(x, op::Set(d, 0x1p-1022))));
    CORVUS_ROUTE(d, "lgamma/zone1",
                 op::And(op::Ge(x, op::Set(d, detail::kLgammaZoneLo)),
                         op::Lt(x, op::Set(d, detail::kLgammaZoneMid))));
    CORVUS_ROUTE(d, "lgamma/zone2",
                 op::And(op::Ge(x, op::Set(d, detail::kLgammaZoneMid)),
                         op::Lt(x, op::Set(d, 0x1.4000000000001p+1))));
    CORVUS_ROUTE(d, "lgamma/table", op::And(over, below_x0));
  }
  Dd<D> res{op::Zero(d), op::Zero(d)};
  if (!op::AllTrue(d, over)) {
    // t = x - centre, exact in every band: x for (0, 1/2), x - 1 by Sterbenz
    // on [1/2, 3/2), x - 2 on [3/2, 5/2]. Table and Stirling lanes run at
    // 5/2, keeping t in the polynomials' range (NaN stays NaN).
    const auto xz = op::IfThenElse(over, zhi, x);
    const auto c1 = op::Lt(xz, op::Set(d, detail::kLgammaZoneMid));
    const auto t = op::IfThenElse(
        small, xz,
        op::IfThenElse(c1, op::Sub(xz, one), op::Sub(xz, op::Set(d, 2.0))));
    // FLOOR: below it only on small lanes, x < ~2^-966 (|B| >= 0.42; other t
    // are >= 2^-53). The <= 2^-1072 loss joins -log x >= 669, ~2^-975 under
    // that sum's own dd rounding -- invisible to lgamma, gamma and beta alike.
    res = DdMulD(d, ZoneBracket(d, t, c1), t);
    if (!op::AllFalse(d, small)) {
      // LogDdAny rather than LogDd because x may be subnormal here. Other
      // lanes log 1, whose table slot carries R = 1, L = 0.
      const auto lg = OutlinedLogDd(d, op::IfThenElse(small, x, one));
      const auto sub = DdAdd(d, res, Dd<D>{op::Neg(lg.hi), op::Neg(lg.lo)});
      res = Dd<D>{op::IfThenElse(small, sub.hi, res.hi),
                  op::IfThenElse(small, sub.lo, res.lo)};
    }
  }
  if (!op::AllFalse(d, over)) {
    CORVUS_ROUTE_SCOPE(d, over);
    const auto x0 = op::Set(d, detail::kLgammaX0);
    const auto tb = LgammaTable(d, op::IfThenElse(over, op::Min(x, x0), x0));
    res = Dd<D>{op::IfThenElse(over, tb.hi, res.hi),
                op::IfThenElse(over, tb.lo, res.lo)};
  }
  return res;
}

// lgamma for x > 0, as a dd. The caller rounds -- the reflection needs the
// full precision, since on the negative axis this value is subtracted from
// three others of comparable size.
template <class D>
CORVUS_OUTLINE Dd<D> LgammaPosDd(D d, op::V<D> x) {
  const auto big = op::Ge(x, op::Set(d, detail::kLgammaX0));
  CORVUS_ROUTE(d, "lgamma/stirling", big);

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
CORVUS_OUTLINE Dd<D> LogSinc(D d, op::V<D> u) {
  // FLOOR: u*u < 2^-968 needs |u| < 2^-484, so x = u (|x| >= 1/2 leaves |u|
  // 0 or >= 2^-53); LogSinc is then < 2^-966 in all, against lgamma(x) ~
  // -log|x| >= 335, so its few-2^-1074 loss is invisible.
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
              DdMul(d, acc, v));  // FLOOR: only with v, argued at v above.
  return DdMul(d, acc, v);        // FLOOR: as acc*v above.
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
    CORVUS_ROUTE_SCOPE(d, neg);
    CORVUS_ROUTE(d, "lgamma/neg", neg);
    // u = x - round(x), exact; |u| <= 1/2. u == 0 marks a pole.
    const auto u = op::Sub(x, op::Round(x));
    const auto au = op::Abs(u);
    // lgamma(x) = -log(|u|*|x|) - logsinc(u) - lgamma(-x). Below |x| = 2^-400
    // u = x, so the product is x^2: both factors scale by 2^600 (exact), and
    // 1200 ln 2 comes back out. Above it the product is 0 at a pole or
    // >= 2^-53 * 1/2 (|x| >= 1/2) or x^2 >= 2^-800.
    const auto tiny = op::Lt(ax, op::Set(d, 0x1p-400));
    const auto sc = op::IfThenElse(tiny, op::Set(d, 0x1p600), op::Set(d, 1.0));
    // FLOOR: at a pole the product is exact zero; otherwise, by the bounds
    // above, it is >= 2^-800, and scaled tiny lanes hold subnormal x * 2^600
    // >= 2^-474, so >= 2^-948.
    const auto pr = TwoProd(d, op::Mul(au, sc), op::Mul(ax, sc));
    auto lp = OutlinedLogDd(d, pr);
    if (!op::AllFalse(d, tiny)) {
      CORVUS_ROUTE(d, "lgamma/neg_tiny", tiny);
      // FLOOR: ln 2 and 1200 are constants, nowhere near it.
      const auto k = DdMulD(d,
                            Dd<D>{op::Set(d, detail::kLogLn2Hi),
                                  op::Set(d, detail::kLogLn2Lo)},
                            op::Set(d, 1200.0));
      const auto lps = DdAdd(d, lp, Dd<D>{op::Neg(k.hi), op::Neg(k.lo)});
      lp = Dd<D>{op::IfThenElse(tiny, lps.hi, lp.hi),
                 op::IfThenElse(tiny, lps.lo, lp.lo)};
    }
    const auto r = DdAdd(d, DdAdd(d, lp, LogSinc(d, u)), g);
    const Dd<D> refl{op::Neg(r.hi), op::Neg(r.lo)};
    const auto pole = op::Eq(u, zero);
    CORVUS_ROUTE(d, "lgamma/special/pole", pole);
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
  // ~2.5599833e305, and the whole negative range below -2.5599833e305, where every
  // double is an integer and therefore a pole (so -inf lands right too). NaN
  // propagates with its payload, matching erf/erfc.
  CORVUS_ROUTE(d, "lgamma/special/zero", op::Eq(x, zero));
  CORVUS_ROUTE(d, "lgamma/special/overflow",
               op::And(op::Gt(ax, op::Set(d, detail::kLgammaMaxArg)),
                       op::Lt(ax, inf)));
  CORVUS_ROUTE(d, "lgamma/special/inf", op::Eq(ax, inf));
  CORVUS_ROUTE(d, "lgamma/special/nan", op::IsNaN(x));
  out = op::IfThenElse(op::Eq(x, zero), inf, out);
  out = op::IfThenElse(op::Gt(ax, op::Set(d, detail::kLgammaMaxArg)), inf, out);
  return op::IfThenElse(op::IsNaN(x), x, out);
}

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#endif  // include guard
