// Double-double natural logarithm: corvus-owned, so no accuracy-critical
// kernel depends on the backend's log. Internal (no public API); the intended
// consumer is lgamma (Stirling's (x-1/2)*log(x), the log of the recurrence
// product, and the reflection term) and after it the incomplete gamma/beta
// prefactor. Per-target include guard (Highway -inl.h idiom).
//
// METHOD
//   x = 2^k * m,  m in [1+53/128 - 1, 1+53/128) halved -- i.e. m centred on 1
//   j = top 7 mantissa bits of x
//   r = R_j*m - 1,  |r| <= 2^-7
//   log(x) = k*ln2 + L_j + log1p(r),   L_j = -log(R_j), tabulated as a dd
//
// WHY r IS EXACT
//   p = fl(R_j*m) lies in [1-2^-7, 1+2^-7], so p - 1 is EXACT by Sterbenz, and
//   ops::ProdLow supplies the exact product residual. Hence r = (p-1) + p_lo
//   is the true R_j*m - 1 with no error at all. A plain fma(R_j, m, -1) would
//   instead carry ~ulp(r)/2 = 2^-62, which is above this kernel's budget --
//   and on non-FMA targets would carry ulp(1)/2 = 2^-53.
//
// WHY THERE IS NO SPECIAL CASE FOR x NEAR 1
//   Two table choices, both made in tools/gen_log_table.py:
//   * m is centred on 1 rather than taken from [1,2). Uncentred, x slightly
//     below 1 gives k = -1 and log(x) = ln2 + log(m): two terms near 0.693
//     cancelling to ~1e-16, which burns ~45 of the representation's ~106 bits
//     exactly where relative accuracy matters most.
//   * the two slots adjacent to m = 1 carry R_j = 1 and L_j = 0 EXACTLY, so on
//     both sides of 1 the result is log1p(m-1) alone, with m-1 exact by
//     Sterbenz and no cancelling term to add. Elsewhere |L_j| >= ~0.0039 while
//     |log1p(r)| <= ~0.0078, so the remaining cancellation is a couple of bits.
//
// ERROR BUDGET (relative), checked numerically by the generator's self-check:
//   argument r            exact
//   log1p truncation      2^-80.6  (series through r^11 over |r| <= 2^-7)
//   log1p tail rounding   ~2^-70   (r^3*q evaluated in plain double)
//   table representation  2^-107.5 dd split of -log(R_j)
//   dd assembly           ~2^-104
// Total ~2^-70 relative, matching exp_dd's ~2^-68.
//
// DOMAIN
//   x must be a POSITIVE NORMAL double. Zero, negatives, Inf and NaN are the
//   caller's responsibility -- the slot index is masked into range so nothing
//   reads out of bounds, but the value in such a lane is unspecified. For
//   subnormal arguments use LogDdAny, which prescales.
#if defined(CORVUS_LOG_DD_INL_H_) == defined(HWY_TARGET_TOGGLE)
#ifdef CORVUS_LOG_DD_INL_H_
#undef CORVUS_LOG_DD_INL_H_
#else
#define CORVUS_LOG_DD_INL_H_
#endif

#include "src/dd-inl.h"
#include "src/log_dd_data.h"
#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {

// log(1 + r) for |r| <= 2^-7, to dd precision RELATIVE to the result.
//
// The first two terms carry the accuracy: r itself is exact, and r^2/2 needs
// its own dd treatment because at |r| = 2^-7 a rounded r^2/2 would land at
// 2^-70 absolute, i.e. 2^-63 relative to log1p(r) -- above budget. Note the
// 2*r.hi*r.lo cross term is NOT droppable for the same reason (it is 2^-69
// absolute, 2^-61 relative). Everything from r^3 down fits in a plain double:
// r^3/3 is already 2^-24, so its own rounding is 2^-77 absolute.
template <class D>
HWY_INLINE Dd<D> Log1pDd(D d, Dd<D> r) {
  const auto rh = r.hi;

  // Series tail q(r) = 1/3 - r/4 + r^2/5 - ... + r^8/11, so that
  // log1p(r) = r - r^2/2 + r^3*q(r) with truncation O(r^12/12).
  auto q = op::Set(d, 1.0 / 11.0);
  q = op::MulAdd(q, rh, op::Set(d, -1.0 / 10.0));
  q = op::MulAdd(q, rh, op::Set(d, 1.0 / 9.0));
  q = op::MulAdd(q, rh, op::Set(d, -1.0 / 8.0));
  q = op::MulAdd(q, rh, op::Set(d, 1.0 / 7.0));
  q = op::MulAdd(q, rh, op::Set(d, -1.0 / 6.0));
  q = op::MulAdd(q, rh, op::Set(d, 1.0 / 5.0));
  q = op::MulAdd(q, rh, op::Set(d, -0.25));
  q = op::MulAdd(q, rh, op::Set(d, 1.0 / 3.0));

  // r^2 in dd (halving by 0.5 is exact), then r - r^2/2.
  const auto sq = TwoProd(d, rh, rh);
  const auto half = op::Set(d, 0.5);
  const Dd<D> half_sq{op::Mul(sq.hi, half),
                      op::MulAdd(rh, r.lo, op::Mul(sq.lo, half))};
  const auto lead = DdAdd(d, r, Dd<D>{op::Neg(half_sq.hi), op::Neg(half_sq.lo)});

  // r^3 * q, plain double: magnitude <= 2^-24.
  const auto tail = op::Mul(op::Mul(sq.hi, rh), q);
  return DdAddD(d, lead, tail);
}

// log(x) for a positive normal double. See the file header.
template <class D>
HWY_INLINE Dd<D> LogDd(D d, op::V<D> x) {
  const op::SignedTag<D> di;

  const auto bits = op::BitCast(di, x);
  // Unbiased exponent of x, and the mantissa re-hosted into [1, 2).
  const auto k0 = op::Sub(op::ShiftRight<52>(bits), op::Set(di, int64_t{1023}));
  const auto m0 = op::BitCast(
      d, op::Or(op::And(bits, op::Set(di, int64_t{0x000FFFFFFFFFFFFF})),
                op::Set(di, int64_t{0x3FF0000000000000})));
  const auto j = op::And(op::ShiftRight<detail::kLogShift>(bits),
                         op::Set(di, detail::kLogN - 1));

  // Centre the mantissa on 1. The threshold is a slot boundary, so this
  // condition is exactly "j >= kLogSplitSlot" and no slot straddles it --
  // which is what lets one R_j serve the whole slot.
  const auto hi_half = op::Ge(m0, op::Set(d, detail::kLogSplit));
  const auto m = op::IfThenElse(hi_half, op::Mul(m0, op::Set(d, 0.5)), m0);
  const auto kf =
      op::Add(op::ConvertToDouble(d, k0),
              op::IfThenElse(hi_half, op::Set(d, 1.0), op::Zero(d)));

  const auto rj = op::GatherIndex(d, detail::kLogTableR, j);
  const Dd<D> lj{op::GatherIndex(d, detail::kLogTableLHi, j),
                 op::GatherIndex(d, detail::kLogTableLLo, j)};

  // r = R_j*m - 1, exact (see the header): p - 1 is exact by Sterbenz and
  // ProdLow supplies the product residual.
  //
  // TwoSum, not a bare pair: the two halves are exact but NOT normalized --
  // |p-1| is ~2^-8 while |p_lo| reaches 2^-53, far above ulp(p-1)/2. Left
  // that way, every term in Log1pDd computed from r.hi alone silently drops
  // its share of r.lo; the r^3 term's share is r^2*r.lo ~ 2^-69, which
  // costs ~2^-63 relative near x = 1. TwoSum rather than
  // Fast2Sum because r passes through zero inside each slot, so no ordering
  // between the two halves can be assumed.
  const auto p = op::Mul(rj, m);
  const auto r = TwoSum(d, op::Sub(p, op::Set(d, 1.0)), op::ProdLow(d, rj, m, p));

  const Dd<D> ln2{op::Set(d, detail::kLogLn2Hi), op::Set(d, detail::kLogLn2Lo)};
  return DdAdd(d, DdAdd(d, DdMulD(d, ln2, kf), lj), Log1pDd(d, r));
}

// log(x_hi + x_lo) for a positive normal dd. The correction is
// log1p(t) with t = x_lo/x_hi <= 2^-53; the quadratic term is kept because
// near x = 1 the whole result can be the same size as t, where dropping it
// would cost ~2^-54 RELATIVE (it is only negligible against an O(1) result).
template <class D>
HWY_INLINE Dd<D> LogDd(D d, Dd<D> x) {
  const auto t = op::Div(x.lo, x.hi);
  const auto c = op::Mul(t, op::MulAdd(t, op::Set(d, -0.5), op::Set(d, 1.0)));
  return DdAddD(d, LogDd(d, x.hi), c);
}

// log(x) for any positive finite dd, SUBNORMALS INCLUDED. LogDd above reads
// the exponent straight out of the bit pattern, which a subnormal does not
// carry, so tiny arguments are first scaled into the normal range by an exact
// power of two and the scaling is taken back out in dd afterwards.
//
// 2^600 lifts the smallest subnormal (2^-1074) to 2^-474, and the threshold is
// low enough that the scaled value cannot overflow. lgamma reaches this
// through two doors: x -> 0+ on the positive axis, and |x - round(x)| for a
// negative subnormal argument.
template <class D>
HWY_INLINE Dd<D> LogDdAny(D d, Dd<D> x) {
  const auto tiny = op::Lt(x.hi, op::Set(d, 0x1p-500));
  const auto s = op::IfThenElse(tiny, op::Set(d, 0x1p600), op::Set(d, 1.0));
  const auto e = op::IfThenElse(tiny, op::Set(d, -600.0), op::Zero(d));
  const Dd<D> xs{op::Mul(x.hi, s), op::Mul(x.lo, s)};  // exact: power of two
  const Dd<D> ln2{op::Set(d, detail::kLogLn2Hi), op::Set(d, detail::kLogLn2Lo)};
  // e is exactly zero on the common path, so this add is a no-op there.
  return DdAdd(d, LogDd(d, xs), DdMulD(d, ln2, e));
}

template <class D>
HWY_INLINE Dd<D> LogDdAny(D d, op::V<D> x) {
  return LogDdAny(d, Dd<D>{x, op::Zero(d)});
}

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#endif  // include guard
