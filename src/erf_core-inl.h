// Shared erf table lookup + local Taylor series, used by both the erf and
// erfc kernels. Returns E_hi and the small part (S*p + E_lo) separately so
// erfc can assemble compensated 1 -/+ erf without losing the low bits.
// Per-target include guard (Highway -inl.h idiom).
#if defined(CORVUS_ERF_CORE_INL_H_) == defined(HWY_TARGET_TOGGLE)
#ifdef CORVUS_ERF_CORE_INL_H_
#undef CORVUS_ERF_CORE_INL_H_
#else
#define CORVUS_ERF_CORE_INL_H_
#endif

#include "src/erf_data.h"
#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {

namespace op = ops;

template <class D>
struct ErfCoreParts {
  op::V<D> e_hi;   // erf(r) high word at the nearest grid point
  op::V<D> small;  // S*p + E_lo: the sub-ulp-of-1 remainder of erf(ac)
};

// erf(ac) = e_hi + small for ac in [0, 6]. Caller must clamp ac into [0, 6]
// and replace NaN lanes (e.g. with 0) -- the gather index is derived from ac.
// Method: table + local Taylor correction (independently derived; see
// tools/gen_erf_table.py header). For r = round(ac*256)/256 (exact) and
// residual d = ac - r (exact by Sterbenz), with w = r^2:
//   erf(ac) = E_hi + (E_lo + S*(d + c2 d^2 + ... + c8 d^8)),
//   c2 = -r, c3 = (2w-1)/3, c4 = r(3-2w)/6, c5 = (4w^2-12w+3)/30,
//   c6 = r(-2w^2/45 + 2w/9 - 1/6), c7 = (8w^3-60w^2+90w-15)/630,
//   c8 = -r(8w^3-84w^2+210w-105)/2520.
// The series runs to d^8, not d^5, because erfc needs the small part to
// RELATIVE precision of erfc(ac) near ac = 6: the c6 truncation term is
// ~3e-30 absolute -- harmless for erf against ulp(1), but up to ~2e-13
// RELATIVE against erfc(6) ~ 2e-17. Worst-case c6 contribution scales as
// ~5e-18 * a^6; terms through c8 push truncation below 1e-18 relative.
template <class D>
HWY_INLINE ErfCoreParts<D> ErfTableCore(D d, op::V<D> ac) {
  const op::SignedTag<D> di;

  const auto rf = op::Round(op::Mul(ac, op::Set(d, 256.0)));  // integral, <= 1536
  const auto ji = op::ConvertToInt(di, rf);                   // exact truncation
  const auto j3 = op::Add(op::ShiftLeft<1>(ji), ji);          // 3*j: field stride

  const auto E  = op::GatherIndex(d, detail::kErfTable + 0, j3);
  const auto S  = op::GatherIndex(d, detail::kErfTable + 1, j3);
  const auto El = op::GatherIndex(d, detail::kErfTable + 2, j3);

  const auto r = op::Mul(rf, op::Set(d, 1.0 / 256.0));  // exact
  const auto dd = op::Sub(ac, r);                       // exact residual
  const auto w = op::Mul(r, r);

  const auto c2 = op::Neg(r);
  const auto c3 = op::MulAdd(w, op::Set(d, 2.0 / 3.0), op::Set(d, -1.0 / 3.0));
  const auto c4 = op::Mul(r, op::MulAdd(w, op::Set(d, -1.0 / 3.0), op::Set(d, 0.5)));
  const auto c5 = op::MulAdd(
      w, op::MulAdd(w, op::Set(d, 2.0 / 15.0), op::Set(d, -2.0 / 5.0)),
      op::Set(d, 0.1));
  const auto c6 = op::Mul(
      r, op::MulAdd(w, op::MulAdd(w, op::Set(d, -2.0 / 45.0), op::Set(d, 2.0 / 9.0)),
                    op::Set(d, -1.0 / 6.0)));
  const auto c7 = op::MulAdd(
      w,
      op::MulAdd(w, op::MulAdd(w, op::Set(d, 4.0 / 315.0), op::Set(d, -2.0 / 21.0)),
                 op::Set(d, 1.0 / 7.0)),
      op::Set(d, -1.0 / 42.0));
  const auto c8 = op::Mul(
      r, op::MulAdd(
             w,
             op::MulAdd(w, op::MulAdd(w, op::Set(d, -1.0 / 315.0), op::Set(d, 1.0 / 30.0)),
                        op::Set(d, -1.0 / 12.0)),
             op::Set(d, 1.0 / 24.0)));

  // t = c2 + d(c3 + d(c4 + ... + d*c8)); p = d + d^2*t
  auto t = op::MulAdd(dd, c8, c7);
  t = op::MulAdd(dd, t, c6);
  t = op::MulAdd(dd, t, c5);
  t = op::MulAdd(dd, t, c4);
  t = op::MulAdd(dd, t, c3);
  t = op::MulAdd(dd, t, c2);
  const auto p = op::MulAdd(op::Mul(dd, dd), t, dd);

  return {E, op::MulAdd(S, p, El)};
}

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#endif  // include guard
