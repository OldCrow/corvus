// Shared erfc assembly, hoisted out of erfc.cpp so erfinv's tail core can
// consume it without either re-deriving it or going through the ROUNDED
// public corvus::erfc/erf (which would floor a Newton/Halley step -- see
// PLAN.md's Phase C erfinv/erfcinv design). Per-target include guard
// (Highway -inl.h idiom).
//
// Two pieces are exposed, both pure refactors of erfc.cpp with NO change to
// its arithmetic -- the existing erfc ULP gates are the regression guard,
// and were re-run after this split to confirm bit-identical output:
//
//   ErfcCoreDd: the |x| <= 6 compensated assembly's dd pair BEFORE the
//   final rounding erfc.cpp used to do inline (hi = 1 - sE by Fast2Sum,
//   lo folds in the exact residual and the erf table's "small" term). This
//   is erfinv's mid-region residual: F = ErfcCoreDd(x0) (-) s stays a dd
//   subtraction, so the ~2^-13-relative cancellation near the root is
//   absorbed with ~2^-77 to spare instead of being floored by a prior
//   rounding to double.
//
//   ErfcTailGFromU: the tail's G(u) polynomial evaluation (u = 1/x, u
//   already computed by the caller), reused by erfinv's far-tail Halley
//   step to write log(erfc(x)) = -x^2 - log(x) + log(G(1/x)) without any
//   exponential (see src/erfinv-inl.h). erfc.cpp's own tail path is
//   unchanged: it still forms u itself (via DdRecip, for the dd-precision
//   1/a it needs elsewhere) and simply calls this helper instead of
//   inlining the same Horner loop.
#if defined(CORVUS_ERFC_CORE_INL_H_) == defined(HWY_TARGET_TOGGLE)
#ifdef CORVUS_ERFC_CORE_INL_H_
#undef CORVUS_ERFC_CORE_INL_H_
#else
#define CORVUS_ERFC_CORE_INL_H_
#endif

#include "src/dd-inl.h"
#include "src/erf_core-inl.h"
#include "src/erfc_tail_data.h"
#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {

// Select one of three per-interval constants by the tail-interval masks
// (m1 innermost, m2 outer -- m1 implies m2). Shared by the tail's own
// scale/shift/coefficient selects and by erfinv's identical 3-interval
// seed-coefficient select.
template <class D, class M>
HWY_INLINE op::V<D> Sel3(D d, M m1, M m2, double v0, double v1, double v2) {
  return op::IfThenElse(m1, op::Set(d, v0),
                        op::IfThenElse(m2, op::Set(d, v1), op::Set(d, v2)));
}

// The erf table's last grid point is r = 6.0 exactly (index kErfTableLastIndex
// = 1536, step 1/256); ErfTableCore's index is round(ac*256), so anything
// with ac*256 < 1536.5 -- i.e. ac < 6 + 1/512 -- still lands on that last
// grid point and reads the table in bounds, extrapolating the local Taylor
// correction by a fraction of the last grid interval. kErfcCoreSafeMax
// leaves a >100x margin below that hard limit for erfinv's mid-region seed
// (see src/erfinv-inl.h), whose ~2^-19 relative error can legitimately place
// x0 a few 2^-17-ish past 6 when the true root sits right at the mid/far
// seam. erfc.cpp itself never needs more than exactly 6.0 (its own domain
// contract), so this is strictly a widening of ErfcCoreDd's safe INPUT
// range, not a change to what erfc.cpp computes -- see the file header.
inline constexpr double kErfcCoreSafeMax = 6.0 + 1.0 / 1024.0;

// erfc(x) for |x| <= kErfcCoreSafeMax (an extrapolation margin past erfc's
// own |x| <= 6 domain -- see kErfcCoreSafeMax), as a dd pair BEFORE the
// final rounding. See the file header. Caller must pass ax already
// NaN-scrubbed to a safe table index (ErfTableCore's precondition) --
// erfc.cpp does this once for the NaN mask (its own callers never exceed 6);
// erfinv's mid-region seed can legitimately land a little past 6, which is
// exactly the margin this function's clamp preserves instead of silently
// mis-evaluating at a clamped-to-6.0 point (see src/erfinv-inl.h).
template <class D>
HWY_INLINE Dd<D> ErfcCoreDd(D d, op::V<D> x, op::V<D> ax_safe) {
  const auto one = op::Set(d, 1.0);
  const auto sgn = op::CopySign(one, x);
  const auto ac = op::Min(ax_safe, op::Set(d, kErfcCoreSafeMax));
  const auto parts = ErfTableCore(d, ac);

  const auto sE = op::Mul(sgn, parts.e_hi);
  const auto hi = op::Sub(one, sE);
  const auto lo = op::Sub(op::Sub(one, hi), sE);  // exact Fast2Sum residual
  return Dd<D>{hi, op::Sub(lo, op::Mul(sgn, parts.small))};
}

// G(u) = a*e^{a^2}*erfc(a) for a = at (tail region, 6 <= at <= 28), u = 1/at
// already computed by the caller (erfc.cpp needs u to dd precision for its
// own 1/a factor and gets it from DdRecip; erfinv's far branch only needs
// double accuracy in u -- see src/erfinv-inl.h -- and computes it more
// cheaply). See tools/gen_erfc_tail_poly.py for the fit.
template <class D>
HWY_INLINE op::V<D> ErfcTailGFromU(D d, op::V<D> at, op::V<D> u) {
  const auto m1 = op::Lt(at, op::Set(d, detail::kErfcTailBound1));
  const auto m2 = op::Lt(at, op::Set(d, detail::kErfcTailBound2));
  const auto scale = Sel3(d, m1, m2, detail::kErfcTailScale[0],
                          detail::kErfcTailScale[1], detail::kErfcTailScale[2]);
  const auto shift = Sel3(d, m1, m2, detail::kErfcTailShift[0],
                          detail::kErfcTailShift[1], detail::kErfcTailShift[2]);
  const auto s = op::MulAdd(u, scale, shift);

  const auto* c = detail::kErfcTailCoef;
  auto poly = Sel3(d, m1, m2, c[0][detail::kErfcTailNCoef - 1],
                   c[1][detail::kErfcTailNCoef - 1],
                   c[2][detail::kErfcTailNCoef - 1]);
  for (int k = detail::kErfcTailNCoef - 2; k >= 0; --k) {
    poly = op::MulAdd(poly, s, Sel3(d, m1, m2, c[0][k], c[1][k], c[2][k]));
  }
  return poly;
}

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#endif  // include guard
