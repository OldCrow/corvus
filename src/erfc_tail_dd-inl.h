// Shared dd erfc-tail factor for the two Temme tails (#69): GammaTailG,
// hoisted out of src/gamma-inl.h so beta's R3 tail (src/beta-inl.h
// BetaR3Temme) evaluates the same G instead of erfc's plain-double one.
// A pure move: gamma's arithmetic is unchanged, and the gamma/gammainv ULP
// gates and a bit dump were re-run to confirm identical output. Per-target
// include guard (Highway -inl.h idiom).
#if defined(CORVUS_ERFC_TAIL_DD_INL_H_) == defined(HWY_TARGET_TOGGLE)
#ifdef CORVUS_ERFC_TAIL_DD_INL_H_
#undef CORVUS_ERFC_TAIL_DD_INL_H_
#else
#define CORVUS_ERFC_TAIL_DD_INL_H_
#endif

#include "src/dd-inl.h"
#include "src/gamma_data.h"
#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {

namespace op = ops;

// Temme tail factor: G(z) = z e^{z^2} erfc(z) as a dd, for z = 1/u in
// [kGammaTailZLo, kGammaTailZHi] = [5, 28.5] (6 until #70), from r = 1/(2z) as a dd
// (the caller's DdRecip, which it also needs for G/(2z)).
//
// WHY A PRIVATE FIT (#64). The tail used erfc's own plain-double G
// (ErfcTailGFromU): ~1.7 x 2^-53 relative (fit 0.81 + Horner rounding 0.89).
// G's relative error reaches the small side at weight (lambda-1)/eta, which
// is 1.28 at lambda = 2 (s_rv < 0 there, so the bracket is smaller than G's
// term), and with a result mantissa near 2 that read 3 ulp at x = 2a - 4 ulp
// (a ~ 1180). The 2-ULP bound was never reachable with a double G. erfc's
// table and ErfcTailGFromU are untouched, so erfc's and erfinv's bits do
// not move.
//
// TWO CALLERS, ONE DOMAIN. gamma-inl.h GammaTemme (a*phi > 25) and
// beta-inl.h BetaR3Temme (cpsi > 25, #69: beta's weight reaches 1.28 too,
// and its plain-double G read 4 ulp against a 3-ULP bound). Both clamp
// their exponent at -800 (kGammaExpFloor, kBetaExpFloor), so on every live
// tail lane z = sqrt(a*phi) or sqrt(cpsi) lies in (5, 28.29], inside the
// fit interval. A third caller must show the same before using it. The
// split was 36 (z = 6) until #70, which moved both callers' core/tail
// split down to 25 and this fit's lower edge with it.
//
// METHOD (tools/gen_gamma_data.py, block "Tail G" and self-check (g)). One
// Chebyshev interpolant in u = 1/z over the whole live range (a*phi in
// (25, 800] puts z in (5, 28.29]; 28.5 is margin), degree 17, re-expanded
// exactly into monomials in the UNSCALED shift s = u - kGammaTailC. s is a
// dd formed by one DdAddD from the dd u = 2r, so no scale constant rounds.
// Horner runs in double on s.hi for coefficients 17..3, then the last three
// steps in dd with dd coefficients 2..0 (the first of them as an exact
// TwoProd with s.hi). The counts are measured, not guessed (on the original
// [6, 28.5] at degree 15): with c0 alone dd the replay reads 2^-58.6, with
// two dd steps 2^-60.9 (the double step-2 rounding, |s q_2|/G ~ 2^-10, is
// the floor), with three 2^-65.3; degree 14 at three steps reads 2^-62.8.
// On [5, 28.5] degree 15 reads 2^-62.5 and degree 17 2^-64.4 (19: no
// better, the dd-step floor). The generator replays this exact op sequence
// and the caller's G/(2z) product under fused AND unfused MulAdd
// (SSE2/SSSE3/SSE4 have no FMA) over uniform and edge-bit-stepped arms, and
// refuses to emit above 2^-63 (measured 2^-64.55 fused, 2^-64.37 unfused,
// worst at z = 5). At weight 1.28 that is ~0.002 ulp of the result.
//
// One interval means no per-lane coefficient selects. Each caller clamps z
// into [5, 28.5] first; that is a no-op on every live tail lane and keeps
// discarded lanes (core lanes, saturated lanes with huge z, NaN-scrubbed
// z = 0) finite. Hot per-vector region inside the outlined Temme cores, so
// CORVUS_OUTLINE (AGENTS.md).
template <class D>
CORVUS_OUTLINE Dd<D> GammaTailG(D d, Dd<D> r) {
  const auto two = op::Set(d, 2.0);
  const Dd<D> u{op::Mul(r.hi, two), op::Mul(r.lo, two)};  // exact
  const auto s = DdAddD(d, u, op::Set(d, -detail::kGammaTailC));
  const double* c = detail::kGammaTailCoef;
  constexpr int kN = detail::kGammaTailNCoef;
  constexpr int kDd = detail::kGammaTailNDd;
  auto q = op::Set(d, c[kN - 1]);
  for (int k = kN - 2; k >= kDd; --k) {
    q = op::MulAdd(q, s.hi, op::Set(d, c[k]));
  }
  Dd<D> g = DdAdd(d,
                  Dd<D>{op::Set(d, c[kDd - 1]),
                        op::Set(d, detail::kGammaTailCoefLo[kDd - 1])},
                  TwoProd(d, s.hi, q));
  for (int k = kDd - 2; k >= 0; --k) {
    g = DdAdd(d,
              Dd<D>{op::Set(d, c[k]), op::Set(d, detail::kGammaTailCoefLo[k])},
              DdMul(d, s, g));
  }
  return g;
}

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#endif  // include guard
