// cos/sin kernels (#32): quadrant reduction x = n*(pi/2) + r with shared
// degree-6 parity cores, in two reduction regions blended per lane.
// Per-target include guard (Highway -inl.h idiom).
//
// PROVENANCE
//   Small region: transcription of the certified clean-room consumer kernel
//   (libstats #95 x86 form <- libhmm #74 <- libstats NEON original; formal
//   divergence audit vs ARM optimized-routines -- see PLAN.md's v0.8.0 S1
//   audit record). Constants: src/trig_data.h (tools/gen_trig_table.py,
//   bit-identical to the donor values by generator self-check).
//   Large region: clean-room Payne-Hanek from published mathematics (Payne
//   & Hanek 1983; Ng 1992; Muller et al., Handbook of FP Arithmetic).
//   Windows: src/trig_ph_data.h (tools/gen_trig_ph_table.py).
//
// STRUCTURE
//   Both exports reduce |x| (never signed x) and share the reduction and the
//   parity cores; sin applies the input's sign bit to its result at the end.
//   That is exact, not an approximation: every term of the sin core is odd
//   in r and IEEE negation commutes with rounding, so sin(-x) computed as
//   -sin(|x|) is bit-identical to reducing the signed argument; cos is even.
//   It also makes sin(+/-0) = +/-0 fall out of construction (the +0 core
//   result carries the sign bit of x), with no special-case blend.
//
// SMALL REGION (|x| <= kTrigDMax = 2^23)
//   n = round(x * 2/pi) via Round (nearest-even); r via the 4-part pi/2
//   split, subtracted in order with each step's residual recovered exactly
//   into rlo. Parts 0-2 carry 30 significant bits, so n*part_k is exact for
//   |n| < 2^23 WITHOUT FMA (30+23 <= 53) -- the reduction's exactness is
//   structural, not FMA-dependent; part 3's product rounding (~2^-124
//   absolute) sits below the reduction's ~2^-60 floor. On non-FMA tiers the
//   MulAdd ladders in the cores fall back to mul+add; the donor gates that
//   exact form at 2 ULP (its SSE2 tier), which is this kernel's expected
//   non-FMA bound too.
//
// LARGE REGION (|x| > 2^23): Payne-Hanek over per-exponent windows
//   W_e = (2^(e-52) * 2/pi) mod 4 as four 53-bit chunks c0..c3 at fixed
//   weights; f = (M * W_e) mod 4 with M = mantissa in [2^52, 2^53), so
//   r = (f - round(f)) * pi/2 and quadrant = round(f) mod 4. Every product
//   M*c_k is an EXACT dd via TwoProd (53-bit chunks need no splitting), and
//   the accumulation below is an EXACT expansion until the final
//   compression:
//     * integer stripping s = hi - 4*round(hi/4) is exact for any double
//       (4*round(hi/4) is representable and shares hi's grid);
//     * every sum in the ladder is a TwoSum (exact by construction);
//     * n = round(s), f1 = s - n is exact (integers on s's grid).
//   The only roundings are the plain adds of the second-order residuals in
//   `lo`, and those SHRINK WITH THE CANCELLATION: each residual is bounded
//   by ulp(its ladder sum), and when f - n is small every ladder sum is
//   small (nothing later in the ladder can cancel an O(1) intermediate), so
//   the absolute error tracks ~2^-105 relative to f - n itself rather than
//   relative to 1. At binary64's deepest reduction cancellation
//   (|f - n| ~ 2^-61.5 at e = 849) the generator-certified window
//   truncation (2^-157) dominates: r is good to ~2^-93 relative there.
//   The reference set carries every exponent's CF worst cases, so the ULP
//   gate certifies this end-to-end.
//   Discarded lanes (small, NaN, inf) still EXECUTE this path
//   (masked-lane doctrine): the row index is clamped into table range, so
//   the gathers stay in bounds and such lanes produce garbage that the
//   blends discard.
//
// SPECIALS
//   cos/sin(+/-inf) = NaN (explicit blend); NaN propagates with payload
//   (explicit IfThenElse(IsNaN(x), x, .) at the end -- the polynomial path
//   also self-propagates, but the blend makes payload preservation a
//   contract, not an accident); sin(+/-0) = +/-0 and cos(+/-0) = 1 exactly
//   by construction (traced in the derivation above).
#if defined(CORVUS_TRIG_INL_H_) == defined(HWY_TARGET_TOGGLE)
#ifdef CORVUS_TRIG_INL_H_
#undef CORVUS_TRIG_INL_H_
#else
#define CORVUS_TRIG_INL_H_
#endif

#include <cstdint>
#include <limits>

#include "src/dd-inl.h"
#include "src/ops-inl.h"
#include "src/trig_data.h"
#include "src/trig_ph_data.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {

// Reduced argument: compensated pair (r, rlo) with |r| <= pi/4 (plus the
// fit pad), and the quadrant integer as a double (small: n up to 5,340,354;
// large: n in [-2, 2]; exact in double either way).
template <class D>
struct TrigRed {
  op::V<D> r;
  op::V<D> rlo;
  op::V<D> nq;
};

// Small-region reduction on ax = |x| <= 2^23. Donor transcription.
template <class D>
HWY_INLINE TrigRed<D> TrigReduceSmall(D d, op::V<D> ax) {
  const auto kf = op::Round(op::Mul(ax, op::Set(d, detail::kTrigTwoOverPi)));
  const auto nkf = op::Neg(kf);

  auto r = op::MulAdd(nkf, op::Set(d, detail::kTrigPio2[0]), ax);  // exact
  auto rlo = op::Zero(d);
  for (int k = 1; k < 4; ++k) {
    const auto pk = op::Set(d, detail::kTrigPio2[k]);
    const auto rk = op::MulAdd(nkf, pk, r);
    // Exact residual of the step: rounding and large cancellation are
    // mutually exclusive here (donor derivation), and n*pk is itself exact
    // for k <= 2.
    const auto e = op::MulAdd(nkf, pk, op::Sub(r, rk));
    rlo = op::Add(rlo, e);
    r = rk;
  }
  return {r, rlo, kf};
}

// Large-region (Payne-Hanek) reduction on ax. Executes for every lane; the
// caller blends. See the file header for the exactness derivation.
template <class D>
HWY_NOINLINE TrigRed<D> TrigReducePh(D d, op::V<D> ax) {
  const op::SignedTag<D> di;

  const auto bits = op::BitCast(di, ax);  // ax >= 0: no sign bit set
  // Row index = biased exponent - (1023 + kTrigPhEMin), clamped into table
  // range IN THE DOUBLE DOMAIN so discarded small/NaN/inf lanes gather in
  // bounds (integer Min/Max would work too; the values are tiny and exact).
  const auto e_raw = op::ShiftRight<52>(bits);
  auto idx_d = op::ConvertToDouble(
      d, op::Sub(e_raw, op::Set(di, int64_t{1023 + detail::kTrigPhEMin})));
  idx_d = op::Min(op::Max(idx_d, op::Zero(d)),
                  op::Set(d, double(detail::kTrigPhRows - 1)));
  const auto idx = op::ConvertToInt(di, idx_d);

  // Mantissa as an integer-valued double M in [2^52, 2^53) (exact: the
  // re-host to [1, 2) is a bit op, the 2^52 scale a power of two).
  const auto m0 = op::BitCast(
      d, op::Or(op::And(bits, op::Set(di, int64_t{0x000FFFFFFFFFFFFF})),
                op::Set(di, int64_t{0x3FF0000000000000})));
  const auto mf = op::Mul(m0, op::Set(d, 0x1p52));

  const auto c0 = op::GatherIndex(d, detail::kTrigPhTab0, idx);
  const auto c1 = op::GatherIndex(d, detail::kTrigPhTab1, idx);
  const auto c2 = op::GatherIndex(d, detail::kTrigPhTab2, idx);
  const auto c3 = op::GatherIndex(d, detail::kTrigPhTab3, idx);

  // Exact dd products (TwoProd; Dekker on non-FMA tiers -- magnitudes are
  // far from the split's overflow ceiling, cf. #12's hazard class).
  const auto p0 = TwoProd(d, mf, c0);  // hi <= 2^55
  const auto p1 = TwoProd(d, mf, c1);  // hi <= 4
  const auto p2 = TwoProd(d, mf, c2);  // hi <= 2^-51
  const auto p3 = TwoProd(d, mf, c3);  // hi <= 2^-104; p3.lo dropped below

  const auto quarter = op::Set(d, 0.25);
  const auto neg4 = op::Set(d, -4.0);

  // Strip multiples of 4 from p0.hi (exact; see header). p0.lo <= 4 already.
  const auto q0 = op::Round(op::Mul(p0.hi, quarter));
  const auto d0 = op::MulAdd(q0, neg4, p0.hi);

  // Exact TwoSum ladder over the O(1) terms, re-stripped once they are
  // summed (the strips drop exact multiples of 4, so quadrant = n mod 4
  // is unaffected).
  const auto a = TwoSum(d, p0.lo, p1.hi);
  const auto b = TwoSum(d, a.hi, d0);
  const auto q1 = op::Round(op::Mul(b.hi, quarter));
  const auto s = op::MulAdd(q1, neg4, b.hi);
  const auto n = op::Round(s);
  const auto f1 = op::Sub(s, n);  // exact: integers on s's grid

  // Second-order ladder: cancellation between f1 and the residuals stays
  // exact; the plain adds in `lo` are the only roundings and their inputs
  // shrink with the cancellation (header derivation).
  const auto g = TwoSum(d, p1.lo, p2.hi);
  const auto h1 = TwoSum(d, f1, b.lo);
  const auto h2 = TwoSum(d, h1.hi, a.lo);
  const auto h3 = TwoSum(d, h2.hi, g.hi);
  const auto lo = op::Add(op::Add(op::Add(h1.lo, h2.lo), op::Add(h3.lo, g.lo)),
                          op::Add(p2.lo, p3.hi));
  const auto fmn = TwoSum(d, h3.hi, lo);  // TwoSum: no ordering assumption

  // r = (f - n) * pi/2 in dd.
  const Dd<D> pio2{op::Set(d, detail::kTrigPio2DdHi),
                   op::Set(d, detail::kTrigPio2DdLo)};
  const auto r = DdMul(d, Dd<D>{fmn.hi, fmn.lo}, pio2);
  return {r.hi, r.lo, n};
}

// Both regions, blended, with a skip guard so vectors with no huge lane
// (the overwhelmingly common case) never pay for the gathers.
template <class D>
HWY_NOINLINE TrigRed<D> TrigReduce(D d, op::V<D> ax) {
  TrigRed<D> red = TrigReduceSmall(d, ax);
  const auto big = op::Gt(ax, op::Set(d, detail::kTrigDMax));  // NaN: false
  if (!op::AllFalse(d, big)) {
    const TrigRed<D> ph = TrigReducePh(d, ax);
    red.r = op::IfThenElse(big, ph.r, red.r);
    red.rlo = op::IfThenElse(big, ph.rlo, red.rlo);
    red.nq = op::IfThenElse(big, ph.nq, red.nq);
  }
  return red;
}

// Shared parity cores on u = r^2 (donor transcription): sin's core is odd
// in r term-by-term (the exactness basis for the |x| symmetry), cos's
// 1 - u/2 head is split into an exact (h, hl) pair -- kTrigCosC[0] == -1/2
// exactly (generator-asserted), and both the split's subtractions are exact
// with or without FMA (Sterbenz + shared grid).
template <class D>
HWY_NOINLINE void TrigCores(D d, op::V<D> r, op::V<D> rlo, op::V<D>& s_core,
                            op::V<D>& c_core) {
  const auto u = op::Mul(r, r);

  auto ps = op::Set(d, detail::kTrigSinC[6]);
  for (int i = 5; i >= 0; --i) {
    ps = op::MulAdd(ps, u, op::Set(d, detail::kTrigSinC[i]));
  }
  s_core = op::Add(r, op::MulAdd(op::Mul(r, u), ps, rlo));

  auto pc = op::Set(d, detail::kTrigCosC[6]);
  for (int i = 5; i >= 1; --i) {
    pc = op::MulAdd(pc, u, op::Set(d, detail::kTrigCosC[i]));
  }
  const auto one = op::Set(d, 1.0);
  const auto half = op::Set(d, 0.5);
  const auto nu = op::Neg(u);
  const auto h = op::MulAdd(nu, half, one);                  // 1 - u/2, exact
  const auto hl = op::MulAdd(nu, half, op::Sub(one, h));     // (1-h) - u/2, exact
  auto mc = op::MulAdd(op::Mul(u, u), pc, hl);
  mc = op::MulAdd(op::Neg(r), rlo, mc);  // first-order effect of rlo on cos
  c_core = op::Add(h, mc);
}

// Quadrant select/sign assembly, shared shape: pick core (bit0 swaps),
// apply the quadrant sign, then the specials blends.
// cos: q=0:+c 1:-s 2:-c 3:+s -> swap on bit0, sign on bit1 XOR bit0.
// sin: q=0:+s 1:+c 2:-s 3:-c -> swap on bit0, sign on bit1; input sign
//      applied afterwards (|x| symmetry, header).
template <bool kIsCos, class D>
HWY_INLINE op::V<D> TrigVec(D d, op::V<D> x) {
  const op::SignedTag<D> di;

  const auto ax = op::Abs(x);
  const TrigRed<D> red = TrigReduce(d, ax);
  op::V<D> s_core, c_core;
  TrigCores(d, red.r, red.rlo, s_core, c_core);

  const auto ki = op::ConvertToInt(di, red.nq);
  const auto one_i = op::Set(di, int64_t{1});
  const auto bit0 = op::And(ki, one_i);
  const auto bit1 = op::And(op::ShiftRight<1>(ki), one_i);

  // All-ones where bit0 = 1 (0 - 1 = all-ones), as a double-domain bitmask.
  const auto swap_v = op::BitCast(d, op::Sub(op::Zero(di), bit0));
  const auto sign_i = kIsCos ? op::Xor(bit1, bit0) : bit1;
  const auto sign_v = op::BitCast(d, op::ShiftLeft<63>(sign_i));

  const auto base = kIsCos ? c_core : s_core;
  const auto other = kIsCos ? s_core : c_core;
  auto res = op::Xor(base, op::And(swap_v, op::Xor(base, other)));
  res = op::Xor(res, sign_v);

  if (!kIsCos) {
    // sin(-x) = -sin(|x|), exact (odd core). Also yields sin(-0) = -0.
    const auto x_sign = op::And(op::BitCast(di, x),
                                op::Set(di, std::numeric_limits<int64_t>::min()));
    res = op::Xor(res, op::BitCast(d, x_sign));
  }

  const auto inf = op::Set(d, std::numeric_limits<double>::infinity());
  const auto qnan = op::Set(d, std::numeric_limits<double>::quiet_NaN());
  res = op::IfThenElse(op::Eq(ax, inf), qnan, res);
  res = op::IfThenElse(op::IsNaN(x), x, res);  // payload-preserving
  return res;
}

template <class D>
HWY_INLINE op::V<D> CosVec(D d, op::V<D> x) {
  return TrigVec<true>(d, x);
}
template <class D>
HWY_INLINE op::V<D> SinVec(D d, op::V<D> x) {
  return TrigVec<false>(d, x);
}

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#endif  // include guard
