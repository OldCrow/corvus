// Public log and log1p kernels (#32): thin assemblies over the corvus-owned
// log_dd core. Per-target include guard (Highway -inl.h idiom).
//
// METHOD
//   log:   LogDdAny(x) -- the centred-mantissa table core (~2^-70 relative,
//          budget certified by tools/gen_log_table.py's self-check) with its
//          2^600 subnormal prescale -- then one final rounding fl(hi + lo).
//          Total <= 0.5 + ~2^-17 ulp.
//   log1p: two regimes, split at |x| <= 2^-30 (the 0x1p-30 literal in the kernel).
//          MAIN: s = TwoSum(1, x) captures 1 + x EXACTLY for every double x
//          (Sterbenz makes 1 + x itself exact on [-1, -1/2]; elsewhere the
//          residual carries what the sum rounded away), then LogDdAny's dd
//          overload: log(s.hi) + t*(1 - t/2) with t = s.lo/s.hi.
//          TINY: the dd overload's correction c = t*fl(1 - t/2) is plain
//          double, and for |x| ~ 2^-53 the result IS c: LogDd(1) is exactly
//          zero, so c's own rounding (the factor quantizes to one part in
//          2^53, and the t^2/2 term drops entirely for |t| <= 2^-53) lands
//          UNATTENUATED in a result of the same size -- up to ~0.66 ulp
//          from truth, misrounding up to 21% of the 2^-53 binade (#35 H1;
//          the v0.8.0 header claimed "full relative accuracy" here, a
//          relative-vs-absolute frame slip). Below the cut, compute the
//          series directly in exact-product arithmetic instead:
//          x - x^2/2 + x^3/3, with x^2 exact by TwoProd, the x^2/2
//          subtraction captured by Fast2Sum (|x^2/2| <= 2^-61|x| << |x|),
//          and the residual + x^2 tail + cubic folded into one low word.
//          Error budget at the cut (worst case both regimes): neglected
//          x^4/4 <= 2^-122 absolute, i.e. 2^-92 relative to the result ~x
//          (it said "relative" until #76, review L-14); series-side
//          rounding ~2^-23 ulp;
//          main-path absolute error ~2^-105 against ulp(2^-30) = 2^-82
//          is ~2^-23 ulp. Total <= 0.5 + ~2^-23 ulp everywhere -- inside
//          the library's CR-except-near-tie standard with margin.
//          NOTE: no representable double x > -1 has a subnormal
//          1 + x (the grid bottoms out at 1 + x = 2^-53 exactly, at
//          nextafter(-1, 0)), so LogDdAny's prescale is never exercised
//          from log1p -- routed through it anyway for one shared shape.
//
// LOG FAST PATH (#43)
//   Same idea as exp's (src/exp-inl.h): LogFast reuses the core's table,
//   centring and exact argument r, but carries log1p(r) to ~2^-64 relative
//   instead of ~2^-70, as an unevaluated y + rho, then applies the same
//   rounding test as exp with c = 1 + 2^-8, which certifies any relative
//   error up to ~2^-62. Lanes that fail it, and
//   every lane that is not a positive normal (zero, negative, subnormal,
//   inf, NaN), take the dd path. Bits are unchanged (both correctly
//   rounded on the 0-ULP gate).
//   Result magnitude: x near 1 lands in slots 0/127 (R = 1, L = 0, k = 0,
//   so the result IS log1p(r) with r = m - 1 exact); every other slot has
//   |log x| >= 2^-8, and k != 0 gives |log x| >= ln2 - 0.35. Budget,
//   relative to the result, on every tier (ProdLow/SquareLow exact
//   without FMA):
//     series truncation after r^9 (|r| <= 2^-7)        2^-65.3
//     r^3*q(r) evaluation, ~3 roundings on <= 2^-22.6  2^-65
//     dropped rl*r^2 and tail-sum roundings             2^-67
//     table and ln2 representation                       2^-99
//   Total < 2^-63.5; the test's 2^-62 leaves a 2.8x margin.
//
// SPECIALS (explicit blends; discarded lanes still execute the core, whose
// slot index is bit-masked into table range, so any input is gather-safe):
//   log:   +-0 -> -inf; x < 0 -> qNaN; +inf -> +inf; NaN propagates.
//   log1p: x == -1 -> -inf; x < -1 -> qNaN; +inf -> +inf; NaN propagates;
//          +-0 -> x (the core would return +0 for -0: RN(+0 + -0) = +0
//          drops the sign, so the sign-of-zero contract is a blend here,
//          unlike sin's by-construction sign).
#if defined(CORVUS_LOG_INL_H_) == defined(HWY_TARGET_TOGGLE)
#ifdef CORVUS_LOG_INL_H_
#undef CORVUS_LOG_INL_H_
#else
#define CORVUS_LOG_INL_H_
#endif

#include <limits>

#include "src/dd-inl.h"
#include "src/dd_special-inl.h"  // Log1pmxDd branches (#45)
#include "src/log_dd-inl.h"
#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {

// The dd path: every input, specials included. Outlined on every compiler
// (rare path behind LogVec's AllTrue skip).
template <class D>
HWY_NOINLINE op::V<D> LogDdVec(D d, op::V<D> x) {
  const Dd<D> r = OutlinedLogDd(d, x);  // shared HWY_NOINLINE wrapper
  auto res = op::Add(r.hi, r.lo);

  const auto zero = op::Zero(d);
  const auto inf = op::Set(d, std::numeric_limits<double>::infinity());
  const auto qnan = op::Set(d, std::numeric_limits<double>::quiet_NaN());
  CORVUS_ROUTE(d, "log/special/zero", op::Eq(x, zero));
  CORVUS_ROUTE(d, "log/special/negative", op::Lt(x, zero));
  CORVUS_ROUTE(d, "log/special/inf", op::Eq(x, inf));
  CORVUS_ROUTE(d, "log/special/nan", op::IsNaN(x));
  res = op::IfThenElse(op::Eq(x, zero), op::Neg(inf), res);  // both zeros
  res = op::IfThenElse(op::Lt(x, zero), qnan, res);
  res = op::IfThenElse(op::Eq(x, inf), inf, res);
  return op::IfThenElse(op::IsNaN(x), x, res);  // payload-preserving
}

// Fast path: log(x) for positive normal x, with a mask of the lanes whose
// result is certified correctly rounded. See the file header.
template <class D>
struct LogFastResult {
  op::V<D> v;
  op::M<D> ok;
};

template <class D>
CORVUS_OUTLINE LogFastResult<D> LogFast(D d, op::V<D> x) {
  const op::SignedTag<D> di;

  // Slot as LogDd (src/log_dd-inl.h); the same centred m and k, by integer
  // arithmetic. c = bits(x) - bits(kLogSplit/2), kLogSplit/2 = 0x1.6ap-1
  // (biased exponent 1022): c >> 52 is (E - 1022), less a borrow exactly
  // when x's mantissa is below kLogSplit's -- so k0 + 1 for the upper
  // part of the binade and k0 below, LogDd's k. Subtracting c's exponent
  // bits from x leaves m = x/2^k in [kLogSplit/2, kLogSplit). Valid for
  // positive normal x (other lanes are discarded).
  static_assert(detail::kLogSplitSlot == 53, "kLogSplit/2 bits below");
  const auto bits = op::BitCast(di, x);
  const auto c = op::Sub(bits, op::Set(di, int64_t{0x3FE6A00000000000}));
  const auto ke = op::And(c, op::Set(di, int64_t(0xFFF0000000000000ull)));
  const auto m = op::BitCast(d, op::Sub(bits, ke));
  const auto kf = op::ConvertToDouble(d, op::ShiftRight<52>(c));
  const auto j = op::And(op::ShiftRight<detail::kLogShift>(bits),
                         op::Set(di, detail::kLogN - 1));

  const auto rj = op::GatherIndex(d, detail::kLogTableR, j);
  const auto lhi = op::GatherIndex(d, detail::kLogTableLHi, j);
  const auto llo = op::GatherIndex(d, detail::kLogTableLLo, j);

  // r = R_j*m - 1 exactly, normalized: p - 1 is exact by Sterbenz and a
  // multiple of ulp(p) >= 2^-53 when nonzero, so its exponent is never
  // below ProdLow's (<= ulp(p)/2) and Fast2Sum is exact.
  const auto p = op::Mul(rj, m);
  const auto a = op::Sub(p, op::Set(d, 1.0));
  // FLOOR: m is centred on 1 (~[0.70, 1.41]) and R_j ~ 1/m, so p ~ 1; lanes
  // outside the positive normal range are discarded (ok below).
  const auto pl = op::ProdLow(d, rj, m, p);
  const auto rh = op::Add(a, pl);
  const auto rl = op::Sub(pl, op::Sub(rh, a));

  // log1p(r) = rh - rh^2/2 + rh^3*q(rh) + rl*(1 - rh), q through r^6/9.
  auto q = op::Set(d, 1.0 / 9.0);
  q = op::MulAdd(q, rh, op::Set(d, -1.0 / 8.0));
  q = op::MulAdd(q, rh, op::Set(d, 1.0 / 7.0));
  q = op::MulAdd(q, rh, op::Set(d, -1.0 / 6.0));
  q = op::MulAdd(q, rh, op::Set(d, 1.0 / 5.0));
  q = op::MulAdd(q, rh, op::Set(d, -0.25));
  q = op::MulAdd(q, rh, op::Set(d, 1.0 / 3.0));
  const auto sq = op::Mul(rh, rh);
  // FLOOR: rh is 0 or a multiple of ulp(R_j)*ulp(m) >= 2^-107 (a is a
  // multiple of ulp(p), pl of that product grid), so rh^2 >= 2^-214.
  const auto sql = op::SquareLow(d, rh, sq);  // exact
  const auto half = op::Set(d, 0.5);
  const auto h = op::Mul(sq, half);  // exact halving

  // k*ln2 + L_hi + rh - h, every step error-free:
  //   A = k*ln2hi exactly as (A, al).
  //   S1 = A + L_hi: Fast2Sum (k != 0: |A| >= 0.69 > |L|; k = 0: A = 0).
  //   S2 = S1 + rh: TwoSum (|L| and |rh| are comparable).
  //   S3 = S2 - h: Fast2Sum (|S2| >= ~2^-8 >> h, or S2 = rh >= h near 1).
  const auto ln2hi = op::Set(d, detail::kLogLn2Hi);
  const auto A = op::Mul(kf, ln2hi);
  // FLOOR: kf is an integer exponent, |kf| <= 1075: k*ln2hi is 0 or >= 0.69.
  const auto al = op::ProdLow(d, kf, ln2hi, A);
  const auto s1 = Fast2Sum(d, A, lhi);
  const auto s2 = TwoSum(d, s1.hi, rh);
  const auto s3 = Fast2Sum(d, s2.hi, op::Neg(h));

  auto t = op::MulAdd(kf, op::Set(d, detail::kLogLn2Lo), llo);
  t = op::Add(t, op::Add(al, op::Add(s1.lo, op::Add(s2.lo, s3.lo))));
  t = op::MulAdd(op::Neg(sql), half, t);
  t = op::MulAdd(op::Neg(rl), rh, op::Add(t, rl));
  t = op::MulAdd(op::Mul(sq, rh), q, t);  // rh^3*q, <= 2^-22.6
  const auto y = op::Add(s3.hi, t);
  const auto rho = op::Sub(t, op::Sub(y, s3.hi));  // Fast2Sum: exact

  // Rounding test (src/exp-inl.h header; c = 1 + 2^-8 certifies a
  // relative error up to ~2^-62) and the domain the fast path covers:
  // positive normal x (NaN compares false).
  const auto ok = op::And(
      op::And(op::Ge(x, op::Set(d, 0x1p-1022)),
              op::Ge(op::Set(d, std::numeric_limits<double>::max()), x)),
      op::Eq(op::Add(y, op::Mul(rho, op::Set(d, 1.0 + 0x1p-8))), y));
  return {y, ok};
}

template <class D>
HWY_INLINE op::V<D> LogVec(D d, op::V<D> x) {
  const LogFastResult<D> f = LogFast(d, x);
  // Routes (census only): the fast path, and the dd fallback by cause
  // (rounding test failed on a positive normal x; subnormal x; specials
  // are routed inside LogDdVec).
  CORVUS_ROUTE(d, "log/fast", f.ok);
  CORVUS_ROUTE(d, "log/dd/round",
               op::And(op::And(op::Ge(x, op::Set(d, 0x1p-1022)),
                               op::Ge(op::Set(d, std::numeric_limits<double>::max()), x)),
                       op::MaskNot(f.ok)));
  CORVUS_ROUTE(d, "log/dd/subnormal",
               op::And(op::Gt(x, op::Zero(d)), op::Gt(op::Set(d, 0x1p-1022), x)));
  if (op::AllTrue(d, f.ok)) return f.v;
  return op::IfThenElse(f.ok, f.v, LogDdVec(d, x));
}

template <class D>
CORVUS_OUTLINE op::V<D> Log1pVec(D d, op::V<D> x) {
  const auto one = op::Set(d, 1.0);
  const auto s = TwoSum(d, one, x);  // exact 1 + x for every double
  const Dd<D> r = OutlinedLogDd(d, Dd<D>{s.hi, s.lo});
  auto res = op::Add(r.hi, r.lo);

  // Tiny regime (see header): series in exact products. sq = x^2 exactly;
  // Fast2Sum's precondition |x| >= |sq.hi/2| holds for every |x| <= 1/2
  // (equality only at 0, where every term is a zero). The x^3 coefficient's
  // sq.lo*(-0.5) is an exact halving and x*sq.hi rounds ~2^-53 relative on
  // a term that is itself <= 2^-60 of the result, so both MulAdds are
  // fusion-indifferent: FMA and mul+add tiers agree bit for bit.
  // |x| <= 2^-30, via Ge (the facade carries no Le); NaN lanes compare
  // false and ride the main path into the payload-preserving blend below.
  const auto tiny = op::Ge(op::Set(d, 0x1p-30), op::Abs(x));
  CORVUS_ROUTE(d, "log1p/tiny", tiny);
  CORVUS_ROUTE(d, "log1p/main", op::MaskNot(tiny));
  if (!op::AllFalse(d, tiny)) {
    CORVUS_ROUTE_SCOPE(d, tiny);
    // FLOOR: |x| < 2^-484; sq.lo is off by <= min(2^-1073, 1.01 x^2), which
    // is below 2^-483 ulp(x) for every such x (result ~ x, subnormals too).
    const auto sq = TwoProd(d, x, x);
    const auto mhalf = op::Set(d, -0.5);
    const auto s1 = Fast2Sum(d, x, op::Mul(sq.hi, mhalf));
    auto lo = op::MulAdd(sq.lo, mhalf, s1.lo);
    // fl(1/3)'s own rounding enters ~2^-60-relative -- far below budget.
    lo = op::MulAdd(op::Mul(x, sq.hi), op::Set(d, 1.0 / 3.0), lo);
    res = op::IfThenElse(tiny, op::Add(s1.hi, lo), res);
  }

  const auto zero = op::Zero(d);
  const auto inf = op::Set(d, std::numeric_limits<double>::infinity());
  const auto qnan = op::Set(d, std::numeric_limits<double>::quiet_NaN());
  const auto neg1 = op::Set(d, -1.0);
  CORVUS_ROUTE(d, "log1p/special/neg1", op::Eq(x, neg1));
  CORVUS_ROUTE(d, "log1p/special/lt_neg1", op::Lt(x, neg1));
  CORVUS_ROUTE(d, "log1p/special/inf", op::Eq(x, inf));
  CORVUS_ROUTE(d, "log1p/special/zero", op::Eq(x, zero));
  CORVUS_ROUTE(d, "log1p/special/nan", op::IsNaN(x));
  res = op::IfThenElse(op::Eq(x, neg1), op::Neg(inf), res);
  res = op::IfThenElse(op::Lt(x, neg1), qnan, res);
  res = op::IfThenElse(op::Eq(x, inf), inf, res);
  res = op::IfThenElse(op::Eq(x, zero), x, res);  // log1p(+-0) = +-0
  return op::IfThenElse(op::IsNaN(x), x, res);  // payload-preserving
}

// ------------------------------------------------------------------------
// Public log1pmx (#45): log1p(x) - x = -phi(x), relative to the result,
// rounded once. phi is src/dd_special-inl.h's Log1pmxDd, run branch by
// branch so each runs only when some lane needs it.
//
// SERIES, 2^-480 <= |x| <= kPhiCut = 1/16: Log1pmxDdSeries, ~2^-74
// relative; u.lo = 0, so its first-order u.lo term is an exact zero.
//
// DIRECT, |x| > 1/16: Log1pmxDdDirect with w = TwoSum(1, x), the EXACT
// 1 + x as a normalized dd for every double x (its low word is zero
// wherever 1 + x is representable, which covers the whole approach to -1,
// so the hazard rule's non-normalized w cannot arise from a double input).
// ~2^-67.6 relative just above the cut, better away from it. x = -1 gives
// w = 0; that lane is overwritten with -inf below.
//
// TINY, |x| < 2^-480: the series' u^2 is formed by TwoProd, whose low word
// loses relative accuracy once u^2 nears the subnormal range, so phi is
// only ABSOLUTELY accurate there (the FLOOR note at Log1pmxDd) -- several
// subnormal ulps of a result that is itself subnormal below |x| = 2^-511.
// Instead x is scaled by 2^600 (exact), phi = x^2/2 is formed exactly as a
// dd in that scale (the dropped x^3/3 is 2^-480 relative), normalized to a
// mantissa in [1, 2), and DdScaleToDouble applies 2^(e - 1200) with ONE
// rounding onto the subnormal grid. Below |x| ~ 2^-537.5 that rounding
// gives zero; the exponent is clamped at -1100, which only ever moves a
// value that already rounds to zero.
//
// Specials: log1pmx(+-0) = +0 (the IEEE value of log1p(x) - x);
// log1pmx(-1) = -inf; log1pmx(x < -1) = NaN; log1pmx(+inf) = -inf; NaN
// propagates (payload preserved). Large finite x is ordinary: the result is
// log1p(x) - x ~ -x, finite down to -DBL_MAX.
template <class D>
CORVUS_OUTLINE op::V<D> Log1pmxVec(D d, op::V<D> x) {
  const op::SignedTag<D> di;
  const auto zero = op::Zero(d);
  const auto one = op::Set(d, 1.0);
  const auto inf = op::Set(d, std::numeric_limits<double>::infinity());
  const auto qnan = op::Set(d, std::numeric_limits<double>::quiet_NaN());
  const auto neg1 = op::Set(d, -1.0);
  const auto ax = op::Abs(x);
  const auto cut = op::Set(d, detail::kPhiCut);
  const auto tcut = op::Set(d, 0x1p-480);
  const auto m_tiny = op::Gt(tcut, ax);  // includes +-0; NaN: false
  const auto m_ser = op::And(op::Ge(ax, tcut), op::Ge(cut, ax));
  // |x| > 1/16 and x > -1, finite: the direct difference. NaN: false.
  const auto m_dir = op::And(
      op::And(op::Gt(ax, cut), op::Gt(x, neg1)),
      op::Gt(inf, x));
  CORVUS_ROUTE(d, "log1pmx/tiny", m_tiny);
  CORVUS_ROUTE(d, "log1pmx/series", m_ser);
  CORVUS_ROUTE(d, "log1pmx/direct", m_dir);

  auto res = zero;
  if (!op::AllFalse(d, m_ser)) {
    CORVUS_ROUTE_SCOPE(d, m_ser);
    const auto xs = op::IfThenElse(m_ser, x, op::Set(d, 0.03125));
    const auto phi = Log1pmxDdSeries(d, Dd<D>{xs, zero});
    res = op::IfThenElse(m_ser, op::Neg(DdToDouble(phi)), res);
  }
  if (!op::AllFalse(d, m_dir)) {
    CORVUS_ROUTE_SCOPE(d, m_dir);
    const auto xd = op::IfThenElse(m_dir, x, one);
    const auto w = TwoSum(d, one, xd);  // exact
    const auto phi = Log1pmxDdDirect(d, Dd<D>{xd, zero}, w);
    res = op::IfThenElse(m_dir, op::Neg(DdToDouble(phi)), res);
  }
  if (!op::AllFalse(d, m_tiny)) {
    CORVUS_ROUTE_SCOPE(d, m_tiny);
    // xt in [2^-474, 2^120) on live lanes; dead lanes run at 1.
    const auto xt =
        op::Mul(op::IfThenElse(op::And(m_tiny, op::Gt(ax, zero)), ax, one),
                op::Set(d, 0x1p600));
    // FLOOR: xt >= 2^-474, so xt^2 >= 2^-948 and its low word is exact.
    const auto sq = TwoProd(d, xt, xt);
    // 2^-eh with eh = exponent of sq.hi (normal, positive): m in [1, 2).
    const auto eh = op::Sub(
        op::And(op::ShiftRight<52>(op::BitCast(di, sq.hi)),
                op::Set(di, int64_t{0x7ff})),
        op::Set(di, int64_t{1023}));
    const auto inv = op::BitCast(
        d, op::ShiftLeft<52>(op::Sub(op::Set(di, int64_t{1023}), eh)));
    const Dd<D> m{op::Mul(sq.hi, inv), op::Mul(sq.lo, inv)};  // exact
    // phi = m * 2^(eh - 1200 - 1): the 2^1200 undoes the scaling of x^2,
    // the extra -1 is the halving.
    const auto e = op::Max(op::Sub(eh, op::Set(di, int64_t{1201})),
                           op::Set(di, int64_t{-1100}));
    const auto phi = DdScaleToDouble(d, m, e);
    res = op::IfThenElse(m_tiny, op::Neg(phi), res);
  }

  CORVUS_ROUTE(d, "log1pmx/special/zero", op::Eq(x, zero));
  CORVUS_ROUTE(d, "log1pmx/special/neg1", op::Eq(x, neg1));
  CORVUS_ROUTE(d, "log1pmx/special/lt_neg1", op::Lt(x, neg1));
  CORVUS_ROUTE(d, "log1pmx/special/inf", op::Eq(x, inf));
  CORVUS_ROUTE(d, "log1pmx/special/nan", op::IsNaN(x));
  res = op::IfThenElse(op::Eq(x, zero), zero, res);  // +-0 -> +0
  res = op::IfThenElse(op::Eq(x, neg1), op::Neg(inf), res);
  res = op::IfThenElse(op::Lt(x, neg1), qnan, res);
  res = op::IfThenElse(op::Eq(x, inf), op::Neg(inf), res);
  return op::IfThenElse(op::IsNaN(x), x, res);  // payload-preserving
}

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#endif  // include guard
