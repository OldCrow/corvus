// Per-component cost and loop-depth evidence for the incomplete gamma/beta
// region cores (#42 lever 1, row C1). Compiles the kernel headers itself
// (kernel-test pattern); the shipped library is untouched. Manual, Release
// only; timings want a quiet machine (tools/quiet_bench.*). Modes:
//   (none)        per-component ns/el by region band (cores, prefactor, CF
//                 depth variants, one #44-style dd quadrature node)
//   depth         per-lane depth each fixed loop needs for a bit-identical dd
//                 result, and the per-vector maximum for 2/4/8 lanes
//   hunt          BetaR2Cf (the precision split and its F gate) against the
//                 pre-#42 all-dd recurrence over R2's routed domain, binned
//                 by F -- the evidence BetaR2Cf's comment cites
//   fstat         how often the F gate fires, per lane and per vector
//   gdump <file>  gamma R2's depth need over a dense (a, x) grid
//   bits <file> [reference files...]  hex beta_p/beta_q/gamma_p/gamma_q over
//                 the reference rows plus 400k random points, to diff builds
// The router replica and the all-dd CF copy below are frozen copies for
// comparison; they do not track the kernel.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <random>
#include <vector>

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "tests/bench_igamma_components.cpp"
#include "hwy/foreach_target.h"  // IWYU pragma: keep

#include "src/beta-inl.h"
#include "src/gamma-inl.h"
#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {

using DTag = op::ScalableTag<double>;
using V = op::V<DTag>;

template <class Body>
static HWY_NOINLINE void Drive3(const double* a, const double* b,
                                const double* x, double* out, size_t n,
                                Body body) {
  const DTag d;
  const size_t N = op::Lanes(d);
  size_t i = 0;
  for (; i + N <= n; i += N) {
    op::Store(body(d, op::Load(d, a + i), op::Load(d, b + i),
                   op::Load(d, x + i)),
              d, out + i);
  }
  if (i < n) {
    const size_t m = n - i;
    op::StoreN(body(d, op::LoadN(d, a + i, m), op::LoadN(d, b + i, m),
                    op::LoadN(d, x + i, m)),
               d, out + i, m);
  }
}

// GammaCfRecip with a runtime depth (same arithmetic as the shipped core).
template <class D>
HWY_NOINLINE Dd<D> CfRecipN(D d, op::V<D> a, op::V<D> x, int nd) {
  const auto ma = op::Neg(a);
  auto kd = op::Sub(op::Add(x, op::Set(d, static_cast<double>(2 * nd + 1))), a);
  for (int j = nd; j > kGammaCfDdSteps; --j) {
    const auto jv = op::Set(d, static_cast<double>(j));
    const auto num = op::Mul(jv, op::Sub(jv, a));
    const auto den =
        op::Sub(op::Add(x, op::Set(d, static_cast<double>(2 * j - 1))), a);
    kd = op::Sub(den, op::Div(num, kd));
  }
  Dd<D> k{kd, op::Zero(d)};
  for (int j = std::min(nd, kGammaCfDdSteps); j >= 1; --j) {
    const auto jv = op::Set(d, static_cast<double>(j));
    const auto num = DdMulD(d, TwoSum(d, jv, ma), jv);
    const auto den = DdAddD(
        d, TwoSum(d, x, op::Set(d, static_cast<double>(2 * j - 1))), ma);
    const auto q = DdMul(d, num, DdRecipDd(d, k));
    k = DdAdd(d, den, Dd<D>{op::Neg(q.hi), op::Neg(q.lo)});
  }
  return DdRecipDd(d, k);
}

// BetaR2Cf (no huge-c prescale: bench bands never reach it) at runtime depth.
template <class D>
HWY_NOINLINE Dd<D> BetaCfN(D d, op::V<D> a, op::V<D> b, Dd<D> c, Dd<D> xi,
                           int nd) {
  const auto one = op::Set(d, 1.0);
  const auto zero = op::Zero(d);
  Dd<D> f{one, zero};
  for (int k = nd; k >= 1; --k) {
    const double mv = static_cast<double>(k / 2);
    const auto m = op::Set(d, mv);
    Dd<D> dk{zero, zero};
    if ((k & 1) == 0) {
      const auto bx = DdSub(d, DdMulD(d, xi, b), DdMulD(d, xi, m));
      const auto r1 = DdDivDd(d, Dd<D>{m, zero},
                              TwoSum(d, a, op::Set(d, 2.0 * mv - 1.0)));
      const auto r2 = DdDivDd(d, bx, TwoSum(d, a, op::Set(d, 2.0 * mv)));
      dk = DdMul(d, r1, r2);
    } else {
      const auto cx = DdAdd(d, DdMul(d, xi, c), DdMulD(d, xi, m));
      const auto r1 =
          DdDivDd(d, TwoSum(d, a, m), TwoSum(d, a, op::Set(d, 2.0 * mv)));
      const auto r2 =
          DdDivDd(d, cx, TwoSum(d, a, op::Set(d, 2.0 * mv + 1.0)));
      dk = DdNeg(d, DdMul(d, r1, r2));
    }
    f = DdAddD(d, DdMul(d, dk, DdRecipDd(d, f)), one);
  }
  return DdRecipDd(d, f);
}

// Beta R2 orientation (pinned rule x < (a+1)/(c+2)); returns alpha, beta, xi.
template <class D>
HWY_INLINE void BetaOrient(D d, op::V<D> a, op::V<D> b, op::V<D> x,
                           op::V<D>& al, op::V<D>& be, Dd<D>& c, Dd<D>& xi) {
  const auto one = op::Set(d, 1.0);
  const auto cc = op::Add(a, b);
  const auto sw = op::Ge(x, op::Div(op::Add(a, one), op::Add(cc, op::Set(d, 2.0))));
  al = op::IfThenElse(sw, b, a);
  be = op::IfThenElse(sw, a, b);
  c = TwoSum(d, a, b);
  const auto y = TwoSum(d, one, op::Neg(x));
  xi = Dd<D>{op::IfThenElse(sw, y.hi, x), op::IfThenElse(sw, y.lo, op::Zero(d))};
}

// Gamma prefactor E = a log x - x - lgamma(a), then e^E as ExpDdFrac.
template <class D>
HWY_INLINE op::V<D> GammaPref(D d, op::V<D> a, op::V<D> x) {
  auto eb = DdMulD(d, LogDdAny(d, x), a);
  eb = DdAddD(d, eb, op::Neg(x));
  const auto lg = LgammaPosDd(d, a);
  eb = DdAdd(d, eb, Dd<D>{op::Neg(lg.hi), op::Neg(lg.lo)});
  const auto ea = GammaClampE(d, eb.hi, eb.lo);
  return ExpDdFrac(d, ea.hi, ea.lo).m.hi;
}


// BetaCfN with the deep steps in plain double and only the last kdd in dd
// (gamma's GammaCfRecip split, applied to beta).
template <class D>
HWY_NOINLINE Dd<D> BetaCfMixed(D d, op::V<D> a, op::V<D> b, Dd<D> c, Dd<D> xi,
                               int nd, int kdd) {
  const auto one = op::Set(d, 1.0);
  const auto zero = op::Zero(d);
  const auto x = xi.hi;
  const auto ch = c.hi;
  auto fd = one;
  int k = nd;
  for (; k > kdd; --k) {
    const double mv = static_cast<double>(k / 2);
    const auto m = op::Set(d, mv);
    op::V<D> dk;
    if ((k & 1) == 0) {
      dk = op::Mul(op::Div(m, op::Add(a, op::Set(d, 2.0 * mv - 1.0))),
                   op::Div(op::Mul(op::Sub(b, m), x), op::Add(a, op::Set(d, 2.0 * mv))));
    } else {
      dk = op::Neg(op::Mul(op::Div(op::Add(a, m), op::Add(a, op::Set(d, 2.0 * mv))),
                           op::Div(op::Mul(op::Add(ch, m), x),
                                   op::Add(a, op::Set(d, 2.0 * mv + 1.0)))));
    }
    fd = op::Add(one, op::Div(dk, fd));
  }
  Dd<D> f{fd, zero};
  for (; k >= 1; --k) {
    const double mv = static_cast<double>(k / 2);
    const auto m = op::Set(d, mv);
    Dd<D> dk{zero, zero};
    if ((k & 1) == 0) {
      const auto bx = DdSub(d, DdMulD(d, xi, b), DdMulD(d, xi, m));
      const auto r1 = DdDivDd(d, Dd<D>{m, zero},
                              TwoSum(d, a, op::Set(d, 2.0 * mv - 1.0)));
      const auto r2 = DdDivDd(d, bx, TwoSum(d, a, op::Set(d, 2.0 * mv)));
      dk = DdMul(d, r1, r2);
    } else {
      const auto cx = DdAdd(d, DdMul(d, xi, c), DdMulD(d, xi, m));
      const auto r1 =
          DdDivDd(d, TwoSum(d, a, m), TwoSum(d, a, op::Set(d, 2.0 * mv)));
      const auto r2 =
          DdDivDd(d, cx, TwoSum(d, a, op::Set(d, 2.0 * mv + 1.0)));
      dk = DdNeg(d, DdMul(d, r1, r2));
    }
    f = DdAddD(d, DdMul(d, dk, DdRecipDd(d, f)), one);
  }
  return DdRecipDd(d, f);
}

// Cost of one dd quadrature node of the #44 kind: t^(a-1) e^(-x t) w at a
// dd node, i.e. a dd log, a dd exp and the dd multiply-accumulate.
template <class D>
HWY_INLINE op::V<D> Quad(D d, op::V<D> a, op::V<D> x, int nodes) {
  Dd<D> s{op::Zero(d), op::Zero(d)};
  const auto am1 = op::Sub(a, op::Set(d, 1.0));
  for (int i = 1; i <= nodes; ++i) {
    const auto t = op::Set(d, (i - 0.5) / nodes);
    const auto lt = LogDdAny(d, t);
    auto e = DdMulD(d, lt, am1);
    e = DdAdd(d, e, DdMulD(d, Dd<D>{t, op::Zero(d)}, op::Neg(x)));
    const auto ex = ExpDdFrac(d, e.hi, e.lo);
    s = DdAdd(d, s, DdMulD(d, ex.m, op::Set(d, 1.0 / nodes)));
  }
  return s.hi;
}

#define ROW(name, expr)                                                    \
  static void name(const double* pa, const double* pb, const double* px,  \
                   double* o, size_t n) {                                  \
    Drive3(pa, pb, px, o, n, [](DTag d, V a, V b, V x) HWY_ATTR {          \
      (void)d; (void)a; (void)b; (void)x;                                  \
      return expr;                                                         \
    });                                                                    \
  }

ROW(FloorImpl, op::Add(a, x))
ROW(GpFullImpl, (GammaVec<true>(d, a, x)))
ROW(GqFullImpl, (GammaVec<false>(d, a, x)))
ROW(GPrefImpl, GammaPref(d, a, x))
ROW(GLgamImpl, LgammaPosDd(d, a).hi)
ROW(GLogImpl, LogDdAny(d, x).hi)
ROW(GExpImpl, ExpDdFrac(d, op::Neg(x), op::Mul(x, op::Set(d, 1e-17))).m.hi)
ROW(GSerImpl, GammaSeriesSum(d, a, x).hi)
ROW(GCf44Impl, GammaCfRecip(d, a, x).hi)
ROW(GCf32Impl, CfRecipN(d, a, x, 32).hi)
ROW(GCf16Impl, CfRecipN(d, a, x, 16).hi)
ROW(GTemImpl, GammaTemme(d, a, x).val.v)
ROW(GR4Impl, GammaSmallQ(d, a, x).m.hi)
ROW(BpFullImpl, (BetaVec<true>(d, a, b, x)))
ROW(BqFullImpl, (BetaVec<false>(d, a, b, x)))
ROW(BSerImpl, BetaR1Series(d, a, b, x).hi)
ROW(BLgdImpl, LgammaDiffDd(d, op::Max(a, b), op::Min(a, b)).hi)
ROW(BLogImpl, OutlinedLogDd(d, x).hi)
static op::V<DTag> BCfRow(DTag d, V a, V b, V x, int nd) {
  // Aggregate-initialized, never default-constructed (BetaPsiCore's note).
  V al = op::Zero(d), be = op::Zero(d);
  Dd<DTag> c{op::Zero(d), op::Zero(d)}, xi{op::Zero(d), op::Zero(d)};
  BetaOrient(d, a, b, x, al, be, c, xi);
  return nd == 0 ? BetaR2Cf(d, al, be, c, xi).hi
                 : BetaCfN(d, al, be, c, xi, nd).hi;
}
ROW(BCf64Impl, BCfRow(d, a, b, x, 0))
ROW(BCf32Impl, BCfRow(d, a, b, x, 32))
ROW(BCf16Impl, BCfRow(d, a, b, x, 16))
static op::V<DTag> BCfMixRow(DTag d, V a, V b, V x, int kdd) {
  // Aggregate-initialized, never default-constructed (BetaPsiCore's note).
  V al = op::Zero(d), be = op::Zero(d);
  Dd<DTag> c{op::Zero(d), op::Zero(d)}, xi{op::Zero(d), op::Zero(d)};
  BetaOrient(d, a, b, x, al, be, c, xi);
  return BetaCfMixed(d, al, be, c, xi, 64, kdd).hi;
}
ROW(BCfM8Impl, BCfMixRow(d, a, b, x, 8))
ROW(BCfM4Impl, BCfMixRow(d, a, b, x, 4))
ROW(Quad24Impl, Quad(d, a, x, 24))
#undef ROW

// ---- depth analysis: one point, broadcast to every lane ------------------
// Smallest depth whose dd result is bit-identical (hi and lo) to the shipped
// depth's. Bit-identical dd => bit-identical output.
static int GammaCfNeed(double a, double x) {
  const DTag d;
  const auto av = op::Set(d, a), xv = op::Set(d, x);
  const auto ref = CfRecipN(d, av, xv, detail::kGammaCfN);
  const double rh = op::GetLane(ref.hi), rl = op::GetLane(ref.lo);
  for (int nd = 1; nd <= detail::kGammaCfN; ++nd) {
    const auto r = CfRecipN(d, av, xv, nd);
    if (op::GetLane(r.hi) == rh && op::GetLane(r.lo) == rl) {
      // Monotone check: confirm every deeper depth agrees too.
      bool ok = true;
      for (int m = nd + 1; m < detail::kGammaCfN && ok; ++m) {
        const auto s = CfRecipN(d, av, xv, m);
        ok = op::GetLane(s.hi) == rh && op::GetLane(s.lo) == rl;
      }
      if (ok) return nd;
    }
  }
  return detail::kGammaCfN;
}
static int BetaCfNeed(double a, double b, double x) {
  const DTag d;
  // Aggregate-initialized, never default-constructed (BetaPsiCore's note).
  V al = op::Zero(d), be = op::Zero(d);
  Dd<DTag> c{op::Zero(d), op::Zero(d)}, xi{op::Zero(d), op::Zero(d)};
  BetaOrient(d, op::Set(d, a), op::Set(d, b), op::Set(d, x), al, be, c, xi);
  const auto ref = BetaCfN(d, al, be, c, xi, detail::kBetaN2);
  const double rh = op::GetLane(ref.hi), rl = op::GetLane(ref.lo);
  for (int nd = 1; nd <= detail::kBetaN2; ++nd) {
    const auto r = BetaCfN(d, al, be, c, xi, nd);
    if (op::GetLane(r.hi) == rh && op::GetLane(r.lo) == rl) {
      bool ok = true;
      for (int m = nd + 1; m < detail::kBetaN2 && ok; ++m) {
        const auto s = BetaCfN(d, al, be, c, xi, m);
        ok = op::GetLane(s.hi) == rh && op::GetLane(s.lo) == rl;
      }
      if (ok) return nd;
    }
  }
  return detail::kBetaN2;
}
// Iteration at which the series' per-lane freeze fires (cap if never).
static int GammaSerNeed(double a, double x) {
  const DTag d;
  const auto one = op::Set(d, 1.0), zero = op::Zero(d);
  const auto av = op::Set(d, a);
  const auto xc = op::Min(op::Set(d, x), op::Add(av, one));
  Dd<DTag> t{one, zero}, s{one, zero};
  for (int n = 1; n <= detail::kGammaSeriesN; ++n) {
    const auto apn = TwoSum(d, av, op::Set(d, static_cast<double>(n)));
    t = DdMul(d, t, DdMulD(d, DdRecipDd(d, apn), xc));
    s = DdAdd(d, s, t);
    if (op::GetLane(t.hi) < op::GetLane(s.hi) * kGammaFreezeEps) return n;
  }
  return detail::kGammaSeriesN;
}
static int BetaSerNeed(double a, double b, double x, double eps) {
  const DTag d;
  const auto one = op::Set(d, 1.0), zero = op::Zero(d);
  const auto av = op::Set(d, a), bv = op::Set(d, b), xv = op::Set(d, x);
  Dd<DTag> t{one, zero}, s{zero, zero};
  for (int n = 1; n <= detail::kBetaN1; ++n) {
    const auto nv = op::Set(d, static_cast<double>(n));
    const auto nb = TwoSum(d, nv, op::Neg(bv));
    const Dd<DTag> rn{op::Set(d, detail::kBetaRecipNHi[n - 1]),
                      op::Set(d, detail::kBetaRecipNLo[n - 1])};
    t = DdMulD(d, DdMul(d, t, DdMul(d, nb, rn)), xv);
    const auto contrib = DdMul(d, t, DdRecipDd(d, TwoSum(d, av, nv)));
    s = DdAdd(d, s, contrib);
    if (std::fabs(op::GetLane(contrib.hi)) <
        std::fabs(op::GetLane(s.hi)) * eps)
      return n;
  }
  return detail::kBetaN1;
}


// Mixed-precision beta CF against the shipped all-dd depth: is the dd F
// bit-identical, and how far apart are they in units of 2^-53 * F?
// Hunt: shipped BetaR2Cf (reference) against the mixed split, one oriented
// point. F = |reference|; err in units of 2^-53 * F.
// Scalar replica of BetaVec's router: true iff (a, b, x) reaches BetaR2Cf
// in the NATURAL orientation (not R4, R1, R3 or the gamma-limit slice).
static int BetaIsR2NatImpl(double a, double b, double x) {
  const double y = 1.0 - x;  // exact enough for routing
  const double tau = std::min(a, b), bmax = std::max(a, b);
  const bool sw4 = b < a;
  const double xt = sw4 ? y : x;
  const double thr = 1.0 / (1.0 + (b + 1.0) / (a + 1.0));
  const bool r2tf = (x < thr) != sw4;
  const bool r4 = detail::kBetaEpsR4 >= tau &&
                  detail::kBetaLn2 >= tau * std::fabs(std::log(xt)) &&
                  (detail::kBetaXi1 >= xt || r2tf) && detail::kBetaB1 >= bmax * xt;
  if (r4) return 0;
  const bool r1n = detail::kBetaXi1 >= x && detail::kBetaB1 >= b * x;
  const bool r1s = detail::kBetaXi1 >= y && detail::kBetaB1 >= a * y;
  if (r1n || r1s) return 0;
  const double nu = 1.0 / (1.0 / a + 1.0 / b);
  const double rat1 = x * (1.0 + b / a), rat2 = y * (1.0 + a / b);
  const bool band = rat1 >= detail::kBetaXiRatioLo && detail::kBetaXiRatioHi >= rat1 &&
                    rat2 >= detail::kBetaXiRatioLo && detail::kBetaXiRatioHi >= rat2;
  const bool glr = bmax >= detail::kBetaGammaLim && tau >= detail::kBetaGlRidgeMin;
  if ((nu >= detail::kBetaTRidge || glr) && band) return 0;
  if (bmax >= detail::kBetaGammaLim && detail::kBetaGammaLim > tau) return 0;
  return x < thr ? 1 : 2;
}

// Hunt on RAW (a, b, x): route exactly, orient as BetaVec does (xi as the
// exact dd 1 - x when swapped), then shipped BetaR2Cf against the mixed split.
// *F < 0 when the point does not reach the CF.
static void BetaHuntImpl(double a0, double b0, double x0, int kdd, double* F,
                         double* err, int* orient) {
  *orient = BetaIsR2NatImpl(a0, b0, x0);
  *F = -1;
  *err = 0;
  if (*orient == 0) return;
  const DTag d;
  const bool sw = *orient == 2;
  const auto one = op::Set(d, 1.0);
  const auto y = TwoSum(d, one, op::Set(d, -x0));
  const auto a = op::Set(d, sw ? b0 : a0), b = op::Set(d, sw ? a0 : b0);
  const auto c = TwoSum(d, a, b);
  const Dd<DTag> xd = sw ? y : Dd<DTag>{op::Set(d, x0), op::Zero(d)};
  const auto r = BetaCfN(d, a, b, c, xd, detail::kBetaN2);  // pre-#42 all-dd
  const auto m = kdd == 0 ? BetaR2Cf(d, a, b, c, xd)
                          : BetaCfMixed(d, a, b, c, xd, detail::kBetaN2, kdd);
  const double rh = op::GetLane(r.hi), rl = op::GetLane(r.lo);
  const double mh = op::GetLane(m.hi), ml = op::GetLane(m.lo);
  *F = std::fabs(rh);
  *err = std::fabs((mh - rh) + (ml - rl)) / (std::fabs(rh) * 0x1p-53);
}

static void BetaMixedCmpImpl(double a, double b, double x, int kdd,
                             int* same, double* rel) {
  const DTag d;
  // Aggregate-initialized, never default-constructed (BetaPsiCore's note).
  V al = op::Zero(d), be = op::Zero(d);
  Dd<DTag> c{op::Zero(d), op::Zero(d)}, xi{op::Zero(d), op::Zero(d)};
  BetaOrient(d, op::Set(d, a), op::Set(d, b), op::Set(d, x), al, be, c, xi);
  const auto r = BetaCfN(d, al, be, c, xi, detail::kBetaN2);
  const auto m = BetaCfMixed(d, al, be, c, xi, detail::kBetaN2, kdd);
  const double rh = op::GetLane(r.hi), rl = op::GetLane(r.lo);
  const double mh = op::GetLane(m.hi), ml = op::GetLane(m.lo);
  *same = (rh == mh && rl == ml);
  *rel = std::fabs((mh - rh) + (ml - rl)) / (std::fabs(rh) * 0x1p-53);
}

static int GammaCfNeedImpl(double a, double x) { return GammaCfNeed(a, x); }
static int BetaCfNeedImpl(double a, double b, double x) {
  return BetaCfNeed(a, b, x);
}
static int GammaSerNeedImpl(double a, double x) { return GammaSerNeed(a, x); }
static int BetaSerNeedImpl(double a, double b, double x, double e) {
  return BetaSerNeed(a, b, x, e);
}
static const char* BenchTargetNameImpl() { return hwy::TargetName(HWY_TARGET); }

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#if HWY_ONCE

#include "corvus/corvus.h"
#include "expect_target.h"

namespace corvus {
using RowFn = void (*)(const double*, const double*, const double*, double*,
                       size_t);
#define EXP(name)                                                          \
  HWY_EXPORT(name##Impl);                                                  \
  static void name(const double* a, const double* b, const double* x,      \
                   double* o, size_t n) {                                  \
    HWY_DYNAMIC_DISPATCH(name##Impl)(a, b, x, o, n);                       \
  }
EXP(Floor)
EXP(GpFull)
EXP(GqFull)
EXP(GPref)
EXP(GLgam)
EXP(GLog)
EXP(GExp)
EXP(GSer)
EXP(GCf44)
EXP(GCf32)
EXP(GCf16)
EXP(GTem)
EXP(GR4)
EXP(BpFull)
EXP(BqFull)
EXP(BSer)
EXP(BLgd)
EXP(BLog)
EXP(BCf64)
EXP(BCf32)
EXP(BCf16)
EXP(BCfM8)
EXP(BCfM4)
EXP(Quad24)
#undef EXP
HWY_EXPORT(GammaCfNeedImpl);
HWY_EXPORT(BetaMixedCmpImpl);
HWY_EXPORT(BetaHuntImpl);
HWY_EXPORT(BetaIsR2NatImpl);
HWY_EXPORT(BetaCfNeedImpl);
HWY_EXPORT(GammaSerNeedImpl);
HWY_EXPORT(BetaSerNeedImpl);
HWY_EXPORT(BenchTargetNameImpl);
static int GammaCfNeedD(double a, double x) {
  return HWY_DYNAMIC_DISPATCH(GammaCfNeedImpl)(a, x);
}
static int BetaCfNeedD(double a, double b, double x) {
  return HWY_DYNAMIC_DISPATCH(BetaCfNeedImpl)(a, b, x);
}
static int GammaSerNeedD(double a, double x) {
  return HWY_DYNAMIC_DISPATCH(GammaSerNeedImpl)(a, x);
}
static int BetaSerNeedD(double a, double b, double x, double e) {
  return HWY_DYNAMIC_DISPATCH(BetaSerNeedImpl)(a, b, x, e);
}
static void BetaMixedCmpD(double a, double b, double x, int k, int* s, double* r) {
  HWY_DYNAMIC_DISPATCH(BetaMixedCmpImpl)(a, b, x, k, s, r);
}
static void BetaHuntD(double a, double b, double x, int k, double* F, double* e,
                      int* o) {
  HWY_DYNAMIC_DISPATCH(BetaHuntImpl)(a, b, x, k, F, e, o);
}
static const char* BenchTargetName() {
  return HWY_DYNAMIC_DISPATCH(BenchTargetNameImpl)();
}
}  // namespace corvus

namespace {

using Clock = std::chrono::steady_clock;
volatile double g_sink;

struct Set3 {
  std::vector<double> a, b, x;
};

using Sampler = void (*)(std::mt19937_64&, double&, double&, double&);

Set3 Make(Sampler s, size_t n) {
  std::mt19937_64 rng(20261029);
  Set3 r{std::vector<double>(n), std::vector<double>(n), std::vector<double>(n)};
  for (size_t i = 0; i < n; ++i) s(rng, r.a[i], r.b[i], r.x[i]);
  return r;
}

double Ns(corvus::RowFn fn, const Set3& s, std::vector<double>& o, int reps) {
  const size_t n = s.a.size();
  fn(s.a.data(), s.b.data(), s.x.data(), o.data(), n);
  std::vector<double> t(static_cast<size_t>(reps));
  for (auto& v : t) {
    const auto t0 = Clock::now();
    fn(s.a.data(), s.b.data(), s.x.data(), o.data(), n);
    const auto t1 = Clock::now();
    g_sink = o[n / 2];
    v = std::chrono::duration<double, std::nano>(t1 - t0).count() /
        static_cast<double>(n);
  }
  std::nth_element(t.begin(), t.begin() + reps / 2, t.end());
  return t[static_cast<size_t>(reps) / 2];
}

struct Row {
  const char* name;
  corvus::RowFn fn;
};

void Band(const char* label, Sampler s, std::initializer_list<Row> rows) {
  const auto set = Make(s, 100000);
  std::vector<double> o(set.a.size());
  std::printf("%s\n", label);
  for (const auto& r : rows) {
    std::printf("  %-10s %9.2f ns/el\n", r.name, Ns(r.fn, set, o, 31));
  }
}

double LogU(std::mt19937_64& g, double lo, double hi) {
  return std::exp(std::uniform_real_distribution<double>(std::log(lo),
                                                         std::log(hi))(g));
}
double U(std::mt19937_64& g, double lo, double hi) {
  return std::uniform_real_distribution<double>(lo, hi)(g);
}
// bench_gamma's bands.
void SG1(std::mt19937_64& g, double& a, double& b, double& x) {
  a = LogU(g, 2, 19); b = 0; x = a * U(g, 0.2, 1.0);
}
void SG2(std::mt19937_64& g, double& a, double& b, double& x) {
  a = LogU(g, 2, 19); b = 0; x = a * U(g, 3.0, 20.0);
}
void SG3(std::mt19937_64& g, double& a, double& b, double& x) {
  a = LogU(g, 50, 1e6); b = 0; x = a * U(g, 0.7, 1.4);
}
void SG4(std::mt19937_64& g, double& a, double& b, double& x) {
  a = LogU(g, 0.05, 1.5); b = 0; x = a * U(g, 0.5, 2.5);
}
// Wider R2 band: a in [0.1, 19], x past a+1 out to 50(a+1) (log-uniform).
void SG2w(std::mt19937_64& g, double& a, double& b, double& x) {
  a = LogU(g, 0.1, 19); b = 0; x = (a + 1) * LogU(g, 1.05, 50);
}
// Near-wall R2: x just past a+1 -- the CF's slowest corner.
void SG2n(std::mt19937_64& g, double& a, double& b, double& x) {
  a = LogU(g, 2, 19); b = 0; x = (a + 1) * U(g, 1.0, 1.3);
}
// bench_beta's bands.
void SB1(std::mt19937_64& g, double& a, double& b, double& x) {
  a = LogU(g, 2, 19); b = LogU(g, 2, 19); x = U(g, 0.01, 0.30);
}
void SB2(std::mt19937_64& g, double& a, double& b, double& x) {
  a = LogU(g, 2, 19); b = LogU(g, 2, 19); x = U(g, 0.46, 0.54);
}
// Wider R2: a, b in [0.5, 30], x uniform in (0.05, 0.95).
void SB2w(std::mt19937_64& g, double& a, double& b, double& x) {
  a = LogU(g, 0.5, 30); b = LogU(g, 0.5, 30); x = U(g, 0.05, 0.95);
}

void Hist(const char* label, std::vector<int> need, int cap) {
  std::vector<int> s = need;
  std::sort(s.begin(), s.end());
  auto pct = [&](double p) {
    return s[std::min(s.size() - 1, static_cast<size_t>(p * s.size()))];
  };
  double mean = 0;
  for (int v : s) mean += v;
  mean /= static_cast<double>(s.size());
  std::printf("%-34s cap %2d | lane mean %5.1f p50 %2d p90 %2d p99 %2d max %2d |",
              label, cap, mean, pct(0.5), pct(0.9), pct(0.99), s.back());
  for (int w : {2, 4, 8}) {
    double vm = 0;
    size_t k = 0;
    for (size_t i = 0; i + w <= need.size(); i += w, ++k) {
      vm += *std::max_element(need.begin() + i, need.begin() + i + w);
    }
    std::printf(" w%d %5.1f", w, vm / static_cast<double>(k));
  }
  std::printf("\n");
}

// Large-F hunt over BetaR2Cf's oriented input domain. xi is taken exact
// (xi.lo = 0); the swapped orientation's inexact xi is a separate question.
int Hunt() {
  std::mt19937_64 g(20261030);
  const int kK[] = {0, 16, 24};
  const int nk = 3;
  const double kBin[] = {0, 4, 10, 20, 30, 40, 64};
  const int nb = 6;
  double mx[2][5][6] = {};
  long cnt[2][6] = {};
  double worst[2][6][3] = {};
  auto one = [&](double a, double b, double x) {
    if (!(x > 0 && x < 1)) return;
    double F, e;
    int o;
    corvus::BetaHuntD(a, b, x, kK[0], &F, &e, &o);
    if (!(F > 0) || !std::isfinite(F)) return;
    const int s = o - 1;
    const double l = std::log2(F);
    int bi = 0;
    while (bi + 1 < nb && l >= kBin[bi + 1]) ++bi;
    ++cnt[s][bi];
    for (int k = 0; k < nk; ++k) {
      if (k) corvus::BetaHuntD(a, b, x, kK[k], &F, &e, &o);
      if (e > mx[s][k][bi]) {
        mx[s][k][bi] = e;
        if (k == 0) { worst[s][bi][0] = a; worst[s][bi][1] = b; worst[s][bi][2] = x; }
      }
    }
  };
  for (int i = 0; i < 600000; ++i) {
    const double a = LogU(g, 0.02, 0x1p47);
    const double b = LogU(g, 0.02, 0x1p47);
    const double p = a / (a + b);
    switch (i % 6) {
      case 0: one(a, b, U(g, 0.0, 1.0)); break;
      case 1: one(a, b, p * (1 + (U(g, -1, 1) > 0 ? 1 : -1) * LogU(g, 1e-14, 0.9))); break;
      case 2: one(a, b, LogU(g, 8, 800) / b); break;          // gamma line
      case 3: one(a, b, 1 - LogU(g, 8, 800) / a); break;      // its mirror
      case 4: one(LogU(g, 0.02, 40), b, LogU(g, 1e-300, 1)); break;
      case 5: one(a, LogU(g, 0.02, 40), 1 - LogU(g, 1e-17, 1)); break;
    }
  }
  one(0.995, 7896.5, 7.4438e-6);  // swapped witness, raw form
  one(1e10, 1, 1 - 8e-8);
  for (int s = 0; s < 2; ++s) {
    std::printf("beta R2 CF, %s orientation: mixed vs shipped, max err in 2^-53 F\n",
                s ? "SWAPPED" : "natural");
    std::printf("%10s %8s", "log2 F", "points");
    for (int k : kK) std::printf("   K=%-2d   ", k);
    std::printf("\n");
    for (int b = 0; b < nb; ++b) {
      std::printf("[%2.0f,%2.0f) %8ld", kBin[b], kBin[b + 1], cnt[s][b]);
      for (int k = 0; k < nk; ++k) std::printf(" %10.3g", mx[s][k][b]);
      std::printf("\n");
    }
    for (int b = 0; b < nb; ++b)
      if (cnt[s][b])
        std::printf("  worst KERNEL bin %d: a=%.17g b=%.17g x=%.17g\n", b,
                    worst[s][b][0], worst[s][b][1], worst[s][b][2]);
  }
  return 0;
}

// Dense dump of gamma R2's bitwise depth need: "a x need" per line over the
// R2 domain {a < 20, x > a+1} u {a >= 20, x >= 2a}, x up to 1000.
int GDump(const char* path) {
  std::ofstream f(path);
  if (!f) return 1;
  char line[96];
  for (int ia = 0; ia < 160; ++ia) {
    const double a = std::exp(std::log(1e-3) + (std::log(3000.0) - std::log(1e-3)) * ia / 159.0);
    const double x0 = a < 20 ? a + 1 : 2 * a;
    for (int ix = 0; ix < 160; ++ix) {
      // log-spaced in x - x0 from 1e-9*x0 to 1000, plus the wall itself
      const double t = ix / 159.0;
      const double x = ix == 0 ? std::nextafter(x0, 2 * x0 + 1)
                               : x0 + std::exp(std::log(1e-6 * (x0 + 1)) +
                                               (std::log(1000.0) - std::log(1e-6 * (x0 + 1))) * t);
      if (x > 1000) continue;
      std::snprintf(line, sizeof line, "%.17g %.17g %d\n", a, x,
                    corvus::GammaCfNeedD(a, x));
      f << line;
    }
  }
  return 0;
}

// Bit dump: beta_p/beta_q/gamma_p/gamma_q over the reference rows named on the
// command line plus 400k wide random points, as hex, to compare two builds.
int Bits(int argc, char** argv) {
  std::vector<double> a, b, x;
  for (int i = 3; i < argc; ++i) {
    std::ifstream in(argv[i]);
    std::string t0, t1, t2, rest;
    while (in >> t0 >> t1 >> t2) {
      std::getline(in, rest);
      if (t0[0] == '#') continue;
      a.push_back(std::strtod(t0.c_str(), nullptr));
      b.push_back(std::strtod(t1.c_str(), nullptr));
      x.push_back(std::strtod(t2.c_str(), nullptr));
    }
  }
  std::mt19937_64 g(7);
  for (int i = 0; i < 400000; ++i) {
    a.push_back(LogU(g, 1e-3, 1e6));
    b.push_back(LogU(g, 1e-3, 1e6));
    x.push_back(i % 3 == 0 ? U(g, 0, 1) : i % 3 == 1 ? LogU(g, 1e-12, 0.5)
                                                    : 1 - LogU(g, 1e-12, 0.5));
  }
  const size_t n = a.size();
  std::vector<double> o1(n), o2(n), o3(n), o4(n), xg(n);
  for (size_t i = 0; i < n; ++i) xg[i] = x[i] * (a[i] + 1) * 2;
  corvus::beta_p(a, b, x, o1);
  corvus::beta_q(a, b, x, o2);
  corvus::gamma_p(a, xg, o3);
  corvus::gamma_q(a, xg, o4);
  std::ofstream f(argv[2]);
  char line[160];
  for (size_t i = 0; i < n; ++i) {
    std::snprintf(line, sizeof line, "%a %a %a %a\n", o1[i], o2[i], o3[i], o4[i]);
    f << line;
  }
  std::printf("%zu points written\n", n);
  return 0;
}

// Gate statistics: fraction of R2-routed lanes, and of w-lane vectors, with
// at least one lane at F >= T, for the bench band and the wide band.
int FStat() {
  for (auto [lbl, smp] : {std::pair<const char*, Sampler>{"bench R2 band", SB2},
                          {"wide (a,b .5..30)", SB2w}}) {
    const auto set = Make(smp, 40000);
    std::vector<double> F;
    for (size_t i = 0; i < set.a.size(); ++i) {
      double f, e;
      int o;
      corvus::BetaHuntD(set.a[i], set.b[i], set.x[i], 16, &f, &e, &o);
      if (f > 0) F.push_back(f);
    }
    std::printf("%s: %zu R2 lanes%s", lbl, F.size(), "\n");
    for (double T : {16.0, 64.0, 256.0, 1024.0}) {
      size_t lane = 0;
      for (double f : F) lane += f >= T;
      std::printf("  T=%6.0f lanes %5.1f%%", T, 100.0 * lane / F.size());
      for (size_t w : {2, 4, 8}) {
        size_t vec = 0, k = 0;
        for (size_t i = 0; i + w <= F.size(); i += w, ++k) {
          bool any = false;
          for (size_t j = 0; j < w; ++j) any |= F[i + j] >= T;
          vec += any;
        }
        std::printf("  w%zu %5.1f%%", w, 100.0 * vec / k);
      }
      std::printf("%s", "\n");
    }
  }
  return 0;
}

int Depth() {
  const size_t n = 4000;
  std::printf("vector-mean columns: mean over consecutive w-lane groups of "
              "the group max (what a per-lane late-start / early-exit loop "
              "would pay)\n");
  {
    for (auto [lbl, s] : {std::pair<const char*, Sampler>{"gamma R2 CF (bench band)", SG2},
                          {"gamma R2 CF (wide)", SG2w},
                          {"gamma R2 CF (near wall)", SG2n}}) {
      const auto set = Make(s, n);
      std::vector<int> need(n);
      for (size_t i = 0; i < n; ++i)
        need[i] = corvus::GammaCfNeedD(set.a[i], set.x[i]);
      Hist(lbl, need, 44);
    }
    for (auto [lbl, s] : {std::pair<const char*, Sampler>{"gamma R1 series (bench band)", SG1}}) {
      const auto set = Make(s, n);
      std::vector<int> need(n);
      for (size_t i = 0; i < n; ++i)
        need[i] = corvus::GammaSerNeedD(set.a[i], set.x[i]);
      Hist(lbl, need, 64);
    }
    for (auto [lbl, s] : {std::pair<const char*, Sampler>{"beta R2 CF (bench band)", SB2},
                          {"beta R2 CF (wide)", SB2w}}) {
      const auto set = Make(s, n / 4);
      std::vector<int> need(n / 4);
      for (size_t i = 0; i < n / 4; ++i)
        need[i] = corvus::BetaCfNeedD(set.a[i], set.b[i], set.x[i]);
      Hist(lbl, need, 64);
    }
    for (auto [lbl, s] : {std::pair<const char*, Sampler>{"bench", SB2},
                          {"wide", SB2w}}) {
      const auto set = Make(s, n);
      for (int kdd : {2, 4, 8, 16}) {
        size_t same = 0;
        double mx = 0, sum = 0;
        for (size_t i = 0; i < n; ++i) {
          int sm; double rel;
          corvus::BetaMixedCmpD(set.a[i], set.b[i], set.x[i], kdd, &sm, &rel);
          same += sm; mx = std::max(mx, rel); sum += rel;
        }
        std::printf("beta R2 CF mixed (%s) last %2d dd: dd-identical %5.1f%%, |dF| mean %.4f max %.4f (units 2^-53 F)\n",
                    lbl, kdd, 100.0 * same / n, sum / n, mx);
      }
    }
    const auto set = Make(SB1, n);
    for (double e : {0x1p-105, 0x1p-60}) {
      std::vector<int> need(n);
      for (size_t i = 0; i < n; ++i)
        need[i] = corvus::BetaSerNeedD(set.a[i], set.b[i], set.x[i], e);
      Hist(e < 1e-20 ? "beta R1 series eps 2^-105" : "beta R1 series eps 2^-60",
           need, 64);
    }
  }
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  if (!corvus_test::ReportAndCheckTarget()) return 2;
  if (std::strcmp(corvus::BenchTargetName(), corvus::active_target()) != 0) {
    std::fprintf(stderr, "FAIL: TU target '%s' != library target '%s'\n",
                 corvus::BenchTargetName(), corvus::active_target());
    return 2;
  }
  if (argc > 1 && std::strcmp(argv[1], "depth") == 0) return Depth();
  if (argc > 1 && std::strcmp(argv[1], "hunt") == 0) return Hunt();
  if (argc > 1 && std::strcmp(argv[1], "fstat") == 0) return FStat();
  if (argc > 2 && std::strcmp(argv[1], "bits") == 0) return Bits(argc, argv);
  if (argc > 2 && std::strcmp(argv[1], "gdump") == 0) return GDump(argv[2]);
  using namespace corvus;
  Band("gamma R1 (a 2..19, x/a .2..1)", SG1,
       {{"floor", Floor}, {"full p", GpFull}, {"pref+exp", GPref},
        {"lgamma", GLgam}, {"log_dd", GLog}, {"exp_dd", GExp},
        {"series", GSer}});
  Band("gamma R2 (a 2..19, x/a 3..20)", SG2,
       {{"floor", Floor}, {"full q", GqFull}, {"pref+exp", GPref},
        {"cf 44", GCf44}, {"cf 32", GCf32}, {"cf 16", GCf16}});
  Band("gamma R3 (a 50..1e6, x/a .7..1.4)", SG3,
       {{"floor", Floor}, {"full p", GpFull}, {"temme", GTem}});
  Band("gamma R4 (a .05..1.5, x/a .5..2.5)", SG4,
       {{"floor", Floor}, {"full q", GqFull}, {"smallq", GR4},
        {"pref+exp", GPref}});
  Band("beta R1 (a,b 2..19, x .01..0.3)", SB1,
       {{"floor", Floor}, {"full p", BpFull}, {"series", BSer},
        {"lgdiff", BLgd}, {"log_dd", BLog}});
  Band("beta R2 (a,b 2..19, x .46..54)", SB2,
       {{"floor", Floor}, {"full q", BqFull}, {"cf 64", BCf64},
        {"cf 32", BCf32}, {"cf 16", BCf16}, {"cfmix 8", BCfM8},
        {"cfmix 4", BCfM4}, {"lgdiff", BLgd}});
  Band("quadrature node cost (gamma R2 band)", SG2,
       {{"floor", Floor}, {"quad 24", Quad24}});
  return 0;
}

#endif  // HWY_ONCE
