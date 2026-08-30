// Per-component cost attribution for the lgamma kernel (#30). Same rules as
// bench_lgamma: manual, quiet machine, Release build only — run it under
// tools/quiet_bench.* with -Targets bench_lgamma_components.
//
// This TU compiles src/lgamma-inl.h itself (kernel-test pattern) and drives
// each component of a band in isolation: the shipped library is untouched, so
// by construction nothing here can move a ULP table.
//
// READING THE NUMBERS. Every driver has the same load/compute/store shape, so
// rows within a band are comparable, but an isolated component still pays the
// full per-element memory traffic the assembled kernel pays once — subtract
// the "floor" row before reasoning about a component's share, and expect the
// component rows to overlap rather than sum to the full-band row (the
// assembled kernel overlaps their dependency chains; a sum-of-parts above the
// whole is normal, a component alone near the whole is the signal #30 wants).
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>
#include <vector>

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "tests/bench_lgamma_components.cpp"
#include "hwy/foreach_target.h"  // IWYU pragma: keep

#include "src/lgamma-inl.h"
#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {

// One masked loop shape for every driver, mirroring src/driver-inl.h. The
// kernel body is a lambda V -> V so each component driver below is one line.
template <class Body>
HWY_NOINLINE void DriveComponent(const double* in, double* out, size_t n,
                                 Body body) {
  const op::ScalableTag<double> d;
  const size_t N = op::Lanes(d);
  size_t i = 0;
  for (; i + N <= n; i += N) {
    op::Store(body(d, op::Load(d, in + i)), d, out + i);
  }
  if (i < n) {
    const size_t m = n - i;
    op::StoreN(body(d, op::LoadN(d, in + i, m)), d, out + i, m);
  }
}

using DTag = op::ScalableTag<double>;
using V = op::V<DTag>;

// --- floor: load + store only, the driver/memory cost every row pays. ------
void FloorImpl(const double* in, double* out, size_t n) {
  DriveComponent(in, out, n, [](DTag, V x) HWY_ATTR { return x; });
}

// --- zone band [0.5, 2.5]: is the Horner chain the zone's cost? -----------
// The shipped all-zone fast path (select, exact t, bracket, final DdMulD).
void ZoneFullImpl(const double* in, double* out, size_t n) {
  DriveComponent(in, out, n, [](DTag d, V x) HWY_ATTR {
    const auto c1 = op::Lt(x, op::Set(d, detail::kLgammaZoneMid));
    const auto t = op::IfThenElse(c1, op::Sub(x, op::Set(d, 1.0)),
                                  op::Sub(x, op::Set(d, 2.0)));
    return DdToDouble(DdMulD(d, ZoneBracket(d, t, c1), t));
  });
}

// The 32-coefficient double Horner exactly as ZoneBracket runs it, per-lane
// Sel2 coefficients included; stops before the dd lead-term ladder.
void ZoneHornerImpl(const double* in, double* out, size_t n) {
  DriveComponent(in, out, n, [](DTag d, V x) HWY_ATTR {
    const auto c1 = op::Lt(x, op::Set(d, detail::kLgammaZoneMid));
    const auto t = op::IfThenElse(c1, op::Sub(x, op::Set(d, 1.0)),
                                  op::Sub(x, op::Set(d, 2.0)));
    const auto* co = detail::kLgammaZoneCoef;
    auto s = Sel2(d, c1, co[0][detail::kLgammaZoneNCoef - 1],
                  co[1][detail::kLgammaZoneNCoef - 1]);
    for (int k = detail::kLgammaZoneNCoef - 2; k >= 0; --k) {
      s = op::MulAdd(s, t, Sel2(d, c1, co[0][k], co[1][k]));
    }
    return s;
  });
}

// The same Horner on ONE coefficient set (no Sel2): the delta against
// zone_horner is what the per-coefficient select machinery costs.
void ZoneHorner1cImpl(const double* in, double* out, size_t n) {
  DriveComponent(in, out, n, [](DTag d, V x) HWY_ATTR {
    const auto t = op::Sub(x, op::Set(d, 1.0));
    const auto* co = detail::kLgammaZoneCoef;
    auto s = op::Set(d, co[0][detail::kLgammaZoneNCoef - 1]);
    for (int k = detail::kLgammaZoneNCoef - 2; k >= 0; --k) {
      s = op::MulAdd(s, t, op::Set(d, co[0][k]));
    }
    return s;
  });
}

// The dd lead-term ladder with a trivial s: the three-level DdAdd/DdMulD
// structure plus the closing DdMulD and rounding, without the Horner.
void ZoneLeadImpl(const double* in, double* out, size_t n) {
  DriveComponent(in, out, n, [](DTag d, V x) HWY_ATTR {
    const auto c1 = op::Lt(x, op::Set(d, detail::kLgammaZoneMid));
    const auto t = op::IfThenElse(c1, op::Sub(x, op::Set(d, 1.0)),
                                  op::Sub(x, op::Set(d, 2.0)));
    const auto* lh = detail::kLgammaZoneLeadHi;
    const auto* ll = detail::kLgammaZoneLeadLo;
    Dd<DTag> acc{Sel2(d, c1, lh[0][2], lh[1][2]),
                 Sel2(d, c1, ll[0][2], ll[1][2])};
    acc = DdAddD(d, acc, op::Mul(t, t));
    for (int k = 1; k >= 0; --k) {
      acc = DdAdd(d, Dd<DTag>{Sel2(d, c1, lh[0][k], lh[1][k]),
                              Sel2(d, c1, ll[0][k], ll[1][k])},
                  DdMulD(d, acc, t));
    }
    return DdToDouble(DdMulD(d, acc, t));
  });
}

// --- recurrence band (2.5, 8): walk-down vs the log on P vs the floor. ----
void RecFullImpl(const double* in, double* out, size_t n) {
  DriveComponent(in, out, n, [](DTag d, V x) HWY_ATTR {
    return DdToDouble(LgammaLow(d, x));
  });
}

// The masked walk-down exactly as LgammaLow runs it, nothing after it.
void RecWalkImpl(const double* in, double* out, size_t n) {
  DriveComponent(in, out, n, [](DTag d, V x) HWY_ATTR {
    const auto one = op::Set(d, 1.0);
    const auto xr = op::Min(x, op::Set(d, detail::kLgammaX0));
    auto y = x;
    Dd<DTag> prod{one, op::Zero(d)};
    for (int k = 1; k <= detail::kLgammaMidSteps; ++k) {
      const auto fire =
          op::Gt(xr, op::Set(d, detail::kLgammaZoneHi + (k - 1)));
      if (op::AllFalse(d, fire)) break;
      const auto step = op::Sub(xr, op::Set(d, static_cast<double>(k)));
      prod = DdMulD(d, prod, op::IfThenElse(fire, step, one));
      y = op::IfThenElse(fire, step, y);
    }
    return op::Add(prod.hi, y);
  });
}

// The outlined dd log alone. Band inputs (2.5, 8) sit inside the magnitude
// range the real P spans on this band, so the table-gather behaviour matches.
void RecLogImpl(const double* in, double* out, size_t n) {
  DriveComponent(in, out, n, [](DTag d, V x) HWY_ATTR {
    return OutlinedLogDd(d, Dd<DTag>{x, op::Zero(d)}).hi;
  });
}

// --- Stirling band [8, 1000]: the log vs the remainder polynomial. --------
void StirFullImpl(const double* in, double* out, size_t n) {
  DriveComponent(in, out, n, [](DTag d, V x) HWY_ATTR {
    return DdToDouble(LgammaStirling(d, x));
  });
}

void StirLogImpl(const double* in, double* out, size_t n) {
  DriveComponent(in, out, n, [](DTag d, V x) HWY_ATTR {
    return LogDd(d, x).hi;
  });
}

// The 1/x^2 remainder Horner plus its divisions, without the log or the dd
// assembly around it.
void StirRemImpl(const double* in, double* out, size_t n) {
  DriveComponent(in, out, n, [](DTag d, V x) HWY_ATTR {
    const auto w = op::Div(op::Set(d, 1.0), op::Mul(x, x));
    auto p = op::Set(d, detail::kLgammaStirCoef[detail::kLgammaStirNCoef - 1]);
    for (int k = detail::kLgammaStirNCoef - 2; k >= 0; --k) {
      p = op::MulAdd(p, w, op::Set(d, detail::kLgammaStirCoef[k]));
    }
    return op::Div(p, x);
  });
}

// --- reflection band [-30, -0.01]: the two extra logs vs LogSinc vs the ---
// positive pipeline it wraps.
void ReflFullImpl(const double* in, double* out, size_t n) {
  DriveComponent(in, out, n, [](DTag d, V x) HWY_ATTR {
    return LgammaVec(d, x);
  });
}

void ReflPosImpl(const double* in, double* out, size_t n) {
  DriveComponent(in, out, n, [](DTag d, V x) HWY_ATTR {
    return DdToDouble(LgammaPosDd(d, op::Abs(x)));
  });
}

// The reflection's two extra dd logs (-log|u| and -log(-x)) and their DdAdd,
// on the exact u the kernel forms.
void ReflLogsImpl(const double* in, double* out, size_t n) {
  DriveComponent(in, out, n, [](DTag d, V x) HWY_ATTR {
    const auto ax = op::Abs(x);
    const auto au = op::Abs(op::Sub(x, op::Round(x)));
    return DdAdd(d, OutlinedLogDd(d, au), OutlinedLogDd(d, ax)).hi;
  });
}

void ReflSincImpl(const double* in, double* out, size_t n) {
  DriveComponent(in, out, n, [](DTag d, V x) HWY_ATTR {
    return LogSinc(d, op::Sub(x, op::Round(x))).hi;
  });
}

// This TU carries its own dispatch table; main() asserts it agrees with the
// library's (test_exp_dd pattern).
const char* BenchTargetNameImpl() { return hwy::TargetName(HWY_TARGET); }

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#if HWY_ONCE

#include "corvus/corvus.h"
#include "expect_target.h"

namespace corvus {
HWY_EXPORT(FloorImpl);
HWY_EXPORT(ZoneFullImpl);
HWY_EXPORT(ZoneHornerImpl);
HWY_EXPORT(ZoneHorner1cImpl);
HWY_EXPORT(ZoneLeadImpl);
HWY_EXPORT(RecFullImpl);
HWY_EXPORT(RecWalkImpl);
HWY_EXPORT(RecLogImpl);
HWY_EXPORT(StirFullImpl);
HWY_EXPORT(StirLogImpl);
HWY_EXPORT(StirRemImpl);
HWY_EXPORT(ReflFullImpl);
HWY_EXPORT(ReflPosImpl);
HWY_EXPORT(ReflSincImpl);
HWY_EXPORT(ReflLogsImpl);
HWY_EXPORT(BenchTargetNameImpl);

// Dispatch inside namespace corvus (single-target cap collapse; see the
// shipped kernels).
using ComponentFn = void (*)(const double*, double*, size_t);
const char* BenchTargetName() {
  return HWY_DYNAMIC_DISPATCH(BenchTargetNameImpl)();
}
#define CORVUS_COMPONENT(name)                                     \
  void name(const double* in, double* out, size_t n) {             \
    HWY_DYNAMIC_DISPATCH(name##Impl)(in, out, n);                  \
  }
CORVUS_COMPONENT(Floor)
CORVUS_COMPONENT(ZoneFull)
CORVUS_COMPONENT(ZoneHorner)
CORVUS_COMPONENT(ZoneHorner1c)
CORVUS_COMPONENT(ZoneLead)
CORVUS_COMPONENT(RecFull)
CORVUS_COMPONENT(RecWalk)
CORVUS_COMPONENT(RecLog)
CORVUS_COMPONENT(StirFull)
CORVUS_COMPONENT(StirLog)
CORVUS_COMPONENT(StirRem)
CORVUS_COMPONENT(ReflFull)
CORVUS_COMPONENT(ReflPos)
CORVUS_COMPONENT(ReflSinc)
CORVUS_COMPONENT(ReflLogs)
#undef CORVUS_COMPONENT
}  // namespace corvus

namespace {

using Clock = std::chrono::steady_clock;

volatile double g_sink;

double NsPerElement(corvus::ComponentFn fn, const std::vector<double>& in,
                    std::vector<double>& out, int reps) {
  const size_t n = in.size();
  fn(in.data(), out.data(), n);
  std::vector<double> times(static_cast<size_t>(reps));
  for (auto& t : times) {
    const auto t0 = Clock::now();
    fn(in.data(), out.data(), n);
    const auto t1 = Clock::now();
    g_sink = out[n / 2];
    t = std::chrono::duration<double, std::nano>(t1 - t0).count() /
        static_cast<double>(n);
  }
  std::nth_element(times.begin(), times.begin() + reps / 2, times.end());
  return times[static_cast<size_t>(reps) / 2];
}

struct Row {
  const char* name;
  corvus::ComponentFn fn;
};

void RunBand(const char* label, double lo, double hi,
             std::initializer_list<Row> rows) {
  std::printf("band [%g, %g] (%s)\n", lo, hi, label);
  std::printf("%16s %14s %14s\n", "component", "10k ns/el", "1M ns/el");
  std::mt19937_64 rng(20260830);
  std::uniform_real_distribution<double> dist(lo, hi);
  std::vector<double> in_small(10000), in_large(1000000);
  for (auto& v : in_small) v = dist(rng);
  for (auto& v : in_large) v = dist(rng);
  std::vector<double> out_small(in_small.size()), out_large(in_large.size());
  for (const auto& row : rows) {
    const double small = NsPerElement(row.fn, in_small, out_small, 51);
    const double large = NsPerElement(row.fn, in_large, out_large, 11);
    std::printf("%16s %14.2f %14.2f\n", row.name, small, large);
  }
  std::printf("\n");
}

}  // namespace

int main() {
  // Numbers attributed to the wrong tier are worse than no numbers.
  if (!corvus_test::ReportAndCheckTarget()) return 2;
  if (std::strcmp(corvus::BenchTargetName(), corvus::active_target()) != 0) {
    std::fprintf(stderr,
                 "FAIL: this TU dispatched '%s' but the library dispatched "
                 "'%s' — the two target sets are not the same build.\n",
                 corvus::BenchTargetName(), corvus::active_target());
    return 2;
  }

  RunBand("zone", 0.5, 2.5,
          {{"floor", corvus::Floor},
           {"zone_full", corvus::ZoneFull},
           {"zone_horner", corvus::ZoneHorner},
           {"zone_horner_1c", corvus::ZoneHorner1c},
           {"zone_lead", corvus::ZoneLead}});
  RunBand("recurrence", 2.5, 8.0,
          {{"floor", corvus::Floor},
           {"rec_full", corvus::RecFull},
           {"rec_walk", corvus::RecWalk},
           {"rec_log", corvus::RecLog},
           {"zone_full", corvus::ZoneFull}});
  RunBand("Stirling", 8.0, 1000.0,
          {{"floor", corvus::Floor},
           {"stir_full", corvus::StirFull},
           {"stir_log", corvus::StirLog},
           {"stir_rem", corvus::StirRem}});
  RunBand("reflection", -30.0, -0.01,
          {{"floor", corvus::Floor},
           {"refl_full", corvus::ReflFull},
           {"refl_pos", corvus::ReflPos},
           {"refl_logs", corvus::ReflLogs},
           {"refl_sinc", corvus::ReflSinc}});
  return 0;
}

#endif  // HWY_ONCE
