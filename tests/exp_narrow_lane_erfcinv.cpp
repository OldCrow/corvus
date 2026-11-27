// EXPERIMENT (#110), not a gate and not registered with ctest.
//
// Question: is a one-value erfcinv cheaper when the SAME kernel runs at one
// or two lanes instead of the target's full vector, and is the result
// bit-identical to the shipped span path?
//
// A span-of-1 call already costs one vector in one region (the driver pads
// the tail with the live element, #42). What is left to remove is the width:
// gathers, and on some CPUs the cost of the widest registers. So this TU
// instantiates ErfcinvVec with CappedTag<double, 1> and <double, 2> inside
// each target and compares them against corvus::erfcinv on a length-1 span.
//
// Usage: exp_narrow_lane_erfcinv tests/data/erfcinv_reference.txt
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <random>
#include <span>
#include <vector>

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "tests/exp_narrow_lane_erfcinv.cpp"
#include "hwy/foreach_target.h"  // IWYU pragma: keep

#include "src/erfinv-inl.h"
#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {

static double ErfcinvN1Impl(double z) {
  const op::CappedTag<double, 1> d;
  return op::GetLane(ErfcinvVec(d, op::Set(d, z)));
}
static double ErfcinvN2Impl(double z) {
  const op::CappedTag<double, 2> d;
  return op::GetLane(ErfcinvVec(d, op::Set(d, z)));
}
// The full vector with the value in every lane: what a scalar entry point in
// the style of DriveScalarBinary would do.
static double ErfcinvFullImpl(double z) {
  const op::ScalableTag<double> d;
  return op::GetLane(ErfcinvVec(d, op::Set(d, z)));
}
// The two halves of the driver's masked tail, one at a time, to see which
// one a span-of-1 call pays for.
static double ErfcinvMaskedLoadImpl(double z) {
  const op::ScalableTag<double> d;
  return op::GetLane(ErfcinvVec(d, op::LoadNOr(op::Set(d, z), d, &z, 1)));
}
static double ErfcinvMaskedStoreImpl(double z) {
  const op::ScalableTag<double> d;
  double r;
  op::StoreN(ErfcinvVec(d, op::Set(d, z)), d, &r, 1);
  return r;
}
// The driver's own shape with a run-time length, as corvus::erfcinv reaches
// it, and the same with a length-1 fast path that skips the masked tail.
static void ErfcinvSpanReplicaImpl(const double* pi, double* po, size_t n) {
  const op::ScalableTag<double> d;
  const size_t N = op::Lanes(d);
  size_t i = 0;
  for (; i + N <= n; i += N) {
    op::Store(ErfcinvVec(d, op::Load(d, pi + i)), d, po + i);
  }
  if (i < n) {
    const size_t m = n - i;
    op::StoreN(ErfcinvVec(d, op::LoadNOr(op::Set(d, pi[i]), d, pi + i, m)), d,
               po + i, m);
  }
}
static void ErfcinvSpanFastImpl(const double* pi, double* po, size_t n) {
  const op::ScalableTag<double> d;
  if (n == 1) {
    po[0] = op::GetLane(ErfcinvVec(d, op::Set(d, pi[0])));
    return;
  }
  ErfcinvSpanReplicaImpl(pi, po, n);
}
static size_t FullLanesImpl() {
  const op::ScalableTag<double> d;
  return op::Lanes(d);
}
static const char* TestTargetNameImpl() { return hwy::TargetName(HWY_TARGET); }

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#if HWY_ONCE

#include "corvus/corvus.h"
#include "expect_target.h"

namespace corvus {
HWY_EXPORT(ErfcinvN1Impl);
HWY_EXPORT(ErfcinvN2Impl);
HWY_EXPORT(ErfcinvFullImpl);
HWY_EXPORT(ErfcinvMaskedLoadImpl);
HWY_EXPORT(ErfcinvMaskedStoreImpl);
HWY_EXPORT(ErfcinvSpanReplicaImpl);
HWY_EXPORT(ErfcinvSpanFastImpl);
HWY_EXPORT(FullLanesImpl);
HWY_EXPORT(TestTargetNameImpl);

static double ErfcinvN1(double z) {
  return HWY_DYNAMIC_DISPATCH(ErfcinvN1Impl)(z);
}
static double ErfcinvN2(double z) {
  return HWY_DYNAMIC_DISPATCH(ErfcinvN2Impl)(z);
}
static double ErfcinvFull(double z) {
  return HWY_DYNAMIC_DISPATCH(ErfcinvFullImpl)(z);
}
static double ErfcinvMaskedLoad(double z) {
  return HWY_DYNAMIC_DISPATCH(ErfcinvMaskedLoadImpl)(z);
}
static double ErfcinvMaskedStore(double z) {
  return HWY_DYNAMIC_DISPATCH(ErfcinvMaskedStoreImpl)(z);
}
static double ErfcinvSpanReplica(double z) {
  double r;
  HWY_DYNAMIC_DISPATCH(ErfcinvSpanReplicaImpl)(&z, &r, 1);
  return r;
}
static double ErfcinvSpanFast(double z) {
  double r;
  HWY_DYNAMIC_DISPATCH(ErfcinvSpanFastImpl)(&z, &r, 1);
  return r;
}
static size_t FullLanes() { return HWY_DYNAMIC_DISPATCH(FullLanesImpl)(); }
static const char* TestTargetName() {
  return HWY_DYNAMIC_DISPATCH(TestTargetNameImpl)();
}
}  // namespace corvus

namespace {

using Clock = std::chrono::steady_clock;
volatile double g_sink;

double Span1(double z) {
  double r;
  corvus::erfcinv(std::span<const double>(&z, 1), std::span<double>(&r, 1));
  return r;
}

bool SameBits(double a, double b) {
  uint64_t ua, ub;
  std::memcpy(&ua, &a, 8);
  std::memcpy(&ub, &b, 8);
  return ua == ub;
}

struct Mismatch {
  size_t n1 = 0, n2 = 0, full = 0, batch = 0;
};

// Every point through five routes; the span-of-1 call is the reference.
void Compare(const std::vector<double>& z, Mismatch* mm, const char* label) {
  std::vector<double> batch(z.size());
  corvus::erfcinv(z, batch);
  Mismatch local;
  for (size_t i = 0; i < z.size(); ++i) {
    const double ref = Span1(z[i]);
    const double a = corvus::ErfcinvN1(z[i]);
    const double b = corvus::ErfcinvN2(z[i]);
    const double c = corvus::ErfcinvFull(z[i]);
    if (!SameBits(ref, a)) {
      if (local.n1++ < 3)
        std::printf("  n1 differs at z=%a: %a vs %a\n", z[i], a, ref);
    }
    if (!SameBits(ref, b)) {
      if (local.n2++ < 3)
        std::printf("  n2 differs at z=%a: %a vs %a\n", z[i], b, ref);
    }
    if (!SameBits(ref, c)) ++local.full;
    if (!SameBits(ref, batch[i])) ++local.batch;
  }
  std::printf("%-22s n=%8zu  mismatches: n1 %zu  n2 %zu  full %zu  batch %zu\n",
              label, z.size(), local.n1, local.n2, local.full, local.batch);
  mm->n1 += local.n1;
  mm->n2 += local.n2;
  mm->full += local.full;
  mm->batch += local.batch;
}

template <class F>
double NsPerCall(F&& one, const std::vector<double>& z) {
  constexpr int kReps = 41;
  std::vector<double> t(kReps);
  for (auto& ti : t) {
    double acc = 0.0;
    const auto t0 = Clock::now();
    for (double zi : z) acc += one(zi);
    const auto t1 = Clock::now();
    g_sink = acc;
    ti = std::chrono::duration<double, std::nano>(t1 - t0).count() /
         static_cast<double>(z.size());
  }
  std::nth_element(t.begin(), t.begin() + kReps / 2, t.end());
  return t[kReps / 2];
}

double NsPerElementBatch(const std::vector<double>& z) {
  constexpr int kReps = 41;
  std::vector<double> out(z.size());
  std::vector<double> t(kReps);
  for (auto& ti : t) {
    const auto t0 = Clock::now();
    corvus::erfcinv(z, out);
    const auto t1 = Clock::now();
    g_sink = out[z.size() / 2];
    ti = std::chrono::duration<double, std::nano>(t1 - t0).count() /
         static_cast<double>(z.size());
  }
  std::nth_element(t.begin(), t.begin() + kReps / 2, t.end());
  return t[kReps / 2];
}

void Time(const char* label, const std::vector<double>& z) {
  // Two passes; the second is reported (the first warms every route).
  double r[9] = {};
  for (int pass = 0; pass < 2; ++pass) {
    r[0] = NsPerElementBatch(z);
    r[1] = NsPerCall(Span1, z);
    r[2] = NsPerCall(corvus::ErfcinvFull, z);
    r[3] = NsPerCall(corvus::ErfcinvN2, z);
    r[4] = NsPerCall(corvus::ErfcinvN1, z);
    r[5] = NsPerCall(corvus::ErfcinvMaskedLoad, z);
    r[6] = NsPerCall(corvus::ErfcinvMaskedStore, z);
    r[7] = NsPerCall(corvus::ErfcinvSpanReplica, z);
    r[8] = NsPerCall(corvus::ErfcinvSpanFast, z);
  }
  std::printf("%-26s %9.1f %9.1f %9.1f %9.1f %9.1f %9.1f %9.1f %9.1f %9.1f\n", label,
              r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8]);
}

}  // namespace

int main(int argc, char** argv) {
  if (!corvus_test::ReportAndCheckTarget()) return 1;
  if (std::strcmp(corvus::TestTargetName(), corvus::active_target()) != 0) {
    std::printf("FAIL: this TU dispatches %s, the library %s\n",
                corvus::TestTargetName(), corvus::active_target());
    return 1;
  }
  std::printf("target %s, full vector %zu lanes\n", corvus::active_target(),
              corvus::FullLanes());

  std::mt19937_64 rng(20261126);
  Mismatch mm;

  // 1. The checked-in reference inputs.
  if (argc > 1) {
    std::vector<double> z;
    if (FILE* f = std::fopen(argv[1], "r")) {
      char line[256];
      while (std::fgets(line, sizeof line, f)) {
        char* end = nullptr;
        const double v = std::strtod(line, &end);
        if (end != line) z.push_back(v);
      }
      std::fclose(f);
    }
    Compare(z, &mm, "reference set");
  }
  // 2. Specials and boundaries.
  {
    const double inf = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> z = {0.0,  -0.0, 2.0,     1.0,     0.5,      1.5,
                             inf,  -inf, nan,     -1.0,    2.5,      0x1p-1074,
                             0x1p-1022, 0x1p-1023, -0x1p-1074, 0x1.fffffffffffffp+0,
                             0x1.fffffffffffffp-2, 0x1.0000000000001p-1,
                             0x1.7ffffffffffffp+0, 0x1.8000000000001p+0,
                             0x1.fffffffffffffp-1, 0x1.0000000000001p+0};
    for (double v : std::vector<double>(z))
      for (int k = 1; k <= 4; ++k) {
        z.push_back(std::nextafter(v, 3.0 * k));
        z.push_back(std::nextafter(v, -3.0 * k));
      }
    Compare(z, &mm, "specials + neighbours");
  }
  // 3. Dense random: uniform on (0, 2), log-uniform down to the subnormals,
  //    and the mirror next to 2.
  {
    std::vector<double> z(1000000);
    std::uniform_real_distribution<double> u(0.0, 2.0);
    for (auto& v : z) v = u(rng);
    Compare(z, &mm, "uniform (0,2)");
    std::uniform_real_distribution<double> e(-1074.0, -1.0);
    for (auto& v : z) v = std::ldexp(1.0 + 0.5 * u(rng), -1) * std::exp2(e(rng));
    Compare(z, &mm, "log-uniform to 2^-1074");
    std::uniform_real_distribution<double> e2(-52.0, -1.0);
    for (auto& v : z) v = 2.0 - std::exp2(e2(rng));
    Compare(z, &mm, "2 - 2^[-52,-1]");
    std::uniform_real_distribution<double> e3(-1074.0, -1.0);
    for (auto& v : z) v = 1.0 + (u(rng) - 1.0) * std::exp2(0.25 * e3(rng));
    Compare(z, &mm, "near 1");
  }
  const bool identical = mm.n1 == 0 && mm.n2 == 0 && mm.full == 0;
  std::printf("%s: narrow-lane results %s the span-of-1 path\n",
              identical ? "IDENTICAL" : "DIFFERENT",
              identical ? "match" : "do not match");

  // Timing, ns per call (batch: ns per element). Indicative unless the
  // machine is quiet.
  std::printf("\n%-26s %9s %9s %9s %9s %9s %9s %9s %9s %9s\n",
              "population (4096 values)", "batch/el", "span-of-1", "full",
              "2 lanes", "1 lane", "m.load", "m.store", "replica", "n1-fast");
  {
    std::vector<double> z(4096);
    std::uniform_real_distribution<double> c(0.5, 1.5);
    for (auto& v : z) v = c(rng);
    Time("central [1/2, 3/2]", z);
    std::uniform_real_distribution<double> lg(-26.0, -1.0);
    for (auto& v : z) v = std::exp2(lg(rng));
    Time("tail mid 2^[-26,-1]", z);
    std::uniform_real_distribution<double> lf(-1000.0, -60.0);
    for (auto& v : z) v = std::exp2(lf(rng));
    Time("tail far 2^[-1000,-60]", z);
    std::uniform_real_distribution<double> u(0.0, 2.0);
    for (auto& v : z) v = u(rng);
    Time("uniform (0,2) (probit)", z);
  }
  return identical ? 0 : 1;
}

#endif  // HWY_ONCE
