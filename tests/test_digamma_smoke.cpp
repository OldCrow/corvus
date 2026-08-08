// Smoke test for corvus::digamma: the full specials table, a handful of
// closed-form values, and lane-mix determinism.
//
// The last one is why this file exists separately from the ULP gate. The
// [2, 8) down-walk stops per lane on a fire mask and lets the whole vector
// break out once no lane fires, so a lane's iteration count is decided by its
// NEIGHBOURS; the region cores additionally run on all lanes with the
// inactive ones scrubbed. If a frozen accumulator were merely having zeros
// added to it rather than being left alone, the extra dd renormalizations
// would show in the last bits -- and only ever in mixed vectors, which a
// single-point test never builds. So every probe here is evaluated alone, at
// several lane offsets, and interleaved with points from other regions, and
// all of it must be bit-identical.
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <vector>

#include "corvus/corvus.h"
#include "expect_target.h"

namespace {

const double kInf = std::numeric_limits<double>::infinity();
const double kNan = std::numeric_limits<double>::quiet_NaN();

int g_fail = 0;

int64_t OrderedBits(double x) {
  int64_t b;
  std::memcpy(&b, &x, sizeof(b));
  return b < 0 ? (INT64_MIN - b) : b;
}
uint64_t UlpDiff(double a, double b) {
  return static_cast<uint64_t>(std::llabs(OrderedBits(a) - OrderedBits(b)));
}
bool SameBits(double a, double b) {
  uint64_t ba, bb;
  std::memcpy(&ba, &a, sizeof(ba));
  std::memcpy(&bb, &b, sizeof(bb));
  return ba == bb;
}

double One(double x) {
  double out = 0.0;
  corvus::digamma(std::span<const double>(&x, 1), std::span<double>(&out, 1));
  return out;
}

// A NaN expectation means "any NaN". An infinite expectation is matched with
// its sign, which is the whole content of the signed-zero pole convention.
void CheckSpecial(const char* what, double x, double want) {
  const double got = One(x);
  bool ok;
  if (std::isnan(want)) {
    ok = std::isnan(got);
  } else {
    ok = got == want && std::signbit(got) == std::signbit(want);
  }
  if (!ok) {
    std::fprintf(stderr, "FAIL: %s: digamma(%.17g) = %.17g, want %.17g\n", what,
                 x, got, want);
    g_fail = 1;
  }
}

void Specials() {
  // psi(+-0) = -+inf: the signed-zero pole convention (scipy parity).
  CheckSpecial("psi(+0)", 0.0, -kInf);
  CheckSpecial("psi(-0)", -0.0, kInf);

  // Every negative integer is a pole -> NaN. The last three carry the rule
  // that every double <= -2^53 is an integer, so the whole far negative axis
  // is poles; -0x1p53 itself is the boundary case.
  const double kNegInts[] = {-1.0,      -2.0,       -3.0,   -17.0,
                             -1.0e6,    -1.0e15,    -0x1p53, -0x1p53 - 2.0,
                             -1.0e300,  -1.7976931348623157e308};
  for (double v : kNegInts) CheckSpecial("negative integer", v, kNan);

  CheckSpecial("psi(+inf)", kInf, kInf);
  CheckSpecial("psi(-inf)", -kInf, kNan);

  // Subnormal arguments: the -1/x limb overflows, so the pole convention
  // carries the sign through.
  CheckSpecial("psi(smallest +subnormal)", 5e-324, -kInf);
  CheckSpecial("psi(+subnormal)", 1e-320, -kInf);
  CheckSpecial("psi(smallest -subnormal)", -5e-324, kInf);
  CheckSpecial("psi(-subnormal)", -1e-320, kInf);

  // NaN propagates WITH its payload, matching every other corvus kernel.
  uint64_t bits = 0x7FF8000000ABCDEFULL;
  double payload;
  std::memcpy(&payload, &bits, sizeof(payload));
  const double got = One(payload);
  if (!SameBits(got, payload)) {
    std::fprintf(stderr, "FAIL: digamma did not propagate the NaN payload\n");
    g_fail = 1;
  }
}

// Closed-form and reference values, one per region. Not an accuracy gate --
// test_digamma_ulp owns that -- just a check that each branch is wired to the
// right formula.
void KnownValues() {
  struct Point {
    double x;
    double want;
    const char* region;
  };
  const Point kPoints[] = {
      {0.25, -0x1.0e8e9943cd7c3p+2, "(0,1)"},
      {0.5, -0x1.f6a897d3214fcp+0, "(0,1)"},
      {1.0, -0x1.2788cfc6fb619p-1, "zone (= -gamma)"},
      {1.5, 0x1.2aed059bd608ap-5, "zone"},
      {2.0, 0x1.b0ee6072093cep-2, "walk"},
      {3.0, 0x1.d8773039049e7p-1, "walk"},
      {8.0, 0x1.02008a3a23e5dp+1, "asymptotic (X0)"},
      {10.0, 0x1.20396dc85cc95p+1, "asymptotic"},
      {100.0, 0x1.26690d4274475p+2, "asymptotic"},
      {1.0e10, 0x1.7069e2aa27361p+4, "asymptotic"},
      {-0.5, 0x1.2aed059bd608ap-5, "reflection"},
      {-0.25, 0x1.750282bca7d92p+1, "reflection"},
      {-1.5, 0x1.680425af12b5ep-1, "reflection"},
      {-2.5, 0x1.1a68793defc15p+0, "reflection"},
  };
  for (const Point& p : kPoints) {
    const double got = One(p.x);
    const uint64_t u = UlpDiff(got, p.want);
    if (u > 2) {
      std::fprintf(stderr,
                   "FAIL: digamma(%.17g) [%s] = %.17g, want %.17g (%llu ULP)\n",
                   p.x, p.region, got, p.want,
                   static_cast<unsigned long long>(u));
      g_fail = 1;
    }
  }
}

// A spread of in-domain points covering every region, both sides of every
// boundary the driver tests, and the walk at every depth 1..6.
std::vector<double> Probes() {
  return {
      // (0, 1), including the up-step's own extremes.
      1e-300, 1e-8, 0.001, 0.1, 0.25, 0.4616321, 0.5, 0.75, 0.9,
      0.99999999999,
      // zone [1, 2), including bit-neighbourhoods of the root.
      1.0, 1.0000000001, 1.2, 1.4616321449683623, 1.4616321449683625, 1.5,
      1.75, 1.9999999999,
      // walk [2, 8) at each depth 1..6, plus both sides of every step wall.
      2.0, 2.5, 2.9999999999, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0,
      7.5, 7.9999999999,
      // asymptotic [8, inf), both sides of X0 and out past the log-only cut.
      7.9999999999999991, 8.0, 8.0000000000000018, 9.0, 12.0, 100.0, 1e6,
      1e17, 1e30, 0x1p85, 0x1p86, 1e150, 1e300,
      // negative axis: between poles, near zeros, near poles, and large.
      -0.5, -0.75, -1.5, -1.9, -2.1, -2.5, -3.5, -4.5, -10.5, -20.5,
      -1.4616321449683625, -3.5442455, -6.6773204, -1e6 - 0.5, -1e15 - 0.5,
      -1.0000000000000002, -0.9999999999999998, -2.0000000000000004,
  };
}

void LaneMix() {
  const auto probes = Probes();
  const size_t n = probes.size();

  std::vector<double> alone(n);
  for (size_t i = 0; i < n; ++i) alone[i] = One(probes[i]);

  // Offsets 0..8 cover every lane position on every tier up to AVX-512, and
  // the leading filler is drawn from the far end of the probe list so each
  // probe sits next to a different region's point every time.
  for (size_t off = 0; off < 9; ++off) {
    std::vector<double> in;
    for (size_t k = 0; k < off; ++k) in.push_back(probes[(n - 1 - k) % n]);
    for (double v : probes) in.push_back(v);
    std::vector<double> got(in.size());
    corvus::digamma(in, got);
    for (size_t i = 0; i < n; ++i) {
      if (!SameBits(got[off + i], alone[i])) {
        std::fprintf(stderr,
                     "FAIL: lane-mix drift at x=%.17g (offset %zu): batched "
                     "%.17g vs alone %.17g\n",
                     probes[i], off, got[off + i], alone[i]);
        g_fail = 1;
      }
    }
  }
}

}  // namespace

int main() {
  if (!corvus_test::ReportAndCheckTarget()) return 2;
  Specials();
  KnownValues();
  LaneMix();
  if (g_fail == 0) std::printf("PASS: digamma smoke\n");
  return g_fail;
}
