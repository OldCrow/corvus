// Measures max ULP deviation of corvus::lgamma against the mpmath-generated
// correctly-rounded reference, split by the kernel's regions.
//
// The negative axis is reported in two bands, and that split is the point
// rather than a convenience. lgamma has zeros there wherever |Gamma(x)| = 1 --
// infinitely many, none with a closed form, so none reproducible by an
// exact-argument form the way x = 1 and x = 2 are. Near them a fixed absolute
// error is an unbounded relative one, so the |result| < 1 band is gated on
// absolute error and the rest on ULP. Reporting one blended number would hide
// which of the two is actually moving.
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#include "corvus/corvus.h"
#include "expect_target.h"

namespace {

// Gates, set from measured values (see PLAN.md / docs/ACCURACY.md) with no
// margin: regressions should trip them. Measured identical on AVX3_ZEN4,
// AVX2, SSE4, SSSE3 and SSE2.
constexpr uint64_t kMaxUlpSmall = 1;   // 0 < x < 1/2   (shift by one log)
constexpr uint64_t kMaxUlpZone = 1;    // 1/2 <= x <= 5/2 (the zeros live here)
constexpr uint64_t kMaxUlpMid = 1;     // 5/2 < x < 8   (product recurrence)
// Stirling is correctly rounded on every reference point, so its gate is 0 --
// the tightest in the project, and deliberately so: this region has no fitted
// approximation carrying error, only log_dd and dd assembly, and a single ULP
// appearing here would mean one of those moved.
constexpr uint64_t kMaxUlpBig = 0;     // x >= 8        (Stirling)
constexpr uint64_t kMaxUlpNeg = 1;     // x < 0, |lgamma| >= 1
constexpr double kMaxAbsNegSmall = 0x1p-53;  // x < 0, |lgamma| < 1

int64_t OrderedBits(double x) {
  int64_t b;
  std::memcpy(&b, &x, sizeof(b));
  return b < 0 ? (INT64_MIN - b) : b;
}

uint64_t UlpDiff(double a, double b) {
  return static_cast<uint64_t>(std::llabs(OrderedBits(a) - OrderedBits(b)));
}

struct Region {
  const char* name;
  uint64_t bound;
  uint64_t max_ulp = 0;
  size_t n = 0;
  size_t miss = 0;  // not correctly rounded
  double worst_x = 0.0;
};

}  // namespace

int main(int argc, char** argv) {
  // Before doing any work: a wrong tier makes every number below meaningless.
  if (!corvus_test::ReportAndCheckTarget()) return 2;

  const char* path = argc > 1 ? argv[1] : "tests/data/lgamma_reference.txt";
  std::ifstream f(path);
  if (!f) {
    std::fprintf(stderr, "cannot open reference file: %s\n", path);
    return 2;
  }

  std::vector<double> in, want;
  std::string sx, sy;
  while (f >> sx >> sy) {
    in.push_back(std::strtod(sx.c_str(), nullptr));
    want.push_back(std::strtod(sy.c_str(), nullptr));
  }
  if (in.size() < 10000) {
    std::fprintf(stderr, "reference file suspiciously small: %zu lines\n",
                 in.size());
    return 2;
  }

  std::vector<double> got(in.size());
  corvus::lgamma(in, got);

  Region regions[] = {
      {"pos small", kMaxUlpSmall}, {"pos zone", kMaxUlpZone},
      {"pos mid", kMaxUlpMid},     {"pos big", kMaxUlpBig},
      {"neg |lg|>=1", kMaxUlpNeg},
  };
  double neg_small_abs = 0.0;
  double neg_small_worst_x = 0.0;
  size_t neg_small_n = 0, neg_small_miss = 0;

  for (size_t i = 0; i < in.size(); ++i) {
    const double x = in[i];
    if (x < 0.0 && std::abs(want[i]) < 1.0) {
      const double e = std::abs(got[i] - want[i]);
      ++neg_small_n;
      if (got[i] != want[i]) ++neg_small_miss;
      if (e > neg_small_abs) {
        neg_small_abs = e;
        neg_small_worst_x = x;
      }
      continue;
    }
    Region& r = x < 0.0        ? regions[4]
                : x < 0.5      ? regions[0]
                : x <= 2.5     ? regions[1]
                : x < 8.0      ? regions[2]
                                : regions[3];
    const uint64_t u = UlpDiff(got[i], want[i]);
    ++r.n;
    if (u > 0) ++r.miss;
    if (u > r.max_ulp) {
      r.max_ulp = u;
      r.worst_x = x;
    }
  }

  int rc = 0;
  for (const Region& r : regions) {
    std::printf(
        "%-13s n=%6zu  max ULP=%3llu (gate %llu)  not-CR: %zu (%.2f%%)  "
        "worst x=%.17g\n",
        r.name, r.n, static_cast<unsigned long long>(r.max_ulp),
        static_cast<unsigned long long>(r.bound), r.miss,
        r.n ? 100.0 * static_cast<double>(r.miss) / static_cast<double>(r.n)
            : 0.0,
        r.worst_x);
    if (r.max_ulp > r.bound) {
      std::fprintf(stderr, "FAIL: %s exceeds gate\n", r.name);
      rc = 1;
    }
  }
  std::printf(
      "%-13s n=%6zu  max abs=%.3e (gate %.1e)  not-CR: %zu (%.2f%%)  "
      "worst x=%.17g\n",
      "neg |lg|<1", neg_small_n, neg_small_abs, kMaxAbsNegSmall, neg_small_miss,
      neg_small_n ? 100.0 * static_cast<double>(neg_small_miss) /
                        static_cast<double>(neg_small_n)
                  : 0.0,
      neg_small_worst_x);
  if (neg_small_abs > kMaxAbsNegSmall) {
    std::fprintf(stderr, "FAIL: neg |lg|<1 exceeds absolute gate\n");
    rc = 1;
  }

  if (rc == 0) std::printf("PASS: all regions within gates\n");
  return rc;
}
