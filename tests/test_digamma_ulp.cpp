// Measures max ULP deviation of corvus::digamma against the mpmath-generated
// correctly-rounded reference, split by the kernel's regions.
//
// The negative axis is reported in two bands, and that split is the point
// rather than a convenience. psi has a zero between every consecutive pair of
// negative poles -- infinitely many, none with a closed form, so none
// reproducible by an exact-argument form the way the positive root at
// x0 ~ 1.4616 is. Near them a fixed absolute error is an unbounded relative
// one, so the |psi| < 1 band is gated on ABSOLUTE error (reported in units of
// 2^-53) and the rest on ULP. That is lgamma's convention verbatim, and it is
// the dual metric PLAN.md's FIRST DESIGN CORRECTION settled on: relative
// where |psi| >= 1, 2^-53-class absolute inside the zero bands. Reporting one
// blended number would hide which of the two is actually moving.
//
// The positive-axis buckets follow the kernel's own region split, so a
// regression names the branch that moved.
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

// Gates PINNED to measured, no margin (G4, 2026-08-08). Identical cells on
// every validated leg: AVX3_ZEN4 native, AVX2/SSE4/SSSE3/SSE2 capped
// (Ryzen), Linux CI sweep, NEON (CI), Windows MSVC — 1 ULP max in all five
// relative buckets, 1.00 x 2^-53 absolute in the negative zero bands (NEON
// matches native down to the not-CR counts).
constexpr uint64_t kMaxUlpRel = 1;
constexpr double kMaxAbsUnits = 1.0;  // in units of 2^-53
constexpr double kAbsUnit = 0x1p-53;

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
  uint64_t max_ulp = 0;
  size_t n = 0;
  size_t miss = 0;  // not correctly rounded
  double worst_x = 0.0;
};

}  // namespace

int main(int argc, char** argv) {
  // Before doing any work: a wrong tier makes every number below meaningless.
  if (!corvus_test::ReportAndCheckTarget()) return 2;

  const char* path = argc > 1 ? argv[1] : "tests/data/digamma_reference.txt";
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
  corvus::digamma(in, got);

  Region regions[] = {
      {"pos (0,1)"},   {"pos zone"},     {"pos walk"},
      {"pos asym"},    {"neg |psi|>=1"},
  };
  double neg_small_abs = 0.0;
  double neg_small_worst_x = 0.0;
  size_t neg_small_n = 0, neg_small_miss = 0;

  for (size_t i = 0; i < in.size(); ++i) {
    const double x = in[i];
    // The lgamma convention: on the negative axis, the metric is chosen by
    // |psi| alone. |psi| < 1 happens only inside the zero bands (|psi| grows
    // past 1 within ~1/8 of the way to either neighbouring pole), so this one
    // test reproduces the per-zero band membership the design describes
    // without carrying a table of zeros the kernel does not have either.
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
    Region& r = x < 0.0   ? regions[4]
                : x < 1.0 ? regions[0]
                : x < 2.0 ? regions[1]
                : x < 8.0 ? regions[2]
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
        static_cast<unsigned long long>(kMaxUlpRel), r.miss,
        r.n ? 100.0 * static_cast<double>(r.miss) / static_cast<double>(r.n)
            : 0.0,
        r.worst_x);
    if (r.max_ulp > kMaxUlpRel) {
      std::fprintf(stderr, "FAIL: %s exceeds gate\n", r.name);
      rc = 1;
    }
  }
  std::printf(
      "%-13s n=%6zu  max abs=%.3e = %.2f x 2^-53 (gate %.1f)  not-CR: %zu "
      "(%.2f%%)  worst x=%.17g\n",
      "neg |psi|<1", neg_small_n, neg_small_abs, neg_small_abs / kAbsUnit,
      kMaxAbsUnits, neg_small_miss,
      neg_small_n ? 100.0 * static_cast<double>(neg_small_miss) /
                        static_cast<double>(neg_small_n)
                  : 0.0,
      neg_small_worst_x);
  if (neg_small_abs > kMaxAbsUnits * kAbsUnit) {
    std::fprintf(stderr, "FAIL: neg |psi|<1 exceeds absolute gate\n");
    rc = 1;
  }

  if (rc == 0) std::printf("PASS: all regions within gates\n");
  return rc;
}
