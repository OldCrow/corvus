// Measures max ULP deviation of corvus::trigamma against the mpmath-generated
// correctly-rounded reference, split by the kernel's regions.
//
// ONE METRIC, EVERYWHERE. This test is structurally simpler than
// test_digamma_ulp, and deliberately so: psi_1(x) = sum over n >= 0 of
// 1/(x + n)^2 is a sum of squares, hence strictly positive wherever finite, so
// it has NO zeros -- not on the positive axis, and not on the negative axis
// either, where its per-interval minima rise monotonically from 8.933 to
// pi^2. digamma and lgamma both need an absolute band because they pass
// through zero at points with no closed form, and a fixed absolute error near
// such a point is an unbounded relative one. Nothing here does. A single
// relative (ULP) metric therefore covers the whole real line, and adding an
// absolute band would report a number about a region that does not exist.
//
// The buckets follow the kernel's own region split, so a regression names the
// branch that moved. The asymptotic bucket spans [X0, inf) as one gate: the
// 1/x-only path above kTrigammaAsymCut is the same branch with its corrections
// selected away, not a separate approximation. Its numbers are also printed
// separately, report-only, because the two sub-paths fail for different
// reasons and the split costs nothing.
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#include "corvus/corvus.h"
#include "expect_target.h"
#include "src/trigamma_data.h"

namespace {

// Gate PINNED to measured, no margin. Identical cells on every validated
// leg — AVX3_ZEN4 native, AVX2/SSE4/SSSE3/SSE2 capped (Ryzen), Linux CI
// sweep, NEON (CI, identical to native including not-CR counts), Windows
// MSVC: (0,1) correctly rounded, all other buckets 1 ULP max. The
// down-walk's predicted 12.36x amplification cost did not consume a bit on
// this reference set — it shows as the walk bucket's elevated not-CR rate
// (2.97% FMA / 4.44% non-FMA) instead.
constexpr uint64_t kMaxUlp = 1;

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

void Accumulate(Region& r, double x, uint64_t u) {
  ++r.n;
  if (u > 0) ++r.miss;
  if (u > r.max_ulp) {
    r.max_ulp = u;
    r.worst_x = x;
  }
}

void Report(const Region& r, bool gated) {
  char gate[24];
  if (gated) {
    std::snprintf(gate, sizeof(gate), "gate %llu",
                  static_cast<unsigned long long>(kMaxUlp));
  } else {
    std::snprintf(gate, sizeof(gate), "report only");
  }
  std::printf(
      "%-20s n=%6zu  max ULP=%3llu (%s)  not-CR: %zu (%.2f%%)  "
      "worst x=%.17g\n",
      r.name, r.n, static_cast<unsigned long long>(r.max_ulp), gate, r.miss,
      r.n ? 100.0 * static_cast<double>(r.miss) / static_cast<double>(r.n)
          : 0.0,
      r.worst_x);
}

}  // namespace

int main(int argc, char** argv) {
  // Before doing any work: a wrong tier makes every number below meaningless.
  if (!corvus_test::ReportAndCheckTarget()) return 2;

  const char* path = argc > 1 ? argv[1] : "tests/data/trigamma_reference.txt";
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
  corvus::trigamma(in, got);

  // Routing mirrors src/trigamma_data.h's own constants, so the test and the
  // kernel cannot drift apart.
  Region regions[] = {
      {"pos (0,1)"}, {"pos zone [1,2)"}, {"pos walk [2,8)"},
      {"pos asym [8,inf)"}, {"negative"},
  };
  // Report-only split of the asymptotic bucket at the 1/x-only cut.
  Region asym_sub[] = {{"  ... [8,2^89)"}, {"  ... [2^89,inf)"}};

  for (size_t i = 0; i < in.size(); ++i) {
    const double x = in[i];
    const uint64_t u = UlpDiff(got[i], want[i]);
    if (x < 0.0) {
      Accumulate(regions[4], x, u);
    } else if (x < corvus::detail::kTrigammaZoneLo) {
      Accumulate(regions[0], x, u);
    } else if (x < corvus::detail::kTrigammaZoneHi) {
      Accumulate(regions[1], x, u);
    } else if (x < corvus::detail::kTrigammaX0) {
      Accumulate(regions[2], x, u);
    } else {
      Accumulate(regions[3], x, u);
      Accumulate(asym_sub[x < corvus::detail::kTrigammaAsymCut ? 0 : 1], x, u);
    }
  }

  int rc = 0;
  for (const Region& r : regions) {
    Report(r, true);
    if (r.max_ulp > kMaxUlp) {
      std::fprintf(stderr, "FAIL: %s exceeds gate\n", r.name);
      rc = 1;
    }
    if (&r == &regions[3]) {
      for (const Region& s : asym_sub) Report(s, false);
    }
  }

  if (rc == 0) std::printf("PASS: all regions within gates\n");
  return rc;
}
