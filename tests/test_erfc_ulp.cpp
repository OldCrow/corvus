// Measures max ULP deviation of corvus::Erfc against the mpmath-generated
// correctly-rounded reference, with a per-region breakdown (the tail is
// gated by the backend Exp accuracy, the core by the table method).
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#include "corvus/corvus.h"

namespace {

// Gates, set from measured values (see PLAN.md erfc section) with no margin:
// regressions should trip them. The tail-normal bound is the backend Exp's
// contribution (see the tail path in src/erfc.cpp) -- a corvus-owned
// compensated exp would tighten it toward the core's 1.
constexpr uint64_t kMaxUlpCore = 1;          // |x| <= 6
constexpr uint64_t kMaxUlpTailNormal = 5;    // |x| > 6, normal results
constexpr uint64_t kMaxUlpTailSubnormal = 1; // |x| > 6, subnormal results

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
  const char* path = argc > 1 ? argv[1] : "tests/data/erfc_reference.txt";
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
    std::fprintf(stderr, "reference file suspiciously small: %zu lines\n", in.size());
    return 2;
  }

  std::vector<double> got(in.size());
  corvus::Erfc(in, got);

  Region regions[] = {
      {"core |x|<=6", kMaxUlpCore},
      {"tail normal", kMaxUlpTailNormal},
      {"tail subnormal", kMaxUlpTailSubnormal},
  };
  constexpr double kMinNormal = 2.2250738585072014e-308;

  for (size_t i = 0; i < in.size(); ++i) {
    Region& r = std::abs(in[i]) <= 6.0 ? regions[0]
                : (want[i] >= kMinNormal ? regions[1] : regions[2]);
    const uint64_t u = UlpDiff(got[i], want[i]);
    ++r.n;
    if (u > 0) {
      ++r.miss;
    }
    if (u > r.max_ulp) {
      r.max_ulp = u;
      r.worst_x = in[i];
    }
  }

  std::printf("corvus active SIMD target: %s\n", corvus::ActiveTarget());
  int rc = 0;
  for (const Region& r : regions) {
    std::printf("%-15s n=%6zu  max ULP=%3llu (gate %llu)  not-CR: %zu (%.2f%%)  worst x=%.17g\n",
                r.name, r.n, static_cast<unsigned long long>(r.max_ulp),
                static_cast<unsigned long long>(r.bound), r.miss,
                r.n ? 100.0 * static_cast<double>(r.miss) / static_cast<double>(r.n) : 0.0,
                r.worst_x);
    if (r.max_ulp > r.bound) {
      std::fprintf(stderr, "FAIL: %s exceeds gate\n", r.name);
      rc = 1;
    }
  }
  if (rc == 0) {
    std::printf("PASS: all regions within gates\n");
  }
  return rc;
}
