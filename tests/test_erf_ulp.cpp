// Measures max ULP deviation of corvus::erf against the mpmath-generated
// correctly-rounded reference (tests/data/erf_reference.txt). Gate: <= 1 ULP.
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "corvus/corvus.h"

namespace {

// Monotonic sign-magnitude mapping: adjacent doubles differ by 1 everywhere,
// including across +/-0.
int64_t OrderedBits(double x) {
  int64_t b;
  std::memcpy(&b, &x, sizeof(b));
  return b < 0 ? (INT64_MIN - b) : b;
}

uint64_t UlpDiff(double a, double b) {
  return static_cast<uint64_t>(std::llabs(OrderedBits(a) - OrderedBits(b)));
}

}  // namespace

int main(int argc, char** argv) {
  const char* path = argc > 1 ? argv[1] : "tests/data/erf_reference.txt";
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
  corvus::erf(in, got);

  uint64_t max_ulp = 0;
  size_t worst = 0;
  size_t over_half = 0;  // count of results differing from correctly-rounded
  for (size_t i = 0; i < in.size(); ++i) {
    const uint64_t u = UlpDiff(got[i], want[i]);
    if (u > 0) {
      ++over_half;
    }
    if (u > max_ulp) {
      max_ulp = u;
      worst = i;
    }
  }

  std::printf("corvus active SIMD target: %s\n", corvus::active_target());
  std::printf("points: %zu   max ULP: %llu   not-correctly-rounded: %zu (%.3f%%)\n",
              in.size(), static_cast<unsigned long long>(max_ulp), over_half,
              100.0 * static_cast<double>(over_half) / static_cast<double>(in.size()));
  if (max_ulp > 0) {
    std::printf("worst: x=%s  got=%s  want=%s\n", std::to_string(in[worst]).c_str(),
                std::to_string(got[worst]).c_str(), std::to_string(want[worst]).c_str());
  }

  if (max_ulp > 1) {
    std::fprintf(stderr, "FAIL: max ULP %llu exceeds bound 1\n",
                 static_cast<unsigned long long>(max_ulp));
    return 1;
  }
  std::printf("PASS: max ULP <= 1\n");
  return 0;
}
