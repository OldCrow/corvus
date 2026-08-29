// Measures max ULP deviation of corvus::erf against the mpmath-generated
// correctly-rounded reference (tests/data/erf_reference.txt). Gate: <= 1 ULP.
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "corvus/corvus.h"
#include "expect_target.h"
#include "ulp_utils.h"

namespace {

using corvus_test::SameBits;
using corvus_test::UlpDiff;

}  // namespace

int main(int argc, char** argv) {
  // Before doing any work: a wrong tier makes every number below meaningless.
  if (!corvus_test::ReportAndCheckTarget()) return 2;

  const char* path = argc > 1 ? argv[1] : "tests/data/erf_reference.txt";
  const auto rows = corvus_test::LoadRef(path, 2, 10000);

  std::vector<double> in, want;
  for (const auto& row : rows) {
    in.push_back(corvus_test::ParseDouble(row.tok[0], path, row.line));
    want.push_back(corvus_test::ParseDouble(row.tok[1], path, row.line));
  }

  std::vector<double> got(in.size());
  {
    // #14 N4: split the whole-set call in two so the masked LoadN/StoreN
    // tail path runs even when the reference set's length happens to be a
    // multiple of every SIMD tier's lane count. The final call's length (3)
    // is below every lane count, and N-3 is odd whenever N is even, so
    // neither call can land on a lane boundary either.
    const std::span<const double> in_s(in);
    const std::span<double> got_s(got);
    const size_t n = in_s.size();
    corvus::erf(in_s.first(n - 3), got_s.first(n - 3));
    corvus::erf(in_s.last(3), got_s.last(3));
  }

  uint64_t max_ulp = 0;
  size_t worst = 0;
  size_t over_half = 0;  // count of results differing from correctly-rounded
  for (size_t i = 0; i < in.size(); ++i) {
    // #14 N6: a reference value of exactly 0.0 must be checked bit-exact --
    // UlpDiff maps +0 and -0 to the same point (distance 0; see
    // ulp_utils.h's policy comment), so it cannot see a sign regression
    // there. This reference set has no exact-zero row today; the check
    // exists so the next regen can't silently drop one through the gate.
    const bool zero_ref = want[i] == 0.0;
    const uint64_t u = zero_ref ? (SameBits(got[i], want[i]) ? 0 : UINT64_MAX)
                                 : UlpDiff(got[i], want[i]);
    if (zero_ref && u != 0) {
      std::fprintf(stderr,
                   "FAIL: signed-zero mismatch at x=%.17g: got=%.17g want=%.17g\n",
                   in[i], got[i], want[i]);
    }
    if (u > 0) {
      ++over_half;
    }
    if (u > max_ulp) {
      max_ulp = u;
      worst = i;
    }
  }

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
