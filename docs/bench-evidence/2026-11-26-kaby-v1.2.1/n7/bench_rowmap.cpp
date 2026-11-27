// Investigation bench (K, v1.2.1): the cost of ONE erf value through a span,
// input by input. Each grid point x is called kCalls times in a row with the
// same x, so every call reads the same table row; the output maps cost
// against x (hence against the row's address). Links libcorvus.a only.
// Output: x median min   (ns per call; 15 repetitions of 2048 calls)
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <span>

#include "corvus/corvus.h"

namespace {
using Clock = std::chrono::steady_clock;
using Fn = void (*)(std::span<const double>, std::span<double>) noexcept;
volatile double g_sink;
__attribute__((noinline)) double Local(Fn f, double x) {
  double r;
  f(std::span<const double>(&x, 1), std::span<double>(&r, 1));
  return r;
}
}  // namespace

int main() {
  std::printf("# target %s\n", corvus::active_target());
  constexpr int kCalls = 2048, kReps = 15;
  for (int g = 0; g <= 1200; ++g) {
    volatile double vx = 0.0025 * g;
    double t[kReps];
    for (double& ti : t) {
      const double x = vx;
      const auto t0 = Clock::now();
      double acc = 0.0;
      for (int i = 0; i < kCalls; ++i) acc += Local(corvus::erf, x);
      const auto t1 = Clock::now();
      g_sink = acc;
      ti = std::chrono::duration<double, std::nano>(t1 - t0).count() / kCalls;
    }
    std::sort(t, t + kReps);
    std::printf("%.4f %7.2f %7.2f\n", 0.0025 * g, t[kReps / 2], t[0]);
  }
  return 0;
}
