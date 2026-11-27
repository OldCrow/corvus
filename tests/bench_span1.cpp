// Benchmark: the cost of ONE value through the span API, in the two shapes a
// caller can give it. Manual, quiet machine, Release only; not registered
// with ctest.
//
//   array:  the value and its result live in arrays (bench_erfinv's scalar
//           walk): out[i] = f(in[i]) through spans of length 1.
//   local:  the value arrives by value and the result is returned, as a
//           scalar wrapper in a consumer does it:
//             double f1(double x) { double r; f({&x, 1}, {&r, 1}); return r; }
//           The span then points at a stack slot written just before the
//           call, and the result is read back just after it.
//
// The drivers' length-1 fast path (v1.2.1, #110) exists for the second shape.
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <random>
#include <span>
#include <vector>

#include "corvus/corvus.h"
#include "expect_target.h"

namespace {

using Clock = std::chrono::steady_clock;
using Fn = void (*)(std::span<const double>, std::span<double>) noexcept;
volatile double g_sink;

#if defined(_MSC_VER)
#define BENCH_NOINLINE __declspec(noinline)
#else
#define BENCH_NOINLINE __attribute__((noinline))
#endif

BENCH_NOINLINE double Local(Fn f, double x) {
  double r;
  f(std::span<const double>(&x, 1), std::span<double>(&r, 1));
  return r;
}

template <class Body>
double Median(Body&& body, size_t n) {
  constexpr int kReps = 41;
  std::vector<double> t(kReps);
  body();
  for (auto& ti : t) {
    const auto t0 = Clock::now();
    body();
    const auto t1 = Clock::now();
    ti = std::chrono::duration<double, std::nano>(t1 - t0).count() /
         static_cast<double>(n);
  }
  std::nth_element(t.begin(), t.begin() + kReps / 2, t.end());
  return t[kReps / 2];
}

void Run(const char* label, Fn f, double lo, double hi) {
  constexpr size_t kN = 4096;
  std::mt19937_64 rng(4096);
  std::uniform_real_distribution<double> u(lo, hi);
  std::vector<double> in(kN), out(kN);
  for (auto& v : in) v = u(rng);

  const double batch = Median(
      [&] {
        f(in, out);
        g_sink = out[kN / 2];
      },
      kN);
  const double array = Median(
      [&] {
        for (size_t i = 0; i < kN; ++i)
          f(std::span<const double>(&in[i], 1), std::span<double>(&out[i], 1));
        g_sink = out[kN / 2];
      },
      kN);
  const double local = Median(
      [&] {
        double acc = 0.0;
        for (size_t i = 0; i < kN; ++i) acc += Local(f, in[i]);
        g_sink = acc;
      },
      kN);
  std::printf("%-28s %9.1f %9.1f %9.1f\n", label, batch, array, local);
}

}  // namespace

int main() {
  if (!corvus_test::ReportAndCheckTarget()) return 1;
  std::printf("%-28s %9s %9s %9s\n", "ns per value", "batch", "array", "local");
  Run("erfinv  [-0.5, 0.5]", corvus::erfinv, -0.5, 0.5);
  Run("erfinv  [-0.999, 0.999]", corvus::erfinv, -0.999, 0.999);
  Run("erfcinv [0.5, 1.5]", corvus::erfcinv, 0.5, 1.5);
  Run("erfcinv (0, 2) (probit)", corvus::erfcinv, 1e-300, 2.0);
  Run("erf     [-3, 3]", corvus::erf, -3.0, 3.0);
  Run("lgamma  [1, 21]", corvus::lgamma, 1.0, 21.0);
  Run("exp     [-20, 20]", corvus::exp, -20.0, 20.0);
  return 0;
}
