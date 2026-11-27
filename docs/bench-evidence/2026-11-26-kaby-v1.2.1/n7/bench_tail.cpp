// Investigation bench (K, v1.2.1): the cost of a short span by caller shape,
// to attribute the driver's tail load and tail store separately. Links
// libcorvus.a and nothing of its internals; built outside any tree:
//   c++ -O3 -DNDEBUG -std=c++20 -fno-fast-math -ffp-contract=off
//       -I<tree>/include bench_tail.cpp <lib>/libcorvus.a
//       -Wl,-rpath,/usr/local/lib /usr/local/lib/libhwy.dylib
//
// Shapes, n values per call:
//   wr+rd  the caller writes n values to a stack array, calls, reads the n
//          results back (bench_span1's "span in place").
//   wr     writes the inputs, calls, does not read the results.
//   rd     inputs come from a large array written long before; results go
//          to a stack array and are read back.
//   none   large arrays both sides, nothing written or read around the call.
// Output: fn shape n median min   (ns per CALL; 41 repetitions of 512 calls)
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <random>
#include <span>
#include <vector>

#include "corvus/corvus.h"

namespace {

using Clock = std::chrono::steady_clock;
using Fn = void (*)(std::span<const double>, std::span<double>) noexcept;
volatile double g_sink;
constexpr size_t kMax = 16;

#define NOINL __attribute__((noinline))

NOINL double WrRd(Fn f, const double* src, double*, size_t n) {
  double a[kMax], r[kMax];
  for (size_t k = 0; k < n; ++k) a[k] = src[k];
  f(std::span<const double>(a, n), std::span<double>(r, n));
  double sum = 0.0;
  for (size_t k = 0; k < n; ++k) sum += r[k];
  return sum;
}
NOINL double Wr(Fn f, const double* src, double* dst, size_t n) {
  double a[kMax];
  for (size_t k = 0; k < n; ++k) a[k] = src[k];
  f(std::span<const double>(a, n), std::span<double>(dst, n));
  return 0.0;
}
NOINL double Rd(Fn f, const double* src, double*, size_t n) {
  double r[kMax];
  f(std::span<const double>(src, n), std::span<double>(r, n));
  double sum = 0.0;
  for (size_t k = 0; k < n; ++k) sum += r[k];
  return sum;
}
NOINL double None(Fn f, const double* src, double* dst, size_t n) {
  f(std::span<const double>(src, n), std::span<double>(dst, n));
  return 0.0;
}

using Shape = double (*)(Fn, const double*, double*, size_t);

void Run(const char* label, Fn f, double lo, double hi, const char* only) {
  if (only && std::strcmp(only, label) != 0) return;
  constexpr size_t kCalls = 512;
  constexpr int kReps = 41;
  std::mt19937_64 rng(4096);
  std::uniform_real_distribution<double> u(lo, hi);
  std::vector<double> in(kCalls * kMax + kMax), out(kCalls * kMax + kMax);
  for (auto& v : in) v = u(rng);
  const struct {
    const char* name;
    Shape s;
  } shapes[] = {{"wr+rd", WrRd}, {"wr", Wr}, {"rd", Rd}, {"none", None}};
  for (const auto& sh : shapes) {
    for (size_t n : {1, 2, 3, 4, 5, 6, 7, 8, 11}) {
      std::vector<double> t(kReps);
      auto body = [&] {
        double acc = 0.0;
        for (size_t i = 0; i < kCalls; ++i)
          acc += sh.s(f, &in[i * kMax], &out[i * kMax], n);
        g_sink = acc + out[kMax];
      };
      body();
      for (auto& ti : t) {
        const auto t0 = Clock::now();
        body();
        const auto t1 = Clock::now();
        ti = std::chrono::duration<double, std::nano>(t1 - t0).count() /
             static_cast<double>(kCalls);
      }
      std::sort(t.begin(), t.end());
      std::printf("%-8s %-6s %2zu %7.1f %7.1f\n", label, sh.name, n,
                  t[kReps / 2], t[0]);
    }
  }
}

}  // namespace

int main(int argc, char** argv) {
  const char* only = argc > 1 ? argv[1] : nullptr;
  std::printf("# target %s\n", corvus::active_target());
  Run("erfinv", corvus::erfinv, -0.5, 0.5, only);
  Run("exp", corvus::exp, -20.0, 20.0, only);
  Run("erf", corvus::erf, -3.0, 3.0, only);
  Run("erfcinv", corvus::erfcinv, 1e-300, 2.0, only);
  return 0;
}
