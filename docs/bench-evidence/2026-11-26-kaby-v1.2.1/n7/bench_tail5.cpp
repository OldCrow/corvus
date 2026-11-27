// Investigation bench 5 (K, v1.2.1): does the cost of a call depend on how
// many DISTINCT calls the timed loop repeats? Each repetition replays the
// same inputs, so a branch predictor can learn a kernel's data-dependent
// branches when the sequence is short. One layout (input and result arrays
// 64-byte aligned, 128 bytes apart), both store widths, n = 1, 4, 7; the
// number of distinct calls per repetition is swept and printed in the `off`
// column as log2. Total work per cell is held near 41 * 16384 calls.
// Original header (bench 2) follows.
// Investigation bench 2 (K, v1.2.1): a short span written, called and read
// back, with the caller's store width and the input array's alignment under
// control. bench_span1's and bench_tail's "write n values" loop compiles to a
// memcpy call on AppleClang 15, so what the driver's loads meet there is
// whatever stores memmove chooses for n*8 bytes, not n 8-byte stores.
// Links libcorvus.a only (build line as bench_tail.cpp).
//
//   store   s8  = n separate 8-byte stores (an empty asm keeps each value in
//                 a register, so the loop cannot become memcpy or vectorise)
//           mc  = one memcpy of n*8 bytes
//   off     the input array starts 8*off bytes past a 64-byte boundary
//           (off 0..7); the result array is 64-byte aligned.
// Output: fn store off n median min   (ns per CALL; 41 reps of 512 calls)
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

#define NOINL __attribute__((noinline))

NOINL double CallS8(Fn f, const double* src, size_t n, size_t off) {
  alignas(64) double buf[32];
  double* a = buf + off;
  double* r = buf + 16;
  for (size_t k = 0; k < n; ++k) {
    double t = src[k];
    __asm__ volatile("" : "+x"(t));
    a[k] = t;
  }
  f(std::span<const double>(a, n), std::span<double>(r, n));
  double sum = 0.0;
  for (size_t k = 0; k < n; ++k) {
    double t = r[k];
    __asm__ volatile("" : "+x"(t));
    sum += t;
  }
  return sum;
}
NOINL double CallMc(Fn f, const double* src, size_t n, size_t off) {
  alignas(64) double buf[32];
  double* a = buf + off;
  double* r = buf + 16;
  std::memcpy(a, src, n * sizeof(double));
  __asm__ volatile("" : : "r"(a) : "memory");
  f(std::span<const double>(a, n), std::span<double>(r, n));
  double sum = 0.0;
  for (size_t k = 0; k < n; ++k) {
    double t = r[k];
    __asm__ volatile("" : "+x"(t));
    sum += t;
  }
  return sum;
}

using Shape = double (*)(Fn, const double*, size_t, size_t);

void Run(const char* label, Fn f, double lo, double hi, const char* only) {
  if (only && std::strcmp(only, label) != 0) return;
  constexpr size_t kMaxCalls = 65536;
  std::mt19937_64 rng(4096);
  std::uniform_real_distribution<double> u(lo, hi);
  std::vector<double> in(kMaxCalls * 8 + 8);
  for (auto& v : in) v = u(rng);
  const struct {
    const char* name;
    Shape s;
  } shapes[] = {{"s8", CallS8}, {"mc", CallMc}};
  volatile size_t voff = 0;
  for (const auto& sh : shapes) {
    for (size_t lg : {6, 9, 12, 14, 16}) {
      const size_t off = lg;
      const size_t kCalls = size_t{1} << lg;
      const int kReps = lg >= 14 ? 21 : 41;
      const size_t kInner = lg >= 14 ? 1 : (size_t{16384} >> lg);
      voff = 0;
      for (size_t n : {1, 4, 7}) {
        std::vector<double> t(kReps);
        auto body = [&] {
          const size_t o = voff;
          double acc = 0.0;
          for (size_t j = 0; j < kInner; ++j)
            for (size_t i = 0; i < kCalls; ++i)
              acc += sh.s(f, &in[i * 8], n, o);
          g_sink = acc;
        };
        body();
        for (auto& ti : t) {
          const auto t0 = Clock::now();
          body();
          const auto t1 = Clock::now();
          ti = std::chrono::duration<double, std::nano>(t1 - t0).count() /
               static_cast<double>(kCalls * kInner);
        }
        std::sort(t.begin(), t.end());
        std::printf("%-8s %-3s %2zu %2zu %7.1f %7.1f\n", label, sh.name, off, n,
                    t[kReps / 2], t[0]);
      }
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
