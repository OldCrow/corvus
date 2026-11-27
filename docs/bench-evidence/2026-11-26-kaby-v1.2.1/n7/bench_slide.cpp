// Investigation bench (K, v1.2.1): does the cost of a short-span call depend
// on WHERE THE STACK IS? The whole measurement runs below an alloca of K
// bytes, so every stack address the caller and the library use moves down by
// K while no instruction changes. A cost that is periodic in K with period
// 4096 is 4K aliasing: a load whose address shares its low 12 bits with an
// earlier store still in flight is held as if it depended on that store.
// Links libcorvus.a only (build line as bench_tail.cpp).
//
// usage: bench_slide [step [limit [fn [span1]]]]   (bytes; default 64, 4096;
//        fn limits the run to one function; a fifth argument selects
//        bench_span1's caller in place of the 8-byte-store caller)
// Output: fn n K median min              (ns per CALL; 21 reps of 512 calls)
// A line "# addr ..." gives the alloca block's address at K = 0.
#include <alloca.h>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <span>
#include <string_view>
#include <vector>

#include "corvus/corvus.h"

namespace {

using Clock = std::chrono::steady_clock;
using Fn = void (*)(std::span<const double>, std::span<double>) noexcept;
volatile double g_sink;
const void* volatile g_addr;

#define NOINL __attribute__((noinline))

// The caller: n separate 8-byte stores, the call, n 8-byte loads.
NOINL double Call(Fn f, const double* src, size_t n) {
  double a[16], r[16];
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

// bench_span1's LocalN, verbatim: 8-element arrays, plain loops (the write
// loop compiles to a memcpy call on AppleClang 15).
NOINL double CallSpan1(Fn f, const double* src, size_t n) {
  double a[8], r[8];
  for (size_t k = 0; k < n; ++k) a[k] = src[k];
  f(std::span<const double>(a, n), std::span<double>(r, n));
  double sum = 0.0;
  for (size_t k = 0; k < n; ++k) sum += r[k];
  return sum;
}
bool g_span1 = false;

struct Band {
  const char* name;
  Fn f;
  double lo, hi;
};

NOINL void Measure(const Band& b, const std::vector<double>& in, size_t n,
                   size_t K) {
  constexpr size_t kCalls = 512;
  constexpr int kReps = 21;
  double t[kReps];
  auto body = [&] {
    double acc = 0.0;
    if (g_span1) {
      for (size_t i = 0; i < kCalls; ++i) acc += CallSpan1(b.f, &in[i * 8], n);
    } else {
      for (size_t i = 0; i < kCalls; ++i) acc += Call(b.f, &in[i * 8], n);
    }
    g_sink = acc;
  };
  body();
  for (double& ti : t) {
    const auto t0 = Clock::now();
    body();
    const auto t1 = Clock::now();
    ti = std::chrono::duration<double, std::nano>(t1 - t0).count() /
         static_cast<double>(kCalls);
  }
  std::sort(t, t + kReps);
  std::printf("%-8s %zu %5zu %7.1f %7.1f\n", b.name, n, K, t[kReps / 2], t[0]);
}

NOINL void AtSlide(const Band& b, const std::vector<double>& in, size_t K) {
  char* pad = static_cast<char*>(alloca(K + 16));
  __asm__ volatile("" : : "r"(pad) : "memory");
  if (K == 0) g_addr = pad;
  for (size_t n : {1, 4, 7}) Measure(b, in, n, K);
}

}  // namespace

int main(int argc, char** argv) {
  const size_t step = argc > 1 ? std::strtoul(argv[1], nullptr, 10) : 64;
  const size_t limit = argc > 2 ? std::strtoul(argv[2], nullptr, 10) : 4096;
  const char* only = argc > 3 ? argv[3] : nullptr;
  g_span1 = argc > 4;
  std::printf("# target %s\n", corvus::active_target());
  const Band bands[] = {{"erf", corvus::erf, -3.0, 3.0},
                        {"erfcinv", corvus::erfcinv, 1e-300, 2.0},
                        {"erfinv", corvus::erfinv, -0.5, 0.5},
                        {"exp", corvus::exp, -20.0, 20.0}};
  for (const Band& b : bands) {
    if (only && std::string_view(only) != b.name) continue;
    std::mt19937_64 rng(4096);
    std::uniform_real_distribution<double> u(b.lo, b.hi);
    std::vector<double> in(512 * 8 + 8);
    for (auto& v : in) v = u(rng);
    for (size_t K = 0; K < limit; K += step) AtSlide(b, in, K);
  }
  std::printf("# addr %p\n", static_cast<const void*>(g_addr));
  return 0;
}
