// Investigation bench (K, v1.2.1): one erf value, the same input every call,
// with the stack slid by K bytes (alloca) for K = 0, 16, ..., 4080. With the
// input fixed, every load the call makes is from a fixed address, so a cost
// that changes with K is a load meeting a stack store 4096*k bytes away.
// Links libcorvus.a only. usage: bench_fixed [x [e|x|i [Klimit [chain]]]]
// (a fifth argument makes each call's input depend on the previous result,
// so the loop runs at the call's latency instead of its throughput)
// Output: K median min  (ns per call)
#include <alloca.h>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <span>

#include "corvus/corvus.h"

namespace {
using Clock = std::chrono::steady_clock;
using Fn = void (*)(std::span<const double>, std::span<double>) noexcept;
volatile double g_sink;
Fn g_fn = corvus::erf;
bool g_chain = false;  // each call's input depends on the previous result
__attribute__((noinline)) double Local(Fn f, double x) {
  double r;
  f(std::span<const double>(&x, 1), std::span<double>(&r, 1));
  return r;
}
__attribute__((noinline)) void AtSlide(size_t K, double x0) {
  char* pad = static_cast<char*>(alloca(K + 16));
  __asm__ volatile("" : : "r"(pad) : "memory");
  constexpr int kCalls = 2048, kReps = 15;
  volatile double vx = x0;
  double t[kReps];
  for (double& ti : t) {
    const double x = vx;
    const auto t0 = Clock::now();
    double acc = 0.0;
    if (g_chain) {
      double xin = x;
      for (int i = 0; i < kCalls; ++i) {
        const double r = Local(g_fn, xin);
        acc += r;
        xin = x + (r - r);  // not folded: -fno-fast-math keeps r - r
      }
    } else {
      for (int i = 0; i < kCalls; ++i) acc += Local(g_fn, x);
    }
    const auto t1 = Clock::now();
    g_sink = acc;
    ti = std::chrono::duration<double, std::nano>(t1 - t0).count() / kCalls;
  }
  std::sort(t, t + kReps);
  std::printf("%5zu %7.2f %7.2f %p\n", K, t[kReps / 2], t[0],
              static_cast<void*>(pad));
}
}  // namespace

int main(int argc, char** argv) {
  const double x0 = argc > 1 ? std::atof(argv[1]) : 0.7;
  if (argc > 2 && argv[2][0] == 'x') g_fn = corvus::exp;
  if (argc > 2 && argv[2][0] == 'i') g_fn = corvus::erfinv;
  std::printf("# target %s x %g\n", corvus::active_target(), x0);
  g_chain = argc > 4;
  const size_t klim = argc > 3 ? std::strtoul(argv[3], nullptr, 10) : 4096;
  for (size_t K = 0; K < klim; K += 16) AtSlide(K, x0);
  return 0;
}
