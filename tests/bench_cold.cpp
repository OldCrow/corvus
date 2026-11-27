// Survey: does any ordinary band of any span function slow down when the
// image's pages are unmapped? Manual, Windows only, not registered with
// ctest.
//
// The method is tests/test_masked_const.cpp's, run over every public span
// function: trim the working set (every page unmapped), time the band (the
// calls map back only what they really read), read one byte of every image
// page, time it again. A band whose two readings differ carries a load that
// executes without reading: a merge-masked load of a rare-path constant with
// an all-false mask (see ErfcinvVec in src/erfinv-inl.h, v1.2.1).
//
// A clean reading says only that THIS link exposes no such constant. A
// consumer's link places constants differently; the static count of such
// loads per kernel is the complement.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <functional>
#include <random>
#include <span>
#include <vector>

#include "corvus/corvus.h"
#include "expect_target.h"

#if defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>

namespace {

using Clock = std::chrono::steady_clock;
using Vec = std::vector<double>;
constexpr size_t kN = 4096;
constexpr int kReps = 30;
int g_flagged = 0;

void TrimWorkingSet() {
  SetProcessWorkingSetSize(GetCurrentProcess(), static_cast<SIZE_T>(-1),
                           static_cast<SIZE_T>(-1));
}

void MapImage() {
  const auto* base =
      reinterpret_cast<const unsigned char*>(GetModuleHandleW(nullptr));
  const auto* dos = reinterpret_cast<const IMAGE_DOS_HEADER*>(base);
  const auto* nt =
      reinterpret_cast<const IMAGE_NT_HEADERS*>(base + dos->e_lfanew);
  const IMAGE_SECTION_HEADER* sec = IMAGE_FIRST_SECTION(nt);
  volatile unsigned char sink = 0;
  for (unsigned i = 0; i < nt->FileHeader.NumberOfSections; ++i) {
    for (DWORD o = 0; o < sec[i].Misc.VirtualSize; o += 4096) {
      sink = base[sec[i].VirtualAddress + o];
    }
  }
  (void)sink;
}

double MinNs(const std::function<void()>& call) {
  double best = 1e300;
  for (int r = 0; r < kReps; ++r) {
    const auto t0 = Clock::now();
    call();
    const auto t1 = Clock::now();
    best = std::min(best,
                    std::chrono::duration<double, std::nano>(t1 - t0).count());
  }
  return best / static_cast<double>(kN);
}

void Survey(const char* label, const std::function<void()>& call) {
  // Twice over: a band that reads the same both times in both states is
  // clean; one pass alone can catch a clock step.
  double ratio[2];
  double trimmed = 0, mapped = 0;
  for (double& r : ratio) {
    TrimWorkingSet();
    trimmed = MinNs(call);
    MapImage();
    mapped = MinNs(call);
    r = trimmed / mapped;
  }
  const bool flag = ratio[0] > 1.15 && ratio[1] > 1.15;
  if (flag) ++g_flagged;
  std::printf("%-30s trimmed %8.2f  mapped %8.2f  ratio %5.2f %5.2f  %s\n",
              label, trimmed, mapped, ratio[0], ratio[1], flag ? "<==" : "");
}

Vec Draw(double lo, double hi, unsigned seed, bool log2_scale = false) {
  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> u(lo, hi);
  Vec v(kN);
  for (auto& x : v) x = log2_scale ? std::exp2(u(rng)) : u(rng);
  return v;
}

using Unary = void (*)(std::span<const double>, std::span<double>) noexcept;
using Binary = void (*)(std::span<const double>, std::span<const double>,
                        std::span<double>) noexcept;
using Ternary = void (*)(std::span<const double>, std::span<const double>,
                         std::span<const double>, std::span<double>) noexcept;

void U(const char* label, Unary f, double lo, double hi,
       bool log2_scale = false) {
  const Vec in = Draw(lo, hi, 4096, log2_scale);
  Vec out(kN);
  Survey(label, [&] { f(in, out); });
}

void B(const char* label, Binary f, double alo, double ahi, double xlo,
       double xhi) {
  const Vec a = Draw(alo, ahi, 11), x = Draw(xlo, xhi, 12);
  Vec out(kN);
  Survey(label, [&] { f(a, x, out); });
}

void T(const char* label, Ternary f, double lo, double hi) {
  const Vec a = Draw(lo, hi, 21), b = Draw(lo, hi, 22),
            x = Draw(0.01, 0.99, 23);
  Vec out(kN);
  Survey(label, [&] { f(a, b, x, out); });
}

}  // namespace

int main() {
  if (!corvus_test::ReportAndCheckTarget()) return 2;
  std::printf(
      "ns per value; '<==' marks a band slower when unmapped in both "
      "passes\n");
  U("erf [-0.5, 0.5]", corvus::erf, -0.5, 0.5);
  U("erf [-3, 3]", corvus::erf, -3.0, 3.0);
  U("erfc [0, 5]", corvus::erfc, 0.0, 5.0);
  U("erfc [5, 25]", corvus::erfc, 5.0, 25.0);
  U("erfinv [-0.5, 0.5]", corvus::erfinv, -0.5, 0.5);
  U("erfinv [-0.999, 0.999]", corvus::erfinv, -0.999, 0.999);
  U("erfcinv [0.5, 1.5]", corvus::erfcinv, 0.5, 1.5);
  U("erfcinv 2^[-26,-1]", corvus::erfcinv, -26.0, -1.0, true);
  U("lgamma [0.01, 1]", corvus::lgamma, 0.01, 1.0);
  U("lgamma [1, 21]", corvus::lgamma, 1.0, 21.0);
  U("lgamma 2^[5, 40]", corvus::lgamma, 5.0, 40.0, true);
  U("lgamma [-20, -0.1]", corvus::lgamma, -20.0, -0.1);
  U("digamma [0.1, 20]", corvus::digamma, 0.1, 20.0);
  U("digamma 2^[5, 40]", corvus::digamma, 5.0, 40.0, true);
  U("digamma [-20, -0.1]", corvus::digamma, -20.0, -0.1);
  U("trigamma [0.1, 20]", corvus::trigamma, 0.1, 20.0);
  U("trigamma 2^[5, 40]", corvus::trigamma, 5.0, 40.0, true);
  U("i0 [0, 10]", corvus::i0, 0.0, 10.0);
  U("i0 [10, 700]", corvus::i0, 10.0, 700.0);
  U("i1 [0, 10]", corvus::i1, 0.0, 10.0);
  U("i1 [10, 700]", corvus::i1, 10.0, 700.0);
  U("i0e [0, 700]", corvus::i0e, 0.0, 700.0);
  U("i1e [0, 700]", corvus::i1e, 0.0, 700.0);
  U("exp [-20, 20]", corvus::exp, -20.0, 20.0);
  U("expm1 [-1, 1]", corvus::expm1, -1.0, 1.0);
  U("expm1 [-20, 20]", corvus::expm1, -20.0, 20.0);
  U("log [0.01, 100]", corvus::log, 0.01, 100.0);
  U("log1p [-0.5, 2]", corvus::log1p, -0.5, 2.0);
  U("log1pmx [-0.5, 2]", corvus::log1pmx, -0.5, 2.0);
  U("cos [-10, 10]", corvus::cos, -10.0, 10.0);
  U("sin [-10, 10]", corvus::sin, -10.0, 10.0);
  U("cos 2^[10, 60]", corvus::cos, 10.0, 60.0, true);
  B("lbeta a,b in [0.5, 50]", corvus::lbeta, 0.5, 50.0, 0.5, 50.0);
  B("gamma_p a[0.5,50] x[0.1,60]", corvus::gamma_p, 0.5, 50.0, 0.1, 60.0);
  B("gamma_q a[0.5,50] x[0.1,60]", corvus::gamma_q, 0.5, 50.0, 0.1, 60.0);
  B("gamma_p_inv a[0.5,50]", corvus::gamma_p_inv, 0.5, 50.0, 0.01, 0.99);
  B("gamma_q_inv a[0.5,50]", corvus::gamma_q_inv, 0.5, 50.0, 0.01, 0.99);
  T("beta_p a,b in [0.5, 50]", corvus::beta_p, 0.5, 50.0);
  T("beta_q a,b in [0.5, 50]", corvus::beta_q, 0.5, 50.0);
  T("beta_p_inv a,b in [0.5, 50]", corvus::beta_p_inv, 0.5, 50.0);
  T("beta_q_inv a,b in [0.5, 50]", corvus::beta_q_inv, 0.5, 50.0);
  std::printf("%d band(s) flagged\n", g_flagged);
  return 0;
}

#else  // !_WIN32

int main() {
  if (!corvus_test::ReportAndCheckTarget()) return 2;
  std::printf("not run: Windows only\n");
  return 0;
}

#endif
