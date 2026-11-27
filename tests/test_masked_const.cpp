// Guard: a kernel's ordinary path must not depend on which of the image's
// pages happen to be mapped.
//
// The defect it holds (v1.2.1): erfcinv ended with unconditional selects of
// its special values, which the compiler turned into merge-masked loads from
// memory. On ordinary inputs the mask is all-false, the load reads nothing
// and so never maps its page, and while that page is unmapped each execution
// takes a slow path: about 90 ns a vector on AVX3_ZEN4, 4x the whole central
// kernel (see ErfcinvVec in src/erfinv-inl.h). It showed only in executables
// whose link left such a constant on a page nothing else read.
//
// The check needs no particular link. It trims the working set, which unmaps
// every page, and times a central band: the calls map back only what they
// really read. Then it reads one byte of every image page and times the same
// band again. A kernel with no such load reads the same in both states; the
// unfixed erfcinv read 4.1x. The limit is 2x: a ratio of two minima in one
// process, so machine speed and ambient load cancel.
//
// Windows only: it leans on trimming the working set and on walking the PE
// section table. Elsewhere it reports that it did not run and passes.
#include <algorithm>
#include <chrono>
#include <cstdio>
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
using Fn = void (*)(std::span<const double>, std::span<double>) noexcept;

constexpr size_t kN = 4096;
constexpr int kReps = 60;
constexpr double kLimit = 2.0;

double MinNsPerValue(Fn f, const std::vector<double>& in,
                     std::vector<double>& out) {
  double best = 1e300;
  for (int r = 0; r < kReps; ++r) {
    const auto t0 = Clock::now();
    f(in, out);
    const auto t1 = Clock::now();
    best = std::min(best,
                    std::chrono::duration<double, std::nano>(t1 - t0).count());
  }
  return best / static_cast<double>(kN);
}

// Unmaps every page of the process; each comes back on its next real access.
void TrimWorkingSet() {
  SetProcessWorkingSetSize(GetCurrentProcess(), static_cast<SIZE_T>(-1),
                           static_cast<SIZE_T>(-1));
}

// Reads one byte from every page of every section of this executable.
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

bool Check(const char* label, Fn f, double lo, double hi) {
  std::mt19937_64 rng(4096);
  std::uniform_real_distribution<double> u(lo, hi);
  std::vector<double> in(kN), out(kN);
  for (auto& v : in) v = u(rng);

  TrimWorkingSet();
  const double trimmed = MinNsPerValue(f, in, out);
  MapImage();
  const double mapped = MinNsPerValue(f, in, out);
  const double ratio = trimmed / mapped;
  const bool ok = ratio < kLimit;
  std::printf("%-22s trimmed %6.2f  mapped %6.2f ns/value  ratio %5.2f  %s\n",
              label, trimmed, mapped, ratio, ok ? "ok" : "FAIL");
  return ok;
}

}  // namespace

int main() {
  if (!corvus_test::ReportAndCheckTarget()) return 2;
  bool ok = true;
  ok &= Check("erfinv  [-0.5, 0.5]", corvus::erfinv, -0.5, 0.5);
  ok &= Check("erfcinv [0.5, 1.5]", corvus::erfcinv, 0.5, 1.5);
  if (!ok) {
    std::printf(
        "FAIL: a band slows down when the image's pages are unmapped "
        "(limit %.1fx)\n",
        kLimit);
    return 1;
  }
  std::printf("PASS\n");
  return 0;
}

#else  // !_WIN32

int main() {
  if (!corvus_test::ReportAndCheckTarget()) return 2;
  std::printf("not run: Windows only\n");
  return 0;
}

#endif
