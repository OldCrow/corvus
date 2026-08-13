// Measures max ULP deviation of corvus::gamma_p / corvus::gamma_q against
// the mpmath / exact-Temme reference, broken down the way the kernel is
// actually built: by region AND by whether the value under test is the side
// the kernel computed DIRECTLY or the 1 (-) direct complement.
//
// The split matters more here than in the other ULP tests. Every region
// computes the smaller of the pair directly and gets the other by
// subtraction from one; the direct side is the one carrying a relative
// accuracy claim down to the subnormal band, while the complement is always
// >= ~0.4 and its bound is a different (easier) statement. Reporting them
// together would average a hard claim with an easy one.
//
// The routing below is re-derived from the same src/gamma_data.h constants
// the kernel reads, so the two cannot drift apart silently.
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#include "corvus/corvus.h"
#include "expect_target.h"
#include "src/gamma_data.h"
#include "src/lgamma_data.h"  // R4's a-bound is lgamma's centre-2 zone, shifted

namespace {

// Gates pinned to the values measured on AVX2 native plus the
// SSE4/SSSE3/SSE2 caps (Kaby Lake, AppleClang) -- max ULP was identical in
// every cell on all four tiers; only two not-CR counts moved by 1. No
// margin, per house rule. Indexed [region][0 = direct, 1 = complement];
// cells that hold no points for a function carry 0 so any routing change
// that populates them trips the gate and forces a re-measure.
constexpr uint64_t kGateP[4][2] = {
    {2, 0},  // R1: direct P max 2
    {0, 1},  // R2: complement P max 1
    {2, 1},  // R3
    {0, 0},  // R4: complement P max 0 (correctly rounded)
};
constexpr uint64_t kGateQ[4][2] = {
    {0, 4},  // R1: complement Q max 4 (worst a=1.5+ulp, x=2.25)
    {2, 0},  // R2: direct Q max 2
    {2, 1},  // R3
    {1, 0},  // R4: direct Q max 1
};

int64_t OrderedBits(double x) {
  int64_t b;
  std::memcpy(&b, &x, sizeof(b));
  return b < 0 ? (INT64_MIN - b) : b;
}

uint64_t UlpDiff(double a, double b) {
  return static_cast<uint64_t>(std::llabs(OrderedBits(a) - OrderedBits(b)));
}

struct Region {
  const char* name;
  uint64_t bound;
  uint64_t max_ulp = 0;
  size_t n = 0;
  size_t miss = 0;
  double worst_a = 0.0;
  double worst_x = 0.0;
};

// Region codes, in the kernel's own order.
enum : int { kR1 = 0, kR2 = 1, kR3 = 2, kR4 = 3 };

// Re-derivation of GammaVec's router. `want_p` selects which function's
// overlap rule applies; `direct_is_p` reports which side that region
// computes directly at this (a, x).
int Route(double a, double x, bool want_p, bool* direct_is_p) {
  const bool small = a < corvus::detail::kGammaAT;
  const bool xle = x <= a + 1.0;
  const bool lo = 2.0 * x <= a;
  const bool hi = x >= 2.0 * a;
  const bool box = a <= corvus::detail::kLgammaZoneHi - 1.0 && x <= 4.0;

  const bool r1 = small ? xle : lo;
  const bool r3 = !small && !lo && !hi;

  const bool r4 = want_p ? (box && !xle) : box;
  if (r4) {
    *direct_is_p = false;
    return kR4;
  }
  if (r1) {
    *direct_is_p = true;
    return kR1;
  }
  if (r3) {
    *direct_is_p = x < a;
    return kR3;
  }
  *direct_is_p = false;
  return kR2;
}

bool LoadReference(const char* path, std::vector<double>* a,
                   std::vector<double>* x, std::vector<double>* p,
                   std::vector<double>* q) {
  std::ifstream f(path);
  if (!f) {
    std::fprintf(stderr, "cannot open reference file: %s\n", path);
    return false;
  }
  std::string sa, sx, sp, sq;
  while (f >> sa >> sx >> sp >> sq) {
    a->push_back(std::strtod(sa.c_str(), nullptr));
    x->push_back(std::strtod(sx.c_str(), nullptr));
    p->push_back(std::strtod(sp.c_str(), nullptr));
    q->push_back(std::strtod(sq.c_str(), nullptr));
  }
  if (a->size() < 10000) {
    std::fprintf(stderr, "reference file suspiciously small: %zu lines\n",
                 a->size());
    return false;
  }
  return true;
}

int ReportRegions(const char* label, const Region* r, int n) {
  int rc = 0;
  for (int i = 0; i < n; ++i) {
    std::printf(
        "%-8s %-14s n=%6zu  max ULP=%3llu (gate %llu)  not-CR: %zu (%.2f%%)  "
        "worst a=%.17g x=%.17g\n",
        label, r[i].name, r[i].n,
        static_cast<unsigned long long>(r[i].max_ulp),
        static_cast<unsigned long long>(r[i].bound), r[i].miss,
        r[i].n ? 100.0 * static_cast<double>(r[i].miss) /
                     static_cast<double>(r[i].n)
               : 0.0,
        r[i].worst_a, r[i].worst_x);
    if (r[i].n > 0 && r[i].max_ulp > r[i].bound) {
      std::fprintf(stderr, "FAIL: %s %s exceeds gate\n", label, r[i].name);
      rc = 1;
    }
  }
  return rc;
}

int Measure(const char* label, bool want_p, const std::vector<double>& a,
            const std::vector<double>& x, const std::vector<double>& got,
            const std::vector<double>& want) {
  // Index [region][0 = direct side, 1 = complement].
  const auto& g = want_p ? kGateP : kGateQ;
  Region reg[8] = {
      {"R1 series dir", g[0][0]},   {"R1 series cmp", g[0][1]},
      {"R2 cf     dir", g[1][0]},   {"R2 cf     cmp", g[1][1]},
      {"R3 temme  dir", g[2][0]},   {"R3 temme  cmp", g[2][1]},
      {"R4 smalla dir", g[3][0]},   {"R4 smalla cmp", g[3][1]},
  };
  for (size_t i = 0; i < a.size(); ++i) {
    bool direct_is_p = false;
    const int code = Route(a[i], x[i], want_p, &direct_is_p);
    const bool is_direct = (direct_is_p == want_p);
    Region& r = reg[2 * code + (is_direct ? 0 : 1)];
    const uint64_t u = UlpDiff(got[i], want[i]);
    ++r.n;
    if (u > 0) ++r.miss;
    if (u > r.max_ulp) {
      r.max_ulp = u;
      r.worst_a = a[i];
      r.worst_x = x[i];
    }
  }
  return ReportRegions(label, reg, 8);
}

}  // namespace

int main(int argc, char** argv) {
  if (!corvus_test::ReportAndCheckTarget()) return 2;

  const char* p_path = argc > 1 ? argv[1] : "tests/data/gamma_p_reference.txt";
  const char* q_path = argc > 2 ? argv[2] : "tests/data/gamma_q_reference.txt";

  int rc = 0;
  {
    std::vector<double> a, x, p, q;
    if (!LoadReference(p_path, &a, &x, &p, &q)) return 2;
    std::vector<double> got(a.size());
    corvus::gamma_p(a, x, got);
    rc |= Measure("gamma_p", true, a, x, got, p);
  }
  {
    std::vector<double> a, x, p, q;
    if (!LoadReference(q_path, &a, &x, &p, &q)) return 2;
    std::vector<double> got(a.size());
    corvus::gamma_q(a, x, got);
    rc |= Measure("gamma_q", false, a, x, got, q);
  }

  if (rc == 0) std::printf("PASS: all regions within gates\n");
  return rc;
}
