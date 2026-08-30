// Measures max ULP deviation of corvus::log and corvus::log1p against the
// mpmath-generated correctly-rounded references (tests/data/log_reference.txt
// and log1p_reference.txt), one binary because the two functions share the
// log TU and the log_dd core. Buckets name the kernels' own hard zones:
//   log:   subnormal-input (the 2^600 prescale seam), near-1 (|x-1| <= 1/16,
//          the centred-mantissa table's relative-accuracy hot zone),
//          general.
//   log1p: corner (x <= -1/2, where 1 + x is Sterbenz-exact and the
//          cancellation is deepest), near-0 (|x| <= 2^-10, where the result
//          degenerates to x - x^2/2), general.
// The log reference carries the exact-zero row log(1) = +0; it is held to
// bit identity (#14 N6 -- UlpDiff cannot see a -0 regression there).
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "corvus/corvus.h"
#include "expect_target.h"
#include "ulp_utils.h"

namespace {

using corvus_test::SameBits;
using corvus_test::UlpDiff;

using Fn = void (*)(std::span<const double>, std::span<double>);

// #14 N4: masked-tail split -- neither subspan is a lane multiple for any
// tier, so the masked LoadN/StoreN path always runs.
void SplitCall(Fn fn, const std::vector<double>& in, std::vector<double>& out) {
  const size_t n = in.size();
  const size_t split = n - 3;
  fn(std::span<const double>(in).subspan(0, split),
     std::span<double>(out).subspan(0, split));
  fn(std::span<const double>(in).subspan(split),
     std::span<double>(out).subspan(split));
}

// Gate PINNED to measured, no margin (docs/ACCURACY.md), 2026-08-30 on
// AVX3_ZEN4 native AND the SSE2 cap: max ULP = 0 in EVERY bucket of both
// functions -- every reference row is correctly rounded (the core's
// ~2^-17 ulp tie window never bites on these sets). A future regen that
// lands a legitimate tie-window row re-pins this with the evidence.
constexpr uint64_t kMaxUlp = 0;

constexpr double kMinNormal = 0x1p-1022;

struct Region {
  const char* name;
  uint64_t max_ulp = 0;
  size_t n = 0;
  size_t miss = 0;
  double worst_x = 0.0;
};

void Accumulate(Region& r, double x, uint64_t u) {
  ++r.n;
  if (u > 0) ++r.miss;
  if (u > r.max_ulp) {
    r.max_ulp = u;
    r.worst_x = x;
  }
}

void Report(const Region& r) {
  std::printf(
      "  %-14s n=%6zu  max ULP=%3llu (gate %llu)  not-CR: %zu (%.2f%%)  "
      "worst x=%.17g\n",
      r.name, r.n, static_cast<unsigned long long>(r.max_ulp),
      static_cast<unsigned long long>(kMaxUlp), r.miss,
      r.n ? 100.0 * static_cast<double>(r.miss) / static_cast<double>(r.n)
          : 0.0,
      r.worst_x);
}

int Gate(const char* label, std::initializer_list<const Region*> regions) {
  int rc = 0;
  for (const Region* r : regions) {
    Report(*r);
    // #14 N10: an empty bucket is a gate that never ran; every bucket is a
    // designed reference stratum.
    if (r->n == 0) {
      std::fprintf(stderr, "FAIL: %s %s bucket is empty (n=0) -- vacuous gate\n",
                   label, r->name);
      rc = 1;
    }
    if (r->max_ulp > kMaxUlp) {
      std::fprintf(stderr, "FAIL: %s %s max ULP %llu exceeds gate %llu\n",
                   label, r->name, static_cast<unsigned long long>(r->max_ulp),
                   static_cast<unsigned long long>(kMaxUlp));
      rc = 1;
    }
  }
  return rc;
}

int RunLog(const char* path) {
  const auto rows = corvus_test::LoadRef(path, 2, 10000);
  std::vector<double> in, want;
  for (const auto& row : rows) {
    in.push_back(corvus_test::ParseDouble(row.tok[0], path, row.line));
    want.push_back(corvus_test::ParseDouble(row.tok[1], path, row.line));
  }
  std::vector<double> got(in.size());
  SplitCall(corvus::log, in, got);

  Region sub{"subnormal-in"};
  Region near1{"near-1"};
  Region general{"general"};
  int rc = 0;
  size_t zero_rows = 0;

  for (size_t i = 0; i < in.size(); ++i) {
    const double x = in[i];
    uint64_t u;
    if (want[i] == 0.0) {
      // log(1) = +0 exactly; bit identity, not ULP distance (#14 N6).
      ++zero_rows;
      const bool ok = SameBits(got[i], want[i]);
      if (!ok) {
        std::fprintf(stderr,
                     "FAIL: log(%.17g) signed-zero mismatch: got=%.17g\n", x,
                     got[i]);
        rc = 1;
      }
      u = ok ? 0 : ~uint64_t{0};
    } else {
      u = UlpDiff(got[i], want[i]);
    }
    if (x < kMinNormal) {
      Accumulate(sub, x, u);
    } else if (std::fabs(x - 1.0) <= 0.0625) {
      Accumulate(near1, x, u);
    } else {
      Accumulate(general, x, u);
    }
  }

  std::printf("log (%s):\n", path);
  rc |= Gate("log", {&sub, &near1, &general});
  if (zero_rows == 0) {
    std::fprintf(stderr,
                 "FAIL: log exact-zero row (x=1) missing from reference -- "
                 "vacuous N6 check\n");
    rc = 1;
  }
  return rc;
}

int RunLog1p(const char* path) {
  const auto rows = corvus_test::LoadRef(path, 2, 10000);
  std::vector<double> in, want;
  for (const auto& row : rows) {
    in.push_back(corvus_test::ParseDouble(row.tok[0], path, row.line));
    want.push_back(corvus_test::ParseDouble(row.tok[1], path, row.line));
  }
  std::vector<double> got(in.size());
  SplitCall(corvus::log1p, in, got);

  Region corner{"corner"};
  Region near0{"near-0"};
  Region general{"general"};
  int rc = 0;

  for (size_t i = 0; i < in.size(); ++i) {
    const double x = in[i];
    uint64_t u;
    if (want[i] == 0.0) {
      // No exact-zero row exists today (x = 0 is excluded); guard the next
      // regen (#14 N6).
      const bool ok = SameBits(got[i], want[i]);
      if (!ok) {
        std::fprintf(stderr,
                     "FAIL: log1p(%.17g) signed-zero mismatch: got=%.17g\n", x,
                     got[i]);
        rc = 1;
      }
      u = ok ? 0 : ~uint64_t{0};
    } else {
      u = UlpDiff(got[i], want[i]);
    }
    if (x <= -0.5) {
      Accumulate(corner, x, u);
    } else if (std::fabs(x) <= 0x1p-10) {
      Accumulate(near0, x, u);
    } else {
      Accumulate(general, x, u);
    }
  }

  std::printf("log1p (%s):\n", path);
  rc |= Gate("log1p", {&corner, &near0, &general});
  return rc;
}

}  // namespace

int main(int argc, char** argv) {
  if (!corvus_test::ReportAndCheckTarget()) return 2;

  const char* log_path = argc > 1 ? argv[1] : "tests/data/log_reference.txt";
  const char* log1p_path =
      argc > 2 ? argv[2] : "tests/data/log1p_reference.txt";

  int rc = 0;
  rc |= RunLog(log_path);
  rc |= RunLog1p(log1p_path);
  if (rc == 0) std::printf("PASS\n");
  return rc;
}
