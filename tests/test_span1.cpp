// Guard for the length-1 fast path in the span drivers (v1.2.1, #110).
//
// A span of one element skips the loop and the masked tail: the driver
// broadcasts the value, runs the kernel once and stores lane 0. This test
// holds every public span function to the property that makes that safe to
// ship as a patch: a value computed alone is the same bits as the same value
// computed inside a long span. The long span has a length that is not a
// multiple of any lane count and mixes regions lane by lane, so a mismatch
// would be either the fast path or a lane leaking into its neighbour.
//
// The functions with scalar entry points are held to the same value too.
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <random>
#include <span>
#include <vector>

#include "corvus/corvus.h"
#include "expect_target.h"

namespace {

using Span = std::span<const double>;
using Out = std::span<double>;

int g_failures = 0;

bool SameBits(double a, double b) {
  if (std::isnan(a) && std::isnan(b)) return true;
  uint64_t ua, ub;
  std::memcpy(&ua, &a, 8);
  std::memcpy(&ub, &b, 8);
  return ua == ub;
}

void Report(const char* name, size_t n, size_t bad) {
  std::printf("%-12s n=%6zu  %s\n", name, n, bad == 0 ? "ok" : "FAIL");
  if (bad != 0) ++g_failures;
}

// Magnitudes from the subnormals to overflow, both signs, the specials, and
// dense draws around the regions most kernels split on.
std::vector<double> UnaryInputs() {
  const double inf = std::numeric_limits<double>::infinity();
  const double nan = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> v = {
      0.0,       -0.0,       1.0,       -1.0,
      2.0,       -2.0,       0.5,       -0.5,
      1.5,       inf,        -inf,      nan,
      0x1p-1074, -0x1p-1074, 0x1p-1022, 0x1.fffffffffffffp+1023};
  for (int e = -1074; e <= 1023; e += 7) {
    v.push_back(std::ldexp(1.37, e));
    v.push_back(-std::ldexp(1.37, e));
  }
  std::mt19937_64 rng(1021);
  std::uniform_real_distribution<double> a(-3.0, 3.0), b(0.0, 2.0),
      c(-40.0, 40.0), d(0.0, 800.0);
  for (int i = 0; i < 300; ++i) {
    v.push_back(a(rng));
    v.push_back(b(rng));
    v.push_back(c(rng));
    v.push_back(d(rng));
  }
  if (v.size() % 2 == 0) v.push_back(0.75);  // odd: never a lane multiple
  return v;
}

template <class F>
void CheckUnary(const char* name, F fn) {
  const std::vector<double> in = UnaryInputs();
  std::vector<double> batch(in.size());
  fn(Span(in), Out(batch));
  size_t bad = 0;
  for (size_t i = 0; i < in.size(); ++i) {
    double one;
    fn(Span(&in[i], 1), Out(&one, 1));
    if (!SameBits(one, batch[i])) {
      if (bad++ < 3)
        std::printf("  %s(%a): alone %a, in span %a\n", name, in[i], one,
                    batch[i]);
    }
  }
  Report(name, in.size(), bad);
}

const double kShapes[] = {1e-3, 0.5, 1.0, 2.5, 10.0, 100.0, 1e4, 1e6};
const double kUnit[] = {0.0, 1e-300, 1e-10, 0.01,        0.3,
                        0.5, 0.7,    0.99,  1.0 - 1e-12, 1.0};
const double kPositive[] = {0.0, 1e-300, 1e-10, 0.01, 0.5, 1.0,
                            3.0, 9.5,    50.0,  1e3,  1e5, 1e7};

template <class F, class S, size_t NX>
void CheckBinary(const char* name, F fn, S scalar, const double (&xs)[NX]) {
  std::vector<double> a, x;
  for (double ai : kShapes)
    for (double xi : xs) {
      a.push_back(ai);
      x.push_back(xi);
    }
  a.push_back(3.0);  // odd length
  x.push_back(xs[NX / 2]);
  std::vector<double> batch(a.size());
  fn(Span(a), Span(x), Out(batch));
  size_t bad = 0;
  for (size_t i = 0; i < a.size(); ++i) {
    double one;
    fn(Span(&a[i], 1), Span(&x[i], 1), Out(&one, 1));
    const bool ok = SameBits(one, batch[i]) &&
                    (scalar == nullptr || SameBits(scalar(a[i], x[i]), one));
    if (!ok && bad++ < 3)
      std::printf("  %s(%a, %a): alone %a, in span %a\n", name, a[i], x[i], one,
                  batch[i]);
  }
  Report(name, a.size(), bad);
}

template <class F, class S>
void CheckTernary(const char* name, F fn, S scalar) {
  std::vector<double> a, b, x;
  for (double ai : kShapes)
    for (double bi : kShapes)
      for (double xi : kUnit) {
        a.push_back(ai);
        b.push_back(bi);
        x.push_back(xi);
      }
  a.push_back(3.0);  // odd length
  b.push_back(4.0);
  x.push_back(0.25);
  std::vector<double> batch(a.size());
  fn(Span(a), Span(b), Span(x), Out(batch));
  size_t bad = 0;
  for (size_t i = 0; i < a.size(); ++i) {
    double one;
    fn(Span(&a[i], 1), Span(&b[i], 1), Span(&x[i], 1), Out(&one, 1));
    const bool ok =
        SameBits(one, batch[i]) && SameBits(scalar(a[i], b[i], x[i]), one);
    if (!ok && bad++ < 3)
      std::printf("  %s(%a, %a, %a): alone %a, in span %a\n", name, a[i], b[i],
                  x[i], one, batch[i]);
  }
  Report(name, a.size(), bad);
}

using Scalar2 = double (*)(double, double);

}  // namespace

int main() {
  if (!corvus_test::ReportAndCheckTarget()) return 1;

#define CORVUS_UNARY(fn) \
  CheckUnary(#fn, [](Span in, Out out) { corvus::fn(in, out); })
  CORVUS_UNARY(erf);
  CORVUS_UNARY(erfc);
  CORVUS_UNARY(erfinv);
  CORVUS_UNARY(erfcinv);
  CORVUS_UNARY(lgamma);
  CORVUS_UNARY(digamma);
  CORVUS_UNARY(trigamma);
  CORVUS_UNARY(i0);
  CORVUS_UNARY(i1);
  CORVUS_UNARY(i0e);
  CORVUS_UNARY(i1e);
  CORVUS_UNARY(exp);
  CORVUS_UNARY(expm1);
  CORVUS_UNARY(log);
  CORVUS_UNARY(log1p);
  CORVUS_UNARY(log1pmx);
  CORVUS_UNARY(cos);
  CORVUS_UNARY(sin);
#undef CORVUS_UNARY

#define CORVUS_BINARY(fn, scalar, xs) \
  CheckBinary(                        \
      #fn, [](Span a, Span x, Out out) { corvus::fn(a, x, out); }, scalar, xs)
  CORVUS_BINARY(lbeta, Scalar2{nullptr}, kPositive);
  CORVUS_BINARY(gamma_p, Scalar2{corvus::gamma_p}, kPositive);
  CORVUS_BINARY(gamma_q, Scalar2{corvus::gamma_q}, kPositive);
  CORVUS_BINARY(gamma_p_inv, Scalar2{corvus::gamma_p_inv}, kUnit);
  CORVUS_BINARY(gamma_q_inv, Scalar2{corvus::gamma_q_inv}, kUnit);
#undef CORVUS_BINARY

#define CORVUS_TERNARY(fn)                                                    \
  CheckTernary(                                                               \
      #fn, [](Span a, Span b, Span x, Out out) { corvus::fn(a, b, x, out); }, \
      [](double a, double b, double x) { return corvus::fn(a, b, x); })
  CORVUS_TERNARY(beta_p);
  CORVUS_TERNARY(beta_q);
  CORVUS_TERNARY(beta_p_inv);
  CORVUS_TERNARY(beta_q_inv);
#undef CORVUS_TERNARY

  if (g_failures != 0) {
    std::printf(
        "FAIL: %d function(s) differ between a length-1 span and a "
        "long one\n",
        g_failures);
    return 1;
  }
  std::printf(
      "PASS: a length-1 span matches a long span on all 27 span "
      "functions\n");
  return 0;
}
