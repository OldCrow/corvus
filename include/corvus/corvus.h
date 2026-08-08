#ifndef CORVUS_CORVUS_H_
#define CORVUS_CORVUS_H_

/// \file corvus.h
/// \brief Public API: SIMD-vectorized statistical special functions.
///
/// This is the only installed header. The SIMD backend is an implementation
/// detail: no Highway types appear here, and consumers need only link
/// libcorvus. All functions dispatch at runtime to the best SIMD tier the
/// CPU supports (SSE2..AVX-512, NEON).
///
/// Common contract for all batch functions:
///  - `in` and `out` must have the same length.
///  - Exact aliasing is allowed (`in.data() == out.data()`); partial overlap
///    is undefined behavior.
///  - No allocation, no exceptions, thread-safe (stateless; the dispatch
///    pointer is resolved on first call).
///  - Accuracy bounds are measured against a correctly-rounded mpmath
///    oracle and enforced per-tier by the test suite; see docs/ACCURACY.md.

#include <cstddef>
#include <span>

namespace corvus {

inline constexpr int kVersionMajor = 0;
inline constexpr int kVersionMinor = 1;
inline constexpr int kVersionPatch = 0;

/// \brief Name of the SIMD target selected by runtime dispatch.
/// \return A static string such as "AVX2", "SSE4", or "NEON". Note Highway
///   names AVX-512 tiers "AVX3" (see README naming note).
const char* active_target();

/// \brief out[i] = erf(in[i]).
///
/// Max 1 ULP over the full domain on validated tiers. Specials: erf(+/-0)
/// = +/-0, erf(+/-inf) = +/-1, NaN propagates (payload preserved).
void erf(std::span<const double> in, std::span<double> out);

/// \brief out[i] = erfc(in[i]).
///
/// Accuracy on validated tiers: max 1 ULP for |x| <= 6 and for subnormal
/// results; max 2 ULP for normal-result x > 6 (bounded by the tail
/// polynomial fit; see docs/ACCURACY.md). Specials: erfc(-inf) = 2,
/// erfc(+inf) = 0, results underflow gradually past x ~ 26.5, NaN
/// propagates (payload preserved).
void erfc(std::span<const double> in, std::span<double> out);

/// \brief out[i] = log |Gamma(in[i])| (C's lgamma, SciPy's gammaln).
///
/// Defined on the whole real axis. No sign output: like SciPy's gammaln,
/// and unlike C's signgam, the sign of Gamma is not reported.
///
/// Accuracy on validated tiers: see docs/ACCURACY.md for the per-region
/// table. Relative accuracy holds on the positive axis, including
/// arbitrarily close to the zeros at x = 1 and x = 2, which are exact. On
/// the negative axis lgamma also passes through zero wherever
/// |Gamma(x)| = 1 -- infinitely many points with no closed form -- and
/// near those the bound is absolute rather than relative; the measured
/// split is documented.
///
/// Specials: lgamma(1) = lgamma(2) = +0; +inf at every pole (x = 0 and the
/// negative integers), on overflow (x above ~2.556e305), and at both
/// infinities; NaN propagates (payload preserved).
void lgamma(std::span<const double> in, std::span<double> out);

/// \brief out[i] = erfinv(in[i]), the inverse error function.
///
/// Accuracy on validated tiers: see docs/ACCURACY.md. Specials:
/// erfinv(+/-0) = +/-0, erfinv(+/-1) = +/-inf, |in[i]| > 1 gives NaN, NaN
/// propagates (payload preserved).
void erfinv(std::span<const double> in, std::span<double> out);

/// \brief out[i] = erfcinv(in[i]), the inverse complementary error function.
///
/// Useful directly as the normal quantile: probit(p) = -sqrt(2)*erfcinv(2p).
/// Accuracy on validated tiers: see docs/ACCURACY.md. Specials:
/// erfcinv(0) = +inf, erfcinv(2) = -inf, in[i] outside [0, 2] gives NaN, NaN
/// propagates (payload preserved).
void erfcinv(std::span<const double> in, std::span<double> out);

/// \brief out[i] = P(a[i], x[i]), the regularized lower incomplete gamma
///   function (SciPy's `gammainc`).
///
/// P(a, x) = 1/Gamma(a) * integral from 0 to x of t^(a-1) e^-t dt, i.e. the
/// CDF of a Gamma(a, 1) variate. "Regularized" is the division by Gamma(a):
/// P is in [0, 1] and P + Q = 1.
///
/// Three spans, all the same length: `a`, `x` and `out`. Whichever of the
/// pair is the smaller at a given (a, x) is the one computed directly, so
/// the accuracy bound is relative on that side however far it underflows;
/// the other side is >= ~0.4 by construction. See docs/ACCURACY.md for the
/// measured per-region table.
///
/// Specials: P(a, 0) = +0 and P(+inf, x) = +0 for finite a, x; P(0, x) = 1
/// and P(a, +inf) = 1 for finite a, x; P(0, 0), P(+inf, +inf) and any
/// negative a or x give NaN; NaN propagates (payload preserved).
void gamma_p(std::span<const double> a, std::span<const double> x,
             std::span<double> out);

/// \brief out[i] = Q(a[i], x[i]) = 1 - P(a[i], x[i]), the regularized upper
///   incomplete gamma function (SciPy's `gammaincc`).
///
/// Q(a, x) = 1/Gamma(a) * integral from x to infinity of t^(a-1) e^-t dt --
/// the survival function of a Gamma(a, 1) variate. Same span contract and
/// same compute-the-smaller-side rule as gamma_p; results underflow
/// gradually and then to zero as x grows, and stay relatively accurate
/// throughout the representable range.
///
/// Specials: Q(a, 0) = 1 and Q(+inf, x) = 1 for finite a, x; Q(0, x) = +0
/// and Q(a, +inf) = +0 for finite a, x; Q(0, 0), Q(+inf, +inf) and any
/// negative a or x give NaN; NaN propagates (payload preserved).
void gamma_q(std::span<const double> a, std::span<const double> x,
             std::span<double> out);

/// \brief out[i] = I_x(a, b), the regularized incomplete beta function
///   (SciPy's `betainc`).
///
/// I_x(a, b) = 1/B(a, b) * integral from 0 to x of t^(a-1) (1-t)^(b-1) dt,
/// i.e. the CDF of a Beta(a, b) variate. "Regularized" is the division by the
/// beta function B(a, b) = Gamma(a)Gamma(b)/Gamma(a+b): I is in [0, 1] and
/// beta_p + beta_q = 1.
///
/// Four spans, all the same length: `a`, `b`, `x` and `out`. Whichever of the
/// pair is the smaller at a given (a, b, x) is the one computed directly, so
/// the accuracy bound is relative on that side however far it underflows; the
/// other side is >= ~0.4 by construction. See docs/ACCURACY.md for the
/// measured per-region table.
///
/// Specials, with `a` and `b` finite and positive unless stated: I_0(a, b) =
/// +0 and I_1(a, b) = 1; a = 0 or b = +inf puts all the mass at 0, so
/// I_x = 1 for x in (0, 1]; b = 0 or a = +inf puts all the mass at 1, so
/// I_x = +0 for x in [0, 1). A degenerate parameter meeting the boundary its
/// own mass sits on gives NaN (a = 0 at x = 0; b = +inf at x = 0; b = 0 at
/// x = 1; a = +inf at x = 1), as does any two-way degeneracy among
/// {a in {0, +inf}, b in {0, +inf}}. Negative a or b, and x outside [0, 1],
/// give NaN; NaN propagates (payload preserved).
void beta_p(std::span<const double> a, std::span<const double> b,
            std::span<const double> x, std::span<double> out);

/// \brief out[i] = 1 - I_x(a, b) = I_{1-x}(b, a), the complementary
///   regularized incomplete beta function.
///
/// The survival function of a Beta(a, b) variate. Same span contract and the
/// same compute-the-smaller-side rule as beta_p, so results stay relatively
/// accurate as they underflow toward either end of the domain.
///
/// Specials are beta_p's with the two limits exchanged: Q_0(a, b) = 1,
/// Q_1(a, b) = +0; a = 0 or b = +inf gives +0 for x in (0, 1]; b = 0 or
/// a = +inf gives 1 for x in [0, 1). The NaN cases are identical to beta_p's.
void beta_q(std::span<const double> a, std::span<const double> b,
            std::span<const double> x, std::span<double> out);

/// \brief out[i] = psi(in[i]), the digamma function (SciPy's `digamma`).
///
/// psi(x) = d/dx log Gamma(x) = Gamma'(x)/Gamma(x), defined on the whole real
/// axis.
///
/// Accuracy on validated tiers: see docs/ACCURACY.md for the per-region
/// table. The bound is RELATIVE on the positive axis, including arbitrarily
/// close to the unique positive zero at x ~ 1.4616321, which the kernel
/// reproduces by construction. On the negative axis psi also has a zero
/// between every consecutive pair of poles -- infinitely many points with no
/// closed form -- so there the bound is relative where |psi| >= 1 and
/// ABSOLUTE, of order 2^-53, inside the bands around those zeros; the
/// measured split is documented.
///
/// Specials: psi(+0) = -inf and psi(-0) = +inf (signed-zero pole
/// convention); NaN at every negative-integer pole, which includes -inf and
/// every double <= -2^53 since all of those are integers; +inf at +inf;
/// arguments small enough that -1/x overflows give -+inf accordingly; NaN
/// propagates (payload preserved).
void digamma(std::span<const double> in, std::span<double> out);

}  // namespace corvus

#endif  // CORVUS_CORVUS_H_
