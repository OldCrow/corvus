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
///  - `in` and `out` must have the same length. Release builds do NOT
///    check lengths — a shorter span anywhere is undefined behavior
///    (out-of-bounds access). Debug builds assert the contract
///    (HWY_DASSERT in the shared driver), so the mistake fails loudly
///    there instead.
///  - Exact aliasing is allowed (`in.data() == out.data()`); partial overlap
///    is undefined behavior. For functions taking several input spans,
///    exactly one input may alias `out` (every input at index i is read
///    before out[i] is written); inputs may alias each other.
///  - No allocation, no exceptions, thread-safe (stateless; the dispatch
///    pointer is resolved on first call).
///  - Accuracy bounds are measured against a correctly-rounded mpmath
///    oracle and enforced per-tier by the test suite; see docs/ACCURACY.md.

#include <cstddef>
#include <span>

namespace corvus {

// Must match CMakeLists.txt's project(VERSION) — enforced at configure
// time (a mismatch is a FATAL_ERROR; added after these constants sat at
// 0.1.0 through the v0.2.0 release unnoticed).
inline constexpr int kVersionMajor = 0;
inline constexpr int kVersionMinor = 5;
inline constexpr int kVersionPatch = 0;

/// \brief Name of the SIMD target selected by runtime dispatch.
/// \return A static string such as "AVX2", "SSE4", or "NEON". Note Highway
///   names AVX-512 tiers "AVX3" (see README naming note).
[[nodiscard]] const char* active_target() noexcept;

/// \brief out[i] = erf(in[i]).
///
/// Max 1 ULP over the full domain on validated tiers. Specials: erf(+/-0)
/// = +/-0, erf(+/-inf) = +/-1, NaN propagates (payload preserved).
void erf(std::span<const double> in, std::span<double> out) noexcept;

/// \brief out[i] = erfc(in[i]).
///
/// Accuracy on validated tiers: max 1 ULP for |x| <= 6 and for subnormal
/// results; max 2 ULP for normal-result x > 6 (bounded by the tail
/// polynomial fit; see docs/ACCURACY.md). Specials: erfc(-inf) = 2,
/// erfc(+inf) = 0, results underflow gradually past x ~ 26.5, NaN
/// propagates (payload preserved).
void erfc(std::span<const double> in, std::span<double> out) noexcept;

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
void lgamma(std::span<const double> in, std::span<double> out) noexcept;

/// \brief out[i] = erfinv(in[i]), the inverse error function.
///
/// Accuracy on validated tiers: see docs/ACCURACY.md. Specials:
/// erfinv(+/-0) = +/-0, erfinv(+/-1) = +/-inf, |in[i]| > 1 gives NaN, NaN
/// propagates (payload preserved).
void erfinv(std::span<const double> in, std::span<double> out) noexcept;

/// \brief out[i] = erfcinv(in[i]), the inverse complementary error function.
///
/// Useful directly as the normal quantile: probit(p) = -sqrt(2)*erfcinv(2p).
/// Accuracy on validated tiers: see docs/ACCURACY.md. Specials:
/// erfcinv(0) = +inf, erfcinv(2) = -inf, in[i] outside [0, 2] gives NaN, NaN
/// propagates (payload preserved).
void erfcinv(std::span<const double> in, std::span<double> out) noexcept;

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
             std::span<double> out) noexcept;

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
             std::span<double> out) noexcept;

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
            std::span<const double> x, std::span<double> out) noexcept;

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
            std::span<const double> x, std::span<double> out) noexcept;

/// \brief out[i] = ln B(a[i], b[i]) = lgamma(a) + lgamma(b) - lgamma(a+b).
///
/// Positive-parameter domain: a > 0 and b > 0, finite; anything else
/// returns NaN. (SciPy's `betaln` accepts non-positive arguments through
/// |Gamma|; corvus deliberately does not -- no statistical consumer needs
/// them.) Computed through the same double-double lgamma-difference
/// machinery as beta_p/beta_q's prefactor, so the a+b cancellation that
/// costs a plain three-lgamma assembly its accuracy for large parameters
/// is removed analytically. Accuracy: CORRECTLY ROUNDED on every measured
/// row of every band and tier — 0 ULP where |ln B| >= 1, half-ulp
/// absolute in the band around ln B's zero curve through (1,1), where the
/// result itself is ill-conditioned (docs/ACCURACY.md). Symmetric in
/// (a, b); saturates to -inf
/// where the true ln B falls below the double range (both parameters
/// huge); NaN propagates as a quiet NaN (the input payload is not preserved).
void lbeta(std::span<const double> a, std::span<const double> b,
           std::span<double> out) noexcept;

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
void digamma(std::span<const double> in, std::span<double> out) noexcept;

/// \brief out[i] = psi_1(in[i]), the trigamma function (SciPy's
///   `polygamma(1, .)`).
///
/// psi_1(x) = d/dx psi(x) = d^2/dx^2 log Gamma(x), defined on the whole real
/// axis.
///
/// Accuracy on validated tiers: see docs/ACCURACY.md for the per-region
/// table. The bound is RELATIVE EVERYWHERE, with no absolute band anywhere --
/// psi_1(x) = sum over n >= 0 of 1/(x + n)^2 is a sum of squares, hence
/// strictly positive wherever it is finite, so unlike digamma and lgamma this
/// function has no zeros for a relative metric to break down at. Its
/// negative-axis minimum is 8.933, at x ~ -0.4957.
///
/// Specials: every pole is a DOUBLE pole and therefore sign-unambiguous, so
/// the answer is +inf at both +0 and -0, at every negative-integer pole
/// (which includes -inf and every double <= -2^53, all of which are
/// integers), and at any argument small enough that 1/x^2 overflows --
/// subnormals of either sign included. psi_1(+inf) = +0; NaN propagates
/// (payload preserved).
void trigamma(std::span<const double> in, std::span<double> out) noexcept;

/// \brief out[i] = x with P(a[i], x) = p[i], the inverse of the regularized
///   lower incomplete gamma function in its second argument (SciPy's
///   `gammaincinv`).
///
/// Equivalently the quantile function of a Gamma(a, 1) variate. Three spans,
/// all the same length: `a`, `p` and `out`.
///
/// Whichever of p and 1 - p is <= 1/2 is the one solved against, the switch
/// being exact, so the accuracy bound is relative on both sides of the
/// median and down to subnormal answers. See docs/ACCURACY.md for the
/// measured per-regime table. Note that for a above ~3e34 the whole
/// transition from P = 0 to P = 1 happens inside one ulp of x, so the answer
/// there is x = a for every p in (0, 1) that is not a hard limit; that is
/// the correctly rounded result, not a shortcut.
///
/// Specials: p = 0 gives +0 and p = 1 gives +inf; p outside [0, 1] gives
/// NaN, as do a <= 0 and a = +inf; NaN propagates (payload preserved).
void gamma_p_inv(std::span<const double> a, std::span<const double> p,
                 std::span<double> out) noexcept;

/// \brief out[i] = x with Q(a[i], x) = q[i], the inverse of the regularized
///   upper incomplete gamma function in its second argument (SciPy's
///   `gammainccinv`).
///
/// The survival-function quantile of a Gamma(a, 1) variate. Same span
/// contract, same exact side switch and the same accuracy statement as
/// gamma_p_inv -- the two are one kernel with one bit of orientation, so a
/// q of 1e-300 is solved as accurately as the corresponding p of
/// 1 - 1e-300 could never be represented.
///
/// Specials: q = 0 gives +inf and q = 1 gives +0; q outside [0, 1] gives
/// NaN, as do a <= 0 and a = +inf; NaN propagates (payload preserved).
void gamma_q_inv(std::span<const double> a, std::span<const double> q,
                 std::span<double> out) noexcept;

/// \brief out[i] = x with I_x(a[i], b[i]) = p[i], the inverse of the
///   regularized incomplete beta function in its third argument (SciPy's
///   `betaincinv`).
///
/// Equivalently the quantile function of a Beta(a, b) variate. Four spans, all
/// the same length: `a`, `b`, `p` and `out`.
///
/// Whichever of p and 1 - p is <= 1/2 is the one solved against, the switch
/// being exact, and the kernel then solves for whichever end of [0, 1] the
/// answer approaches -- so the accuracy bound is relative on both sides of the
/// median and down to subnormal answers. See docs/ACCURACY.md for the measured
/// per-regime table.
///
/// THE SWAP IDENTITY IS THE LOSSLESS ROUTE FOR AN ANSWER NEAR 1. A double
/// cannot hold 1 - 1e-30 to better than 1e-16 absolute, so an x near the top
/// of the range is resolution-limited no matter how it is computed. The
/// identity I_x(a, b) = 1 - I_{1-x}(b, a) turns that into a question with a
/// representable answer: `1 - x` AT FULL RELATIVE PRECISION is
/// `beta_p_inv(b, a, q)` (equivalently `beta_q_inv(b, a, p)`), because the two
/// calls solve the SAME equation with the parameters and the probability side
/// exchanged. Callers who need the distance from 1 rather than x itself should
/// use that form; it is exactly what the kernel does internally.
///
/// Two named regimes carry a weaker guarantee, both by the nature of the
/// problem rather than by the method. Where BOTH a and b are below ~1e-16 the
/// distribution's interior density is ~4 min(a, b), so an interior x is not
/// resolvable to 1 ULP by any double-precision inverse; there the returned x
/// satisfies a BACKWARD-error contract instead (its forward value is within a
/// few ulp of the requested probability). Where the shape-side parameter
/// exceeds ~1e32 the entire transition from I = 0 to I = 1 happens inside one
/// or two ulp of x, so every interior probability has effectively the same
/// answer; that is the correctly rounded result, not a shortcut.
///
/// Specials, with `a` and `b` finite and positive unless stated: p = 0 gives
/// +0 and p = 1 gives 1. a = 0 or b = +inf puts all the mass at 0, so every
/// quantile is +0; b = 0 or a = +inf puts all the mass at 1, so every quantile
/// is 1. Any two-way degeneracy among {a in {0, +inf}, b in {0, +inf}} gives
/// NaN, as do negative a or b and p outside [0, 1]; NaN propagates (payload
/// preserved).
void beta_p_inv(std::span<const double> a, std::span<const double> b,
                std::span<const double> p, std::span<double> out) noexcept;

/// \brief out[i] = x with 1 - I_x(a[i], b[i]) = q[i], the inverse of the
///   complementary regularized incomplete beta function in its third argument.
///
/// The survival-function quantile of a Beta(a, b) variate. Same span contract,
/// same exact side switch and the same accuracy statement as beta_p_inv -- the
/// two are one kernel with one bit of orientation, so a q of 1e-300 is solved
/// as accurately as the corresponding p of 1 - 1e-300 could never be
/// represented. The swap identity documented on beta_p_inv applies here too:
/// `1 - x` at full relative precision is `beta_q_inv(b, a, p)`.
///
/// Specials: q = 0 gives 1 and q = 1 gives +0; the parameter degeneracies and
/// the NaN cases are identical to beta_p_inv's.
void beta_q_inv(std::span<const double> a, std::span<const double> b,
                std::span<const double> q, std::span<double> out) noexcept;

/// \brief out[i] = I0(in[i]), the modified Bessel function of the first
///   kind, order 0 (SciPy's `iv(0, .)` / `i0`).
///
/// Method: a truncated power series in q = x^2/4 for |x| <= 8 (all terms
/// positive, well conditioned), and a Chebyshev fit in 1/x times
/// exp(|x|)/sqrt(2*pi*|x|) for |x| > 8. Accuracy: max 1 ULP over the full
/// real axis, measured and gate-pinned on every SIMD tier
/// (docs/ACCURACY.md).
///
/// Specials: I0(0) = 1; I0(+-inf) = +inf; NaN propagates (payload
/// preserved). I0 saturates to +inf for |x| above the last finite double
/// (~713.99, the exact boundary is documented in docs/ACCURACY.md).
void i0(std::span<const double> in, std::span<double> out) noexcept;

/// \brief out[i] = I1(in[i]), the modified Bessel function of the first
///   kind, order 1 (SciPy's `iv(1, .)` / `i1`).
///
/// Same method as i0, on the ODD series/fit pair: I1 is odd, so the sign of
/// a nonzero result matches the sign of the input. Accuracy: max 1 ULP
/// over the full real axis, measured and gate-pinned on every SIMD tier
/// (docs/ACCURACY.md).
///
/// Specials: I1(+0) = +0, I1(-0) = -0; I1(+inf) = +inf, I1(-inf) = -inf; NaN
/// propagates (payload preserved). I1 saturates to +-inf for |x| above the
/// last finite double (~713.99, close to but not identical to i0's
/// boundary; see docs/ACCURACY.md).
void i1(std::span<const double> in, std::span<double> out) noexcept;

/// \brief out[i] = exp(-|in[i]|) * I0(in[i]), the exponentially-scaled
///   modified Bessel function of the first kind, order 0 (SciPy's
///   `ive(0, .)`).
///
/// Useful wherever I0 itself would overflow (I0 overflows past |x| ~ 714):
/// this stays finite over the whole real axis. A(kappa) = i1e(kappa) /
/// i0e(kappa) is the exact mean resultant length of a von Mises(kappa)
/// distribution, computed without ever forming the unscaled I0/I1.
/// Accuracy: max 1 ULP over the full real axis, measured and gate-pinned
/// on every SIMD tier (docs/ACCURACY.md).
///
/// Specials: i0e(0) = 1; i0e(+-inf) = +0; NaN propagates (payload
/// preserved). Never underflows on a finite double (minimum ~3e-155 at
/// DBL_MAX).
void i0e(std::span<const double> in, std::span<double> out) noexcept;

/// \brief out[i] = sign(in[i]) * exp(-|in[i]|) * I1(|in[i]|), the
///   exponentially-scaled modified Bessel function of the first kind,
///   order 1 (SciPy's `ive(1, .)`).
///
/// Same scaling rationale as i0e, on the odd pair; see i0e's doc for the von
/// Mises mean-resultant-length identity. Accuracy: max 1 ULP over the full
/// real axis, measured and gate-pinned on every SIMD tier
/// (docs/ACCURACY.md).
///
/// Specials: i1e(+0) = +0, i1e(-0) = -0; i1e(+inf) = +0, i1e(-inf) = -0
/// (sign follows the input at both zero and infinity); NaN propagates
/// (payload preserved).
void i1e(std::span<const double> in, std::span<double> out) noexcept;

/// \brief out[i] = cos(in[i]), over the FULL double range.
///
/// Quadrant reduction with an exact 4-part pi/2 split for |x| <= 2^23 and a
/// Payne-Hanek reduction beyond it (certified at each exponent's worst
/// reduction cancellation, including binary64's global worst case) -- there
/// is no domain cutoff and no libm fallback on any tier. Accuracy: measured
/// and gate-pinned per SIMD tier (docs/ACCURACY.md).
///
/// Specials: cos(+-0) = 1 exactly; cos(+-inf) = NaN; NaN propagates
/// (payload preserved), evaluated in-vector in every lane position.
void cos(std::span<const double> in, std::span<double> out) noexcept;

/// \brief out[i] = sin(in[i]), over the FULL double range.
///
/// Same reduction and cores as cos (sin is computed from its own quadrant
/// table, not as cos(x - pi/2)). sin is odd by construction: sin(-x) is the
/// exact negation of sin(x) in every lane. Accuracy: measured and
/// gate-pinned per SIMD tier (docs/ACCURACY.md).
///
/// Specials: sin(+-0) = +-0 with the sign preserved exactly;
/// sin(+-inf) = NaN; NaN propagates (payload preserved), evaluated
/// in-vector in every lane position.
void sin(std::span<const double> in, std::span<double> out) noexcept;

}  // namespace corvus

#endif  // CORVUS_CORVUS_H_
