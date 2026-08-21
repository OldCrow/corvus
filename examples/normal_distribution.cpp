/// @file normal_distribution.cpp
/// @brief The normal CDF and its quantile, built from corvus::erfc and
///        corvus::erfcinv — including the far tail, where the obvious
///        formulation fails and the right one does not.
///
/// Two identities carry the whole example:
///
///     Phi(z)      = 0.5 * erfc(-z / sqrt(2))        lower tail
///     1 - Phi(z)  = 0.5 * erfc( z / sqrt(2))        upper tail
///     Phi^-1(p)   = -sqrt(2) * erfcinv(2p)          quantile / probit
///
/// The pair of tail forms is the point. Both are exact in real arithmetic, so
/// a reader can be forgiven for computing the upper tail as `1 - Phi(z)` — but
/// in floating point that subtraction destroys every significant digit once
/// Phi(z) rounds to 1, which happens around z = 8.3. Choosing the form that
/// computes the SMALL quantity directly is what keeps relative accuracy in the
/// tail, and it costs nothing.

#include <corvus/corvus.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <span>

namespace {

// sqrt(2) to double precision. Spelled out rather than computed so the example
// has no dependency on the compiler folding std::sqrt at compile time.
constexpr double kSqrt2 = 1.4142135623730951;

/// @brief Standard normal CDF for a whole batch: Phi(z) = 0.5 * erfc(-z/sqrt2).
///
/// Note the shape — one corvus call for the entire span, not a call per
/// element. That is the intended idiom: the runtime-dispatched kernel picks the
/// widest SIMD target the CPU supports and processes the array in vector-width
/// chunks, with a masked tail for the remainder. Calling it on one element at a
/// time works, but throws away the reason it exists.
void normal_cdf(std::span<const double> z, std::span<double> out) {
    for (std::size_t i = 0; i < z.size(); ++i) {
        out[i] = -z[i] / kSqrt2;
    }
    corvus::erfc(out, out);  // in-place is fine: same span, elementwise kernel
    for (double& v : out) {
        v *= 0.5;
    }
}

/// @brief Upper-tail probability Q(z) = 1 - Phi(z), computed directly.
///
/// Same function as `1 - normal_cdf(z)` in exact arithmetic, and a completely
/// different function in double. See the tail section in main().
void normal_sf(std::span<const double> z, std::span<double> out) {
    for (std::size_t i = 0; i < z.size(); ++i) {
        out[i] = z[i] / kSqrt2;
    }
    corvus::erfc(out, out);
    for (double& v : out) {
        v *= 0.5;
    }
}

/// @brief Standard normal quantile: Phi^-1(p) = -sqrt(2) * erfcinv(2p).
///
/// corvus's erfcinv is max 1 ULP over its whole domain, including results deep
/// in the subnormal range, so this is a genuine inverse rather than an
/// approximation that degrades where it is most often needed.
void normal_ppf(std::span<const double> p, std::span<double> out) {
    for (std::size_t i = 0; i < p.size(); ++i) {
        out[i] = 2.0 * p[i];
    }
    corvus::erfcinv(out, out);
    for (double& v : out) {
        v *= -kSqrt2;
    }
}

bool close_enough(double got, double want, double rel_tol) {
    if (got == want) {
        return true;
    }
    const double scale = std::abs(want);
    return std::abs(got - want) <= rel_tol * (scale > 0.0 ? scale : 1.0);
}

}  // namespace

int main() {
    std::cout << "corvus — normal distribution\n";
    std::cout << "SIMD target: " << corvus::active_target() << "\n\n";
    std::cout << std::scientific << std::setprecision(16);

    bool ok = true;

    // -----------------------------------------------------------------------
    // 1. CDF and quantile over an everyday range, batched in one call each.
    // -----------------------------------------------------------------------
    constexpr std::array<double, 7> kZ{-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0};
    std::array<double, kZ.size()> p{};
    std::array<double, kZ.size()> back{};

    normal_cdf(kZ, p);
    normal_ppf(p, back);

    std::cout << "  z          Phi(z)                    Phi^-1(Phi(z))\n";
    for (std::size_t i = 0; i < kZ.size(); ++i) {
        std::cout << std::fixed << std::setprecision(1) << std::setw(6) << kZ[i] << "  "
                  << std::scientific << std::setprecision(16) << p[i] << "  " << back[i] << "\n";
        // Round-tripping through both directions is the sharpest cheap check
        // available: it exercises the CDF and the quantile against each other
        // with no reference table involved.
        if (!close_enough(back[i], kZ[i], 1e-14)) {
            std::cerr << "round-trip failed at z = " << kZ[i] << "\n";
            ok = false;
        }
    }

    // Phi(0) = 1/2 exactly — erfc(0) = 1 is exact and halving is exact.
    if (p[3] != 0.5) {
        std::cerr << "Phi(0) should be exactly 0.5, got " << p[3] << "\n";
        ok = false;
    }

    std::cout << "\n  Two details in that table are worth reading, not skipping.\n"
                 "\n"
                 "  The z = 0 row round-trips to NEGATIVE zero. erfcinv(1) is +0, and\n"
                 "  the -sqrt(2) factor carries the sign onto it. -0.0 == 0.0 compares\n"
                 "  true and behaves identically in arithmetic, so this is correct, not\n"
                 "  a defect — but it will show up in printed output and in any code\n"
                 "  that inspects the sign bit.\n"
                 "\n"
                 "  The z = +3 row does NOT round-trip exactly, while z = -3 does. The\n"
                 "  functions are not at fault: Phi(3) = 0.9986 stores its information\n"
                 "  in the digits just below 1, so the double simply does not have the\n"
                 "  resolution left to name z again. The negative side keeps its value\n"
                 "  small and loses nothing. That asymmetry is the same effect the tail\n"
                 "  section below is about, showing up early and mildly.\n";

    // -----------------------------------------------------------------------
    // 2. The far tail, where the formulation decides the answer.
    // -----------------------------------------------------------------------
    std::cout << "\nupper tail — Q(z) = 1 - Phi(z), two ways\n\n";
    std::cout << "  z     0.5*erfc(z/sqrt2)         1 - Phi(z)\n";

    constexpr std::array<double, 5> kTail{2.0, 5.0, 8.0, 10.0, 20.0};
    std::array<double, kTail.size()> direct{};
    std::array<double, kTail.size()> cdf{};

    normal_sf(kTail, direct);
    normal_cdf(kTail, cdf);

    for (std::size_t i = 0; i < kTail.size(); ++i) {
        const double naive = 1.0 - cdf[i];
        std::cout << std::fixed << std::setprecision(1) << std::setw(6) << kTail[i] << "  "
                  << std::scientific << std::setprecision(16) << direct[i] << "  " << naive << "\n";
    }

    std::cout << "\n  Past z = 8.3 or so, Phi(z) rounds to 1.0 and the naive column\n"
                 "  collapses to exactly zero. The direct column is still carrying\n"
                 "  full relative accuracy — Q(20) is around 2.75e-89, and corvus's\n"
                 "  erfc holds its bound down through the subnormals.\n";

    // The naive route must actually have failed, or this example is not
    // demonstrating what its text claims.
    if ((1.0 - cdf[4]) != 0.0) {
        std::cerr << "expected 1 - Phi(20) to underflow to zero\n";
        ok = false;
    }
    if (!(direct[4] > 0.0)) {
        std::cerr << "expected a positive Q(20) from the direct form\n";
        ok = false;
    }

    // -----------------------------------------------------------------------
    // 3. The quantile follows the tail all the way down.
    // -----------------------------------------------------------------------
    std::array<double, kTail.size()> tail_z{};
    normal_ppf(direct, tail_z);

    std::cout << "\nquantile of those tail probabilities (expect -z)\n\n";
    for (std::size_t i = 0; i < kTail.size(); ++i) {
        std::cout << std::scientific << std::setprecision(3) << "  p = " << direct[i]
                  << "   ->  " << std::fixed << std::setprecision(12) << std::setw(18) << tail_z[i]
                  << "\n";
        if (!close_enough(tail_z[i], -kTail[i], 1e-13)) {
            std::cerr << "tail quantile failed at z = " << kTail[i] << "\n";
            ok = false;
        }
    }

    std::cout << "\n" << (ok ? "All checks passed.\n" : "FAILURES — see above.\n");
    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
