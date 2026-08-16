/// @file log_space_counting.cpp
/// @brief Binomial coefficients and Beta normalisers in log space, with
///        corvus::lgamma and corvus::lbeta.
///
/// Combinatorial quantities leave the double range almost immediately —
/// 171! already overflows — so anything counting at scale works in logs. The
/// standard assemblies are
///
///     log C(n,k)   = lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1)
///     log B(a,b)   = lgamma(a) + lgamma(b) - lgamma(a+b)
///
/// and both are differences of large nearly-equal numbers, which by now is a
/// familiar warning sign in these examples. The other four showed the fix as
/// "call the member of the pair that computes your quantity directly". This
/// one shows the other fix: a DEDICATED function whose derivation removes the
/// cancellation before any floating point happens.
///
/// corvus::lbeta is that function. It is not lgamma(a) + lgamma(b) -
/// lgamma(a+b) evaluated carefully; the a+b cancellation is taken out
/// analytically and the result assembled through the beta family's
/// double-double lgamma-difference machinery, with one final rounding. On
/// every measured row it is CORRECTLY ROUNDED — 0 ULP wherever |ln B| >= 1.
///
/// The two connect through an identity worth knowing, since it gives a second
/// route to the binomial coefficient that never forms a large difference:
///
///     C(n,k) = 1 / ( (n+1) * B(k+1, n-k+1) )
///     log C(n,k) = -log(n+1) - lbeta(k+1, n-k+1)
///
/// At n = 1e15 the lgamma route has to subtract two numbers near 3.35e16 to
/// produce an answer near 135. One ulp at that magnitude is 4. The lbeta route
/// never goes near 1e16 at all.

#include "corvus/corvus.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <span>
#include <vector>

namespace {

/// @brief log C(n,k) the usual way — three lgammas and two subtractions.
void log_binom_via_lgamma(std::span<const double> n, std::span<const double> k,
                          std::span<double> out) {
    const std::size_t m = n.size();
    std::vector<double> arg(m);
    std::vector<double> lg_n(m);
    std::vector<double> lg_k(m);
    std::vector<double> lg_nk(m);

    for (std::size_t i = 0; i < m; ++i) {
        arg[i] = n[i] + 1.0;
    }
    corvus::lgamma(arg, lg_n);
    for (std::size_t i = 0; i < m; ++i) {
        arg[i] = k[i] + 1.0;
    }
    corvus::lgamma(arg, lg_k);
    for (std::size_t i = 0; i < m; ++i) {
        arg[i] = n[i] - k[i] + 1.0;
    }
    corvus::lgamma(arg, lg_nk);

    for (std::size_t i = 0; i < m; ++i) {
        out[i] = lg_n[i] - lg_k[i] - lg_nk[i];
    }
}

/// @brief log C(n,k) through the Beta identity — no large difference formed.
void log_binom_via_lbeta(std::span<const double> n, std::span<const double> k,
                         std::span<double> out) {
    const std::size_t m = n.size();
    std::vector<double> a(m);
    std::vector<double> b(m);
    for (std::size_t i = 0; i < m; ++i) {
        a[i] = k[i] + 1.0;
        b[i] = n[i] - k[i] + 1.0;
    }
    corvus::lbeta(a, b, out);
    for (std::size_t i = 0; i < m; ++i) {
        out[i] = -std::log(n[i] + 1.0) - out[i];
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
    std::cout << "corvus — counting in log space\n";
    std::cout << "SIMD target: " << corvus::active_target() << "\n\n";

    bool ok = true;

    // -----------------------------------------------------------------------
    // 1. Exact anchors. B(1,1) = 1 so its log is zero; B(1/2,1/2) = pi.
    // -----------------------------------------------------------------------
    constexpr std::array<double, 2> kAnchorA{1.0, 0.5};
    constexpr std::array<double, 2> kAnchorB{1.0, 0.5};
    std::array<double, 2> anchor{};
    corvus::lbeta(kAnchorA, kAnchorB, anchor);

    std::cout << std::scientific << std::setprecision(17);
    std::cout << "  lbeta(1, 1)     = " << anchor[0] << "   (B = 1, so ln B = 0)\n";
    std::cout << "  lbeta(1/2, 1/2) = " << anchor[1] << "   (B = pi, log pi = "
              << 1.1447298858494002 << ")\n";

    // The first anchor is NOT exactly zero, and should not be expected to be.
    // ln B has a zero curve through (1,1), and a result of zero has no relative
    // sense — every nonzero double is infinitely far from it in relative terms.
    // So corvus documents an ABSOLUTE bound in that band, 0.5*2^-53, and switches
    // to the relative bound only where |ln B| >= 1. This is the same shape as
    // lgamma's negative axis and digamma's zeros: an ill-conditioned band gets
    // an absolute guarantee, and a library that claimed 1 ULP there would be
    // claiming something that cannot be measured, never mind met.
    constexpr double kHalfUlpAbs = 0.5 * 1.1102230246251565e-16;  // 0.5 * 2^-53
    std::cout << "    ln B has a zero here, so the bound is absolute (0.5*2^-53 = "
              << kHalfUlpAbs << "),\n    not relative — and " << std::abs(anchor[0])
              << " is well inside it.\n";
    if (std::abs(anchor[0]) > kHalfUlpAbs) {
        std::cerr << "lbeta(1,1) outside its documented absolute bound: " << anchor[0] << "\n";
        ok = false;
    }
    if (!close_enough(anchor[1], 1.1447298858494002, 2.3e-16)) {
        std::cerr << "lbeta(1/2,1/2) missed log(pi)\n";
        ok = false;
    }

    // -----------------------------------------------------------------------
    // 2. lbeta against the naive assembly.
    // -----------------------------------------------------------------------
    constexpr std::array<double, 5> kA{500.0, 1.0e6, 1.0e12, 1.0e15, 300.0};
    constexpr std::array<double, 5> kB{300.0, 1.0e6, 3.0, 1.0e-5, 300.0};
    // mpmath, 50 digits.
    constexpr std::array<double, 5> kRefLB{-530.94820113822811, -1386300.0033629211,
                                           -82.199916167228693, 11.512574305131876,
                                           -417.47427078333686};

    std::array<double, 5> direct{};
    corvus::lbeta(kA, kB, direct);

    // The naive route, assembled from three lgammas exactly as a user would.
    std::array<double, 5> lg_a{};
    std::array<double, 5> lg_b{};
    std::array<double, 5> lg_ab{};
    std::array<double, 5> sum_ab{};
    for (std::size_t i = 0; i < kA.size(); ++i) {
        sum_ab[i] = kA[i] + kB[i];
    }
    corvus::lgamma(kA, lg_a);
    corvus::lgamma(kB, lg_b);
    corvus::lgamma(sum_ab, lg_ab);

    std::cout << "\n  log B(a,b): the dedicated function vs the three-lgamma assembly\n\n";
    std::cout << "         a            b            lbeta              naive           "
                 "largest term\n";
    for (std::size_t i = 0; i < kA.size(); ++i) {
        const double naive = lg_a[i] + lg_b[i] - lg_ab[i];
        std::cout << std::scientific << std::setprecision(2) << std::setw(10) << kA[i] << "  "
                  << std::setw(11) << kB[i] << "  " << std::fixed << std::setprecision(9)
                  << std::setw(19) << direct[i] << "  " << std::setw(19) << naive << "  "
                  << std::scientific << std::setprecision(2) << std::setw(10) << lg_ab[i] << "\n";
        if (!close_enough(direct[i], kRefLB[i], 1e-15)) {
            std::cerr << "lbeta wrong at a = " << kA[i] << ", b = " << kB[i] << "\n";
            ok = false;
        }
    }

    std::cout << "\n  Read the last column with the two before it. Where the terms being\n"
                 "  subtracted are of the same size as the answer, the routes agree. Where\n"
                 "  the terms are 1e13 and the answer is -82, the naive column has quietly\n"
                 "  lost four digits. Where the terms are 3.35e16 and the answer is 11.5,\n"
                 "  one ulp of the terms is 4 — the answer is smaller than the noise, and\n"
                 "  the naive column is not approximately right, it is unrelated.\n";

    // The extreme row is the fourth: a = 1e15, b = 1e-5.
    const double naive_extreme = lg_a[3] + lg_b[3] - lg_ab[3];
    if (close_enough(naive_extreme, kRefLB[3], 1e-3)) {
        std::cerr << "expected the naive assembly to fail badly at a = 1e15, b = 1e-5\n";
        ok = false;
    }

    // -----------------------------------------------------------------------
    // 3. Binomial coefficients, both routes.
    // -----------------------------------------------------------------------
    constexpr std::array<double, 5> kN{10.0, 1000.0, 100000.0, 1.0e9, 1.0e15};
    constexpr std::array<double, 5> kK{5.0, 500.0, 3.0, 4.0, 4.0};
    constexpr std::array<double, 5> kRefLC{5.5294290875114234, 689.46726156785121,
                                           32.74698692543263, 79.7150095114377,
                                           134.97705174929479};

    std::array<double, 5> via_g{};
    std::array<double, 5> via_b{};
    log_binom_via_lgamma(kN, kK, via_g);
    log_binom_via_lbeta(kN, kK, via_b);

    std::cout << "\n  log C(n,k): three lgammas vs the Beta identity\n\n";
    std::cout << "            n     k      via lgamma        via lbeta         reference\n";
    for (std::size_t i = 0; i < kN.size(); ++i) {
        std::cout << std::scientific << std::setprecision(2) << std::setw(13) << kN[i] << "  "
                  << std::fixed << std::setprecision(0) << std::setw(4) << kK[i] << "  "
                  << std::setprecision(9) << std::setw(15) << via_g[i] << "  " << std::setw(15)
                  << via_b[i] << "  " << std::setw(15) << kRefLC[i] << "\n";
        // Only the lbeta route is held to the full bound.
        if (!close_enough(via_b[i], kRefLC[i], 1e-14)) {
            std::cerr << "lbeta route wrong at n = " << kN[i] << "\n";
            ok = false;
        }
    }

    std::cout << "\n  C(10,5) = 252 either way. At n = 1e15 the lgamma route is wrong in\n"
                 "  the first decimal place — it is subtracting numbers near 3.35e16 —\n"
                 "  while the identity route never forms anything larger than 170.\n";

    if (close_enough(via_g[4], kRefLC[4], 1e-9)) {
        std::cerr << "expected the lgamma route to degrade visibly at n = 1e15\n";
        ok = false;
    }

    // -----------------------------------------------------------------------
    // 4. What it is for: a binomial log-probability at a scale that has no
    //    non-log formulation at all.
    // -----------------------------------------------------------------------
    constexpr std::array<double, 1> kBigN{1.0e6};
    constexpr std::array<double, 1> kBigK{500000.0};
    std::array<double, 1> log_c{};
    log_binom_via_lbeta(kBigN, kBigK, log_c);

    const double log_pmf = log_c[0] + 1.0e6 * std::log(0.5);
    std::cout << "\n  log P(X = 500000) for X ~ Binomial(1e6, 1/2)\n\n";
    std::cout << "    log C(1e6, 5e5) = " << std::fixed << std::setprecision(8) << log_c[0] << "\n";
    std::cout << "    log P           = " << std::setprecision(12) << log_pmf << "\n";
    std::cout << "    P               = " << std::scientific << std::setprecision(9)
              << std::exp(log_pmf) << "\n";
    std::cout << "    sqrt(2/(pi n))  = " << std::sqrt(2.0 / (3.141592653589793 * 1.0e6))
              << "   <- the Stirling limit it should approach\n";

    if (!close_enough(log_c[0], 693140.04701306368, 1e-14)) {
        std::cerr << "log C(1e6, 5e5) wrong\n";
        ok = false;
    }

    // log_c is good to about 1e-14 RELATIVE, which at 693140 is 1e-8 absolute —
    // and the final line subtracts 693147 from 693140 to reach -7.13. Roughly
    // seventeen bits go in that step, so the log-probability cannot be better
    // than about 1e-11 relative however exact lbeta was. Assert the floor, not
    // a number the arithmetic cannot deliver.
    constexpr double kPmfFloor = 1e-10;
    if (!close_enough(log_pmf, -7.1335468816268648, kPmfFloor)) {
        std::cerr << "log PMF outside even its cancellation floor: " << log_pmf << "\n";
        ok = false;
    }

    std::cout << "\n  C(1e6, 5e5) itself is about e^693140, so there is no version of this\n"
                 "  calculation that touches the coefficient directly. Logs are not an\n"
                 "  optimisation here; they are the only representation available.\n"
                 "\n"
                 "  And the last line pays the usual toll one final time: log P comes from\n"
                 "  subtracting 693147 from 693140, so about seventeen bits go, and the\n"
                 "  answer is good to ~1e-11 rather than ~1e-16. lbeta being correctly\n"
                 "  rounded does not survive that step, and nothing could make it — the\n"
                 "  probability genuinely is a small difference of two large logs. Knowing\n"
                 "  the floor is the whole of the remedy.\n";

    std::cout << "\n" << (ok ? "All checks passed.\n" : "FAILURES — see above.\n");
    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
