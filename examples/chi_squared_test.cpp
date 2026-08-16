/// @file chi_squared_test.cpp
/// @brief Chi-squared p-values and critical values from corvus::gamma_p,
///        corvus::gamma_q and corvus::gamma_p_inv.
///
/// Chi-squared with k degrees of freedom is a Gamma, so the regularized
/// incomplete gamma pair covers the whole distribution:
///
///     CDF        F(x)     = P(k/2, x/2)          -> gamma_p
///     p-value    1 - F(x) = Q(k/2, x/2)          -> gamma_q
///     critical   F^-1(p)  = 2 * P^-1(k/2, p)     -> gamma_p_inv
///
/// A significance test wants the UPPER tail, and `gamma_q` computes it
/// directly rather than as `1 - gamma_p`. That matters more here than for most
/// distributions: interesting p-values are small by definition, and a p-value
/// obtained by subtracting a number near 1 from 1 has no significant digits
/// left. corvus's routing always evaluates whichever of P/Q is the smaller,
/// so its bound stays RELATIVE all the way into the subnormals — a p-value of
/// 1e-107 is as accurate as one of 0.05.
///
/// Note the batch shape below: gamma_p takes spans for the shape AND the
/// argument. Evaluating many x at one fixed dof means materialising the shape
/// as an array of equal values. That is the API being honest about what it
/// does — it evaluates elementwise over parameter arrays, and does not
/// silently broadcast a scalar.

#include "corvus/corvus.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <span>
#include <vector>

namespace {

/// @brief Fill `out` with the chi-squared CDF at each (x, dof) pair.
void chi2_cdf(std::span<const double> x, std::span<const double> dof, std::span<double> out) {
    std::vector<double> a(x.size());
    std::vector<double> half_x(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) {
        a[i] = 0.5 * dof[i];
        half_x[i] = 0.5 * x[i];
    }
    corvus::gamma_p(a, half_x, out);
}

/// @brief Upper-tail probability — the p-value — at each (x, dof) pair.
void chi2_sf(std::span<const double> x, std::span<const double> dof, std::span<double> out) {
    std::vector<double> a(x.size());
    std::vector<double> half_x(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) {
        a[i] = 0.5 * dof[i];
        half_x[i] = 0.5 * x[i];
    }
    corvus::gamma_q(a, half_x, out);
}

/// @brief Critical value: the x with CDF(x) = p, at each (p, dof) pair.
void chi2_ppf(std::span<const double> p, std::span<const double> dof, std::span<double> out) {
    std::vector<double> a(p.size());
    for (std::size_t i = 0; i < p.size(); ++i) {
        a[i] = 0.5 * dof[i];
    }
    corvus::gamma_p_inv(a, p, out);
    for (double& v : out) {
        v *= 2.0;
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
    std::cout << "corvus — chi-squared test statistics\n";
    std::cout << "SIMD target: " << corvus::active_target() << "\n\n";

    bool ok = true;

    // -----------------------------------------------------------------------
    // 1. The textbook 5% critical values, recovered from the quantile.
    // -----------------------------------------------------------------------
    constexpr std::array<double, 6> kDof{1.0, 2.0, 3.0, 5.0, 10.0, 30.0};
    std::array<double, kDof.size()> p95{};
    p95.fill(0.95);

    std::array<double, kDof.size()> crit{};
    std::array<double, kDof.size()> back{};
    chi2_ppf(p95, kDof, crit);
    chi2_cdf(crit, kDof, back);

    std::cout << "  5% critical values (upper tail)\n\n";
    std::cout << "   dof    critical x       CDF at that x\n";
    for (std::size_t i = 0; i < kDof.size(); ++i) {
        std::cout << std::fixed << std::setprecision(0) << std::setw(6) << kDof[i] << "  "
                  << std::setprecision(6) << std::setw(12) << crit[i] << "  " << std::setprecision(16)
                  << back[i] << "\n";
        // Round-trip: the quantile and the CDF have to agree, and this checks
        // both against each other without a reference table.
        if (!close_enough(back[i], 0.95, 1e-14)) {
            std::cerr << "round-trip failed at dof = " << kDof[i] << "\n";
            ok = false;
        }
    }

    // Two values every statistics text lists: 3.841459 at 1 dof, 18.307038 at
    // 10. Worth asserting explicitly — a reader recognises them, and it pins
    // the example to something outside itself.
    if (!close_enough(crit[0], 3.841458820694124, 1e-12)) {
        std::cerr << "chi2(1) 5% critical value wrong: " << crit[0] << "\n";
        ok = false;
    }
    if (!close_enough(crit[4], 18.307038053275146, 1e-12)) {
        std::cerr << "chi2(10) 5% critical value wrong: " << crit[4] << "\n";
        ok = false;
    }

    // -----------------------------------------------------------------------
    // 2. p-values far into the tail, where subtraction has nothing left.
    // -----------------------------------------------------------------------
    constexpr std::array<double, 6> kStatX{3.841459, 20.0, 50.0, 100.0, 200.0, 500.0};
    constexpr std::array<double, 6> kStatDof{1.0, 1.0, 1.0, 1.0, 10.0, 4.0};

    std::array<double, kStatX.size()> pv{};
    std::array<double, kStatX.size()> cdf{};
    chi2_sf(kStatX, kStatDof, pv);
    chi2_cdf(kStatX, kStatDof, cdf);

    std::cout << "\n  p-values, computed directly vs by subtraction\n\n";
    std::cout << "   dof         x     gamma_q (direct)        1 - gamma_p\n";
    for (std::size_t i = 0; i < kStatX.size(); ++i) {
        std::cout << std::fixed << std::setprecision(0) << std::setw(6) << kStatDof[i] << "  "
                  << std::setprecision(3) << std::setw(9) << kStatX[i] << "  " << std::scientific
                  << std::setprecision(12) << pv[i] << "  " << (1.0 - cdf[i]) << "\n";

        // Where both routes still carry digits they must agree; the point of
        // the example is where they stop agreeing, not a disagreement here.
        if (pv[i] > 1e-8 && !close_enough(1.0 - cdf[i], pv[i], 1e-6)) {
            std::cerr << "the two routes disagree while both are still viable, at x = "
                      << kStatX[i] << "\n";
            ok = false;
        }
    }

    std::cout << "\n  The right-hand column reaches zero and stays there. The left-hand\n"
                 "  column keeps going: chi2(4) at x = 500 has p around 6.7e-107, which\n"
                 "  is a real number a test can report, and 1 - F(x) cannot represent it\n"
                 "  at all because F(x) rounded to 1.0 long before.\n"
                 "\n"
                 "  The row that should worry you is x = 50, though, not x = 500. There\n"
                 "  the two columns still LOOK equally healthy and they differ in the\n"
                 "  fifth significant digit — the subtraction has already thrown away\n"
                 "  most of the answer while still returning something plausible. A\n"
                 "  p-value that collapses to zero announces itself; one that is quietly\n"
                 "  wrong in the fourth digit gets published.\n";

    if ((1.0 - cdf[5]) != 0.0) {
        std::cerr << "expected 1 - F(500) to underflow to zero at 4 dof\n";
        ok = false;
    }
    if (!(pv[5] > 0.0 && pv[5] < 1e-100)) {
        std::cerr << "expected a tiny positive p-value from gamma_q, got " << pv[5] << "\n";
        ok = false;
    }

    // -----------------------------------------------------------------------
    // 3. P and Q are complements — check the identity where it is meaningful.
    // -----------------------------------------------------------------------
    std::cout << "\n  P + Q = 1 across the moderate range:\n\n";
    for (std::size_t i = 0; i < 3; ++i) {
        const double sum = cdf[i] + pv[i];
        std::cout << std::fixed << std::setprecision(0) << "   dof " << std::setw(2) << kStatDof[i]
                  << ", x = " << std::setprecision(3) << std::setw(8) << kStatX[i]
                  << "   P + Q - 1 = " << std::scientific << std::setprecision(3) << (sum - 1.0)
                  << "\n";
        if (std::abs(sum - 1.0) > 4e-16) {
            std::cerr << "P + Q strayed from 1 at x = " << kStatX[i] << "\n";
            ok = false;
        }
    }

    std::cout << "\n" << (ok ? "All checks passed.\n" : "FAILURES — see above.\n");
    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
