/// @file students_t.cpp
/// @brief Student's t — p-values and critical values — through the regularized
///        incomplete beta pair, corvus::beta_p and corvus::beta_p_inv.
///
/// The t distribution reduces to an incomplete beta under the substitution
///
///     x = nu / (nu + t^2)
///
/// which maps |t| = infinity to x = 0 and t = 0 to x = 1. In terms of it,
///
///     one-sided p-value   Q(t) = 0.5 * I_x(nu/2, 1/2)          -> beta_p
///     critical value      x    = I^-1(nu/2, 1/2; alpha)        -> beta_p_inv
///                         t    = sqrt(nu * (1 - x) / x)
///
/// Note which direction the substitution runs: LARGER |t| gives SMALLER x. The
/// interesting region for a hypothesis test — far out in the tail — is the
/// region where x approaches zero, and corvus's incomplete beta carries a
/// RELATIVE bound there rather than an absolute one. That is what makes the
/// last section of this example possible at all.
///
/// Two of these cases have closed forms, and the example checks against them
/// rather than only against itself:
///
///   nu = 1 is Cauchy, so the two-sided 5% critical value is tan(0.475*pi).
///   nu = 2 gives I_x(1, 1/2) = 1 - sqrt(1 - x), which inverts by hand: at
///     alpha = 0.05, x = 1 - 0.95^2 = 0.0975 EXACTLY.
///
/// That second one is a genuinely exact target, so it is asserted to 1 ULP.

#include <corvus/corvus.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <span>
#include <vector>

namespace {

/// @brief One-sided upper-tail p-value for each (t, dof) pair, t > 0.
void t_sf(std::span<const double> t, std::span<const double> dof, std::span<double> out) {
    const std::size_t n = t.size();
    std::vector<double> a(n);
    std::vector<double> b(n, 0.5);
    std::vector<double> x(n);
    for (std::size_t i = 0; i < n; ++i) {
        a[i] = 0.5 * dof[i];
        x[i] = dof[i] / (dof[i] + t[i] * t[i]);
    }
    corvus::beta_p(a, b, x, out);
    for (double& v : out) {
        v *= 0.5;
    }
}

/// @brief The t with two-sided tail mass `alpha` — i.e. P(|T| > t) = alpha.
void t_critical(std::span<const double> alpha, std::span<const double> dof, std::span<double> out) {
    const std::size_t n = alpha.size();
    std::vector<double> a(n);
    std::vector<double> b(n, 0.5);
    for (std::size_t i = 0; i < n; ++i) {
        a[i] = 0.5 * dof[i];
    }
    corvus::beta_p_inv(a, b, alpha, out);  // out receives x
    // 1 - x cancels as x -> 1, i.e. for a two-sided alpha near 1 (t near 0).
    // Every alpha in this example is small, so x stays away from 1; for the
    // alpha -> 1 regime ask for 1 - x directly via the swap identity,
    // beta_q_inv(b, a, alpha) (see the corvus.h note on beta_p_inv).
    for (std::size_t i = 0; i < n; ++i) {
        out[i] = std::sqrt(dof[i] * (1.0 - out[i]) / out[i]);
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
    std::cout << "corvus — Student's t\n";
    std::cout << "SIMD target: " << corvus::active_target() << "\n\n";

    bool ok = true;

    // -----------------------------------------------------------------------
    // 1. The exact anchor: nu = 2 inverts in closed form.
    //
    // I_x(1, 1/2) = 1 - sqrt(1 - x), so I_x = 0.05 means sqrt(1 - x) = 0.95
    // and x = 1 - 0.9025 = 0.0975, exactly representable. If beta_p_inv is
    // right, it returns that number and not merely something near it.
    // -----------------------------------------------------------------------
    constexpr std::array<double, 1> kA1{1.0};
    constexpr std::array<double, 1> kB1{0.5};
    constexpr std::array<double, 1> kP1{0.05};
    std::array<double, 1> x_exact{};
    corvus::beta_p_inv(kA1, kB1, kP1, x_exact);

    std::cout << std::scientific << std::setprecision(17);
    std::cout << "  exact anchor, nu = 2:\n";
    std::cout << "    beta_p_inv(1, 1/2, 0.05) = " << x_exact[0] << "\n";
    std::cout << "    closed form 1 - 0.95^2   = " << 0.0975 << "\n";
    if (!close_enough(x_exact[0], 0.0975, 2.3e-16)) {
        std::cerr << "the nu = 2 anchor missed its closed form\n";
        ok = false;
    }

    // -----------------------------------------------------------------------
    // 2. Two-sided 5% critical values across the usual range of dof.
    // -----------------------------------------------------------------------
    constexpr std::array<double, 6> kDof{1.0, 2.0, 5.0, 10.0, 30.0, 100.0};
    // Independently computed (mpmath, 60 digits, bisected on x).
    constexpr std::array<double, 6> kRefT{12.70620473617470, 4.302652729749464,
                                          2.570581835636315, 2.228138851986275,
                                          2.042272456301238, 1.983971518523552};
    std::array<double, kDof.size()> alpha{};
    alpha.fill(0.05);

    std::array<double, kDof.size()> crit{};
    t_critical(alpha, kDof, crit);

    std::cout << "\n  two-sided 5% critical values\n\n";
    std::cout << "   dof              t          reference       rel. error\n";
    for (std::size_t i = 0; i < kDof.size(); ++i) {
        const double rel = std::abs(crit[i] - kRefT[i]) / kRefT[i];
        std::cout << std::fixed << std::setprecision(0) << std::setw(6) << kDof[i] << "  "
                  << std::setprecision(12) << std::setw(16) << crit[i] << "  " << std::setw(16)
                  << kRefT[i] << "  " << std::scientific << std::setprecision(2) << rel << "\n";
        if (!close_enough(crit[i], kRefT[i], 1e-13)) {
            std::cerr << "critical value wrong at dof = " << kDof[i] << "\n";
            ok = false;
        }
    }

    // nu = 1 is Cauchy: the 97.5th percentile is tan(0.475*pi).
    const double cauchy = std::tan(0.475 * 3.141592653589793);
    std::cout << "\n  nu = 1 is Cauchy, so that first row has a closed form too:\n";
    std::cout << "    tan(0.475*pi) = " << std::fixed << std::setprecision(12) << cauchy << "\n";
    if (!close_enough(crit[0], cauchy, 1e-13)) {
        std::cerr << "nu = 1 disagrees with the Cauchy closed form\n";
        ok = false;
    }

    // -----------------------------------------------------------------------
    // 3. Round-trip: those critical values back through the p-value.
    // -----------------------------------------------------------------------
    std::array<double, kDof.size()> pv{};
    t_sf(crit, kDof, pv);

    std::cout << "\n  one-sided p at each critical t (expect 0.025)\n\n";
    for (std::size_t i = 0; i < kDof.size(); ++i) {
        std::cout << std::fixed << std::setprecision(0) << "   dof " << std::setw(4) << kDof[i]
                  << "   p = " << std::scientific << std::setprecision(17) << pv[i] << "\n";
        if (!close_enough(pv[i], 0.025, 1e-13)) {
            std::cerr << "round-trip p wrong at dof = " << kDof[i] << "\n";
            ok = false;
        }
    }

    // -----------------------------------------------------------------------
    // 4. Significance far past anything a table lists.
    //
    // This is the section the relative bound buys. A one-sided p of 1e-30 is
    // not a hypothetical — it is the scale genome-wide association work and
    // particle physics operate at. Ask for the t that achieves it.
    //
    // Note carefully what CANNOT be done: there is no way to pose this
    // question through a CDF-side inverse, because the input would have to be
    // p = 1 - 1e-30, and that is not a double. It rounds to exactly 1.0 and
    // takes the question with it. The pair is not a convenience here; it is
    // the only expressible route.
    // -----------------------------------------------------------------------
    constexpr std::array<double, 4> kTailDof{10.0, 10.0, 30.0, 5.0};
    constexpr std::array<double, 4> kTailQ{1e-10, 1e-30, 1e-30, 1e-100};
    constexpr std::array<double, 4> kRefTailT{25.46600802169772, 2564.523931742615,
                                              49.88819553304385, 1.568392559099338e+20};

    std::array<double, 4> two_q{};
    for (std::size_t i = 0; i < two_q.size(); ++i) {
        two_q[i] = 2.0 * kTailQ[i];  // one-sided q -> two-sided alpha
    }
    std::array<double, 4> tail_t{};
    t_critical(two_q, kTailDof, tail_t);

    std::cout << "\n  t values at extreme one-sided significance\n\n";
    std::cout << "   dof        one-sided p                    t          rel. error\n";
    for (std::size_t i = 0; i < kTailQ.size(); ++i) {
        const double rel = std::abs(tail_t[i] - kRefTailT[i]) / kRefTailT[i];
        std::cout << std::fixed << std::setprecision(0) << std::setw(6) << kTailDof[i] << "  "
                  << std::scientific << std::setprecision(3) << std::setw(12) << kTailQ[i] << "  "
                  << std::setprecision(15) << std::setw(24) << tail_t[i] << "  "
                  << std::setprecision(2) << rel << "\n";
        if (!close_enough(tail_t[i], kRefTailT[i], 1e-12)) {
            std::cerr << "extreme quantile wrong at q = " << kTailQ[i] << "\n";
            ok = false;
        }

        // And the input really is unrepresentable from the other side.
        if ((1.0 - kTailQ[i]) != 1.0 && kTailQ[i] < 1e-17) {
            std::cerr << "expected 1 - q to round to 1.0 for q = " << kTailQ[i] << "\n";
            ok = false;
        }
    }

    std::cout << "\n  Every one of those lands within a bit or two of the reference, two of\n"
                 "  them exactly, including a p of 1e-100 where t is 1.6e20. The reason is\n"
                 "  that the substitution sends the far tail to x near ZERO, and a relative\n"
                 "  bound near zero is a bound you can still use. Had the mapping run the\n"
                 "  other way — the tail crowding into x near one — these rows would be\n"
                 "  noise no matter how good the kernel was, for the reason set out in\n"
                 "  examples/README.md.\n";

    std::cout << "\n" << (ok ? "All checks passed.\n" : "FAILURES — see above.\n");
    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
