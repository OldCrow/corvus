/// @file gamma_mle_fit.cpp
/// @brief Maximum-likelihood fit of a Gamma shape parameter — many fits at
///        once — using corvus::digamma and corvus::trigamma.
///
/// This is what digamma and trigamma are FOR. The Gamma log-likelihood's
/// derivative in the shape k contains digamma, so the MLE condition is
///
///     log(k) - digamma(k) = log(mean(x)) - mean(log(x)) = s
///
/// where s is the dataset's whole contribution: two sums, computed once. The
/// scale then follows in closed form, theta = mean(x) / k, so the entire fit
/// reduces to solving one scalar equation.
///
/// Newton needs the derivative of the left side, which is 1/k - trigamma(k) —
/// so the pair is exactly what the algorithm consumes, and having both as exact
/// functions rather than finite differences is what makes the iteration
/// converge in a handful of steps from a decent start.
///
/// It does NOT converge to the last bit, and the reason is instructive rather
/// than a defect: the residual subtracts two nearly equal numbers, so its
/// accuracy floor is set by that cancellation and not by digamma's 1 ULP. The
/// program works the arithmetic out at the end.
///
/// THE BATCH DIMENSION HERE IS FITS, NOT DATA POINTS. Each dataset collapses to
/// its own s before any special function is touched; the vector being processed
/// is the vector of current Newton iterates, one per fit. Six fits over five
/// steps means ten corvus calls in total, not ten per dataset. That shape —
/// vectorize across problems, not within one — is often the one that pays in
/// statistical work, and it is easy to miss while looking for a long array.
///
/// Note also what corvus does NOT provide: log() and the like. Basic
/// transcendentals belong to the SIMD backend's own math contrib, and the
/// scalar work below just uses <cmath>. corvus's scope is the special
/// functions that gate statistical work.

#include <corvus/corvus.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <span>

namespace {

/// Minka's starting approximation for the Gamma shape MLE. Good to about four
/// digits over the whole useful range, which is why the iteration below needs
/// so few steps.
double minka_initial_shape(double s) {
    return (3.0 - s + std::sqrt((s - 3.0) * (s - 3.0) + 24.0 * s)) / (12.0 * s);
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
    std::cout << "corvus — Gamma shape MLE, six fits at once\n";
    std::cout << "SIMD target: " << corvus::active_target() << "\n\n";

    bool ok = true;

    // -----------------------------------------------------------------------
    // 0. Two exact anchors, so the rest is not checking corvus against itself.
    //
    // digamma(1) = -gamma (Euler-Mascheroni) and trigamma(1) = pi^2/6 are
    // known in closed form. The fit below uses digamma on both the way in and
    // the way out, so without an outside reference a systematic error could
    // cancel and leave the example looking healthy.
    // -----------------------------------------------------------------------
    constexpr std::array<double, 1> kOne{1.0};
    std::array<double, 1> anchor{};

    corvus::digamma(kOne, anchor);
    if (!close_enough(anchor[0], -0.5772156649015329, 1e-15)) {
        std::cerr << "digamma(1) should be -0.5772156649015329, got " << anchor[0] << "\n";
        ok = false;
    }
    std::cout << std::scientific << std::setprecision(16);
    std::cout << "  digamma(1)  = " << anchor[0] << "   (-Euler-Mascheroni)\n";

    corvus::trigamma(kOne, anchor);
    if (!close_enough(anchor[0], 1.6449340668482264, 1e-15)) {
        std::cerr << "trigamma(1) should be 1.6449340668482264, got " << anchor[0] << "\n";
        ok = false;
    }
    std::cout << "  trigamma(1) = " << anchor[0] << "   (pi^2 / 6)\n";

    // -----------------------------------------------------------------------
    // 1. Six datasets, given by their sufficient statistics.
    //
    // These s values were produced independently (mpmath at 40 digits) from
    // the shapes listed beside them, so recovering those shapes is a genuine
    // end-to-end check rather than a round trip through corvus and back.
    // -----------------------------------------------------------------------
    constexpr std::array<double, 6> kTrueShape{0.5, 1.0, 2.5, 7.0, 40.0, 300.0};
    constexpr std::array<double, 6> kS{
        1.2703628454614782,   0.57721566490153287, 0.21313409122891189,
        0.07312581395684617,  0.012552080079093177, 0.0016675925915637915};

    std::array<double, 6> k{};
    for (std::size_t i = 0; i < k.size(); ++i) {
        k[i] = minka_initial_shape(kS[i]);
    }

    std::cout << "\n  Newton on  f(k) = log(k) - digamma(k) - s,   f'(k) = 1/k - trigamma(k)\n";
    std::cout << "\n  iter    max |step / k|\n";

    std::array<double, 6> psi0{};
    std::array<double, 6> psi1{};

    // A batched Newton runs every lane for the same number of steps: there is
    // no per-lane early exit, so the slowest fit sets the count and the others
    // spend their last iterations taking steps that do nothing. That is the
    // trade — a few wasted lanes in exchange for evaluating the special
    // functions across the whole batch at once.
    //
    // Stop when the step size stops shrinking rather than at a fixed
    // tolerance. Newton here does not converge to zero step, and the reason is
    // the third lesson of this example — see the note printed below.
    constexpr int kMaxIters = 12;
    int iters_used = 0;
    double prev_worst = std::numeric_limits<double>::infinity();
    for (int iter = 0; iter < kMaxIters; ++iter) {
        corvus::digamma(k, psi0);
        corvus::trigamma(k, psi1);

        double worst = 0.0;
        for (std::size_t i = 0; i < k.size(); ++i) {
            const double f = std::log(k[i]) - psi0[i] - kS[i];
            const double fp = 1.0 / k[i] - psi1[i];
            const double step = f / fp;
            k[i] -= step;
            worst = std::max(worst, std::abs(step / k[i]));
        }
        ++iters_used;
        std::cout << std::setw(6) << iter + 1 << "    " << std::scientific << std::setprecision(3)
                  << worst << "\n";
        if (worst == 0.0 || worst >= prev_worst) {
            break;
        }
        prev_worst = worst;
    }

    // -----------------------------------------------------------------------
    // 2. Did it land on the shapes the data came from?
    // -----------------------------------------------------------------------
    std::cout << "\n  true shape        recovered            relative error\n";
    for (std::size_t i = 0; i < k.size(); ++i) {
        const double rel = std::abs(k[i] - kTrueShape[i]) / kTrueShape[i];
        std::cout << std::fixed << std::setprecision(1) << std::setw(11) << kTrueShape[i] << "  "
                  << std::setprecision(15) << std::setw(20) << k[i] << "    " << std::scientific
                  << std::setprecision(2) << rel << "\n";
        if (!close_enough(k[i], kTrueShape[i], 1e-12)) {
            std::cerr << "failed to recover shape " << kTrueShape[i] << "\n";
            ok = false;
        }
    }

    std::cout << "\n  " << iters_used << " iterations for all six, so " << 2 * iters_used
              << " corvus calls in total —\n"
                 "  not that many per fit. The shapes span three orders of magnitude and\n"
                 "  converge together, because Minka's start is uniformly good and the\n"
                 "  derivative is exact rather than a finite difference.\n";

    std::cout << "\n  Two things above are worth accounting for, and they are the same\n"
                 "  thing: the step size stops shrinking around 1e-13 instead of\n"
                 "  reaching zero, and the recovered error grows with the shape — exact\n"
                 "  at the small shapes, a few parts in 1e13 at k = 300.\n"
                 "\n"
                 "  Neither is corvus losing accuracy. digamma is 1 ULP at every one of\n"
                 "  these points. The loss is in OUR residual: at k = 300, log(k) is\n"
                 "  5.7038 and digamma(k) is 5.7021, and f subtracts them to get\n"
                 "  1.6676e-3. Twelve bits cancel. Half an ulp of 5.7 is 4.4e-16, which\n"
                 "  against an answer of 1.7e-3 is a relative error near 3e-13 — and\n"
                 "  that is the floor the iteration hits, and the error it leaves in k.\n"
                 "\n"
                 "  Same shape as the tail problem in the other two examples: an input\n"
                 "  good to 1 ULP, an expression that cancels most of it away, and a\n"
                 "  result whose accuracy is set by the expression rather than by the\n"
                 "  function. Here it is unavoidable — the MLE condition IS that\n"
                 "  difference — so the honest response is to know the floor rather than\n"
                 "  to iterate below it and believe the digits.\n";

    std::cout << "\n" << (ok ? "All checks passed.\n" : "FAILURES — see above.\n");
    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
