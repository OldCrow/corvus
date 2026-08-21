/// @file von_mises_density.cpp
/// @brief The von Mises log-density and mean resultant length, built on the
///        SCALED Bessel functions corvus::i0e and corvus::i1e.
///
/// The von Mises distribution on the circle is
///
///     f(theta) = exp(kappa * cos(theta - mu)) / (2*pi*I0(kappa))
///
/// so everything hinges on I0, and I0 overflows a double at kappa ~ 713.99.
/// A concentration of 714 is not exotic — it is an angular standard deviation
/// of about two degrees — so any implementation resting on the unscaled I0
/// simply stops working on ordinary data.
///
/// The scaled forms exist for this. i0e(x) = exp(-|x|)*I0(x) and likewise for
/// i1e, and they stay finite across the whole double range (their minimum is
/// around 3e-155 at DBL_MAX) and never underflow. Two compositions follow:
///
///     log I0(kappa) = log(i0e(kappa)) + kappa
///     A(kappa)      = I1/I0 = i1e/i0e          <- the scalings cancel EXACTLY
///
/// The second is worth pausing on. exp(-kappa) divides out of numerator and
/// denominator identically, so the ratio is not an approximation of I1/I0 that
/// happens to avoid overflow — it is the same number, obtained from operands
/// that stay in range. A(kappa) is the mean resultant length, which is what
/// the MLE for kappa is solved against, so this is the composition a fitting
/// routine actually needs.
///
/// The first carries a caveat worth stating: relative error under 1 ULP for
/// kappa above about 2, and absolute error at most 3.3e-16 everywhere. Near
/// kappa = 0, log I0 approaches zero while i0e does not, so the sum loses the
/// relative sense — sufficient for log-density work, and the docs say so
/// rather than claiming a relative bound that does not hold.

#include <corvus/corvus.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <span>
#include <vector>

namespace {

constexpr double kTwoPi = 6.283185307179586;

/// @brief log I0(kappa) for a batch, valid past the point I0 itself overflows.
void log_i0(std::span<const double> kappa, std::span<double> out) {
    corvus::i0e(kappa, out);
    for (std::size_t i = 0; i < kappa.size(); ++i) {
        out[i] = std::log(out[i]) + kappa[i];
    }
}

/// @brief von Mises log-density at each (theta, mu, kappa).
void von_mises_logpdf(std::span<const double> theta, double mu, std::span<const double> kappa,
                      std::span<double> out) {
    std::vector<double> log_norm(kappa.size());
    log_i0(kappa, log_norm);
    for (std::size_t i = 0; i < theta.size(); ++i) {
        out[i] = kappa[i] * std::cos(theta[i] - mu) - std::log(kTwoPi) - log_norm[i];
    }
}

/// @brief A(kappa) = I1(kappa)/I0(kappa), the mean resultant length.
void mean_resultant(std::span<const double> kappa, std::span<double> out) {
    std::vector<double> num(kappa.size());
    corvus::i1e(kappa, num);
    corvus::i0e(kappa, out);
    for (std::size_t i = 0; i < kappa.size(); ++i) {
        out[i] = num[i] / out[i];
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
    std::cout << "corvus — von Mises density and concentration\n";
    std::cout << "SIMD target: " << corvus::active_target() << "\n\n";

    bool ok = true;

    // -----------------------------------------------------------------------
    // 1. Where the unscaled form stops, and the scaled form does not.
    // -----------------------------------------------------------------------
    constexpr std::array<double, 5> kK{700.0, 713.0, 714.0, 800.0, 1.0e6};
    std::array<double, kK.size()> plain{};
    std::array<double, kK.size()> scaled{};
    corvus::i0(kK, plain);
    corvus::i0e(kK, scaled);

    std::cout << "  kappa            I0(kappa)              i0e(kappa)\n";
    for (std::size_t i = 0; i < kK.size(); ++i) {
        std::cout << std::scientific << std::setprecision(3) << std::setw(11) << kK[i] << "  "
                  << std::setw(22) << plain[i] << "  " << std::setprecision(10) << std::setw(18)
                  << scaled[i] << "\n";
    }
    std::cout << "\n  I0 saturates to +inf at the true overflow boundary (|x| ~ 713.99),\n"
                 "  exactly, rather than drifting into wrong finite values first. i0e is\n"
                 "  still returning ordinary numbers at kappa = 1e6.\n";

    if (!std::isinf(plain[3]) || !std::isinf(plain[4])) {
        std::cerr << "expected I0 to overflow past the boundary\n";
        ok = false;
    }
    if (!std::isfinite(scaled[4]) || !(scaled[4] > 0.0)) {
        std::cerr << "expected i0e to stay finite and positive at kappa = 1e6\n";
        ok = false;
    }

    // -----------------------------------------------------------------------
    // 2. log I0 through the composition, against independent references.
    // -----------------------------------------------------------------------
    constexpr std::array<double, 7> kLogK{0.5, 2.0, 10.0, 100.0, 700.0, 5000.0, 1.0e6};
    // mpmath, 40 digits.
    constexpr std::array<double, 7> kRefLogI0{0.061549719185481307, 0.82399354148295634,
                                              7.9429720831186952,   96.779732689942577,
                                              695.8056999984434,    4994.822489873588,
                                              999992.17330631276};
    std::array<double, kLogK.size()> got{};
    log_i0(kLogK, got);

    std::cout << "\n  log I0(kappa) = log(i0e(kappa)) + kappa\n\n";
    std::cout << "      kappa            computed             reference     rel. error\n";
    for (std::size_t i = 0; i < kLogK.size(); ++i) {
        const double rel = std::abs(got[i] - kRefLogI0[i]) / std::abs(kRefLogI0[i]);
        std::cout << std::scientific << std::setprecision(3) << std::setw(11) << kLogK[i] << "  "
                  << std::fixed << std::setprecision(10) << std::setw(18) << got[i] << "  "
                  << std::setw(18) << kRefLogI0[i] << "  " << std::scientific
                  << std::setprecision(2) << rel << "\n";
        if (!close_enough(got[i], kRefLogI0[i], 1e-14)) {
            std::cerr << "log I0 wrong at kappa = " << kLogK[i] << "\n";
            ok = false;
        }
    }
    std::cout << "\n  The last three rows are past the point I0 itself exists as a double.\n";

    // -----------------------------------------------------------------------
    // 3. A(kappa) = i1e/i0e — the scalings cancel exactly.
    // -----------------------------------------------------------------------
    constexpr std::array<double, 7> kRefA{0.24249961258080194, 0.69777465796400795,
                                          0.94859982595484593, 0.99498737300516882,
                                          0.99928545881842612, 0.99989999499899973,
                                          0.99999949999987503};
    std::array<double, kLogK.size()> a_scaled{};
    mean_resultant(kLogK, a_scaled);

    // The same ratio from the UNSCALED pair, for the range where it survives.
    std::array<double, kLogK.size()> i0_plain{};
    std::array<double, kLogK.size()> i1_plain{};
    corvus::i0(kLogK, i0_plain);
    corvus::i1(kLogK, i1_plain);

    std::cout << "\n  A(kappa) = I1/I0, from the scaled pair and the unscaled pair\n\n";
    std::cout << "      kappa       i1e/i0e            i1/i0           rel. error vs ref\n";
    for (std::size_t i = 0; i < kLogK.size(); ++i) {
        const double unscaled = i1_plain[i] / i0_plain[i];
        const double rel = std::abs(a_scaled[i] - kRefA[i]) / kRefA[i];
        std::cout << std::scientific << std::setprecision(3) << std::setw(11) << kLogK[i] << "  "
                  << std::fixed << std::setprecision(15) << std::setw(17) << a_scaled[i] << "  "
                  << std::setw(17) << unscaled << "  " << std::scientific << std::setprecision(2)
                  << rel << "\n";
        if (!close_enough(a_scaled[i], kRefA[i], 1e-15)) {
            std::cerr << "A(kappa) wrong at kappa = " << kLogK[i] << "\n";
            ok = false;
        }
    }
    std::cout << "\n  Past kappa ~ 714 the middle column is inf/inf — a NaN where the\n"
                 "  scaled route is still exact. Below that the two agree, because they\n"
                 "  are the same quantity: exp(-kappa) divides out identically.\n";

    if (!std::isnan(i1_plain[6] / i0_plain[6])) {
        std::cerr << "expected the unscaled ratio to go NaN at kappa = 1e6\n";
        ok = false;
    }

    // -----------------------------------------------------------------------
    // 4. A log-density, and the limitation worth knowing about.
    // -----------------------------------------------------------------------
    constexpr std::array<double, 4> kTheta{0.0, 0.1, 0.5, 3.0};
    std::array<double, kTheta.size()> conc{};
    conc.fill(800.0);  // past where I0 exists
    std::array<double, kTheta.size()> lpdf{};
    von_mises_logpdf(kTheta, 0.0, conc, lpdf);

    std::cout << "\n  log-density at kappa = 800, mu = 0 (I0 unavailable at this kappa)\n\n";
    for (std::size_t i = 0; i < kTheta.size(); ++i) {
        std::cout << std::fixed << std::setprecision(1) << "   theta = " << std::setw(5)
                  << kTheta[i] << "   log f = " << std::setprecision(10) << std::setw(18)
                  << lpdf[i] << "\n";
        if (!std::isfinite(lpdf[i])) {
            std::cerr << "log-density not finite at theta = " << kTheta[i] << "\n";
            ok = false;
        }
    }

    // The circular variance is 1 - A(kappa), and it is the one thing here that
    // the scaled forms do NOT rescue. Worth saying plainly in an example that
    // has otherwise been showing them rescue everything.
    std::cout << "\n  One limitation, stated rather than skipped. The circular variance is\n"
                 "  1 - A(kappa), and A approaches 1 as kappa grows — so forming that\n"
                 "  complement in double sheds about log2(2*kappa) bits NO MATTER HOW\n"
                 "  ACCURATE A IS. At kappa = 1e6 that is around 21 bits gone, and the\n"
                 "  exactness of i1e/i0e does not help, because the loss happens after\n"
                 "  corvus has returned. It is the same rule as everywhere else in these\n"
                 "  examples: a subtraction of near-equal quantities sets the accuracy,\n"
                 "  not the function feeding it.\n\n";

    std::array<double, 3> big{{1.0e3, 1.0e5, 1.0e6}};
    std::array<double, 3> abig{};
    mean_resultant(big, abig);
    std::cout << "      kappa      1 - A(kappa)     bits left of a 53-bit double\n";
    for (std::size_t i = 0; i < big.size(); ++i) {
        const double comp = 1.0 - abig[i];
        const double bits_lost = std::log2(2.0 * big[i]);
        std::cout << std::scientific << std::setprecision(3) << std::setw(11) << big[i] << "  "
                  << std::setw(14) << comp << "         " << std::fixed << std::setprecision(0)
                  << std::setw(3) << (53.0 - bits_lost) << "\n";
    }

    std::cout << "\n" << (ok ? "All checks passed.\n" : "FAILURES — see above.\n");
    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
