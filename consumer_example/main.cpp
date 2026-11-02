/// @file main.cpp
/// @brief Minimal consumer example for corvus via find_package.
///
/// Verifies that an installed corvus can be found, linked, and used
/// by an external project.

#include <corvus/corvus.h>

#include <array>
#include <cmath>
#include <iostream>

int main() {
    const std::array<double, 5> in{-2.0, -1.0, 0.0, 1.0, 2.0};
    std::array<double, 5> out{};

    corvus::erf(in, out);

    std::cout << "corvus consumer example\n";
    std::cout << "========================\n";
    std::cout << "SIMD target: " << corvus::active_target() << "\n\n";
    for (std::size_t i = 0; i < in.size(); ++i) {
        std::cout << "  erf(" << in[i] << ") = " << out[i] << "\n";
    }

    // Verify known values bit-exactly. erf(1) = 0.84270079294971486934...
    // lies 0.22 ulp from its correctly rounded double, so every tier must
    // return exactly that double; erf is odd and erf(+0) = +0.
    const double erf1 = 0x1.af767a741088bp-1;
    if (out[3] != erf1 || out[1] != -erf1 || out[2] != 0.0 || std::signbit(out[2])) {
        std::cerr.precision(17);
        std::cerr << "Verification failed: erf(1) = " << out[3] << ", erf(-1) = " << out[1]
                  << ", erf(0) = " << out[2] << "; expected " << erf1 << ", " << -erf1
                  << ", +0\n";
        return 1;
    }

    std::cout << "\nVerification passed.\n";
    return 0;
}
