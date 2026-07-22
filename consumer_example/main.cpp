/// @file main.cpp
/// @brief Minimal consumer example for corvus via find_package.
///
/// Verifies that an installed corvus can be found, linked, and used
/// by an external project.

#include "corvus/corvus.h"

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

    // Verify a known value: erf(0) = 0.
    if (std::abs(out[2] - 0.0) > 1e-12) {
        std::cerr << "Verification failed: erf(0) = " << out[2] << ", expected 0\n";
        return 1;
    }

    std::cout << "\nVerification passed.\n";
    return 0;
}
