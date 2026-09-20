/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2026, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: filter_design.hpp
** -----
** File Created: Friday, 9th January 2026 14:58:10
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Friday, 6th March 2026 09:48:34
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_FILTER_DESIGN_HPP
#define MSL_FILTER_DESIGN_HPP

#include <vector>

namespace msl::signal {

// ============================================================================
// Filter specifications
// ============================================================================

/**
 * @brief Filter type based on frequency response
 *
 * - lowpass: Pass frequencies below cutoff, attenuate above cutoff
 * - highpass: Pass frequencies above cutoff, attenuate below cutoff
 * - bandpass: Pass frequencies between low and high cutoff, attenuate outside
 * - bandstop: Attenuate frequencies between low and high cutoff, pass outside
 */
enum class FilterType {
    lowpass,  // Low-pass filter
    highpass, // High-pass filter
    bandpass, // Band-pass filter
    bandstop  // Band-stop (notch) filter
};

/**
 * @brief Filter design method
 *
 * - butterworth: Maximally flat passband, monotonic stopband
 * - chebyshev1: Equiripple passband, monotonic stopband (future)
 * - chebyshev2: Monotonic passband, equiripple stopband (future)
 * - elliptic: Equiripple in both bands (future)
 */
enum class FilterMethod {
    butterworth, // Maximally flat passband
    chebyshev1,  // Equiripple passband, monotonic stopband (future)
    chebyshev2,  // Monotonic passband, equiripple stopband (future)
    elliptic     // Equiripple in both bands (future)
};

/**
 * @brief Filter coefficients (numerator and denominator)
 *
 * Represents digital filter as transfer function:
 * H(z) = (b[0] + b[1]*z^-1 + ... + b[M]*z^-M) / (a[0] + a[1]*z^-1 + ... +
 * a[N]*z^-N)
 *
 * Where a[0] is typically normalized to 1.0
 */
struct FilterCoefficients {
    // Numerator coefficients (b[0] corresponds to z^0 term, b[1] to z^-1, etc.)
    std::vector<double> b{};
    // Denominator coefficients (a[0] = 1.0, a[1] corresponds to z^-1 term,etc.)
    std::vector<double> a{};

    FilterCoefficients() = default;
    FilterCoefficients(std::vector<double> num, std::vector<double> den)
        : b(std::move(num)), a(std::move(den)) {}

    // Get filter order (max of numerator and denominator order)
    [[nodiscard]] size_t order() const {
        return std::max(a.size(), b.size()) - 1;
    }
};

} // namespace msl::signal

#endif // MSL_FILTER_DESIGN_HPP