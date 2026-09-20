/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2026, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: filter.hpp
** -----
** File Created: Friday, 9th January 2026 14:58:10
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Friday, 6th March 2026 09:44:27
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_FILTER_HPP
#define MSL_FILTER_HPP

#include <algorithm>
#include <span>
#include <stdexcept>
#include <vector>

#include "filter_design.hpp"
#include "matrix/real_matrix_owned.hpp"

namespace msl::signal {

// ============================================================================
// Direct Form II Transposed (most numerically stable)
// ============================================================================

/**
 * @brief Apply IIR filter to signal using Direct Form II Transposed
 *
 * This is the most numerically stable form for IIR filters.
 *
 * @param signal Input signal
 * @param result Output signal (must be same size as input)
 * @param coeffs Filter coefficients (a[0] must be 1.0)
 */
inline void filter(std::span<const double> signal,
                   std::span<double> result,
                   const FilterCoefficients &coeffs) {
    if (signal.size() != result.size()) {
        throw std::invalid_argument(
            "Filter: input and output spans must have the same size");
    }

    if (coeffs.a.empty() || coeffs.b.empty()) {
        throw std::invalid_argument("Filter: coefficients cannot be empty");
    }

    if (std::abs(coeffs.a[0] - 1.0) > 1e-10) {
        throw std::invalid_argument(
            "Filter: first denominator coefficient must be 1.0");
    }

    const size_t n = signal.size();
    const size_t na = coeffs.a.size();
    const size_t nb = coeffs.b.size();
    const size_t nz = std::max(na, nb) - 1; // Number of delay states

    std::vector<double> z(nz, 0.0); // State vector (delay line)

    // Direct Form II Transposed implementation
    for (size_t i = 0; i < n; ++i) {
        double x = signal[i];

        // Output = b[0]*x + z[0]
        result[i] = (nb > 0 ? coeffs.b[0] : 0.0) * x + (nz > 0 ? z[0] : 0.0);

        // Update state vector
        for (size_t j = 0; j < nz - 1; ++j) {
            double bj = (j + 1 < nb) ? coeffs.b[j + 1] : 0.0;
            double aj = (j + 1 < na) ? coeffs.a[j + 1] : 0.0;
            z[j] = bj * x - aj * result[i] + z[j + 1];
        }

        // Last state
        if (nz > 0) {
            double bn = (nz < nb) ? coeffs.b[nz] : 0.0;
            double an = (nz < na) ? coeffs.a[nz] : 0.0;
            z[nz - 1] = bn * x - an * result[i];
        }
    }
}

/**
 * @brief Apply IIR filter to signal using Direct Form II Transposed
 *
 * This is the most numerically stable form for IIR filters.
 *
 * @param signal Input signal
 * @param coeffs Filter coefficients (a[0] must be 1.0)
 * @return Filtered signal
 *
 * @note Uses zero initial conditions
 */
inline std::vector<double> filter(std::span<const double> signal,
                                  const FilterCoefficients &coeffs) {

    std::vector<double> output(signal.size());
    filter(signal, output, coeffs);
    return output;
}

/**
 * @brief Apply filter to each column of a matrix
 *
 * Convenience function to apply filter to each column of a matrix
 * Equivalent to applying filter to each column separately
 *
 * @param signals Matrix where each column is a signal
 * @param coeffs Filter coefficients
 * @return Filtered matrix
 */
inline matrix::matrixd filter_columns(const matrix::real_matrix_base &signals,
                                      const FilterCoefficients &coeffs) {
    matrix::matrixd output(signals.rows(), signals.cols());

    for (size_t j = 0; j < signals.cols(); ++j) {
        auto col_span = signals.column(j);
        auto filtered = filter(col_span, coeffs);

        for (size_t i = 0; i < signals.rows(); ++i) {
            output(i, j) = filtered[i];
        }
    }

    return output;
}

} // namespace msl::signal

#endif // MSL_FILTER_APPLY_HPP