/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2026, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: detrend.hpp
** -----
** File Created: Friday, 9th January 2026 14:58:10
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Thursday, 5th March 2026 22:52:30
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_DETREND_HPP
#define MSL_DETREND_HPP

#include <cstddef>
#include <span>
#include <stdexcept>
#include <vector>

#include "matrix/real_matrix_base.hpp"
#include "polynomial/polynomial.hpp"

namespace msl::signal {
/**
 * @brief Detrend data by removing polynomial trend of degree n
 *
 * @param data Input data to detrend
 * @param result Output buffer for detrended data (must have same size as input)
 * @param n Degree of polynomial trend to remove (n=0 removes mean, n=1 removes
 * linear trend, etc.)
 */
inline void
detrend(std::span<const double> data, std::span<double> result, std::size_t n) {
    if (data.size() != result.size()) {
        throw std::invalid_argument(
            "Detrend: data and result must have same size");
    }

    if (data.size() <= n) {
        throw std::invalid_argument("Detrend: data size must be greater than "
                                    "polynomial degree for detrending");
    }

    if (n == 0) {
        auto mean_val = std::accumulate(data.begin(), data.end(), 0.0)
                        / static_cast<double>(data.size());
        for (std::size_t i = 0; i < data.size(); ++i) {
            result[i] = data[i] - mean_val;
        }
    }

    // Generate x values
    std::vector<double> x(data.size());
    std::iota(x.begin(), x.end(), 0.0);

    // Fit polynomial to data
    auto poly = msl::polynomial::Polynomial::from_fit(x, data, n);
    // Subtract trend
    for (std::size_t i = 0; i < data.size(); ++i) {
        result[i] = data[i] - poly(x[i]);
    }
}

/**
 * @brief Detrend data by removing polynomial trend of degree n
 *
 * @param data Input data to detrend
 * @param n Degree of polynomial trend to remove (n=0 removes mean, n=1 removes
 * linear trend, etc.)
 * @return std::vector<double> Detrended data
 */
inline std::vector<double> detrend(std::span<const double> data,
                                   std::size_t n = 1) {
    std::vector<double> detrended(data.size());
    detrend(data, detrended, n);
    return detrended;
}

/**
 * @brief Detrend data by removing polynomial trend of degree n
 *
 * @param x X-values for the data
 * @param data Input data to detrend
 * @param n Degree of polynomial trend to remove (n=0 removes mean, n=1 removes
 * linear trend, etc.)
 * @return std::vector<double> Detrended data
 */
inline void detrend(std::span<const double> x,
                    std::span<const double> data,
                    std::span<double> result,
                    std::size_t n = 1) {
    if (x.size() != data.size() || x.size() != result.size()) {
        throw std::invalid_argument(
            "Detrend: x, data, and result must have same size");
    }

    if (data.size() <= n) {
        throw std::invalid_argument("Detrend: data size must be greater than "
                                    "polynomial degree for detrending");
    }

    if (n == 0) {
        auto mean_val = std::accumulate(data.begin(), data.end(), 0.0)
                        / static_cast<double>(data.size());
        for (std::size_t i = 0; i < data.size(); ++i) {
            result[i] = data[i] - mean_val;
        }
        return;
    }

    // Fit polynomial to data
    auto poly = msl::polynomial::Polynomial::from_fit(x, data, n);

    // Subtract trend
    for (std::size_t i = 0; i < data.size(); ++i) {
        result[i] = data[i] - poly(x[i]);
    }
}

/**
 * @brief Detrend data by removing polynomial trend of degree n
 *
 * @param x X-values for the data
 * @param data Input data to detrend
 * @param n Degree of polynomial trend to remove (n=0 removes mean, n=1 removes
 * linear trend, etc.)
 * @return std::vector<double> Detrended data
 */
inline std::vector<double> detrend(std::span<const double> x,
                                   std::span<const double> data,
                                   std::size_t n = 1) {
    std::vector<double> detrended(data.size());
    detrend(x, data, detrended, n);
    return detrended;
}

/**
 * @brief Detrend data by removing polynomial trend of degree n
 *
 * @param data Input data to detrend
 * @param n Degree of polynomial trend to remove (n=0 removes mean, n=1 removes
 * linear trend, etc.)
 * @return matrix::matrixd Detrended data
 */
inline matrix::matrixd detrend(const matrix::real_matrix_base &data,
                               std::size_t n = 1) {
    matrix::matrixd detrended(data.rows(), data.cols());
    for (std::size_t col = 0; col < data.cols(); ++col) {
        // Extract column data
        auto col_data = data.column(col);

        // Detrend column
        auto detrended_col = detrended.column(col);
        detrend(col_data, detrended_col, n);
    }
    return detrended;
}

} // namespace msl::signal

#endif // MSL_DETREND_HPP