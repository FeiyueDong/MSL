/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2026, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: cross_covariance.hpp
** -----
** Author: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_CROSS_COVARIANCE_HPP
#define MSL_CROSS_COVARIANCE_HPP

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <numeric>
#include <span>
#include <stdexcept>
#include <vector>

namespace msl::signal {

/**
 * @brief Compute an unbiased cross-covariance sequence for real signals.
 *
 * The result contains lags `-max_lag, ..., 0, ..., max_lag`. For lag `k`,
 * the convention is `sum(x[n + k] * y[n]) / (N - abs(k))`, matching MATLAB
 * `xcov(x, y, max_lag, "unbiased")`. Both input means are removed before the
 * lag products are formed.
 */
inline std::vector<double> xcov_unbiased(std::span<const double> x,
                                         std::span<const double> y,
                                         size_t max_lag) {
    if (x.size() != y.size()) {
        throw std::invalid_argument(
            "Unbiased cross-covariance: signals must have the same length");
    }
    if (x.empty()) {
        throw std::invalid_argument(
            "Unbiased cross-covariance: signals must not be empty");
    }
    if (max_lag >= x.size()) {
        throw std::invalid_argument(
            "Unbiased cross-covariance: max_lag must be less than signal "
            "length");
    }

    const double x_mean = std::accumulate(x.begin(), x.end(), 0.0) / x.size();
    const double y_mean = std::accumulate(y.begin(), y.end(), 0.0) / y.size();
    std::vector<double> result(2 * max_lag + 1, 0.0);

    for (std::ptrdiff_t lag = -static_cast<std::ptrdiff_t>(max_lag);
         lag <= static_cast<std::ptrdiff_t>(max_lag);
         ++lag) {
        const size_t x_start = lag > 0 ? static_cast<size_t>(lag) : 0;
        const size_t y_start = lag < 0 ? static_cast<size_t>(-lag) : 0;
        const size_t count = x.size() - static_cast<size_t>(std::abs(lag));
        double sum = 0.0;
        for (size_t i = 0; i < count; ++i) {
            sum += (x[x_start + i] - x_mean) * (y[y_start + i] - y_mean);
        }
        result[static_cast<size_t>(lag
                                   + static_cast<std::ptrdiff_t>(max_lag))] =
            sum / static_cast<double>(count);
    }
    return result;
}

/**
 * @brief Compute all lags from `-(N-1)` through `N-1`.
 */
inline std::vector<double> xcov_unbiased(std::span<const double> x,
                                         std::span<const double> y) {
    if (x.empty()) {
        throw std::invalid_argument(
            "Unbiased cross-covariance: signals must not be empty");
    }
    return xcov_unbiased(x, y, x.size() - 1);
}

} // namespace msl::signal

#endif // MSL_CROSS_COVARIANCE_HPP
