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
#include <cmath>
#include <limits>
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

/**
 * Overlap-local Pearson correlation at lags [-max_lag, +max_lag].
 * Positive k compares x[n+k] with y[n]. Inputs must be equal-length,
 * finite, N >= 2, max_lag <= N-2. Zero overlap variance returns NaN.
 * Direct O(N * max_lag) evaluation; no FFT or interpolation.
 */
inline std::vector<double> xcorr_normalized(std::span<const double> x,
                                           std::span<const double> y,
                                           std::size_t max_lag) {
    if (x.size() != y.size() || x.size() < 2 || max_lag > x.size() - 2
        || max_lag > (std::numeric_limits<std::size_t>::max() - 1) / 2
        || x.size() > static_cast<std::size_t>(
                           std::numeric_limits<std::ptrdiff_t>::max()))
        throw std::invalid_argument("Normalized cross-correlation: invalid lengths/lag");
    for (std::size_t i = 0; i < x.size(); ++i)
        if (!std::isfinite(x[i]) || !std::isfinite(y[i]))
            throw std::invalid_argument("Normalized cross-correlation: finite inputs required");
    // Same-sign subtraction preserves small changes on a large DC offset.
    // Opposite-sign values are scaled first to avoid subtraction overflow.
    const auto delta = [](double value, double anchor, double scale) {
        if (scale == 0.) return 0.;
        return std::signbit(value) == std::signbit(anchor)
                   ? (value-anchor)/scale : value/scale-anchor/scale;
    };
    std::vector<double> result(2 * max_lag + 1,
                               std::numeric_limits<double>::quiet_NaN());
    const auto limit = static_cast<std::ptrdiff_t>(max_lag);
    for (auto k = -limit; k <= limit; ++k) {
        const auto ax = k > 0 ? static_cast<std::size_t>(k) : 0;
        const auto by = k < 0 ? static_cast<std::size_t>(-k) : 0;
        const auto count = x.size() - static_cast<std::size_t>(std::abs(k));
        double sx = 0., sy = 0.;
        for (std::size_t n = 0; n < count; ++n) {
            sx = std::max(sx, std::abs(x[ax+n]));
            sy = std::max(sy, std::abs(y[by+n]));
        }
        // Scale each overlap separately: extreme samples outside it must
        // not erase its variation. Anchor each overlap at its first value. Constant overlaps stay
        // exactly zero even when summing a non-binary-representable value.
        double mx = 0., my = 0.;
        for (std::size_t n = 0; n < count; ++n) {
            mx += delta(x[ax+n], x[ax], sx);
            my += delta(y[by+n], y[by], sy);
        }
        mx /= static_cast<double>(count);
        my /= static_cast<double>(count);
        double xx = 0., yy = 0., xy = 0.;
        for (std::size_t n = 0; n < count; ++n) {
            const double dx = delta(x[ax+n], x[ax], sx) - mx;
            const double dy = delta(y[by+n], y[by], sy) - my;
            xx += dx*dx;
            yy += dy*dy;
            xy += dx*dy;
        }
        if (xx > 0. && yy > 0.)
            result[static_cast<std::size_t>(k+limit)] =
                std::clamp((xy / std::sqrt(xx)) / std::sqrt(yy), -1., 1.);
    }
    return result;
}

} // namespace msl::signal

#endif // MSL_CROSS_COVARIANCE_HPP
