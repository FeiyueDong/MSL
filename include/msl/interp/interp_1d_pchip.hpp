/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2026, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: interp_1d_pchip.hpp
** -----
** File Created: Sunday, 19th April 2026 20:46:34
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Sunday, 19th April 2026 21:01:49
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_INTERP_1D_PCHIP_HPP
#define MSL_INTERP_1D_PCHIP_HPP

#include <cmath>
#include <stdexcept>
#include <vector>

#include "interp_1d_base.hpp"

namespace msl::interp {

/**
 * @brief Piecewise Cubic Hermite Interpolating Polynomial (PCHIP)
 * * This implementation strictly follows MATLAB's 'pchip' algorithm,
 * ensuring shape-preserving (monotonic) behavior.
 */
class PchipSpline : public InterpolatorBase {
public:
    PchipSpline() = default;

    /**
     * @brief Create PCHIP interpolator from data points.
     *
     * @param x Independent variable samples (strictly increasing)
     * @param y Dependent variable samples
     * @param extrap_mode Extrapolation mode
     * @return PchipSpline interpolator
     */
    [[nodiscard]] static PchipSpline
    from_data(std::span<const double> x,
              std::span<const double> y,
              ExtrapolationMode extrap_mode = ExtrapolationMode::Polynomial) {
        PchipSpline spline;
        spline.extrap_mode_ = extrap_mode;
        spline.set_data(x, y);
        return spline;
    }

    /**
     * @brief Set data points
     *
     * @param x Independent variable values
     * @param y Dependent variable values (must have same size as x)
     */
    void set_data(std::span<const double> x,
                  std::span<const double> y) override {
        x_.assign(x.begin(), x.end());
        y_.assign(y.begin(), y.end());
        validate_input();
        if (x_.size() < 2) {
            throw std::runtime_error(
                "At least 2 points are required for PCHIP interpolation.");
        }
        compute_coefficients();
    }

    /**
     * @brief Interpolate at a single point
     *
     * @param x Query point
     * @return Interpolated value at @p x
     */
    double interpolate(double x) const override {
        size_t i = find_interval(x);
        double dx = x - x_[i];
        return y_[i] + b_[i] * dx + c_[i] * dx * dx + d_[i] * dx * dx * dx;
    }

private:
    std::vector<double> b_; // 1st order coefficients
    std::vector<double> c_; // 2nd order coefficients
    std::vector<double> d_; // 3rd order coefficients

    void compute_coefficients() {
        size_t n = x_.size();
        if (n < 2) {
            throw std::runtime_error(
                "At least 2 points are required for PCHIP interpolation.");
        }

        // Calculate step sizes h[i] and divided differences delta[i]
        std::vector<double> h(n - 1);
        std::vector<double> delta(n - 1);
        for (size_t i = 0; i < n - 1; ++i) {
            h[i] = x_[i + 1] - x_[i];
            delta[i] = (y_[i + 1] - y_[i]) / h[i];
        }

        // Array to store the derivatives (slopes) at each point
        std::vector<double> m(n, 0.0);

        if (n == 2) {
            // Degenerates to linear interpolation for 2 points
            m[0] = delta[0];
            m[1] = delta[0];
        } else {
            // Compute internal slopes (MATLAB's weighted harmonic mean)
            for (size_t k = 1; k < n - 1; ++k) {
                // If the signs of adjacent divided differences are opposite,
                // it's a local extremum. Slope must be 0 to prevent overshoot.
                if (delta[k - 1] * delta[k] > 0.0) {
                    double w1 = 2.0 * h[k] + h[k - 1];
                    double w2 = h[k] + 2.0 * h[k - 1];
                    m[k] = (w1 + w2) / (w1 / delta[k - 1] + w2 / delta[k]);
                } else {
                    m[k] = 0.0;
                }
            }

            // Compute boundary slopes (MATLAB's specific shape-preserving
            // formulas) Left boundary (k = 0)
            m[0] = ((2.0 * h[0] + h[1]) * delta[0] - h[0] * delta[1])
                   / (h[0] + h[1]);
            if (m[0] * delta[0] <= 0.0) {
                m[0] = 0.0;
            } else if (delta[0] * delta[1] <= 0.0
                       && std::abs(m[0]) > std::abs(3.0 * delta[0])) {
                m[0] = 3.0 * delta[0];
            }

            // Right boundary (k = n - 1)
            m[n - 1] = ((2.0 * h[n - 2] + h[n - 3]) * delta[n - 2]
                        - h[n - 2] * delta[n - 3])
                       / (h[n - 2] + h[n - 3]);
            if (m[n - 1] * delta[n - 2] <= 0.0) {
                m[n - 1] = 0.0;
            } else if (delta[n - 2] * delta[n - 3] <= 0.0
                       && std::abs(m[n - 1]) > std::abs(3.0 * delta[n - 2])) {
                m[n - 1] = 3.0 * delta[n - 2];
            }
        }

        // Compute final polynomial coefficients (b, c, d)
        // a_i is implicitly y_[i]
        b_.resize(n - 1);
        c_.resize(n - 1);
        d_.resize(n - 1);

        for (size_t i = 0; i < n - 1; ++i) {
            b_[i] = m[i];
            c_[i] = (3.0 * delta[i] - 2.0 * m[i] - m[i + 1]) / h[i];
            d_[i] = (m[i] + m[i + 1] - 2.0 * delta[i]) / (h[i] * h[i]);
        }
    }
};

/**
 * @brief PCHIP interpolation (zero-copy output).
 *
 * @param x Independent variable samples
 * @param y Dependent variable samples
 * @param x_new Query points
 * @param result Output buffer (must have same size as @p x_new)
 */
inline void interp1_pchip(std::span<const double> x,
                          std::span<const double> y,
                          std::span<const double> x_new,
                          std::span<double> result) {
    if (x_new.size() != result.size()) {
        throw std::invalid_argument(
            "interp1_pchip: x_new and result spans must have same size");
    }
    auto pchip = PchipSpline::from_data(x, y);
    pchip(x_new, result);
}

inline std::vector<double> interp1_pchip(std::span<const double> x,
                                         std::span<const double> y,
                                         std::span<const double> x_new) {
    auto interp = PchipSpline::from_data(x, y);
    return interp(x_new);
}

} // namespace msl::interp

#endif // MSL_INTERP_1D_PCHIP_HPP