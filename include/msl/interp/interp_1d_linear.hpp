/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: interp_1d_linear.hpp
** -----
** File Created: Wednesday, 15th October 2025 08:57:54
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Sunday, 14th December 2025 17:03:44
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_INTERP_1D_LINEAR
#define MSL_INTERP_1D_LINEAR

#include "interp_1d_base.hpp"

namespace msl::interp {
/**
 * @brief Piecewise linear interpolation.
 *
 * Uses the two neighboring samples around the query point to perform
 * first-order interpolation.
 */
class Linear : public InterpolatorBase {
public:
    Linear() = default;

    /**
     * @brief Create linear interpolator from data points.
     *
     * @param x Independent variable samples (strictly increasing)
     * @param y Dependent variable samples
     * @return Linear interpolator
     */
    [[nodiscard]] static Linear
    from_data(std::span<const double> x, std::span<const double> y) {
        Linear interp;
        interp.set_data(x, y);
        return interp;
    }

    /**
     * @brief Set interpolation data.
     *
     * @param x Independent variable samples (strictly increasing)
     * @param y Dependent variable samples
     */
    void set_data(std::span<const double> x,
                  std::span<const double> y) override {
        x_.assign(x.begin(), x.end());
        y_.assign(y.begin(), y.end());
        validate_input();
    }

    /**
     * @brief Interpolate at one point.
     *
     * @param x Query point
     * @return Interpolated value at @p x
     */
    double interpolate(double x) const override {
        size_t i = find_interval(x);

        // Linear interpolation: y = y0 + (y1-y0)/(x1-x0) * (x-x0)
        double t = (x - x_[i]) / (x_[i + 1] - x_[i]);
        return y_[i] + t * (y_[i + 1] - y_[i]);
    }

    /**
     * @brief Evaluate first derivative (segment slope) at one point.
     *
     * @param x Query point
     * @return Slope of the active linear segment
     */
    [[nodiscard]] double derivative(double x) const {
        size_t i = find_interval(x);
        return (y_[i + 1] - y_[i]) / (x_[i + 1] - x_[i]);
    }
};

/**
 * @brief Piecewise linear interpolation (zero-copy output).
 *
 * @param x Independent variable samples
 * @param y Dependent variable samples
 * @param x_new Query points
 * @param result Output buffer (must have same size as @p x_new)
 */
inline void interp1_linear(std::span<const double> x,
                           std::span<const double> y,
                           std::span<const double> x_new,
                           std::span<double> result) {
    if (x_new.size() != result.size()) {
        throw std::invalid_argument(
            "interp1_linear: x_new and result spans must have same size");
    }

    auto interp = Linear::from_data(x, y);
    interp(x_new, result);
}

/**
 * @brief Piecewise linear interpolation.
 *
 * @param x Independent variable samples
 * @param y Dependent variable samples
 * @param x_new Query points
 * @return Interpolated results at @p x_new
 */
inline std::vector<double> interp1_linear(std::span<const double> x,
                                          std::span<const double> y,
                                          std::span<const double> x_new) {
    return Linear::from_data(x, y)(x_new);
}

} // namespace msl::interp

#endif // MSL_INTERP_1D_LINEAR