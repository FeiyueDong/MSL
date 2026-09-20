/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2026, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: interp_1d_near.hpp
** -----
** File Created: Sunday, 19th April 2026 21:01:16
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Sunday, 19th April 2026 21:02:07
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_INTERP_1D_NEAR_HPP
#define MSL_INTERP_1D_NEAR_HPP

#include "interp_1d_base.hpp"

namespace msl::interp {
/**
 * @brief Interpolation using near neighbor (piecewise constant)
 */
class Near : public InterpolatorBase {
public:
    // Type of near neighbor selection
    enum class NearType {
        Previous, // Use previous neighbor (floor)
        Next,     // Use next neighbor (ceil)
        Nearest   // Use nearest neighbor (round)
    };

    Near() = default;

    /**
     * @brief Create nearest-neighbor interpolator from data points.
     *
     * @param x Independent variable samples (strictly increasing)
     * @param y Dependent variable samples
     * @param type Near neighbor selection type
     * @return Near interpolator
     */
    [[nodiscard]] static Near
    from_data(std::span<const double> x,
              std::span<const double> y,
              NearType type = NearType::Nearest) {
        Near interp;
        interp.set_data(x, y);
        interp.near_type_ = type;
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
    void set_near_type(NearType type) { near_type_ = type; }

    /**
     * @brief Interpolate at one point.
     *
     * @param x Query point
     * @return Interpolated value at @p x
     */
    double interpolate(double x) const override {
        size_t i = find_interval(x);
        switch (near_type_) {
            case NearType::Previous:
                return y_[i];
            case NearType::Next:
                return y_[i + 1];
            case NearType::Nearest:
                if (x - x_[i] < x_[i + 1] - x) {
                    return y_[i];
                } else {
                    return y_[i + 1];
                }
            default:
                return y_[i];
        }
    }

private:
    NearType near_type_ = NearType::Nearest;
};

/**
 * @brief Interpolation using near neighbor (piecewise constant) (zero-copy
 * output).
 *
 * @param x Independent variable samples
 * @param y Dependent variable samples
 * @param x_new Query points
 * @param result Output buffer (must have same size as @p x_new)
 */
inline void interp1_near(std::span<const double> x,
                         std::span<const double> y,
                         std::span<const double> x_new,
                         std::span<double> result) {
    if (x_new.size() != result.size()) {
        throw std::invalid_argument(
            "interp1_near: x_new and result spans must have same size");
    }

    auto interp = Near::from_data(x, y);
    interp(x_new, result);
}

/**
 * @brief Interpolation using near neighbor (piecewise constant).
 *
 * @param x Independent variable samples
 * @param y Dependent variable samples
 * @param x_new Query points
 * @return Interpolated results at @p x_new
 */
inline std::vector<double> interp1_near(std::span<const double> x,
                                        std::span<const double> y,
                                        std::span<const double> x_new) {
    return Near::from_data(x, y)(x_new);
}

} // namespace msl::interp

#endif // MSL_INTERP_1D_NEAR_HPP