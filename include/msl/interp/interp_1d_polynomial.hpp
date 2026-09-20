/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: interp_1d_polynomial.hpp
** -----
** File Created: Wednesday, 15th October 2025 16:51:18
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Sunday, 14th December 2025 17:03:49
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_INTERP_1D_POLYNOMIAL_HPP
#define MSL_INTERP_1D_POLYNOMIAL_HPP

#include <span>
#include <vector>

#include "interp_1d_base.hpp"

namespace msl::interp {

/**
 * @brief Polynomial interpolation using divided differences (Newton form)
 *
 * More stable than Lagrange form, allows incremental updates.
 * Warning: High degree polynomials can oscillate (Runge phenomenon)
 *
 * Recommended: Use for <= 10 points, or consider splines instead
 */
class Polynomial : public InterpolatorBase {
public:
    Polynomial() = default;

    /**
     * @brief Create polynomial interpolator from data points.
     *
     * @param x Independent variable samples (strictly increasing)
     * @param y Dependent variable samples
     * @return Polynomial interpolator
     */
    [[nodiscard]] static Polynomial from_data(std::span<const double> x,
                                              std::span<const double> y) {
        Polynomial interp;
        interp.set_data(x, y);
        return interp;
    }

    /**
     * @brief Set interpolation data and recompute Newton coefficients.
     *
     * @param x Independent variable samples (strictly increasing)
     * @param y Dependent variable samples
     */
    void set_data(std::span<const double> x,
                  std::span<const double> y) override {
        x_.assign(x.begin(), x.end());
        y_.assign(y.begin(), y.end());
        validate_input();
        compute_divided_differences();
    }

    /**
     * @brief Interpolate at a single point.
     *
     * @param x Query point
     * @return Interpolated value at @p x
     */
    double interpolate(double x) const override {
        // Newton's form: P(x) = c0 + c1(x-x0) + c2(x-x0)(x-x1) + ...
        size_t n = x_.size();
        double result = coeffs_[n - 1];

        // Horner's method (reverse order for stability)
        for (int i = n - 2; i >= 0; --i) {
            result = result * (x - x_[i]) + coeffs_[i];
        }

        return result;
    }

    /**
     * @brief Get polynomial degree
     */
    [[nodiscard]] size_t degree() const { return x_.size() - 1; }

    /**
     * @brief Evaluate derivative at point
     *
     * Computes derivative together with polynomial evaluation using nested
     * multiplication in Newton form.
     */
    [[nodiscard]] double derivative(double x) const {
        size_t n = x_.size();
        if (n < 2) {
            return 0.0;
        }

        double p = coeffs_[n - 1];
        double dp = 0.0;

        for (size_t i = n - 1; i-- > 0;) {
            dp = dp * (x - x_[i]) + p;
            p = p * (x - x_[i]) + coeffs_[i];
        }

        return dp;
    }

private:
    std::vector<double> coeffs_; // Divided difference coefficients

    /**
     * @brief Build Newton divided-difference coefficients.
     */
    void compute_divided_differences() {
        size_t n = x_.size();
        coeffs_.resize(n);

        // In-place divided differences: coeffs_[k] becomes k-th Newton
        // coefficient after each order update.
        for (size_t i = 0; i < n; ++i) {
            coeffs_[i] = y_[i];
        }

        for (size_t order = 1; order < n; ++order) {
            for (size_t i = n - 1; i >= order; --i) {
                coeffs_[i] =
                    (coeffs_[i] - coeffs_[i - 1]) / (x_[i] - x_[i - order]);

                if (i == order) {
                    break;
                }
            }
        }
    }
};

/**
 * @brief Polynomial interpolation (zero-copy output).
 *
 * @param x Independent variable samples
 * @param y Dependent variable samples
 * @param x_new Query points
 * @param result Output buffer (must have same size as @p x_new)
 */
inline void interp1_polynomial(std::span<const double> x,
                               std::span<const double> y,
                               std::span<const double> x_new,
                               std::span<double> result) {
    if (x_new.size() != result.size()) {
        throw std::invalid_argument(
            "interp1_polynomial: x_new and result spans must have same size");
    }

    auto interp = Polynomial::from_data(x, y);
    interp(x_new, result);
}

/**
 * @brief Convenience function for polynomial interpolation
 */
inline std::vector<double> interp1_polynomial(std::span<const double> x,
                                              std::span<const double> y,
                                              std::span<const double> x_new) {
    return Polynomial::from_data(x, y)(x_new);
}

} // namespace msl::interp

#endif // MSL_INTERP_1D_POLYNOMIAL_HPP