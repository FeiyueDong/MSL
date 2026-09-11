/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: interp_1d_akima.hpp
** -----
** File Created: Wednesday, 15th October 2025 14:12:29
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Sunday, 14th December 2025 17:03:33
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_INTERP_1D_AKIMA
#define MSL_INTERP_1D_AKIMA

#include "interp_1d_base.hpp"

namespace msl::interp {
/**
 * @brief Akima spline interpolation (smoother than cubic, limited overshoot).
 *
 * Uses modified Akima weights by default, matching MATLAB `makima`.
 * Set `modified_akima = false` for the original Akima scheme.
 */
class AkimaSpline : public InterpolatorBase {
public:
    AkimaSpline() = default;

    /**
     * @brief Create Akima spline interpolator from data points.
     *
     * @param x Independent variable samples (strictly increasing)
     * @param y Dependent variable samples
     * @param modified_akima Use modified Akima weights (MATLAB makima)
     * @return AkimaSpline interpolator
     */
    [[nodiscard]] static AkimaSpline from_data(std::span<const double> x,
                                               std::span<const double> y,
                                               bool modified_akima = true) {
        AkimaSpline spline;
        spline.modified_akima_ = modified_akima;
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
        if (x_.size() < 4) {
            throw std::invalid_argument(
                "interp1_akima: need at least 4 points");
        }
        compute_coefficients();
    }

    void set_modified_akima(bool modified) {
        if (modified != modified_akima_) {
            modified_akima_ = modified;
            compute_coefficients();
        }
    }

    /**
     * @brief Interpolate at a single point
     *
     * @param x Evaluation point
     * @return Interpolated value at x
     */
    double interpolate(double x) const override {
        size_t i = find_interval(x);
        double dx = x - x_[i];
        return y_[i] + b_[i] * dx + c_[i] * dx * dx + d_[i] * dx * dx * dx;
    }

private:
    std::vector<double> b_; // First derivative coefficients
    std::vector<double> c_; // Second derivative coefficients
    std::vector<double> d_; // Third derivative coefficients

    bool modified_akima_ =
        true; // Use modified Akima weights to match MATLAB makima

    /**
     * @brief Compute polynomial coefficients
     */
    void compute_coefficients() {
        size_t n = x_.size();

        // Compute slopes
        std::vector<double> m(n + 3);

        // Interior slopes
        for (size_t i = 2; i < n + 1; ++i) {
            m[i] = (y_[i - 1] - y_[i - 2]) / (x_[i - 1] - x_[i - 2]);
        }

        // Extrapolate for boundary
        m[0] = 3.0 * m[2] - 2.0 * m[3];
        m[1] = 2.0 * m[2] - m[3];
        m[n + 1] = 2.0 * m[n] - m[n - 1];
        m[n + 2] = 3.0 * m[n] - 2.0 * m[n - 1];

        // Compute Akima weights
        std::vector<double> t(n);
        for (size_t i = 0; i < n; ++i) {
            if (modified_akima_) {
                // Modified Akima weights: |d_{i+1} - d_i| + |d_{i+1} + d_i|/2
                double w1 = std::abs(m[i + 3] - m[i + 2])
                            + std::abs(m[i + 3] + m[i + 2]) / 2.0;
                double w2 =
                    std::abs(m[i + 1] - m[i]) + std::abs(m[i + 1] + m[i]) / 2.0;

                if (w1 + w2 < 1e-10) {
                    t[i] = 0.5 * (m[i + 1] + m[i + 2]);
                } else {
                    t[i] = (w1 * m[i + 1] + w2 * m[i + 2]) / (w1 + w2);
                }
            } else {
                // Original Akima weights: |d_{i+1} - d_i|
                double w1 = std::abs(m[i + 3] - m[i + 2]);
                double w2 = std::abs(m[i + 1] - m[i]);

                if (w1 + w2 < 1e-10) {
                    t[i] = 0.5 * (m[i + 1] + m[i + 2]);
                } else {
                    t[i] = (w1 * m[i + 1] + w2 * m[i + 2]) / (w1 + w2);
                }
            }
        }

        // Compute polynomial coefficients
        b_.resize(n - 1);
        c_.resize(n - 1);
        d_.resize(n - 1);

        for (size_t i = 0; i < n - 1; ++i) {
            double h = x_[i + 1] - x_[i];
            b_[i] = t[i];
            c_[i] = (3.0 * m[i + 2] - 2.0 * t[i] - t[i + 1]) / h;
            d_[i] = (t[i] + t[i + 1] - 2.0 * m[i + 2]) / (h * h);
        }
    }
};

/**
 * @brief Akima spline interpolation(Zero-copy version)
 *
 * @param x Independent variable values
 * @param y Dependent variable values (must have same size as x)
 * @param x_new Evaluation points
 * @param result Output buffer for interpolated values (must have same size as
 * x_new)
 */
inline void interp1_akima(std::span<const double> x,
                          std::span<const double> y,
                          std::span<const double> x_new,
                          std::span<double> result) {
    if (x_new.size() != result.size()) {
        throw std::invalid_argument(
            "interp1_akima: x_new and result spans must have same size");
    }

    auto spline = AkimaSpline::from_data(x, y);
    spline(x_new, result);
}

inline std::vector<double> interp1_akima(std::span<const double> x,
                                         std::span<const double> y,
                                         std::span<const double> x_new) {
    return AkimaSpline::from_data(x, y)(x_new);
}

} // namespace msl::interp
#endif // MSL_INTERP_1D_AKIMA