/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: interp_1d_cubic.hpp
** -----
** File Created: Wednesday, 15th October 2025 14:11:26
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Sunday, 14th December 2025 17:03:56
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_INTERP_1D_CUBIC
#define MSL_INTERP_1D_CUBIC

#include "interp_1d_base.hpp"

namespace msl::interp {
/**
 * @brief Cubic spline interpolation.
 *
 * The default boundary condition is Not-a-knot, which matches MATLAB
 * `interp1(..., 'spline')`. Natural and clamped boundary conditions are
 * available through `BoundaryCondition`.
 */
class CubicSpline : public InterpolatorBase {
public:
    // Boundary condition types
    enum class BoundaryCondition {
        Natural,  // Second derivative is zero at endpoints
        NotAKnot, // Third derivative is continuous at second and penultimate
                  // points
        Clamped   // First derivative specified at endpoints
    };

    CubicSpline() = default;

    /**
     * @brief Create cubic spline interpolator from data points.
     *
     * @param x Independent variable samples (strictly increasing)
     * @param y Dependent variable samples
     * @param bc_type Boundary condition type
     * @return CubicSpline interpolator
     */
    [[nodiscard]] static CubicSpline
    from_data(std::span<const double> x,
              std::span<const double> y,
              BoundaryCondition bc_type = BoundaryCondition::NotAKnot) {
        CubicSpline spline;
        spline.bc_type_ = bc_type;
        spline.set_data(x, y);
        return spline;
    }

    /**
     * @brief Set data points and precompute spline coefficients.
     *
     * @param x Independent variable samples (strictly increasing)
     * @param y Dependent variable samples
     */
    void set_data(std::span<const double> x,
                  std::span<const double> y) override {
        x_.assign(x.begin(), x.end());
        y_.assign(y.begin(), y.end());
        validate_input();
        compute_coefficients();
    }

    /**
     * @brief Set boundary condition type
     *
     * @param bc_type Boundary condition type to set
     * @param left_slope Left endpoint slope (only used for Clamped BC)
     * @param right_slope Right endpoint slope (only used for Clamped BC)
     */
    void set_boundary_condition(BoundaryCondition bc_type,
                                double left_slope = 0.0,
                                double right_slope = 0.0) {
        if (bc_type != bc_type_) {
            bc_type_ = bc_type;
            if (bc_type == BoundaryCondition::Clamped) {
                left_slope_ = left_slope;
                right_slope_ = right_slope;
            }
            compute_coefficients();
        }
    }

    /**
     * @brief Interpolate at a single point.
     *
     * @param x Query point
     * @return Interpolated value at @p x
     */
    double interpolate(double x) const override {
        size_t i = find_interval(x);

        // Cubic spline: S_i(x) = a_i + b_i*(x-x_i) + c_i*(x-x_i)^2 +
        // d_i*(x-x_i)^3
        double dx = x - x_[i];
        return y_[i] + b_[i] * dx + c_[i] * dx * dx + d_[i] * dx * dx * dx;
    }

    /**
     * @brief Evaluate first derivative at a point.
     *
     * @param x Query point
     * @return First derivative value
     */
    [[nodiscard]] double derivative(double x) const {
        size_t i = find_interval(x);
        double dx = x - x_[i];
        return b_[i] + 2.0 * c_[i] * dx + 3.0 * d_[i] * dx * dx;
    }

    /**
     * @brief Evaluate second derivative at a point.
     *
     * @param x Query point
     * @return Second derivative value
     */
    [[nodiscard]] double second_derivative(double x) const {
        size_t i = find_interval(x);
        double dx = x - x_[i];
        return 2.0 * c_[i] + 6.0 * d_[i] * dx;
    }

private:
    std::vector<double> b_; // First derivative coefficients
    std::vector<double> c_; // Second derivative coefficients
    std::vector<double> d_; // Third derivative coefficients

    BoundaryCondition bc_type_ = BoundaryCondition::NotAKnot;
    double left_slope_ = 0.0;  // For clamped BC
    double right_slope_ = 0.0; // For clamped BC

    void compute_coefficients() {
        size_t n = x_.size();
        if (n < 3) {
            throw std::runtime_error(
                "At least 3 points are required for cubic spline.");
        }

        // Calculate step sizes h[i] = x[i+1] - x[i]
        std::vector<double> h(n - 1);
        for (size_t i = 0; i < n - 1; ++i) {
            h[i] = x_[i + 1] - x_[i];
        }

        // Calculate standard right-hand side differences (alpha)
        std::vector<double> alpha(n, 0.0);
        for (size_t i = 1; i < n - 1; ++i) {
            alpha[i] =
                3.0
                * ((y_[i + 1] - y_[i]) / h[i] - (y_[i] - y_[i - 1]) / h[i - 1]);
        }

        // Allocate memory for coefficient c
        c_.assign(n, 0.0);

        // Dispatch to the corresponding solver based on boundary condition
        switch (bc_type_) {
            case BoundaryCondition::Natural:
                compute_natural(h, alpha);
                break;
            case BoundaryCondition::NotAKnot:
                compute_not_a_knot(h, alpha);
                break;
            case BoundaryCondition::Clamped:
                compute_clamped(h, alpha);
                break;
        }

        // Calculate final b_i and d_i coefficients
        b_.resize(n - 1);
        d_.resize(n - 1);
        for (size_t i = 0; i < n - 1; ++i) {
            b_[i] = (y_[i + 1] - y_[i]) / h[i]
                    - h[i] * (c_[i + 1] + 2.0 * c_[i]) / 3.0;
            d_[i] = (c_[i + 1] - c_[i]) / (3.0 * h[i]);
        }
    }

    void compute_natural(const std::vector<double> &h,
                         const std::vector<double> &alpha) {
        size_t n = x_.size();
        std::vector<double> L(n, 0.0), D(n, 1.0), U(n, 0.0), B(n, 0.0);

        for (size_t i = 1; i < n - 1; ++i) {
            L[i] = h[i - 1];
            D[i] = 2.0 * (h[i - 1] + h[i]);
            U[i] = h[i];
            B[i] = alpha[i];
        }

        c_ = solve_tridiagonal(L, D, U, B);
    }

    void compute_clamped(const std::vector<double> &h,
                         const std::vector<double> &alpha) {
        size_t n = x_.size();
        std::vector<double> L(n, 0.0), D(n, 1.0), U(n, 0.0), B(n, 0.0);

        for (size_t i = 1; i < n - 1; ++i) {
            L[i] = h[i - 1];
            D[i] = 2.0 * (h[i - 1] + h[i]);
            U[i] = h[i];
            B[i] = alpha[i];
        }

        // Apply clamped boundary specific formulas to the first and last row
        D[0] = 2.0 * h[0];
        U[0] = h[0];
        B[0] = 3.0 * ((y_[1] - y_[0]) / h[0] - left_slope_);

        D[n - 1] = 2.0 * h[n - 2];
        L[n - 1] = h[n - 2];
        B[n - 1] = 3.0 * (right_slope_ - (y_[n - 1] - y_[n - 2]) / h[n - 2]);

        c_ = solve_tridiagonal(L, D, U, B);
    }

    void compute_not_a_knot(const std::vector<double> &h,
                            const std::vector<double> &alpha) {
        size_t n = x_.size();

        // Not-a-knot strictly requires at least 4 points to form the (N-2)
        // subsystem. Fallback to Natural boundary if only 3 points are
        // provided.
        if (n == 3) {
            compute_natural(h, alpha);
            return;
        }

        size_t N_reduced = n - 2;
        std::vector<double> L(N_reduced, 0.0), D(N_reduced, 0.0),
            U(N_reduced, 0.0), B(N_reduced, 0.0);

        // Construct the strictly tridiagonal reduced system
        D[0] = (h[0] + h[1]) * (h[0] + 2.0 * h[1]);
        U[0] = h[1] * h[1] - h[0] * h[0];
        B[0] = h[1] * alpha[1];

        for (size_t i = 1; i < N_reduced - 1; ++i) {
            L[i] = h[i];
            D[i] = 2.0 * (h[i] + h[i + 1]);
            U[i] = h[i + 1];
            B[i] = alpha[i + 1];
        }

        L[N_reduced - 1] = h[n - 3] * h[n - 3] - h[n - 2] * h[n - 2];
        D[N_reduced - 1] = (h[n - 3] + h[n - 2]) * (h[n - 2] + 2.0 * h[n - 3]);
        B[N_reduced - 1] = h[n - 3] * alpha[n - 2];

        // Solve the reduced (N-2) x (N-2) system
        std::vector<double> c_sol = solve_tridiagonal(L, D, U, B);

        // Back-substitute to reconstruct the full c_ array of size N
        c_[0] = (h[0] + h[1]) / h[1] * c_sol[0] - h[0] / h[1] * c_sol[1];
        for (size_t i = 0; i < N_reduced; ++i) {
            c_[i + 1] = c_sol[i];
        }
        c_[n - 1] = (h[n - 3] + h[n - 2]) / h[n - 3] * c_sol[N_reduced - 1]
                    - h[n - 2] / h[n - 3] * c_sol[N_reduced - 2];
    }

    std::vector<double> solve_tridiagonal(const std::vector<double> &L,
                                          const std::vector<double> &D,
                                          const std::vector<double> &U,
                                          const std::vector<double> &B) const {
        size_t n = D.size();
        std::vector<double> c_sol(n);
        std::vector<double> cp(n); // c-prime
        std::vector<double> bp(n); // B-prime

        if (std::abs(D[0]) < 1e-14) {
            throw std::runtime_error(
                "Zero pivot encountered in Thomas algorithm.");
        }

        // Forward elimination
        cp[0] = U[0] / D[0];
        bp[0] = B[0] / D[0];

        for (size_t i = 1; i < n; i++) {
            double denom = D[i] - L[i] * cp[i - 1];
            if (std::abs(denom) < 1e-14) {
                throw std::runtime_error(
                    "Division by zero in tridiagonal solver.");
            }
            if (i < n - 1)
                cp[i] = U[i] / denom;
            bp[i] = (B[i] - L[i] * bp[i - 1]) / denom;
        }

        // Backward substitution
        c_sol[n - 1] = bp[n - 1];
        for (int i = static_cast<int>(n) - 2; i >= 0; i--) {
            c_sol[i] = bp[i] - cp[i] * c_sol[i + 1];
        }
        return c_sol;
    }
};

/**
 * @brief Not-a-knot cubic spline interpolation (zero-copy output).
 *
 * @param x Independent variable samples
 * @param y Dependent variable samples
 * @param x_new Query points
 * @param result Output buffer (must have same size as @p x_new)
 */
inline void interp1_cubic(std::span<const double> x,
                          std::span<const double> y,
                          std::span<const double> x_new,
                          std::span<double> result) {
    if (x_new.size() != result.size()) {
        throw std::invalid_argument(
            "interp1_cubic: x_new and result spans must have same size");
    }

    auto spline = CubicSpline::from_data(x, y);
    spline(x_new, result);
}

/**
 * @brief Not-a-knot cubic spline interpolation.
 *
 * @param x Independent variable samples
 * @param y Dependent variable samples
 * @param x_new Query points
 * @return Interpolated results at @p x_new
 */
inline std::vector<double> interp1_cubic(std::span<const double> x,
                                         std::span<const double> y,
                                         std::span<const double> x_new) {
    return CubicSpline::from_data(x, y)(x_new);
}

} // namespace msl::interp

#endif // MSL_INTERP_1D_CUBIC