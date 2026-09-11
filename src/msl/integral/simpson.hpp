/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2026, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: simpson.hpp
** -----
** File Created: Friday, 9th January 2026 14:58:10
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Thursday, 5th March 2026 14:58:00
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_SIMPSON_HPP
#define MSL_SIMPSON_HPP

#include <cstddef>
#include <functional>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "matrix/real_matrix_base.hpp"

namespace msl::integral {

namespace detail {
inline double simpson_nonuniform_segment(double x0,
                                         double x1,
                                         double x2,
                                         double y0,
                                         double y1,
                                         double y2) {
    double h0 = x1 - x0;
    double h1 = x2 - x1;
    if (h0 <= 0.0 || h1 <= 0.0) {
        throw std::invalid_argument(
            "Simpson: x values must be strictly increasing");
    }

    double hsum = h0 + h1;
    return hsum / 6.0
           * ((2.0 - h1 / h0) * y0 + (hsum * hsum / (h0 * h1)) * y1
              + (2.0 - h0 / h1) * y2);
}

/**
 * @brief Integral of the quadratic through (x0,y0), (x1,y1), (x2,y2) over the
 * first half interval [x0, x1].
 *
 * Used for cumulative Simpson integration at odd sample indices, so that the
 * cumulative value at every point is obtained from the same local quadratic
 * interpolant.
 */
inline double simpson_nonuniform_half_segment(double x0,
                                              double x1,
                                              double x2,
                                              double y0,
                                              double y1,
                                              double y2) {
    double h0 = x1 - x0;
    double h1 = x2 - x1;
    if (h0 <= 0.0 || h1 <= 0.0) {
        throw std::invalid_argument(
            "Simpson: x values must be strictly increasing");
    }

    const double slope01 = (y1 - y0) / h0;
    const double slope12 = (y2 - y1) / h1;
    const double second_divided_difference = (slope12 - slope01) / (x2 - x0);

    return y0 * h0 + 0.5 * h0 * (y1 - y0)
           - (h0 * h0 * h0 / 6.0) * second_divided_difference;
}
} // namespace detail

// ============================================================================
// Total Integral (Simpson's Rule Integration)
// ============================================================================

/**
 * @brief Simpson's 1/3 rule (requires odd number of points)
 *
 * More accurate than trapezoidal for smooth functions
 *
 * @param y Function values (size must be odd)
 * @param dx Uniform spacing
 * @return Integral value
 */
inline double simpson(std::span<const double> y, double dx) {
    if (y.size() < 3) {
        throw std::invalid_argument("Simpson: needs at least 3 points");
    }

    size_t n = y.size();

    // If even number of points, use trapezoidal for last segment
    bool use_trap_last = (n % 2 == 0);
    size_t n_simp = use_trap_last ? n - 1 : n;

    double sum = y[0] + y[n_simp - 1];

    // Odd indices (weight = 4)
    for (size_t i = 1; i < n_simp - 1; i += 2) {
        sum += 4.0 * y[i];
    }

    // Even indices (weight = 2)
    for (size_t i = 2; i < n_simp - 1; i += 2) {
        sum += 2.0 * y[i];
    }

    double result = sum * dx / 3.0;

    // Add trapezoidal correction for last segment if needed
    if (use_trap_last) {
        result += 0.5 * (y[n - 2] + y[n - 1]) * dx;
    }

    return result;
}

/**
 * @brief Simpson's rule for non-uniformly spaced samples
 *
 * Uses three-point quadratic integration for each pair of intervals. If the
 * number of points is even, the last interval is integrated with trapezoidal
 * rule.
 *
 * @param x Independent variable values, strictly increasing
 * @param y Function values at x points
 * @return Integral value
 */
inline double simpson(std::span<const double> x, std::span<const double> y) {
    if (x.size() != y.size()) {
        throw std::invalid_argument("Simpson: x and y must have same size");
    }
    if (x.size() < 3) {
        throw std::invalid_argument("Simpson: needs at least 3 points");
    }

    double result = 0.0;
    size_t i = 2;
    for (; i < y.size(); i += 2) {
        result += detail::simpson_nonuniform_segment(
            x[i - 2], x[i - 1], x[i], y[i - 2], y[i - 1], y[i]);
    }

    if (y.size() % 2 == 0) {
        size_t n = y.size();
        double dx = x[n - 1] - x[n - 2];
        if (dx <= 0.0) {
            throw std::invalid_argument(
                "Simpson: x values must be strictly increasing");
        }
        result += 0.5 * (y[n - 2] + y[n - 1]) * dx;
    }

    return result;
}

/**
 * @brief Simpson's rule for function on uniform grid
 *
 * @param f Function to integrate
 * @param a Lower limit
 * @param b Upper limit
 * @param intervals Number of subintervals. If odd, the last interval uses
 * trapezoidal rule.
 * @return Integral value
 */
template <typename Func>
    requires std::is_invocable_r_v<double, Func, double>
inline double simpson(Func &&f, double a, double b, size_t intervals) {
    if (a >= b) {
        throw std::invalid_argument(
            "Simpson: lower limit must be less than upper limit");
    }
    if (intervals < 2) {
        throw std::invalid_argument(
            "Simpson: needs at least 2 intervals for integration");
    }

    double dx = (b - a) / static_cast<double>(intervals);
    bool use_trap_last = (intervals % 2 != 0);
    size_t simp_intervals = use_trap_last ? intervals - 1 : intervals;

    double sum = std::invoke(f, a)
                 + std::invoke(f, a + static_cast<double>(simp_intervals) * dx);
    for (size_t i = 1; i < simp_intervals; i += 2) {
        sum += 4.0 * std::invoke(f, a + static_cast<double>(i) * dx);
    }
    for (size_t i = 2; i < simp_intervals; i += 2) {
        sum += 2.0 * std::invoke(f, a + static_cast<double>(i) * dx);
    }

    double result = sum * dx / 3.0;
    if (use_trap_last) {
        double x0 = a + static_cast<double>(intervals - 1) * dx;
        result += 0.5 * (std::invoke(f, x0) + std::invoke(f, b)) * dx;
    }

    return result;
}

/**
 * @brief Cumulative Simpson's rule for vector (uniform spacing, with output
 * buffer)
 *
 * Computes the cumulative integral `result[i] = ∫[x0 to xi] y dx`. For even
 * indices, uses Simpson's rule over the two preceding intervals. For odd
 * indices, integrates the same three-point interpolating parabola over the
 * single preceding interval:
 *
 *   result[i-1] = result[i-2] + dx * (5*y[i-2] + 8*y[i-1] - y[i]) / 12
 *   result[i]   = result[i-2] + dx * (y[i-2] + 4*y[i-1] + y[i]) / 3
 *
 * Both values are therefore exact for quadratics. If the number of points is
 * even, the last interval is integrated with the trapezoidal rule.
 *
 * Zero-copy operation: directly fills the provided result span without internal
 * allocation.
 *
 * @param y Function values at equally spaced points
 * @param result Output buffer for cumulative integral values (must have same
 * size as y)
 * @param dx Spacing between points
 */
inline void
cumsimpson(std::span<const double> y, std::span<double> result, double dx) {
    if (y.size() < 3) {
        throw std::invalid_argument("Simpson: needs at least 3 points");
    }

    if (result.size() != y.size()) {
        throw std::invalid_argument(
            "Simpson: result span must have same size as input");
    }

    // First point
    result[0] = 0.0;

    // Integrate each pair of intervals with the three-point parabola
    for (size_t i = 2; i < y.size(); i += 2) {
        result[i - 1] = result[i - 2]
                        + (5.0 * y[i - 2] + 8.0 * y[i - 1] - y[i]) * dx / 12.0;
        result[i] =
            result[i - 2] + (y[i - 2] + 4.0 * y[i - 1] + y[i]) * dx / 3.0;
    }

    // Handle last point if even number of points
    if (y.size() % 2 == 0) {
        size_t n = y.size();
        result[n - 1] = result[n - 2] + 0.5 * (y[n - 2] + y[n - 1]) * dx;
    }
}

/**
 * @brief Cumulative Simpson's rule for non-uniformly spaced samples
 *
 * Even indices use three-point quadratic integration. Odd indices use the
 * integral of the same three-point quadratic over the preceding single
 * interval, so cumulative values are exact for quadratics at every point. If
 * the number of points is even, the last interval is integrated with the
 * trapezoidal rule.
 *
 * @param x Independent variable values, strictly increasing
 * @param y Function values at x points
 * @param result Output buffer for cumulative integral values
 */
inline void cumsimpson(std::span<const double> x,
                       std::span<const double> y,
                       std::span<double> result) {
    if (x.size() != y.size()) {
        throw std::invalid_argument("Simpson: x and y must have same size");
    }
    if (x.size() < 3) {
        throw std::invalid_argument("Simpson: needs at least 3 points");
    }
    if (result.size() != y.size()) {
        throw std::invalid_argument(
            "Simpson: result span must have same size as input");
    }

    result[0] = 0.0;
    for (size_t i = 2; i < y.size(); i += 2) {
        result[i - 1] = result[i - 2]
                        + detail::simpson_nonuniform_half_segment(
                            x[i - 2], x[i - 1], x[i], y[i - 2], y[i - 1], y[i]);
        result[i] = result[i - 2]
                    + detail::simpson_nonuniform_segment(
                        x[i - 2], x[i - 1], x[i], y[i - 2], y[i - 1], y[i]);
    }

    if (y.size() % 2 == 0) {
        size_t n = y.size();
        double dx = x[n - 1] - x[n - 2];
        if (dx <= 0.0) {
            throw std::invalid_argument(
                "Simpson: x values must be strictly increasing");
        }
        result[n - 1] = result[n - 2] + 0.5 * (y[n - 2] + y[n - 1]) * dx;
    }
}

/**
 * @brief Cumulative Simpson's rule for vector (uniform spacing)
 *
 * Computes cumulative integral: output[i] = ∫[0 to i] y dx using Simpson's
 * rule. Even indices use Simpson's rule over the two preceding intervals; odd
 * indices integrate the same interpolating parabola over the single preceding
 * interval. If the number of points is even, the last interval is integrated
 * with the trapezoidal rule.
 *
 * @param y Function values at equally spaced points
 * @param dx Spacing between points
 * @return Cumulative integral values
 */
inline std::vector<double> cumsimpson(std::span<const double> y, double dx) {
    std::vector<double> result(y.size());
    cumsimpson(y, result, dx);
    return result;
}

/**
 * @brief Cumulative Simpson's rule for non-uniformly spaced samples
 *
 * @param x Independent variable values
 * @param y Function values at x points
 * @return Cumulative integral values
 */
inline std::vector<double> cumsimpson(std::span<const double> x,
                                      std::span<const double> y) {
    std::vector<double> result(y.size());
    cumsimpson(x, y, result);
    return result;
}

/**
 * @brief Simpson's rule for matrix (column-wise, uniform spacing)
 *
 * Integrates each column independently using Simpson's rule.
 *
 * @param mat Input matrix (each column is a function)
 * @param dx Spacing between rows
 * @return Vector of integral values for each column
 */
inline void simpson(const matrix::real_matrix_base &mat,
                    std::span<double> result,
                    double dx) {
    if (mat.rows() < 3) {
        throw std::invalid_argument("Simpson: needs at least 3 rows");
    }
    if (result.size() != mat.cols()) {
        throw std::invalid_argument(
            "Simpson: result span must have same size as number of columns");
    }

    for (size_t j = 0; j < mat.cols(); ++j) {
        std::vector<double> col(mat.rows());
        for (size_t i = 0; i < mat.rows(); ++i) {
            col[i] = mat(i, j);
        }
        result[j] = simpson(col, dx);
    }
}

/**
 * @brief Simpson's rule for matrix columns with non-uniform row coordinates
 *
 * @param x Independent variable values, one per matrix row
 * @param mat Input matrix (each column is a function)
 * @param result Output buffer for integral values
 */
inline void simpson(std::span<const double> x,
                    const matrix::real_matrix_base &mat,
                    std::span<double> result) {
    if (x.size() != mat.rows()) {
        throw std::invalid_argument(
            "Simpson: x size must match number of rows");
    }
    if (mat.rows() < 3) {
        throw std::invalid_argument("Simpson: needs at least 3 rows");
    }
    if (result.size() != mat.cols()) {
        throw std::invalid_argument(
            "Simpson: result span must have same size as number of columns");
    }

    for (size_t j = 0; j < mat.cols(); ++j) {
        std::vector<double> col(mat.rows());
        for (size_t i = 0; i < mat.rows(); ++i) {
            col[i] = mat(i, j);
        }
        result[j] = simpson(x, col);
    }
}

/**
 * @brief Simpson's rule for matrix (column-wise, uniform spacing)
 *
 * Integrates each column independently using Simpson's rule.
 *
 * @param mat Input matrix (each column is a function)
 * @param dx Spacing between rows
 * @return Vector of integral values for each column
 */
inline std::vector<double> simpson(const matrix::real_matrix_base &mat,
                                   double dx) {
    std::vector<double> result(mat.cols());
    simpson(mat, result, dx);
    return result;
}

/**
 * @brief Simpson's rule for matrix columns with non-uniform row coordinates
 *
 * @param x Independent variable values, one per matrix row
 * @param mat Input matrix
 * @return Vector of integral values for each column
 */
inline std::vector<double> simpson(std::span<const double> x,
                                   const matrix::real_matrix_base &mat) {
    std::vector<double> result(mat.cols());
    simpson(x, mat, result);
    return result;
}

} // namespace msl::integral

#endif // MSL_SIMPSON_HPP
