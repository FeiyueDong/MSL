/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2026, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: romberg.hpp
** -----
** File Created: Friday, 9th January 2026 14:58:10
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Thursday, 5th March 2026 14:57:54
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_ROMBERG_HPP
#define MSL_ROMBERG_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "matrix/real_matrix_base.hpp"
#include "trapz.hpp"

namespace msl::integral {

namespace detail {
inline size_t romberg_levels_from_sample_count(size_t n) {
    if (n < 3) {
        throw std::invalid_argument("Romberg: needs at least 3 points");
    }

    size_t intervals = n - 1;
    size_t levels = 0;
    size_t check = 1;
    while (check < intervals) {
        check *= 2;
        ++levels;
    }
    if (check != intervals) {
        throw std::invalid_argument("Romberg: needs 2^k + 1 points");
    }

    return levels;
}

inline void validate_uniform_spacing(std::span<const double> x) {
    if (x.size() < 2) {
        throw std::invalid_argument("Romberg: needs at least 2 x values");
    }

    double dx = x[1] - x[0];
    if (dx <= 0.0) {
        throw std::invalid_argument(
            "Romberg: x values must be strictly increasing");
    }

    double scale = std::max(1.0, std::abs(dx));
    double tol = 64.0 * std::numeric_limits<double>::epsilon() * scale;
    for (size_t i = 2; i < x.size(); ++i) {
        double current_dx = x[i] - x[i - 1];
        if (current_dx <= 0.0 || std::abs(current_dx - dx) > tol) {
            throw std::invalid_argument(
                "Romberg: x values must be uniformly spaced");
        }
    }
}
} // namespace detail

// ============================================================================
// Total Integral (Romberg Integration)
// ============================================================================

/**
 * @brief Romberg integration for uniformly sampled data
 *
 * Uses Richardson extrapolation on a nested trapezoidal table. The samples must
 * contain 2^k + 1 equally spaced points.
 *
 * @param y Function values (size must be 2^k + 1)
 * @param dx Uniform spacing
 * @param tol Tolerance for convergence (default: 1e-10)
 * @return Integral value
 */
inline double
romberg(std::span<const double> y, double dx, double tol = 1e-10) {
    if (dx <= 0.0) {
        throw std::invalid_argument("Romberg: dx must be positive");
    }
    if (tol <= 0.0) {
        throw std::invalid_argument("Romberg: tolerance must be positive");
    }

    size_t n = y.size();
    size_t k = detail::romberg_levels_from_sample_count(n);
    size_t intervals = n - 1;

    // Romberg table
    std::vector<std::vector<double>> R(k + 1);
    for (size_t i = 0; i <= k; ++i) {
        R[i].resize(i + 1);
    }

    // R[0,0] is the one-panel trapezoidal rule over the whole interval.
    R[0][0] =
        0.5 * (y.front() + y.back()) * dx * static_cast<double>(intervals);

    for (size_t i = 1; i <= k; ++i) {
        size_t panels = size_t{1} << i;
        size_t stride = intervals / panels;
        double h = dx * static_cast<double>(stride);
        double sum = 0.5 * (y.front() + y.back());
        for (size_t j = stride; j < intervals; j += stride) {
            sum += y[j];
        }
        R[i][0] = sum * h;

        double power = 4.0;
        for (size_t j = 1; j <= i; ++j) {
            R[i][j] = (power * R[i][j - 1] - R[i - 1][j - 1]) / (power - 1.0);
            power *= 4.0;
        }

        if (std::abs(R[i][i] - R[i - 1][i - 1]) < tol) {
            return R[i][i];
        }
    }

    return R[k][k];
}

/**
 * @brief Romberg integration for uniformly spaced x/y samples
 *
 * @param x Independent variable values, uniformly spaced
 * @param y Function values at x points
 * @param tol Tolerance for convergence
 * @return Integral value
 */
inline double romberg(std::span<const double> x,
                      std::span<const double> y,
                      double tol = 1e-10) {
    if (x.size() != y.size()) {
        throw std::invalid_argument("Romberg: x and y must have same size");
    }
    detail::validate_uniform_spacing(x);
    return romberg(y, x[1] - x[0], tol);
}

/**
 * @brief Romberg integration for function
 *
 * Uses the standard nested composite trapezoidal sequence and Richardson
 * extrapolation.
 *
 * @param f Function to integrate
 * @param a Lower limit
 * @param b Upper limit
 * @param tol Tolerance for convergence
 * @param max_iter Maximum number of Romberg levels
 * @return Integral value
 */
template <typename Func>
    requires std::is_invocable_r_v<double, Func, double>
inline double romberg(Func &&f,
                      double a,
                      double b,
                      double tol = 1e-10,
                      size_t max_iter = 20) {
    if (a >= b) {
        throw std::invalid_argument(
            "Romberg: lower limit must be less than upper limit");
    }
    if (tol <= 0.0) {
        throw std::invalid_argument("Romberg: tolerance must be positive");
    }
    if (max_iter == 0) {
        throw std::invalid_argument("Romberg: max_iter must be positive");
    }

    std::vector<std::vector<double>> R(max_iter);
    for (size_t i = 0; i < max_iter; ++i) {
        R[i].resize(i + 1);
    }

    double h = b - a;
    R[0][0] = 0.5 * h * (std::invoke(f, a) + std::invoke(f, b));

    for (size_t i = 1; i < max_iter; ++i) {
        h *= 0.5;
        size_t new_points = size_t{1} << (i - 1);
        double sum = 0.0;
        for (size_t j = 0; j < new_points; ++j) {
            double x = a + static_cast<double>(2 * j + 1) * h;
            sum += std::invoke(f, x);
        }

        R[i][0] = 0.5 * R[i - 1][0] + h * sum;

        double power = 4.0;
        for (size_t j = 1; j <= i; ++j) {
            R[i][j] = (power * R[i][j - 1] - R[i - 1][j - 1]) / (power - 1.0);
            power *= 4.0;
        }

        if (std::abs(R[i][i] - R[i - 1][i - 1]) < tol) {
            return R[i][i];
        }
    }

    return R[max_iter - 1][max_iter - 1];
}

/**
 * @brief Romberg integration for matrix (column-wise, uniform spacing)
 *
 * Zero-copy operation: directly fills the provided result span without
 * internal allocation. Integrates each column independently using Romberg's
 * method.
 *
 * @param mat Input matrix (each column is a function, rows must be 2^k + 1)
 * @param result Output buffer for integral values (must have same size as
 * number of columns)
 * @param dx Spacing between rows
 * @param tol Tolerance for convergence (default: 1e-10)
 */
inline void romberg(const matrix::real_matrix_base &mat,
                    std::span<double> result,
                    double dx,
                    double tol = 1e-10) {
    if (mat.rows() < 3) {
        throw std::invalid_argument("Romberg: needs at least 3 rows");
    }

    if (result.size() != mat.cols()) {
        throw std::invalid_argument(
            "Romberg: result span must have same size as number of columns");
    }

    for (size_t j = 0; j < mat.cols(); ++j) {
        std::vector<double> col(mat.rows());
        for (size_t i = 0; i < mat.rows(); ++i) {
            col[i] = mat(i, j);
        }
        result[j] = romberg(col, dx, tol);
    }
}

/**
 * @brief Romberg integration for matrix columns with uniform row coordinates
 *
 * @param x Independent variable values, uniformly spaced
 * @param mat Input matrix
 * @param result Output buffer for integral values
 * @param tol Tolerance for convergence
 */
inline void romberg(std::span<const double> x,
                    const matrix::real_matrix_base &mat,
                    std::span<double> result,
                    double tol = 1e-10) {
    if (x.size() != mat.rows()) {
        throw std::invalid_argument(
            "Romberg: x size must match number of rows");
    }
    if (result.size() != mat.cols()) {
        throw std::invalid_argument(
            "Romberg: result span must have same size as number of columns");
    }

    detail::validate_uniform_spacing(x);
    for (size_t j = 0; j < mat.cols(); ++j) {
        std::vector<double> col(mat.rows());
        for (size_t i = 0; i < mat.rows(); ++i) {
            col[i] = mat(i, j);
        }
        result[j] = romberg(x, col, tol);
    }
}

/**
 * @brief Romberg integration for matrix (column-wise, uniform spacing)
 *
 * Convenience wrapper that allocates and returns a vector. For zero-copy
 * operations, use the void version with span parameter.
 *
 * Integrates each column independently using Romberg's method.
 *
 * @param mat Input matrix (each column is a function, rows must be 2^k + 1)
 * @param dx Spacing between rows
 * @param tol Tolerance for convergence (default: 1e-10)
 * @return Vector of integral values for each column
 */
inline std::vector<double>
romberg(const matrix::real_matrix_base &mat, double dx, double tol = 1e-10) {
    std::vector<double> result(mat.cols());
    romberg(mat, result, dx, tol);
    return result;
}

/**
 * @brief Romberg integration for matrix columns with uniform row coordinates
 *
 * @param x Independent variable values, uniformly spaced
 * @param mat Input matrix
 * @param tol Tolerance for convergence
 * @return Vector of integral values for each column
 */
inline std::vector<double> romberg(std::span<const double> x,
                                   const matrix::real_matrix_base &mat,
                                   double tol = 1e-10) {
    std::vector<double> result(mat.cols());
    romberg(x, mat, result, tol);
    return result;
}

} // namespace msl::integral

#endif // MSL_ROMBERG_HPP
