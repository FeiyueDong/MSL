/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2026, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: trapz.hpp
** -----
** File Created: Friday, 9th January 2026 14:58:10
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Thursday, 5th March 2026 14:57:47
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_TRAPZ_HPP
#define MSL_TRAPZ_HPP

#include <cstddef>
#include <functional>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "matrix/real_matrix_base.hpp"
#include "matrix/real_matrix_owned.hpp"

namespace msl::integral {

// ============================================================================
// Total Integral (Trapezoidal Rule)
// ============================================================================

/**
 * @brief Compute total integral using trapezoidal rule (uniform spacing)
 *
 * @param y Function values
 * @param dx Spacing
 * @return Total integral value
 */
inline double trapz(std::span<const double> y, double dx) {
    if (y.size() < 2) {
        throw std::invalid_argument(
            "Trapz: need at least 2 points for integration");
    }

    double sum = 0.5 * (y.front() + y.back());
    for (size_t i = 1; i < y.size() - 1; ++i) {
        sum += y[i];
    }

    return sum * dx;
}

/**
 * @brief Compute total integral using trapezoidal rule (non-uniform spacing)
 *
 * @param x Independent variable values
 * @param y Function values at x points
 * @return Total integral value
 */
inline double trapz(std::span<const double> x, std::span<const double> y) {
    if (x.size() != y.size()) {
        throw std::invalid_argument("Trapz: x and y must have same size");
    }
    if (x.size() < 2) {
        throw std::invalid_argument(
            "Trapz: need at least 2 points for integration");
    }

    double sum = 0.0;
    for (size_t i = 1; i < y.size(); ++i) {
        sum += 0.5 * (y[i] + y[i - 1]) * (x[i] - x[i - 1]);
    }

    return sum;
}

/**
 * @brief Compute total integral using composite trapezoidal rule for function
 *
 * @param f Function to integrate
 * @param a Lower limit
 * @param b Upper limit
 * @param intervals Number of subintervals
 * @return Total integral value
 */
template <typename Func>
    requires std::is_invocable_r_v<double, Func, double>
inline double trapz(Func &&f, double a, double b, size_t intervals) {
    if (a >= b) {
        throw std::invalid_argument(
            "Trapz: lower limit must be less than upper limit");
    }
    if (intervals == 0) {
        throw std::invalid_argument(
            "Trapz: need at least 1 interval for integration");
    }

    double dx = (b - a) / static_cast<double>(intervals);
    double sum = 0.5 * (std::invoke(f, a) + std::invoke(f, b));
    for (size_t i = 1; i < intervals; ++i) {
        sum += std::invoke(f, a + static_cast<double>(i) * dx);
    }

    return sum * dx;
}

/**
 * @brief Compute total integral using trapezoidal rule (uniform spacing)
 *
 * Convenience wrapper that allocates and returns a vector.
 * For zero-copy operations, use the void version with span parameter.
 *
 * @param y Function values
 * @param dx Spacing
 * @return Total integral value
 */
inline void trapz(const matrix::real_matrix_base &mat,
                  std::span<double> result,
                  double dx) {
    if (mat.rows() < 2) {
        throw std::invalid_argument(
            "Trapz: need at least 2 rows for integration");
    }
    if (result.size() != mat.cols()) {
        throw std::invalid_argument(
            "Trapz: result span must have same size as number of columns");
    }

    for (size_t j = 0; j < mat.cols(); ++j) {
        double sum = 0.5 * (mat(0, j) + mat(mat.rows() - 1, j));
        for (size_t i = 1; i < mat.rows() - 1; ++i) {
            sum += mat(i, j);
        }
        result[j] = sum * dx;
    }
}

/**
 * @brief Compute total integral for matrix columns (non-uniform spacing)
 *
 * @param x Independent variable values, one per matrix row
 * @param mat Input matrix (each column is a function)
 * @param result Output buffer for integral values
 */
inline void trapz(std::span<const double> x,
                  const matrix::real_matrix_base &mat,
                  std::span<double> result) {
    if (x.size() != mat.rows()) {
        throw std::invalid_argument("Trapz: x size must match number of rows");
    }
    if (mat.rows() < 2) {
        throw std::invalid_argument(
            "Trapz: need at least 2 rows for integration");
    }
    if (result.size() != mat.cols()) {
        throw std::invalid_argument(
            "Trapz: result span must have same size as number of columns");
    }

    for (size_t j = 0; j < mat.cols(); ++j) {
        double sum = 0.0;
        for (size_t i = 1; i < mat.rows(); ++i) {
            sum += 0.5 * (mat(i, j) + mat(i - 1, j)) * (x[i] - x[i - 1]);
        }
        result[j] = sum;
    }
}

/**
 * @brief Compute total integral using trapezoidal rule (uniform spacing)
 *
 * Convenience wrapper that allocates and returns a vector. For zero-copy
 * operations, use the void version with span parameter.
 *
 * @param mat Input matrix (each column is a function)
 * @param dx Spacing between rows
 * @return Vector of integral values for each column
 */
inline std::vector<double> trapz(const matrix::real_matrix_base &mat,
                                 double dx) {
    std::vector<double> result(mat.cols());
    trapz(mat, result, dx);
    return result;
}

/**
 * @brief Compute total integral for matrix columns (non-uniform spacing)
 *
 * @param x Independent variable values, one per matrix row
 * @param mat Input matrix
 * @return Vector of integral values for each column
 */
inline std::vector<double> trapz(std::span<const double> x,
                                 const matrix::real_matrix_base &mat) {
    std::vector<double> result(mat.cols());
    trapz(x, mat, result);
    return result;
}

// ============================================================================
// Cumulative Integral (Trapezoidal Rule)
// ============================================================================

/**
 * @brief Cumulative trapezoidal integration for vector (uniform spacing, with
 * output buffer)
 *
 * Computes cumulative integral: result[i] = integral from 0 to i of y dx.
 * Zero-copy: directly fills the provided result span without internal
 * allocation.
 *
 * @param y Function values at equally spaced points
 * @param result Output buffer for cumulative integral values
 * @param dx Spacing between points
 */
inline void
cumtrapz(std::span<const double> y, std::span<double> result, double dx) {
    if (y.size() < 2) {
        throw std::invalid_argument(
            "Cumtrapz: needs at least 2 points for integration");
    }
    if (result.size() != y.size()) {
        throw std::invalid_argument(
            "Cumtrapz: result span must have same size as y");
    }

    result[0] = 0.0;
    for (size_t i = 1; i < y.size(); ++i) {
        result[i] = result[i - 1] + 0.5 * (y[i] + y[i - 1]) * dx;
    }
}

/**
 * @brief Cumulative trapezoidal integration for vector (uniform spacing)
 *
 * @param y Function values at equally spaced points
 * @param dx Spacing between points
 * @return Cumulative integral values
 */
inline std::vector<double> cumtrapz(std::span<const double> y, double dx) {
    std::vector<double> result(y.size());
    cumtrapz(y, result, dx);
    return result;
}

/**
 * @brief Cumulative trapezoidal integration for vector (non-uniform spacing,
 * with output buffer)
 *
 * @param x Independent variable values
 * @param y Function values at x points
 * @param result Output buffer for cumulative integral values
 */
inline void cumtrapz(std::span<const double> x,
                     std::span<const double> y,
                     std::span<double> result) {
    if (x.size() != y.size()) {
        throw std::invalid_argument("Cumtrapz: x and y must have same size");
    }
    if (x.size() < 2) {
        throw std::invalid_argument(
            "Cumtrapz: needs at least 2 points for integration");
    }
    if (result.size() != y.size()) {
        throw std::invalid_argument(
            "Cumtrapz: result span must have same size as y");
    }

    result[0] = 0.0;
    for (size_t i = 1; i < y.size(); ++i) {
        double dx = x[i] - x[i - 1];
        result[i] = result[i - 1] + 0.5 * (y[i] + y[i - 1]) * dx;
    }
}

/**
 * @brief Cumulative trapezoidal integration for vector (non-uniform spacing)
 *
 * @param x Independent variable values
 * @param y Function values at x points
 * @return Cumulative integral values
 */
inline std::vector<double> cumtrapz(std::span<const double> x,
                                    std::span<const double> y) {
    std::vector<double> result(y.size());
    cumtrapz(x, y, result);
    return result;
}

/**
 * @brief Cumulative trapezoidal integration for matrix columns (uniform
 * spacing)
 *
 * @param mat Input matrix (each column is a function)
 * @param dx Spacing between rows
 * @return Matrix of cumulative integrals
 */
inline matrix::matrixd cumtrapz(const matrix::real_matrix_base &mat,
                                double dx) {
    if (mat.rows() < 2) {
        throw std::invalid_argument(
            "Cumtrapz: needs at least 2 rows for integration");
    }

    matrix::matrixd result(mat.rows(), mat.cols(), 0.0);

    for (size_t j = 0; j < mat.cols(); ++j) {
        for (size_t i = 1; i < mat.rows(); ++i) {
            result(i, j) =
                result(i - 1, j) + 0.5 * (mat(i, j) + mat(i - 1, j)) * dx;
        }
    }

    return result;
}

/**
 * @brief Cumulative trapezoidal integration for matrix columns (non-uniform
 * spacing)
 *
 * @param x Independent variable values, one per matrix row
 * @param mat Input matrix
 * @return Matrix of cumulative integrals
 */
inline matrix::matrixd cumtrapz(std::span<const double> x,
                                const matrix::real_matrix_base &mat) {
    if (x.size() != mat.rows()) {
        throw std::invalid_argument(
            "Cumtrapz: x size must match number of rows");
    }
    if (mat.rows() < 2) {
        throw std::invalid_argument(
            "Cumtrapz: needs at least 2 rows for integration");
    }

    matrix::matrixd result(mat.rows(), mat.cols(), 0.0);

    for (size_t j = 0; j < mat.cols(); ++j) {
        for (size_t i = 1; i < mat.rows(); ++i) {
            double dx = x[i] - x[i - 1];
            result(i, j) =
                result(i - 1, j) + 0.5 * (mat(i, j) + mat(i - 1, j)) * dx;
        }
    }

    return result;
}

} // namespace msl::integral

#endif // MSL_TRAPZ_HPP
