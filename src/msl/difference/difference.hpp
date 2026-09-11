/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2026, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: difference.hpp
** -----
** File Created: Friday, 9th January 2026 14:58:10
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Thursday, 5th March 2026 14:58:13
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_DIFFERENCE_HPP
#define MSL_DIFFERENCE_HPP

#include <algorithm>
#include <span>
#include <stdexcept>
#include <vector>

#include "matrix/martrix_decompose.hpp"
#include "matrix/real_matrix_base.hpp"
#include "matrix/real_matrix_owned.hpp"

namespace msl::difference {

namespace detail {
/**
 * @brief Reject non-positive uniform spacing.
 */
inline void validate_spacing(double dx) {
    if (!(dx > 0.0)) {
        throw std::invalid_argument("Difference: spacing dx must be positive");
    }
}

/**
 * @brief Reject empty or non-increasing coordinate arrays.
 */
inline void validate_coordinates(std::span<const double> x) {
    if (x.size() < 2) {
        throw std::invalid_argument(
            "Difference: coordinate array must contain at least 2 points");
    }
    for (size_t i = 1; i < x.size(); ++i) {
        if (!(x[i] > x[i - 1])) {
            throw std::invalid_argument(
                "Difference: coordinates must be strictly increasing");
        }
    }
}

/**
 * @brief Three-point central difference for non-uniform coordinates.
 *
 * The derivative at the center point is the exact derivative of the
 * quadratic through the three samples.
 */
inline double nonuniform_central_difference(double x_left,
                                            double x_center,
                                            double x_right,
                                            double y_left,
                                            double y_center,
                                            double y_right) {
    const double dx_left = x_center - x_left;
    const double dx_right = x_right - x_center;
    return (dx_left * (y_right - y_center) / dx_right
            + dx_right * (y_center - y_left) / dx_left)
           / (dx_left + dx_right);
}
} // namespace detail

// ============================================================================
// Difference (First-order)
// ============================================================================

/**
 * @brief First-order difference: diff[i] = y[i+1] - y[i]
 *
 * Zero-copy operation: directly fills the provided result span.
 * Output size is n-1
 *
 * @param y Input values
 * @param result Output buffer for differences (must have size n-1)
 */
inline void diff(std::span<const double> y, std::span<double> result) {
    if (y.size() < 2) {
        throw std::invalid_argument("Diff: need at least 2 points for diff");
    }

    if (result.size() != y.size() - 1) {
        throw std::invalid_argument("Diff: result span must have size n-1");
    }

    for (size_t i = 0; i < result.size(); ++i) {
        result[i] = y[i + 1] - y[i];
    }
}

/**
 * @brief First-order difference: diff[i] = y[i+1] - y[i]
 *
 * Convenience wrapper that allocates and returns a vector.
 * For zero-copy operations, use the void version with span parameter.
 *
 * @param y Input values (can be vector, array, or span)
 * @return Vector of differences (size n-1)
 */
inline std::vector<double> diff(std::span<const double> y) {
    if (y.size() < 2) {
        throw std::invalid_argument("Diff: need at least 2 points for diff");
    }

    std::vector<double> result(y.size() - 1);

    diff(y, result);

    return result;
}

/**
 * @brief First-order difference for matrix (along specified axis)
 *
 * @param mat Input matrix
 * @param axis 0 = row-wise (vertical diff), 1 = column-wise (horizontal diff)
 * @return Difference matrix
 */
inline matrix::matrixd diff(const matrix::real_matrix_base &mat, int axis = 0) {
    if (axis == 0) {
        // Row-wise: diff along rows (vertical)
        if (mat.rows() < 2) {
            throw std::invalid_argument(
                "Diff: need at least 2 rows for row-wise diff");
        }

        matrix::matrixd result(mat.rows() - 1, mat.cols());
        for (size_t j = 0; j < mat.cols(); ++j) {
            for (size_t i = 0; i < result.rows(); ++i) {
                result(i, j) = mat(i + 1, j) - mat(i, j);
            }
        }
        return result;
    } else if (axis == 1) {
        // Column-wise: diff along columns (horizontal)
        if (mat.cols() < 2) {
            throw std::invalid_argument(
                "Diff: need at least 2 cols for column-wise diff");
        }

        matrix::matrixd result(mat.rows(), mat.cols() - 1);
        for (size_t i = 0; i < mat.rows(); ++i) {
            for (size_t j = 0; j < result.cols(); ++j) {
                result(i, j) = mat(i, j + 1) - mat(i, j);
            }
        }
        return result;
    } else {
        throw std::invalid_argument("Diff: axis must be 0 or 1");
    }
}

// ============================================================================
// Gradient (Derivative Approximation)
// ============================================================================

/**
 * @brief Numerical gradient using forward differences (uniform spacing)
 *
 * Zero-copy operation: directly fills the provided gradient span.
 * Uses forward difference approximation.
 * Output size is n-1
 *
 * @param y Function values
 * @param grad Output buffer for gradient values (must have size n-1)
 * @param dx Spacing between points (default: 1.0)
 */
inline void forward_gradient(std::span<const double> y,
                             std::span<double> grad,
                             double dx = 1.0) {
    if (y.size() < 2) {
        throw std::invalid_argument(
            "Gradient: need at least 2 points for gradient");
    }

    if (grad.size() != y.size() - 1) {
        throw std::invalid_argument("Gradient: span must have size n-1");
    }

    detail::validate_spacing(dx);

    // Forward difference, the first point is 0
    for (size_t i = 0; i < y.size() - 1; ++i) {
        grad[i] = (y[i + 1] - y[i]) / dx;
    }
}

/**
 * @brief Numerical gradient using forward differences (uniform spacing)
 *
 * Convenience wrapper that allocates and returns a vector.
 * For zero-copy operations, use the void version with span parameter.
 *
 * @param y Function values
 * @param dx Spacing between points (default: 1.0)
 * @return Vector of gradient values (size n-1)
 */
inline std::vector<double> forward_gradient(std::span<const double> y,
                                            double dx = 1.0) {
    if (y.size() < 2) {
        throw std::invalid_argument(
            "Gradient: need at least 2 points for gradient");
    }

    std::vector<double> grad(y.size() - 1);

    forward_gradient(y, grad, dx);

    return grad;
}

/**
 * @brief Numerical gradient using forward differences (non-uniform spacing)
 *
 * Zero-copy operation: directly fills the provided gradient span.
 * Uses forward difference with non-uniform spacing.
 * Output size is n-1
 *
 * @param y Function values
 * @param grad Output buffer for gradient values (must have size n-1)
 * @param x Independent variable values (must have same size as y)
 */
inline void forward_gradient(std::span<const double> y,
                             std::span<double> grad,
                             std::span<const double> x) {
    if (x.size() != y.size()) {
        throw std::invalid_argument("Gradient: x and y must have same size");
    }
    if (x.size() < 2) {
        throw std::invalid_argument(
            "Gradient: need at least 2 points for gradient");
    }
    if (grad.size() != y.size() - 1) {
        throw std::invalid_argument("Gradient: span must have size n-1");
    }

    detail::validate_coordinates(x);

    // Forward difference, the first point is 0
    for (size_t i = 0; i < y.size() - 1; ++i) {
        double dx_local = x[i + 1] - x[i];
        grad[i] = (y[i + 1] - y[i]) / dx_local;
    }
}

/**
 * @brief Numerical gradient using forward differences (non-uniform spacing)
 *
 * Convenience wrapper that allocates and returns a vector.
 * For zero-copy operations, use the void version with span parameter.
 *
 * @param y Function values
 * @param x Independent variable values (must have same size as y)
 * @return Vector of gradient values (size n-1)
 */
inline std::vector<double> forward_gradient(std::span<const double> y,
                                            std::span<const double> x) {
    if (y.size() < 2) {
        throw std::invalid_argument(
            "Gradient: need at least 2 points for gradient");
    }

    std::vector<double> grad(y.size() - 1);

    forward_gradient(y, grad, x);

    return grad;
}

/**
 * @brief Gradient for matrix (along specified axis), uses forward differences
 *
 * @param mat Input matrix
 * @param dx Spacing
 * @param axis 0 = gradient along rows (vertical, size n-1 x m), 1 = gradient
 * along columns (horizontal, size n x m-1)
 * @return Gradient matrix
 */
inline matrix::matrixd forward_gradient(const matrix::real_matrix_base &mat,
                                        double dx = 1.0,
                                        int axis = 0) {
    detail::validate_spacing(dx);

    if (axis == 0) {
        // Gradient along rows (vertical direction)
        if (mat.rows() < 2) {
            throw std::invalid_argument("Gradient: need at least 2 rows");
        }

        matrix::matrixd grad(mat.rows() - 1, mat.cols());

        // First row gradient is zero
        for (size_t j = 0; j < mat.cols(); ++j) {
            // Forward difference
            for (size_t i = 0; i < mat.rows() - 1; ++i) {
                grad(i, j) = (mat(i + 1, j) - mat(i, j)) / dx;
            }
        }

        return grad;
    } else if (axis == 1) {
        // Gradient along columns (horizontal direction)
        if (mat.cols() < 2) {
            throw std::invalid_argument("Gradient: need at least 2 columns");
        }

        matrix::matrixd grad(mat.rows(), mat.cols() - 1);

        // First column gradient is zero
        for (size_t i = 0; i < mat.rows(); ++i) {
            // Central at interior columns
            for (size_t j = 0; j < mat.cols() - 1; ++j) {
                grad(i, j) = (mat(i, j + 1) - mat(i, j)) / dx;
            }
        }

        return grad;
    } else {
        throw std::invalid_argument("axis must be 0 or 1");
    }
}

/**
 * @brief Gradient for matrix with non-uniform spacing (along specified axis),
 * uses forward differences
 *
 * @param mat Input matrix
 * @param x Spacing values
 * @param axis 0 = gradient along rows (vertical, size n-1 x m), 1 = gradient
 * along columns (horizontal, size n x m-1)
 * @return Gradient matrix
 */
inline matrix::matrixd forward_gradient(const matrix::real_matrix_base &mat,
                                        std::span<const double> x,
                                        int axis = 0) {
    if (axis == 0) {
        // Gradient along rows (vertical direction)
        if (mat.rows() < 2) {
            throw std::invalid_argument("Gradient: need at least 2 rows");
        }

        if (x.size() != mat.rows()) {
            throw std::invalid_argument(
                "Gradient: x size must be number of rows for axis=0");
        }

        detail::validate_coordinates(x);

        matrix::matrixd grad(mat.rows() - 1, mat.cols());

        for (size_t j = 0; j < mat.cols(); ++j) {
            // Forward difference
            for (size_t i = 0; i < mat.rows() - 1; ++i) {
                double dx_local = x[i + 1] - x[i];
                grad(i, j) = (mat(i + 1, j) - mat(i, j)) / dx_local;
            }
        }

        return grad;
    } else if (axis == 1) {
        // Gradient along columns (horizontal direction)
        if (mat.cols() < 2) {
            throw std::invalid_argument("Gradient: need at least 2 columns");
        }

        if (x.size() != mat.cols()) {
            throw std::invalid_argument(
                "Gradient: x size must be number of cols for axis=1");
        }

        detail::validate_coordinates(x);

        matrix::matrixd grad(mat.rows(), mat.cols() - 1);

        for (size_t i = 0; i < mat.rows(); ++i) {
            // Forward difference
            for (size_t j = 0; j < mat.cols() - 1; ++j) {
                double dx_local = x[j + 1] - x[j];
                grad(i, j) = (mat(i, j + 1) - mat(i, j)) / dx_local;
            }
        }

        return grad;
    } else {
        throw std::invalid_argument("axis must be 0 or 1");
    }
}

/**
 * @brief Numerical gradient using central differences (uniform spacing)
 *
 * Zero-copy operation: directly fills the provided gradient span.
 * Uses:
 * - Forward difference at first point
 * - Central difference at interior points
 * - Backward difference at last point
 * Output size equals input size
 *
 * @param y Function values
 * @param grad Output buffer for gradient values (must have same size as y)
 * @param dx Spacing between points (default: 1.0)
 */
inline void central_gradient(std::span<const double> y,
                             std::span<double> grad,
                             double dx = 1.0) {
    if (y.size() < 2) {
        throw std::invalid_argument(
            "Gradient: need at least 2 points for gradient");
    }

    if (grad.size() != y.size()) {
        throw std::invalid_argument(
            "Gradient: span must have same size as input");
    }

    detail::validate_spacing(dx);

    // Forward difference at first point
    grad[0] = (y[1] - y[0]) / dx;

    // Central difference at interior points
    for (size_t i = 1; i < y.size() - 1; ++i) {
        grad[i] = (y[i + 1] - y[i - 1]) / (2.0 * dx);
    }

    // Backward difference at last point
    grad[y.size() - 1] = (y[y.size() - 1] - y[y.size() - 2]) / dx;
}

/**
 * @brief Numerical gradient using central differences (uniform spacing)
 *
 * Convenience wrapper that allocates and returns a vector.
 * For zero-copy operations, use the void version with span parameter.
 *
 * @param y Function values
 * @param dx Spacing between points (default: 1.0)
 * @return Vector of gradient values (same size as input)
 */
inline std::vector<double> central_gradient(std::span<const double> y,
                                            double dx = 1.0) {
    if (y.size() < 2) {
        throw std::invalid_argument(
            "Gradient: need at least 2 points for gradient");
    }

    std::vector<double> grad(y.size());

    central_gradient(y, grad, dx);

    return grad;
}

/**
 * @brief Numerical gradient with non-uniform spacing
 *
 * Zero-copy operation: directly fills the provided gradient span.
 * Uses weighted central differences to handle non-uniform spacing.
 * Output size equals input size
 *
 * @param x Independent variable values
 * @param y Dependent variable values (must have same size as x)
 * @param grad Output buffer for gradient values (must have same size as y)
 */
inline void central_gradient(std::span<const double> x,
                             std::span<const double> y,
                             std::span<double> grad) {
    if (x.size() != y.size()) {
        throw std::invalid_argument("Gradient: x and y must have same size");
    }
    if (x.size() < 2) {
        throw std::invalid_argument(
            "Gradient: need at least 2 points for gradient");
    }
    if (grad.size() != y.size()) {
        throw std::invalid_argument(
            "Gradient: span must have same size as input");
    }

    detail::validate_coordinates(x);

    // Forward difference at first point
    double dx0 = x[1] - x[0];
    grad[0] = (y[1] - y[0]) / dx0;

    // Central difference at interior points
    for (size_t i = 1; i < y.size() - 1; ++i) {
        grad[i] = detail::nonuniform_central_difference(
            x[i - 1], x[i], x[i + 1], y[i - 1], y[i], y[i + 1]);
    }

    // Backward difference at last point
    size_t n = y.size();
    double dx_last = x[n - 1] - x[n - 2];
    grad[n - 1] = (y[n - 1] - y[n - 2]) / dx_last;
}

/**
 * @brief Numerical gradient with non-uniform spacing
 *
 * Convenience wrapper that allocates and returns a vector.
 * For zero-copy operations, use the void version with span parameter.
 *
 * @param x Independent variable values
 * @param y Dependent variable values (must have same size as x)
 * @return Vector of gradient values (same size as input)
 */
inline std::vector<double> central_gradient(std::span<const double> x,
                                            std::span<const double> y) {
    if (y.size() < 2) {
        throw std::invalid_argument(
            "Gradient: need at least 2 points for gradient");
    }

    std::vector<double> grad(y.size());

    central_gradient(x, y, grad);

    return grad;
}

/**
 * @brief Gradient for matrix (along specified axis)
 *
 * @param mat Input matrix
 * @param dx Spacing
 * @param axis 0 = gradient along rows (vertical, size n x m), 1 = gradient
 * along columns (horizontal, size n x m)
 * @return Gradient matrix (same size as input)
 */
inline matrix::matrixd central_gradient(const matrix::real_matrix_base &mat,
                                        double dx = 1.0,
                                        int axis = 0) {
    detail::validate_spacing(dx);

    if (axis == 0) {
        // Gradient along rows (vertical direction)
        if (mat.rows() < 2) {
            throw std::invalid_argument("Gradient: need at least 2 rows");
        }

        matrix::matrixd grad(mat.rows(), mat.cols());

        for (size_t j = 0; j < mat.cols(); ++j) {
            // Forward at first row
            grad(0, j) = (mat(1, j) - mat(0, j)) / dx;

            // Central at interior rows
            for (size_t i = 1; i < mat.rows() - 1; ++i) {
                grad(i, j) = (mat(i + 1, j) - mat(i - 1, j)) / (2.0 * dx);
            }

            // Backward at last row
            size_t n = mat.rows() - 1;
            grad(n, j) = (mat(n, j) - mat(n - 1, j)) / dx;
        }

        return grad;
    } else if (axis == 1) {
        // Gradient along columns (horizontal direction)
        if (mat.cols() < 2) {
            throw std::invalid_argument("Gradient: need at least 2 columns");
        }

        matrix::matrixd grad(mat.rows(), mat.cols());

        for (size_t i = 0; i < mat.rows(); ++i) {
            // Forward at first column
            grad(i, 0) = (mat(i, 1) - mat(i, 0)) / dx;

            // Central at interior columns
            for (size_t j = 1; j < mat.cols() - 1; ++j) {
                grad(i, j) = (mat(i, j + 1) - mat(i, j - 1)) / (2.0 * dx);
            }

            // Backward at last column
            size_t n = mat.cols() - 1;
            grad(i, n) = (mat(i, n) - mat(i, n - 1)) / dx;
        }

        return grad;
    } else {
        throw std::invalid_argument("Gradient: axis must be 0 or 1");
    }
}

/**
 * @brief Gradient for matrix with non-uniform spacing (along specified axis)
 *
 * @param mat Input matrix
 * @param x Spacing values
 * @param axis 0 = gradient along rows (vertical, size n x m), 1 = gradient
 * along columns (horizontal, size n x m)
 * @return Gradient matrix (same size as input)
 */
inline matrix::matrixd central_gradient(const matrix::real_matrix_base &mat,
                                        std::span<const double> x,
                                        int axis = 0) {
    if (axis == 0) {
        // Gradient along rows (vertical direction)
        if (mat.rows() < 2) {
            throw std::invalid_argument("Gradient: need at least 2 rows");
        }

        if (x.size() != mat.rows()) {
            throw std::invalid_argument(
                "Gradient: x size must be number of rows for axis=0");
        }

        detail::validate_coordinates(x);

        matrix::matrixd grad(mat.rows(), mat.cols());

        for (size_t j = 0; j < mat.cols(); ++j) {
            // Forward at first row
            grad(0, j) = (mat(1, j) - mat(0, j)) / (x[1] - x[0]);

            // Central at interior rows
            for (size_t i = 1; i < mat.rows() - 1; ++i) {
                grad(i, j) =
                    detail::nonuniform_central_difference(x[i - 1],
                                                          x[i],
                                                          x[i + 1],
                                                          mat(i - 1, j),
                                                          mat(i, j),
                                                          mat(i + 1, j));
            }

            // Backward at last row
            size_t n = mat.rows() - 1;
            grad(n, j) = (mat(n, j) - mat(n - 1, j)) / (x[n] - x[n - 1]);
        }

        return grad;
    } else if (axis == 1) {
        // Gradient along columns (horizontal direction)
        if (mat.cols() < 2) {
            throw std::invalid_argument("Gradient: need at least 2 columns");
        }

        if (x.size() != mat.cols()) {
            throw std::invalid_argument(
                "Gradient: x size must be number of cols for axis=1");
        }

        detail::validate_coordinates(x);

        matrix::matrixd grad(mat.rows(), mat.cols());

        for (size_t i = 0; i < mat.rows(); ++i) {
            // Forward at first column
            grad(i, 0) = (mat(i, 1) - mat(i, 0)) / (x[1] - x[0]);

            // Central at interior columns
            for (size_t j = 1; j < mat.cols() - 1; ++j) {
                grad(i, j) =
                    detail::nonuniform_central_difference(x[j - 1],
                                                          x[j],
                                                          x[j + 1],
                                                          mat(i, j - 1),
                                                          mat(i, j),
                                                          mat(i, j + 1));
            }

            // Backward at last column
            size_t n = mat.cols() - 1;
            grad(i, n) = (mat(i, n) - mat(i, n - 1)) / (x[n] - x[n - 1]);
        }

        return grad;
    } else {
        throw std::invalid_argument("Gradient: axis must be 0 or 1");
    }
}

// ============================================================================
// 2D Gradient (returns both directions)
// ============================================================================

/**
 * @brief Compute 2D gradient (both directions)
 *
 * Returns {dy/dx (horizontal), dy/dy (vertical)}
 *
 * @param mat Input matrix
 * @param dx Horizontal spacing
 * @param dy Vertical spacing
 * @return Pair of gradient matrices {grad_x, grad_y}
 */
inline std::pair<matrix::matrixd, matrix::matrixd>
central_gradient2d(const matrix::real_matrix_base &mat,
                   double dx = 1.0,
                   double dy = 1.0) {
    if (mat.rows() < 2 || mat.cols() < 2) {
        throw std::invalid_argument(
            "Gradient: need at least 2x2 matrix for 2D gradient");
    }

    matrix::matrixd grad_x =
        central_gradient(mat, dx, 1); // Horizontal gradient
    matrix::matrixd grad_y = central_gradient(mat, dy, 0); // Vertical gradient

    return {grad_x, grad_y};
}

// ============================================================================
// Higher-order Derivatives
// ============================================================================

/**
 * @brief Second derivative using central differences
 *
 * Zero-copy operation: directly fills the provided gradient span.
 * d²y/dx² ≈ (y[i+1] - 2*y[i] + y[i-1]) / dx²
 * Output size equals input size
 *
 * @param y Function values
 * @param grad2 Output buffer for second derivative values (must have same size
 * as y)
 * @param dx Spacing between points (default: 1.0)
 */
inline void central_gradient2(std::span<const double> y,
                              std::span<double> grad2,
                              double dx = 1.0) {
    if (y.size() < 3) {
        throw std::invalid_argument(
            "Gradient: need at least 3 points for second derivative");
    }

    if (grad2.size() != y.size()) {
        throw std::invalid_argument(
            "Gradient: second derivative span must have same size as input");
    }

    detail::validate_spacing(dx);

    double dx2 = dx * dx;

    // Forward difference at first point (less accurate)
    grad2[0] = (y[2] - 2.0 * y[1] + y[0]) / dx2;

    // Central difference at interior points
    for (size_t i = 1; i < y.size() - 1; ++i) {
        grad2[i] = (y[i + 1] - 2.0 * y[i] + y[i - 1]) / dx2;
    }

    // Backward difference at last point
    size_t n = y.size() - 1;
    grad2[n] = (y[n] - 2.0 * y[n - 1] + y[n - 2]) / dx2;
}

/**
 * @brief Second derivative using central differences
 *
 * Convenience wrapper that allocates and returns a vector.
 * For zero-copy operations, use the void version with span parameter.
 *
 * @param y Function values
 * @param dx Spacing between points (default: 1.0)
 * @return Vector of second derivative values (same size as input)
 */
inline std::vector<double> central_gradient2(std::span<const double> y,
                                             double dx = 1.0) {
    std::vector<double> grad2(y.size());
    central_gradient2(y, grad2, dx);
    return grad2;
}

// ============================================================================
// Laplacian (for 2D data)
// ============================================================================

/**
 * @brief Compute Laplacian (∇²f = ∂²f/∂x² + ∂²f/∂y²)
 *
 * @param mat Input matrix
 * @param dx Horizontal spacing
 * @param dy Vertical spacing
 * @return Laplacian matrix (same size, boundaries set to 0)
 */
inline matrix::matrixd laplacian(const matrix::real_matrix_base &mat,
                                 double dx = 1.0,
                                 double dy = 1.0) {
    if (mat.rows() < 3 || mat.cols() < 3) {
        throw std::invalid_argument(
            "Gradient: need at least 3x3 matrix for Laplacian");
    }

    detail::validate_spacing(dx);
    detail::validate_spacing(dy);

    matrix::matrixd result(mat.rows(), mat.cols(), 0.0);

    double dx2 = dx * dx;
    double dy2 = dy * dy;

    // Interior points only
    for (size_t i = 1; i < mat.rows() - 1; ++i) {
        for (size_t j = 1; j < mat.cols() - 1; ++j) {
            // ∂²f/∂x²
            double d2_dx2 =
                (mat(i, j + 1) - 2.0 * mat(i, j) + mat(i, j - 1)) / dx2;

            // ∂²f/∂y²
            double d2_dy2 =
                (mat(i + 1, j) - 2.0 * mat(i, j) + mat(i - 1, j)) / dy2;

            result(i, j) = d2_dx2 + d2_dy2;
        }
    }

    return result;
}

// ============================================================================
// Divergence and Curl (for vector fields)
// ============================================================================

/**
 * @brief Compute divergence of 2D vector field
 *
 * div(F) = ∂Fx/∂x + ∂Fy/∂y
 *
 * @param Fx X-component of vector field
 * @param Fy Y-component of vector field
 * @param dx Horizontal spacing
 * @param dy Vertical spacing
 * @return Divergence (scalar field)
 */
inline matrix::matrixd divergence(const matrix::real_matrix_base &Fx,
                                  const matrix::real_matrix_base &Fy,
                                  double dx = 1.0,
                                  double dy = 1.0) {
    if (Fx.rows() != Fy.rows() || Fx.cols() != Fy.cols()) {
        throw std::invalid_argument("Gradient: Fx and Fy must have same size");
    }

    auto dFx_dx = central_gradient(Fx, dx, 1); // ∂Fx/∂x
    auto dFy_dy = central_gradient(Fy, dy, 0); // ∂Fy/∂y

    matrix::matrixd div(Fx.rows(), Fx.cols());
    for (size_t i = 0; i < div.rows(); ++i) {
        for (size_t j = 0; j < div.cols(); ++j) {
            div(i, j) = dFx_dx(i, j) + dFy_dy(i, j);
        }
    }

    return div;
}

/**
 * @brief Compute curl of 2D vector field
 *
 * curl(F) = ∂Fy/∂x - ∂Fx/∂y
 *
 * @param Fx X-component of vector field
 * @param Fy Y-component of vector field
 * @param dx Horizontal spacing
 * @param dy Vertical spacing
 * @return Curl (scalar field, z-component)
 */
inline matrix::matrixd curl(const matrix::real_matrix_base &Fx,
                            const matrix::real_matrix_base &Fy,
                            double dx = 1.0,
                            double dy = 1.0) {
    if (Fx.rows() != Fy.rows() || Fx.cols() != Fy.cols()) {
        throw std::invalid_argument("Gradient: Fx and Fy must have same size");
    }

    auto dFy_dx = central_gradient(Fy, dx, 1); // ∂Fy/∂x
    auto dFx_dy = central_gradient(Fx, dy, 0); // ∂Fx/∂y

    matrix::matrixd curl_z(Fx.rows(), Fx.cols());
    for (size_t i = 0; i < curl_z.rows(); ++i) {
        for (size_t j = 0; j < curl_z.cols(); ++j) {
            curl_z(i, j) = dFy_dx(i, j) - dFx_dy(i, j);
        }
    }

    return curl_z;
}

// ============================================================================
// Savitzky-Golay Filter for Smoothed Derivatives
// ============================================================================

/**
 * @brief Savitzky-Golay derivative (smoothed, for noisy data)
 *
 * Zero-copy operation: directly fills the provided gradient span.
 * Fits a polynomial of degree `poly_order` to a moving window of
 * `window_size` samples and evaluates its derivative. Interior points use the
 * classic symmetric Savitzky-Golay weights; boundary points use a
 * least-squares fit over the first/last window, matching MATLAB
 * `sgolayfilt` boundary handling. Output size equals input size.
 *
 * @param y Function values
 * @param grad Output buffer for smoothed gradient values (must have same size
 * as y)
 * @param window_size Window size (must be odd, >= 5)
 * @param poly_order Polynomial order (< window_size)
 * @param dx Spacing
 */
inline void savgol_gradient(std::span<const double> y,
                            std::span<double> grad,
                            int window_size = 5,
                            int poly_order = 2,
                            double dx = 1.0) {
    if (window_size % 2 == 0) {
        throw std::invalid_argument("Gradient: window_size must be odd");
    }
    if (poly_order >= window_size) {
        throw std::invalid_argument(
            "Gradient: poly_order must be < window_size");
    }
    if (y.size() < static_cast<size_t>(window_size)) {
        throw std::invalid_argument(
            "Gradient: signal too short for window_size");
    }
    if (grad.size() != y.size()) {
        throw std::invalid_argument(
            "Gradient: span must have same size as input");
    }

    const int half_window = window_size / 2;
    const int n = static_cast<int>(y.size());

    if (poly_order == 0) {
        std::fill(grad.begin(), grad.end(), 0.0);
        return;
    }

    // Build the least-squares derivative weights for a given set of local
    // offsets: the derivative at offset 0 of the polynomial fit is the
    // second row of the pseudoinverse of the Vandermonde matrix.
    const auto derivative_weights =
        [poly_order](std::span<const double> offsets) {
            const int points = static_cast<int>(offsets.size());
            const int order = std::min(poly_order, points - 1);
            matrix::matrixd vandermonde(points, order + 1);
            for (int i = 0; i < points; ++i) {
                double power = 1.0;
                for (int j = 0; j <= order; ++j) {
                    vandermonde(i, j) = power;
                    power *= offsets[static_cast<size_t>(i)];
                }
            }
            const auto pseudo_inverse = matrix::pinv(vandermonde);
            std::vector<double> weights(static_cast<size_t>(points));
            for (int i = 0; i < points; ++i) {
                weights[static_cast<size_t>(i)] = pseudo_inverse(1, i);
            }
            return weights;
        };

    // Interior points: one symmetric window template reused for all samples
    std::vector<double> interior_offsets(static_cast<size_t>(window_size));
    for (int k = -half_window; k <= half_window; ++k) {
        interior_offsets[static_cast<size_t>(k + half_window)] =
            static_cast<double>(k);
    }
    const auto interior_weights = derivative_weights(interior_offsets);

    // Boundary points: first/last full window evaluated at each sample
    std::vector<std::vector<double>> left_weights(
        static_cast<size_t>(half_window));
    std::vector<std::vector<double>> right_weights(
        static_cast<size_t>(half_window));
    for (int i = 0; i < half_window; ++i) {
        std::vector<double> offsets(static_cast<size_t>(window_size));
        for (int j = 0; j < window_size; ++j) {
            offsets[static_cast<size_t>(j)] = static_cast<double>(j - i);
        }
        left_weights[static_cast<size_t>(i)] = derivative_weights(offsets);
        for (int j = 0; j < window_size; ++j) {
            offsets[static_cast<size_t>(j)] =
                static_cast<double>(j - (window_size - 1) + i);
        }
        right_weights[static_cast<size_t>(i)] = derivative_weights(offsets);
    }

    for (int i = 0; i < n; ++i) {
        double value = 0.0;
        if (i >= half_window && i + half_window < n) {
            for (int k = -half_window; k <= half_window; ++k) {
                value += interior_weights[static_cast<size_t>(k + half_window)]
                         * y[static_cast<size_t>(i + k)];
            }
        } else if (i < half_window) {
            for (int j = 0; j < window_size; ++j) {
                value +=
                    left_weights[static_cast<size_t>(i)][static_cast<size_t>(j)]
                    * y[static_cast<size_t>(j)];
            }
        } else {
            const int start = n - window_size;
            for (int j = 0; j < window_size; ++j) {
                value += right_weights[static_cast<size_t>(n - 1 - i)]
                                      [static_cast<size_t>(j)]
                         * y[static_cast<size_t>(start + j)];
            }
        }
        grad[static_cast<size_t>(i)] = value / dx;
    }
}

/**
 * @brief Savitzky-Golay derivative (smoothed, for noisy data)
 *
 * Convenience wrapper that allocates and returns a vector.
 * For zero-copy operations, use the void version with span parameter.
 *
 * @param y Function values
 * @param window_size Window size (must be odd, >= 5)
 * @param poly_order Polynomial order (< window_size)
 * @param dx Spacing
 * @return Vector of smoothed gradient values (same size as input)
 */
inline std::vector<double> savgol_gradient(std::span<const double> y,
                                           int window_size = 5,
                                           int poly_order = 2,
                                           double dx = 1.0) {
    std::vector<double> grad(y.size());

    savgol_gradient(y, grad, window_size, poly_order, dx);

    return grad;
}

} // namespace msl::difference

#endif // MSL_DIFFERENCE_HPP
