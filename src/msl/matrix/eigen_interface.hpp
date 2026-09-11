/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: eigen_interface.hpp
** -----
** File Created: Saturday, 11th October 2025 23:14:26
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Sunday, 14th December 2025 17:04:27
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_EIGEN_INTERFACE_HPP
#define MSL_EIGEN_INTERFACE_HPP

#include <complex>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <stdexcept>

#include "matrix/complex_matrix_base.hpp"
#include "matrix/complex_matrix_owned.hpp"
#include "matrix/real_matrix_owned.hpp"
#include "matrix/real_matrix_view.hpp"

namespace msl::matrix::eigen_interface {

// ============================================================================
// 1. MSL -> Eigen (Zero-copy views using Eigen::Map)
// ============================================================================

/**
 * @brief Create Eigen::Map from matrixd (non-owning view)
 * @note Zero-copy operation. The Eigen::Map is only valid while mat exists.
 * @warning mat must be in contiguous column-major layout
 */
inline Eigen::Map<Eigen::MatrixXd> as_eigen(matrixd &mat) {
    return Eigen::Map<Eigen::MatrixXd>(mat.data(), mat.rows(), mat.cols());
}

/**
 * @brief Create const Eigen::Map from const matrixd
 */
inline Eigen::Map<const Eigen::MatrixXd> as_eigen(const matrixd &mat) {
    return Eigen::Map<const Eigen::MatrixXd>(
        mat.data(), mat.rows(), mat.cols());
}

/**
 * @brief Create Eigen::Map from matrixd_view
 * @warning The view must point to contiguous column-major data
 */
inline Eigen::Map<Eigen::MatrixXd> as_eigen(matrixd_view &view) {
    return Eigen::Map<Eigen::MatrixXd>(view.data(), view.rows(), view.cols());
}

/**
 * @brief Create const Eigen::Map from const matrixd_view
 */
inline Eigen::Map<const Eigen::MatrixXd> as_eigen(const matrixd_view &view) {
    return Eigen::Map<const Eigen::MatrixXd>(
        view.data(), view.rows(), view.cols());
}

/**
 * @brief Create const Eigen::Map from const_matrixd_view
 */
inline Eigen::Map<const Eigen::MatrixXd>
as_eigen(const const_matrixd_view &view) {
    return Eigen::Map<const Eigen::MatrixXd>(
        view.data(), view.rows(), view.cols());
}

/**
 * @brief Create Eigen::Map from real_matrix_base
 * @warning The view must point to contiguous column-major data
 */
inline Eigen::Map<const Eigen::MatrixXd> as_eigen(real_matrix_base &base) {
    return Eigen::Map<const Eigen::MatrixXd>(
        base.data(), base.rows(), base.cols());
}

/**
 * @brief Create const Eigen::Map from const real_matrix_base
 * @warning The view must point to contiguous column-major data
 */
inline Eigen::Map<const Eigen::MatrixXd>
as_eigen(const real_matrix_base &base) {
    return Eigen::Map<const Eigen::MatrixXd>(
        base.data(), base.rows(), base.cols());
}

/**
 * @brief Create Eigen::Map from complex_matrix_base
 * @warning The view must point to contiguous column-major data
 */
inline Eigen::Map<
    const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>>
as_eigen(complex_matrix_base &base) {
    return Eigen::Map<const Eigen::Matrix<std::complex<double>,
                                          Eigen::Dynamic,
                                          Eigen::Dynamic>>(
        base.data(), base.rows(), base.cols());
}

/**
 * @brief Create Eigen::Map from complex_matrix_base
 * @warning The view must point to contiguous column-major data
 */
inline Eigen::Map<
    const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>>
as_eigen(const complex_matrix_base &base) {
    return Eigen::Map<const Eigen::Matrix<std::complex<double>,
                                          Eigen::Dynamic,
                                          Eigen::Dynamic>>(
        base.data(), base.rows(), base.cols());
}

// ============================================================================
// 2. Eigen -> MSL (Deep copy - safe)
// ============================================================================

/**
 * @brief Create matrixd from any Eigen matrix expression (deep copy)
 * @param eig Any Eigen matrix type or expression
 * @return New matrixd owning the data
 * @note This always performs a deep copy, ensuring safety
 */
template <typename Derived>
inline matrixd from_eigen(const Eigen::MatrixBase<Derived> &eig) {
    // Force evaluation of expression templates
    Eigen::MatrixXd evaluated = eig;
    matrixd result(evaluated.rows(), evaluated.cols());
    std::copy(
        evaluated.data(), evaluated.data() + evaluated.size(), result.data());
    return result;
}

/**
 * @brief Specialization for Eigen::MatrixXd (may avoid one copy)
 */
inline matrixd from_eigen(const Eigen::MatrixXd &eig) {
    matrixd result(eig.rows(), eig.cols());
    std::copy(eig.data(), eig.data() + eig.size(), result.data());
    return result;
}

/**
 * @brief Specialization for Eigen::MatrixXd (may avoid one copy)
 */
inline matrixc from_eigen(
    const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>
        &eig) {
    matrixc result(eig.rows(), eig.cols());
    std::copy(eig.data(), eig.data() + eig.size(), result.data());
    return result;
}

/**
 * @brief Create matrixc from any Eigen complex matrix expression (deep copy)
 */
inline matrixd from_eigen(
    const Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>
        &eig) {
    matrixd result(eig.rows(), eig.cols());
    std::copy(eig.data(), eig.data() + eig.size(), result.data());
    return result;
}

/**
 * @brief Create matrixc from any Eigen complex matrix expression (deep copy)
 */
inline matrixc
from_eigen(const Eigen::Map<
           Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>>
               &eig) {
    matrixc result(eig.rows(), eig.cols());
    std::copy(eig.data(), eig.data() + eig.size(), result.data());
    return result;
}

// ============================================================================
// 3. Eigen -> MSL View (Use with extreme caution!)
// ============================================================================

/**
 * @brief Create matrixd_view from Eigen::MatrixXd
 *
 * @warning DANGEROUS: The view is only valid while the Eigen matrix exists!
 * @warning Only use for named lvalue Eigen matrices, never temporaries
 *
 * Safe usage:
 *   Eigen::MatrixXd A(10, 10);
 *   auto view = view_from_eigen(A);  // OK: A is a named variable
 *
 * UNSAFE usage:
 *   auto view = view_from_eigen(get_eigen_matrix());  // DELETED by
 * overload auto view = view_from_eigen(A + B);               // DELETED by
 * overload
 *
 * @param eig Named Eigen matrix (lvalue reference only)
 * @return Non-owning view into Eigen's data
 */
inline matrixd_view view_from_eigen(Eigen::MatrixXd &eig) {
    return matrixd_view(eig.data(), eig.rows(), eig.cols());
}

/**
 * @brief Create const view from const Eigen::MatrixXd
 */
inline const_matrixd_view view_from_eigen(const Eigen::MatrixXd &eig) {
    return const_matrixd_view(eig.data(), eig.rows(), eig.cols());
}

// Delete rvalue overloads to prevent dangling references from temporaries
matrixd_view view_from_eigen(Eigen::MatrixXd &&) = delete;
const_matrixd_view view_from_eigen(const Eigen::MatrixXd &&) = delete;

// Also delete for expression templates
template <typename Derived>
matrixd_view view_from_eigen(Eigen::MatrixBase<Derived> &&) = delete;

template <typename Derived>
const_matrixd_view
view_from_eigen(const Eigen::MatrixBase<Derived> &&) = delete;

// ============================================================================
// 4. In-place copy operations
// ============================================================================

/**
 * @brief Copy from any Eigen expression to existing MSL matrix
 * @param dest Destination MSL matrix (will be resized if needed)
 * @param eig Source Eigen matrix or expression
 */
template <typename Derived>
inline void copy_from_eigen(matrixd &dest,
                            const Eigen::MatrixBase<Derived> &eig) {
    if (dest.rows() != static_cast<size_t>(eig.rows())
        || dest.cols() != static_cast<size_t>(eig.cols())) {
        dest.resize(eig.rows(), eig.cols());
    }

    // Force evaluation for expression templates
    Eigen::MatrixXd evaluated = eig;
    std::copy(
        evaluated.data(), evaluated.data() + evaluated.size(), dest.data());
}

/**
 * @brief Copy from MSL matrix to Eigen matrix
 * @param src Source MSL matrix
 * @param dest Destination Eigen matrix (will be resized if needed)
 */
inline void copy_to_eigen(const matrixd &src, Eigen::MatrixXd &dest) {
    dest.resize(src.rows(), src.cols());
    std::copy(src.data(), src.data() + src.size(), dest.data());
}

/**
 * @brief Copy from MSL view to Eigen matrix
 */
inline void copy_to_eigen(const matrixd_view &src, Eigen::MatrixXd &dest) {
    dest.resize(src.rows(), src.cols());
    std::copy(src.data(), src.data() + src.size(), dest.data());
}

/**
 * @brief Copy from const MSL view to Eigen matrix
 */
inline void copy_to_eigen(const const_matrixd_view &src,
                          Eigen::MatrixXd &dest) {
    dest.resize(src.rows(), src.cols());
    std::copy(src.data(), src.data() + src.size(), dest.data());
}

// ============================================================================
// 5. Convenience wrappers for common operations
// ============================================================================

/**
 * @brief Perform Eigen operation on MSL matrix and return result
 *
 * Example:
 *   matrixd A(10, 10);
 *   auto B = with_eigen(A, [](auto& eig) {
 *       return eig.transpose() * eig;  // Eigen operations
 *   });
 *
 * @param mat Input MSL matrix
 * @param func Lambda taking Eigen::Map reference
 * @return New matrixd with result
 */
template <typename Func>
inline matrixd with_eigen(matrixd &mat, Func &&func) {
    auto eig_map = as_eigen(mat);
    auto result = func(eig_map);
    return from_eigen(result);
}

/**
 * @brief Const version
 */
template <typename Func>
inline matrixd with_eigen(const matrixd &mat, Func &&func) {
    auto eig_map = as_eigen(mat);
    auto result = func(eig_map);
    return from_eigen(result);
}

/**
 * @brief In-place Eigen operation on MSL matrix
 *
 * Example:
 *   matrixd A(10, 10);
 *   inplace_eigen(A, [](auto& eig) {
 *       eig.array() *= 2.0;  // Modify in-place
 *   });
 *
 * @param mat MSL matrix to modify
 * @param func Lambda taking Eigen::Map reference
 */
template <typename Func>
inline void inplace_eigen(matrixd &mat, Func &&func) {
    auto eig_map = as_eigen(mat);
    func(eig_map);
}

// ============================================================================
// 6. Common operations (convenience functions)
// ============================================================================

/**
 * @brief Matrix-matrix multiplication using Eigen
 */
inline matrixd matmul(const matrixd &A, const matrixd &B) {
    if (A.cols() != B.rows()) {
        throw std::invalid_argument("matmul: inner dimensions must agree");
    }
    auto eig_A = as_eigen(A);
    auto eig_B = as_eigen(B);
    return from_eigen(eig_A * eig_B);
}

/**
 * @brief Matrix-vector multiplication
 */
inline std::vector<double> matvec(const matrixd &A,
                                  const std::vector<double> &x) {
    if (A.cols() != x.size()) {
        throw std::invalid_argument(
            "matvec: matrix columns must match vector size");
    }

    auto eig_A = as_eigen(A);
    Eigen::Map<const Eigen::VectorXd> eig_x(x.data(), x.size());
    Eigen::VectorXd result = eig_A * eig_x;

    return std::vector<double>(result.data(), result.data() + result.size());
}

/**
 * @brief Solve linear system Ax = b
 */
inline std::vector<double> solve(const matrixd &A,
                                 const std::vector<double> &b) {
    if (A.rows() != A.cols()) {
        throw std::invalid_argument("solve: matrix must be square");
    }
    if (A.rows() != b.size()) {
        throw std::invalid_argument(
            "solve: right-hand side size must match matrix rows");
    }

    auto eig_A = as_eigen(A);
    Eigen::Map<const Eigen::VectorXd> eig_b(b.data(), b.size());

    Eigen::VectorXd x = eig_A.lu().solve(eig_b);

    return std::vector<double>(x.data(), x.data() + x.size());
}

} // namespace msl::matrix::eigen_interface

#endif // MSL_EIGEN_INTERFACE_HPP