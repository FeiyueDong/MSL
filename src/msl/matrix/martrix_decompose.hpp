/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: martrix_decompose.hpp
** -----
** File Created: Friday, 7th November 2025 21:22:35
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Sunday, 14th December 2025 17:04:31
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_MATRIX_DECOMPOSE_HPP
#define MSL_MATRIX_DECOMPOSE_HPP

#include <array>
#include <complex>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/SVD>
#include <limits>
#include <stdexcept>
#include <vector>

#include "complex_matrix_base.hpp"
#include "complex_matrix_owned.hpp"
#include "eigen_interface.hpp"
#include "real_matrix_base.hpp"
#include "real_matrix_owned.hpp"

namespace msl::matrix {
// ============================================================================
// 1. Matrix SVD Functions
// ============================================================================
struct truncated_svd_result {
    real_matrix_owned U;
    std::vector<double> singular_values;
    real_matrix_owned V;
};

// Decompose real matrix A into U, S, V such that A = U * S * V^T
// U and V are orthogonal matrices, S is diagonal matrix of singular values
inline std::array<real_matrix_owned, 3> svd(const real_matrix_base &A) {
    // Placeholder implementation (to be replaced with actual SVD algorithm)
    size_t m = A.rows();
    size_t n = A.cols();
    std::array<real_matrix_owned, 3> result;
    auto &U = result[0] = real_matrix_owned(m, m); // Orthogonal matrix
    auto &S = result[1] = real_matrix_owned(m, n); // Diagonal matrix
    auto &V = result[2] = real_matrix_owned(n, n); // Orthogonal matrix

    auto eig_A = eigen_interface::as_eigen(A);
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(
        eig_A, Eigen::ComputeFullU | Eigen::ComputeFullV);

    Eigen::MatrixXd eig_U = svd.matrixU();
    Eigen::MatrixXd eig_V = svd.matrixV();
    Eigen::MatrixXd eig_S = svd.singularValues();

    // Copy
    std::copy(eig_U.data(), eig_U.data() + eig_U.size(), U.data());
    std::copy(eig_V.data(), eig_V.data() + eig_V.size(), V.data());
    std::fill(S.data(), S.data() + S.size(), 0.0);
    for (size_t i = 0; i < std::min(m, n); ++i) {
        S(i, i) = eig_S(i);
    }

    return result;
}

// Decompose complex matrix A into U, S, V such that A = U * S * V^T
// U and V are orthogonal matrices, S is diagonal matrix of singular values
inline std::array<complex_matrix_owned, 3> svd(const complex_matrix_base &A) {
    // Placeholder implementation (to be replaced with actual SVD algorithm)
    size_t m = A.rows();
    size_t n = A.cols();
    std::array<complex_matrix_owned, 3> result;
    auto &U = result[0] = complex_matrix_owned(m, m); // Orthogonal matrix
    auto &S = result[1] = complex_matrix_owned(m, n); // Diagonal matrix
    auto &V = result[2] = complex_matrix_owned(n, n); // Orthogonal matrix

    auto eig_A = eigen_interface::as_eigen(A);
    Eigen::JacobiSVD<
        Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>>
        svd(eig_A, Eigen::ComputeFullU | Eigen::ComputeFullV);

    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> eig_U =
        svd.matrixU();
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> eig_V =
        svd.matrixV();
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> eig_S =
        svd.singularValues();

    // Copy
    std::copy(eig_U.data(), eig_U.data() + eig_U.size(), U.data());
    std::copy(eig_V.data(), eig_V.data() + eig_V.size(), V.data());
    std::fill(S.data(), S.data() + S.size(), 0.0);
    for (size_t i = 0; i < std::min(m, n); ++i) {
        S(i, i) = eig_S(i);
    }

    return result;
}

/**
 * @brief Compute the leading singular triplets of a real matrix.
 *
 * The returned matrices satisfy `A ~= U * diag(singular_values) * V^T`.
 * `V` is empty when `compute_right_vectors` is false.
 */
inline truncated_svd_result truncated_svd(const real_matrix_base &A,
                                          size_t rank,
                                          bool compute_right_vectors = true) {
    const size_t max_rank = std::min(A.rows(), A.cols());
    if (rank == 0 || rank > max_rank) {
        throw std::invalid_argument(
            "Truncated SVD: rank must be in [1, min(rows, cols)]");
    }

    const auto eig_A = eigen_interface::as_eigen(A);
    unsigned int options = Eigen::ComputeThinU;
    if (compute_right_vectors) {
        options |= Eigen::ComputeThinV;
    }
    Eigen::JacobiSVD<Eigen::MatrixXd> decomposition(eig_A, options);

    truncated_svd_result result;
    result.U = eigen_interface::from_eigen(
        decomposition.matrixU().leftCols(static_cast<Eigen::Index>(rank)));
    result.singular_values.resize(rank);
    for (size_t i = 0; i < rank; ++i) {
        result.singular_values[i] =
            decomposition.singularValues()(static_cast<Eigen::Index>(i));
    }
    if (compute_right_vectors) {
        result.V = eigen_interface::from_eigen(
            decomposition.matrixV().leftCols(static_cast<Eigen::Index>(rank)));
    }
    return result;
}

namespace detail {
inline real_matrix_owned pinv_impl(const real_matrix_base &A,
                                   double tolerance,
                                   bool use_default_tolerance) {
    const auto eig_A = eigen_interface::as_eigen(A);
    Eigen::JacobiSVD<Eigen::MatrixXd> decomposition(
        eig_A, Eigen::ComputeThinU | Eigen::ComputeThinV);

    const auto &singular_values = decomposition.singularValues();
    if (use_default_tolerance) {
        const double largest =
            singular_values.size() == 0 ? 0.0 : singular_values(0);
        tolerance = static_cast<double>(std::max(A.rows(), A.cols()))
                    * std::numeric_limits<double>::epsilon() * largest;
    }

    Eigen::VectorXd inverse_values = singular_values;
    for (Eigen::Index i = 0; i < inverse_values.size(); ++i) {
        inverse_values(i) =
            singular_values(i) > tolerance ? 1.0 / singular_values(i) : 0.0;
    }

    return eigen_interface::from_eigen(decomposition.matrixV()
                                       * inverse_values.asDiagonal()
                                       * decomposition.matrixU().adjoint());
}
} // namespace detail

/**
 * @brief Compute the Moore-Penrose pseudoinverse of a real matrix.
 *
 * Singular values less than or equal to
 * `max(rows, cols) * epsilon * largest_singular_value` are discarded.
 */
inline real_matrix_owned pinv(const real_matrix_base &A) {
    return detail::pinv_impl(A, 0.0, true);
}

/**
 * @brief Compute the Moore-Penrose pseudoinverse with an explicit tolerance.
 */
inline real_matrix_owned pinv(const real_matrix_base &A, double tolerance) {
    if (tolerance < 0.0) {
        throw std::invalid_argument(
            "Pseudoinverse: tolerance must be non-negative");
    }
    return detail::pinv_impl(A, tolerance, false);
}

// ============================================================================
// 2. Matrix Eigenvalue Functions
// ============================================================================
// Decompose real matrix A into V, D such that A * V = V * D;
// V is matrix of eigenvectors, D is diagonal matrix of eigenvalues
inline std::array<msl::matrix::matrixc, 2> eig(const real_matrix_base &A) {
    if (A.rows() != A.cols()) {
        throw std::invalid_argument(
            "Eigenvalue decomposition requires a square matrix");
    }

    // Placeholder implementation (to be replaced with actual Eigenvalue
    // algorithm)
    size_t n = A.rows();
    std::array<msl::matrix::matrixc, 2> result;
    auto &V = result[0] = msl::matrix::matrixc(n, n); // Eigenvector matrix
    auto &D = result[1] =
        msl::matrix::matrixc(n, n); // Diagonal eigenvalue matrix

    auto eig_A = eigen_interface::as_eigen(A);
    Eigen::EigenSolver<Eigen::MatrixXd> eig(eig_A);

    auto eig_V = eig.eigenvectors();
    auto &eig_D = eig.eigenvalues();

    // Copy
    std::copy(eig_V.data(), eig_V.data() + eig_V.size(), V.data());
    std::fill(D.data(), D.data() + D.size(), 0.0);
    for (size_t i = 0; i < n; ++i) {
        D(i, i) = eig_D(i);
    }

    return result;
}

// Decompose real matrix A and B into V, D such that A * V = B * V * D;
// V is matrix of eigenvectors, D is diagonal matrix of eigenvalues
inline std::array<msl::matrix::matrixc, 2> eig(const real_matrix_base &A,
                                               const real_matrix_base &B) {
    if (A.rows() != A.cols() || B.rows() != B.cols() || A.rows() != B.rows()) {
        throw std::invalid_argument(
            "Generalized eigenvalue decomposition requires square matrices of "
            "the same size");
    }

    // Placeholder implementation (to be replaced with actual Eigenvalue
    // algorithm)
    size_t n = A.rows();
    std::array<msl::matrix::matrixc, 2> result;
    auto &V = result[0] = msl::matrix::matrixc(n, n); // Eigenvector matrix
    auto &D = result[1] =
        msl::matrix::matrixc(n, n); // Diagonal eigenvalue matrix

    auto eig_A = eigen_interface::as_eigen(A);
    auto eig_B = eigen_interface::as_eigen(B);
    Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> eig(eig_A, eig_B);

    auto eig_V = eig.eigenvectors();
    auto eig_D = eig.eigenvalues();

    // Copy
    std::copy(eig_V.data(), eig_V.data() + eig_V.size(), V.data());
    std::fill(D.data(), D.data() + D.size(), 0.0);
    for (size_t i = 0; i < n; ++i) {
        D(i, i) = eig_D(i);
    }

    return result;
}

// ============================================================================
// 3. Matrix LU Decomposition Functions
// ============================================================================
// Decompose real matrix A into L and U such that A = L * U
// L is lower triangular matrix, U is upper triangular matrix
inline std::array<real_matrix_owned, 2> lu(const real_matrix_base &A) {
    if (A.rows() != A.cols()) {
        throw std::invalid_argument(
            "LU decomposition requires a square matrix");
    }

    size_t n = A.rows();
    std::array<real_matrix_owned, 2> result;
    auto &L = result[0] = real_matrix_owned(n, n); // Lower triangular matrix
    auto &U = result[1] = real_matrix_owned(n, n); // Upper triangular matrix

    auto eig_A = eigen_interface::as_eigen(A);
    Eigen::PartialPivLU<Eigen::MatrixXd> lu(eig_A);

    Eigen::MatrixXd eig_L = lu.matrixLU().triangularView<Eigen::UnitLower>();
    Eigen::MatrixXd eig_U = lu.matrixLU().triangularView<Eigen::Upper>();

    // Copy L
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            L(i, j) = eig_L(i, j);
        }
    }

    // Copy U
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i; j < n; ++j) {
            U(i, j) = eig_U(i, j);
        }
    }

    return result;
}

// Decompose complex matrix A into L and U such that A = L * U
// L is lower triangular matrix, U is upper triangular matrix
inline std::array<complex_matrix_owned, 2> lu(const complex_matrix_base &A) {
    if (A.rows() != A.cols()) {
        throw std::invalid_argument(
            "LU decomposition requires a square matrix");
    }

    size_t n = A.rows();
    std::array<complex_matrix_owned, 2> result;
    auto &L = result[0] = complex_matrix_owned(n, n); // Lower triangular matrix
    auto &U = result[1] = complex_matrix_owned(n, n); // Upper triangular matrix

    auto eig_A = eigen_interface::as_eigen(A);
    Eigen::PartialPivLU<
        Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>>
        lu(eig_A);

    const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>
        &matLU = lu.matrixLU();
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> eig_L =
        matLU.triangularView<Eigen::UnitLower>();
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> eig_U =
        matLU.triangularView<Eigen::Upper>();

    // Copy L
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            L(i, j) = eig_L(i, j);
        }
    }

    // Copy U
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i; j < n; ++j) {
            U(i, j) = eig_U(i, j);
        }
    }

    return result;
}

// ============================================================================
// 4. Matrix QR Decomposition Functions
// ============================================================================
// Decompose real matrix A into Q and R such that A = Q * R
// Q is orthogonal matrix, R is upper triangular matrix
inline std::array<real_matrix_owned, 2> qr(const real_matrix_base &A,
                                           bool full = false) {
    size_t m = A.rows();
    size_t n = A.cols();
    std::array<real_matrix_owned, 2> result;
    auto &Q = result[0];
    auto &R = result[1];

    auto eig_A = eigen_interface::as_eigen(A);
    Eigen::HouseholderQR<Eigen::MatrixXd> qr(eig_A);

    Eigen::MatrixXd eig_R = qr.matrixQR().triangularView<Eigen::Upper>();

    if (full) {
        // full QR: Q is m x m
        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(m, m);
        Q = eigen_interface::from_eigen(qr.householderQ() * I);
        R = eigen_interface::from_eigen(eig_R);
    } else {
        // thin QR: Q is m x n
        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(m, n);
        Q = eigen_interface::from_eigen(qr.householderQ() * I);
        R = eigen_interface::from_eigen(eig_R.topRows(n));
    }

    return result;
}

// Decompose complex matrix A into Q and R such that A = Q * R
// Q is orthogonal matrix, R is upper triangular matrix
inline std::array<complex_matrix_owned, 2> qr(const complex_matrix_base &A) {
    size_t m = A.rows();
    size_t n = A.cols();
    std::array<complex_matrix_owned, 2> result;
    auto &Q = result[0] = complex_matrix_owned(m, m); // Orthogonal matrix
    auto &R = result[1] = complex_matrix_owned(m, n); // Upper triangular matrix

    auto eig_A = eigen_interface::as_eigen(A);
    Eigen::HouseholderQR<
        Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>>
        qr(eig_A);

    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> eig_Q =
        qr.householderQ() * Eigen::MatrixXcd::Identity(m, m);
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> eig_R =
        qr.matrixQR().triangularView<Eigen::Upper>();

    // Copy Q
    // Note: eig_Q is complete m x m matrix
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < m; ++j) {
            Q(i, j) = eig_Q(i, j);
        }
    }

    // Copy R
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            R(i, j) = eig_R(i, j);
        }
    }

    return result;
}

} // namespace msl::matrix

#endif // MSL_MATRIX_DECOMPOSE_HPP
