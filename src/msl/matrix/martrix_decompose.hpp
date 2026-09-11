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

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/SVD>
#include <limits>
#include <numeric>
#include <random>
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

struct truncated_svd_options {
    size_t oversampling = 10;
    size_t power_iterations = 2;
    std::uint64_t seed = 0x5eedULL;
    bool compute_right_vectors = true;
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

namespace detail {
inline Eigen::MatrixXd thin_q(const Eigen::MatrixXd &matrix) {
    Eigen::HouseholderQR<Eigen::MatrixXd> qr(matrix);
    return qr.householderQ()
           * Eigen::MatrixXd::Identity(matrix.rows(), matrix.cols());
}

inline std::vector<Eigen::Index>
descending_singular_value_order(const Eigen::VectorXd &singular_values) {
    std::vector<Eigen::Index> order(
        static_cast<size_t>(singular_values.size()));
    std::iota(order.begin(), order.end(), Eigen::Index{0});
    std::stable_sort(order.begin(),
                     order.end(),
                     [&singular_values](Eigen::Index lhs, Eigen::Index rhs) {
                         return singular_values(lhs) > singular_values(rhs);
                     });
    return order;
}

inline void validate_truncated_svd_input(const real_matrix_base &A,
                                         size_t rank) {
    const size_t max_rank = std::min(A.rows(), A.cols());
    if (rank == 0 || rank > max_rank) {
        throw std::invalid_argument(
            "Truncated SVD: rank must be in [1, min(rows, cols)]");
    }
    if (!std::all_of(A.data(), A.data() + A.size(), [](double value) {
            return std::isfinite(value);
        })) {
        throw std::invalid_argument(
            "Truncated SVD: input must contain only finite values");
    }
}

inline truncated_svd_result
exact_truncated_svd(const Eigen::Map<const Eigen::MatrixXd> &A,
                    size_t rank,
                    bool compute_right_vectors) {
    unsigned int decomposition_options = Eigen::ComputeThinU;
    if (compute_right_vectors) {
        decomposition_options |= Eigen::ComputeThinV;
    }
    Eigen::JacobiSVD<Eigen::MatrixXd> decomposition(A, decomposition_options);
    const auto order =
        descending_singular_value_order(decomposition.singularValues());

    const Eigen::Index result_rank = static_cast<Eigen::Index>(rank);
    Eigen::MatrixXd selected_U(A.rows(), result_rank);
    Eigen::MatrixXd selected_V;
    if (compute_right_vectors) {
        selected_V.resize(A.cols(), result_rank);
    }

    truncated_svd_result result;
    result.singular_values.resize(rank);
    for (size_t i = 0; i < rank; ++i) {
        const auto source = order[i];
        selected_U.col(static_cast<Eigen::Index>(i)) =
            decomposition.matrixU().col(source);
        result.singular_values[i] = decomposition.singularValues()(source);
        if (compute_right_vectors) {
            selected_V.col(static_cast<Eigen::Index>(i)) =
                decomposition.matrixV().col(source);
        }
    }
    result.U = eigen_interface::from_eigen(selected_U);
    if (compute_right_vectors) {
        result.V = eigen_interface::from_eigen(selected_V);
    }
    return result;
}
} // namespace detail

/**
 * @brief Compute the leading singular triplets using deterministic randomized
 * SVD.
 *
 * The returned matrices satisfy `A ~= U * diag(singular_values) * V^T`.
 * `V` is empty when `options.compute_right_vectors` is false.
 */
inline truncated_svd_result
truncated_svd(const real_matrix_base &A,
              size_t rank,
              const truncated_svd_options &options) {
    detail::validate_truncated_svd_input(A, rank);
    const size_t max_rank = std::min(A.rows(), A.cols());
    const size_t sample_size = options.oversampling >= max_rank - rank
                                   ? max_rank
                                   : rank + options.oversampling;
    const auto eig_A = eigen_interface::as_eigen(A);

    if (sample_size == max_rank) {
        return detail::exact_truncated_svd(
            eig_A, rank, options.compute_right_vectors);
    }

    const Eigen::Index samples = static_cast<Eigen::Index>(sample_size);
    Eigen::MatrixXd omega(A.cols(), samples);
    std::mt19937_64 generator(options.seed);
    std::normal_distribution<double> distribution(0.0, 1.0);
    for (Eigen::Index j = 0; j < omega.cols(); ++j) {
        for (Eigen::Index i = 0; i < omega.rows(); ++i) {
            omega(i, j) = distribution(generator);
        }
    }

    Eigen::MatrixXd Q = detail::thin_q(eig_A * omega);
    for (size_t iteration = 0; iteration < options.power_iterations;
         ++iteration) {
        Eigen::MatrixXd Z = detail::thin_q(eig_A.transpose() * Q);
        Q = detail::thin_q(eig_A * Z);
    }

    Eigen::MatrixXd projected = Q.transpose() * eig_A;
    unsigned int decomposition_options = Eigen::ComputeThinU;
    if (options.compute_right_vectors) {
        decomposition_options |= Eigen::ComputeThinV;
    }
    Eigen::JacobiSVD<Eigen::MatrixXd> decomposition(projected,
                                                    decomposition_options);
    const auto order =
        detail::descending_singular_value_order(decomposition.singularValues());

    const Eigen::Index result_rank = static_cast<Eigen::Index>(rank);
    Eigen::MatrixXd selected_Ub(samples, result_rank);
    Eigen::MatrixXd selected_V;
    if (options.compute_right_vectors) {
        selected_V.resize(A.cols(), result_rank);
    }

    truncated_svd_result result;
    result.singular_values.resize(rank);
    for (size_t i = 0; i < rank; ++i) {
        const auto source = order[i];
        selected_Ub.col(static_cast<Eigen::Index>(i)) =
            decomposition.matrixU().col(source);
        result.singular_values[i] = decomposition.singularValues()(source);
        if (options.compute_right_vectors) {
            selected_V.col(static_cast<Eigen::Index>(i)) =
                decomposition.matrixV().col(source);
        }
    }
    result.U = eigen_interface::from_eigen(Q * selected_Ub);
    if (options.compute_right_vectors) {
        result.V = eigen_interface::from_eigen(selected_V);
    }
    return result;
}

/**
 * @brief Compute the leading singular triplets with default randomized SVD
 * options.
 */
inline truncated_svd_result truncated_svd(const real_matrix_base &A,
                                          size_t rank,
                                          bool compute_right_vectors = true) {
    truncated_svd_options options;
    options.compute_right_vectors = compute_right_vectors;
    return truncated_svd(A, rank, options);
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
/**
 * @brief Result of a real LU decomposition with partial pivoting.
 *
 * The decomposition satisfies `P * A = L * U`, where `P` is the row
 * permutation defined by `permutation`: row `i` of `P * A` is row
 * `permutation[i]` of `A`.
 */
struct lu_result {
    std::vector<size_t> permutation; // Row permutation of P * A = L * U
    real_matrix_owned L;             // Unit lower triangular matrix
    real_matrix_owned U;             // Upper triangular matrix
};

/**
 * @brief Result of a complex LU decomposition with partial pivoting.
 *
 * The decomposition satisfies `P * A = L * U`, where `P` is the row
 * permutation defined by `permutation`: row `i` of `P * A` is row
 * `permutation[i]` of `A`.
 */
struct complex_lu_result {
    std::vector<size_t> permutation; // Row permutation of P * A = L * U
    complex_matrix_owned L;          // Unit lower triangular matrix
    complex_matrix_owned U;          // Upper triangular matrix
};

// Decompose real matrix A with partial pivoting such that P * A = L * U
// L is unit lower triangular, U is upper triangular, P is a row permutation
inline lu_result lu(const real_matrix_base &A) {
    if (A.rows() != A.cols()) {
        throw std::invalid_argument(
            "LU decomposition requires a square matrix");
    }

    size_t n = A.rows();
    lu_result result;
    result.L = real_matrix_owned(n, n);
    result.U = real_matrix_owned(n, n);

    auto eig_A = eigen_interface::as_eigen(A);
    Eigen::PartialPivLU<Eigen::MatrixXd> lu(eig_A);

    Eigen::MatrixXd eig_L = lu.matrixLU().triangularView<Eigen::UnitLower>();
    Eigen::MatrixXd eig_U = lu.matrixLU().triangularView<Eigen::Upper>();

    // Copy L
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            result.L(i, j) = eig_L(i, j);
        }
    }

    // Copy U
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i; j < n; ++j) {
            result.U(i, j) = eig_U(i, j);
        }
    }

    // Row permutation: (P * A)(i, j) = A(permutation[i], j)
    result.permutation.resize(n);
    const auto &indices = lu.permutationP().indices();
    for (size_t i = 0; i < n; ++i) {
        const auto permuted_row =
            static_cast<size_t>(indices[static_cast<Eigen::Index>(i)]);
        result.permutation[permuted_row] = i;
    }

    return result;
}

// Decompose complex matrix A with partial pivoting such that P * A = L * U
// L is unit lower triangular, U is upper triangular, P is a row permutation
inline complex_lu_result lu(const complex_matrix_base &A) {
    if (A.rows() != A.cols()) {
        throw std::invalid_argument(
            "LU decomposition requires a square matrix");
    }

    size_t n = A.rows();
    complex_lu_result result;
    result.L = complex_matrix_owned(n, n);
    result.U = complex_matrix_owned(n, n);

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
            result.L(i, j) = eig_L(i, j);
        }
    }

    // Copy U
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i; j < n; ++j) {
            result.U(i, j) = eig_U(i, j);
        }
    }

    // Row permutation: (P * A)(i, j) = A(permutation[i], j)
    result.permutation.resize(n);
    const auto &indices = lu.permutationP().indices();
    for (size_t i = 0; i < n; ++i) {
        const auto permuted_row =
            static_cast<size_t>(indices[static_cast<Eigen::Index>(i)]);
        result.permutation[permuted_row] = i;
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
