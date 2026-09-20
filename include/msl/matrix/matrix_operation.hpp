/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: matrix_operation.hpp
** -----
** File Created: Sunday, 14th December 2025 16:38:53
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Sunday, 14th December 2025 17:04:40
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_MATRIX_OPERATION_HPP
#define MSL_MATRIX_OPERATION_HPP

#include "complex_matrix_base.hpp"
#include "complex_matrix_owned.hpp"
#include "eigen_interface.hpp"
#include "real_matrix_base.hpp"
#include "real_matrix_owned.hpp"

namespace msl::matrix {
// - Matrix transpose
inline real_matrix_owned transpose(const real_matrix_base &A) {
    size_t m = A.rows();
    size_t n = A.cols();
    real_matrix_owned At(n, m);
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            At(j, i) = A(i, j);
        }
    }
    return At;
}

// Matrix transpose
inline complex_matrix_owned transpose(const complex_matrix_base &A) {
    size_t m = A.rows();
    size_t n = A.cols();
    complex_matrix_owned At(n, m);
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            At(j, i) = A(i, j);
        }
    }
    return At;
}

// Matrix conjugate transpose
inline complex_matrix_owned conjugate_transpose(const complex_matrix_base &A) {
    matrixc result = transpose(A);
    for (auto &val : result) {
        val = std::conj(val);
    }

    return result;
}

// - Matrix trace
inline double trace(const real_matrix_base &A) {
    if (A.rows() != A.cols()) {
        throw std::invalid_argument("Trace requires a square matrix");
    }

    double tr = 0.0;
    size_t n = A.rows();
    for (size_t i = 0; i < n; ++i) {
        tr += A(i, i);
    }
    return tr;
}

inline std::complex<double> trace(const complex_matrix_base &A) {
    if (A.rows() != A.cols()) {
        throw std::invalid_argument("Trace requires a square matrix");
    }

    std::complex<double> tr = 0.0;
    size_t n = A.rows();
    for (size_t i = 0; i < n; ++i) {
        tr += A(i, i);
    }
    return tr;
}

// - Inverse of matrix
inline real_matrix_owned inverse(const real_matrix_base &A) {
    if (A.rows() != A.cols()) {
        throw std::invalid_argument("Inverse requires a square matrix");
    }

    const auto &eig_A = eigen_interface::as_eigen(A);
    Eigen::MatrixXd eig_A_inv = eig_A.inverse();

    return eigen_interface::from_eigen(eig_A_inv);
}

// Inverse of complex matrix
inline complex_matrix_owned inverse(const complex_matrix_base &A) {
    if (A.rows() != A.cols()) {
        throw std::invalid_argument("Inverse requires a square matrix");
    }

    size_t n = A.rows();
    complex_matrix_owned A_inv(n, n);

    auto eig_A = eigen_interface::as_eigen(A);
    Eigen::MatrixXcd eig_A_inv = eig_A.inverse();

    // Copy back
    std::copy(
        eig_A_inv.data(), eig_A_inv.data() + eig_A_inv.size(), A_inv.data());

    return A_inv;
}

// - determinant
inline double determinant(const real_matrix_base &A) {
    if (A.rows() != A.cols()) {
        throw std::invalid_argument("Determinant requires a square matrix");
    }

    auto eig_A = eigen_interface::as_eigen(A);
    return eig_A.determinant();
}

inline std::complex<double> determinant(const complex_matrix_base &A) {
    if (A.rows() != A.cols()) {
        throw std::invalid_argument("Determinant requires a square matrix");
    }

    auto eig_A = eigen_interface::as_eigen(A);
    return eig_A.determinant();
}

} // namespace msl::matrix

#endif // MSL_MATRIX_OPERATION_HPP