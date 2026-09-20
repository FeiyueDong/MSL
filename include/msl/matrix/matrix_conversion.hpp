/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: matrix_conversion.hpp
** -----
** File Created: Wednesday, 15th October 2025 17:51:18
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Sunday, 14th December 2025 17:04:36
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_REAL_MATRIX_CONVERSION_HPP
#define MSL_REAL_MATRIX_CONVERSION_HPP

#include "complex_matrix_owned.hpp"
#include "complex_matrix_view.hpp"
#include "real_matrix_owned.hpp"
#include "real_matrix_view.hpp"

namespace msl::matrix {
// - Function of real_matrix_owned
// Construct from real_matrix_view (deep copy)
inline real_matrix_owned::real_matrix_owned(const real_matrix_view &view)
    : Base() {
    this->rows_ = view.rows();
    this->cols_ = view.cols();
    storage_.resize(view.size());
    update_span();
    // Manual copy required as view might not be contiguous column-major
    for (size_t j = 0; j < view.cols(); ++j) {
        for (size_t i = 0; i < view.rows(); ++i) {
            (*this)(i, j) = view(i, j);
        }
    }
}

// Copy assignment from real_matrix_view (deep copy)
inline real_matrix_owned &
real_matrix_owned::operator=(const real_matrix_view &view) {
    this->resize(view.rows(), view.cols());
    // Note: std::copy cannot be used directly if view is not contiguous
    // column-major. A manual loop is safer.
    for (size_t j = 0; j < view.cols(); ++j) {
        for (size_t i = 0; i < view.rows(); ++i) {
            (*this)(i, j) = view(i, j);
        }
    }
    update_span();
    return *this;
}

// - Function of real_matrix_view
// Construct from real_matrix_owned (most common use case)
inline real_matrix_view::real_matrix_view(real_matrix_owned &owner)
    : Base(owner.rows(), owner.cols(), owner.data_span()) {}

// Convert to owned matrix (deep copy)
inline real_matrix_owned real_matrix_view::to_owned() const {
    return real_matrix_owned(this->rows_, this->cols_, this->data_);
}

// - Function of const_real_matrix_view
// Construct from real_matrix_owned (most common use case)
inline const_real_matrix_view::const_real_matrix_view(
    const real_matrix_owned &owner)
    : data_(owner.data_span()), rows_(owner.rows()), cols_(owner.cols()) {}

inline real_matrix_owned const_real_matrix_view::to_owned() const {
    return real_matrix_owned(
        this->rows_, this->cols_, std::span<const double>(this->data_));
}

// - Function of complex_matrix_owned
// Construct from complex_matrix_view (deep copy)
inline complex_matrix_owned::complex_matrix_owned(
    const complex_matrix_view &view)
    : Base() {
    this->rows_ = view.rows();
    this->cols_ = view.cols();
    storage_.resize(view.size());
    // Manual copy required as view might not be contiguous column-major
    for (size_t j = 0; j < view.cols(); ++j) {
        for (size_t i = 0; i < view.rows(); ++i) {
            (*this)(i, j) = view(i, j);
        }
    }
    update_span();
}

// Copy assignment from complex_matrix_view (deep copy)
inline complex_matrix_owned &
complex_matrix_owned::operator=(const complex_matrix_view &view) {
    this->resize(view.rows(), view.cols());
    // Note: std::copy cannot be used directly if view is not contiguous
    // column-major. A manual loop is safer.
    for (size_t j = 0; j < view.cols(); ++j) {
        for (size_t i = 0; i < view.rows(); ++i) {
            (*this)(i, j) = view(i, j);
        }
    }
    update_span();
    return *this;
}

// - Function of complex_matrix_view
// Construct from complex_matrix_owned (most common use case)
inline complex_matrix_view::complex_matrix_view(complex_matrix_owned &owner)
    : Base(owner.rows(), owner.cols(), owner.data_span()) {}

// Convert to owned matrix (deep copy)
inline complex_matrix_owned complex_matrix_view::to_owned() const {
    return complex_matrix_owned(this->rows_, this->cols_, this->data_);
}

// - Function of const_complex_matrix_view
// Construct from complex_matrix_owned (most common use case)
inline const_complex_matrix_view::const_complex_matrix_view(
    const complex_matrix_owned &owner)
    : data_(owner.data_span()), rows_(owner.rows()), cols_(owner.cols()) {}

// Convert to owned matrix
inline complex_matrix_owned const_complex_matrix_view::to_owned() const {
    return complex_matrix_owned(
        this->rows_,
        this->cols_,
        std::span<const std::complex<double>>(this->data_));
}

} // namespace msl::matrix

#endif // MSL_REAL_MATRIX_CONVERSION_HPP