/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: complex_matrix_owned.hpp
** -----
** File Created: Saturday, 11th October 2025 21:20:57
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Sunday, 14th December 2025 17:04:16
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_COMPLEX_MATRIX_OWNED_HPP
#define MSL_COMPLEX_MATRIX_OWNED_HPP

#include <algorithm>
#include <initializer_list>
#include <vector>

#include "complex_matrix_base.hpp"
#include "real_matrix_owned.hpp"

namespace msl::matrix {
class complex_matrix_view; // Forward declaration

class complex_matrix_owned : public complex_matrix_base {
private:
    using Base = complex_matrix_base;
    std::vector<std::complex<double>> storage_; // Actual data storage

    // Update base span when storage changes
    void update_span() {
        this->data_ = std::span<std::complex<double>>(storage_);
    }

public:
    // --- Constructors ---
    // Default constructor
    complex_matrix_owned() : Base() {}

    // Size constructor (zero-initialized)
    complex_matrix_owned(size_t rows, size_t cols)
        : Base(), storage_(rows * cols) {
        this->rows_ = rows;
        this->cols_ = cols;
        update_span();
    }

    // Size + value constructor
    complex_matrix_owned(size_t rows,
                         size_t cols,
                         const std::complex<double> &value)
        : Base(), storage_(rows * cols, value) {
        this->rows_ = rows;
        this->cols_ = cols;
        update_span();
    }

    // Construct from span (deep copy)
    complex_matrix_owned(size_t rows,
                         size_t cols,
                         std::span<const std::complex<double>> data)
        : Base(), storage_(data.begin(), data.end()) {
        if (data.size() != rows * cols) {
            throw std::invalid_argument(
                "Matrix: data span size must equal rows * cols");
        }
        this->rows_ = rows;
        this->cols_ = cols;
        update_span();
    }

    // Construct from initializer list (row-major input)
    complex_matrix_owned(size_t rows,
                         size_t cols,
                         std::initializer_list<std::complex<double>> values)
        : Base(), storage_(rows * cols) {
        if (values.size() != rows * cols) {
            throw std::invalid_argument(
                "Matrix: initializer list size must equal rows * cols");
        }
        this->rows_ = rows;
        this->cols_ = cols;

        // Convert row-major to column-major
        auto it = values.begin();
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                storage_[j * rows + i] = *it++;
            }
        }
        update_span();
    }
    // --- Copy Constructors ---
    complex_matrix_owned(const complex_matrix_owned &other)
        : Base(), storage_(other.storage_) {
        this->rows_ = other.rows_;
        this->cols_ = other.cols_;
        update_span();
    }
    explicit complex_matrix_owned(const complex_matrix_view &view);

    // --- Copy operations (deep copy) ---
    complex_matrix_owned &operator=(const complex_matrix_owned &other) {
        if (this != &other) {
            storage_ = other.storage_;
            this->rows_ = other.rows_;
            this->cols_ = other.cols_;
            update_span();
        }
        return *this;
    }
    complex_matrix_owned &operator=(const complex_matrix_view &view);

    // --- Move operations (no-throw guarantee) ---
    complex_matrix_owned(complex_matrix_owned &&other) noexcept
        : Base(), storage_(std::move(other.storage_)) {
        this->rows_ = other.rows_;
        this->cols_ = other.cols_;
        update_span();

        other.rows_ = 0;
        other.cols_ = 0;
        other.data_ = {};
    }
    complex_matrix_owned &operator=(complex_matrix_owned &&other) noexcept {
        if (this != &other) {
            storage_ = std::move(other.storage_);
            this->rows_ = other.rows_;
            this->cols_ = other.cols_;
            update_span();

            other.rows_ = 0;
            other.cols_ = 0;
            other.data_ = {};
        }
        return *this;
    }

    // --- Destructor ---
    ~complex_matrix_owned() = default;

    // --- Size management ---
    void resize(size_t rows,
                size_t cols,
                const std::complex<double> &value = std::complex<double>()) {
        this->rows_ = rows;
        this->cols_ = cols;
        storage_.resize(rows * cols, value);
        update_span();
    }

    void clear() {
        this->rows_ = 0;
        this->cols_ = 0;
        storage_.clear();
        update_span();
    }

    // Reserve capacity (useful for avoiding reallocations)
    void reserve(size_t capacity) { storage_.reserve(capacity); }

    [[nodiscard]] size_t capacity() const noexcept {
        return storage_.capacity();
    }

    void shrink_to_fit() { storage_.shrink_to_fit(); }

    // --- Fill operations ---
    void fill(const std::complex<double> &value) {
        std::fill(storage_.begin(), storage_.end(), value);
    }

    void fill_zeros() {
        std::fill(storage_.begin(), storage_.end(), std::complex<double>{});
    }

    // --- Row/Column extraction (returns owning copies) ---
    [[nodiscard]] complex_matrix_owned row_copy(size_t i) const {
        if (i >= this->rows_) {
            throw std::out_of_range("Matrix: row index out of range");
        }
        complex_matrix_owned row(1, this->cols_);
        for (size_t j = 0; j < this->cols_; ++j) {
            row(0, j) = (*this)(i, j);
        }
        return row;
    }
    [[nodiscard]] complex_matrix_owned column_copy(size_t j) const {
        if (j >= this->cols_) {
            throw std::out_of_range("Matrix: column index out of range");
        }
        complex_matrix_owned col(this->rows_, 1);
        const auto &col_span = this->column(j);
        std::copy(col_span.begin(), col_span.end(), col.storage_.begin());
        return col;
    }

    // --- Submatrix extraction (returns an owning copy) ---
    [[nodiscard]] complex_matrix_owned submatrix_copy(size_t row_start,
                                                      size_t row_end,
                                                      size_t col_start,
                                                      size_t col_end) const {
        if (row_start >= row_end || row_end > this->rows_) {
            throw std::out_of_range("Matrix: row range is out of bounds");
        }
        if (col_start >= col_end || col_end > this->cols_) {
            throw std::out_of_range("Matrix: column range is out of bounds");
        }

        size_t sub_rows = row_end - row_start;
        size_t sub_cols = col_end - col_start;
        complex_matrix_owned result(sub_rows, sub_cols);

        for (size_t j = 0; j < sub_cols; ++j) {
            for (size_t i = 0; i < sub_rows; ++i) {
                result(i, j) = (*this)(row_start + i, col_start + j);
            }
        }
        return result;
    }

    // --- In-place operations ---
    complex_matrix_owned &operator+=(const complex_matrix_base &other) {
        if (this->rows_ != other.rows() || this->cols_ != other.cols()) {
            throw std::invalid_argument(
                "Matrix: operands must have the same shape");
        }
        for (size_t i = 0; i < storage_.size(); ++i) {
            storage_[i] += other.data()[i];
        }
        return *this;
    }
    complex_matrix_owned &operator-=(const complex_matrix_base &other) {
        if (this->rows_ != other.rows() || this->cols_ != other.cols()) {
            throw std::invalid_argument(
                "Matrix: operands must have the same shape");
        }
        for (size_t i = 0; i < storage_.size(); ++i) {
            storage_[i] -= other.data()[i];
        }
        return *this;
    }
    complex_matrix_owned &operator*=(const std::complex<double> &scalar) {
        for (auto &val : storage_) {
            val *= scalar;
        }
        return *this;
    }
    complex_matrix_owned &operator/=(const std::complex<double> &scalar) {
        for (auto &val : storage_) {
            val /= scalar;
        }
        return *this;
    }

    // --- Swap ---
    void swap(complex_matrix_owned &other) noexcept {
        using std::swap;
        swap(storage_, other.storage_);
        swap(this->rows_, other.rows_);
        swap(this->cols_, other.cols_);
        update_span();
        other.update_span();
    }

    // --- Factory methods ---
    // Create diagonal matrix from vector
    [[nodiscard]] static complex_matrix_owned
    diagonal(std::span<const std::complex<double>> diag) {
        size_t n = diag.size();
        complex_matrix_owned result(n, n, 0.0);
        for (size_t i = 0; i < n; ++i) {
            result(i, i) = diag[i];
        }
        return result;
    }

    // --- complex to real
    // Extract real part as real_matrix_owned
    [[nodiscard]] real_matrix_owned real() const {
        real_matrix_owned real_mat(this->rows_, this->cols_);
        for (size_t i = 0; i < this->rows_; ++i) {
            for (size_t j = 0; j < this->cols_; ++j) {
                real_mat(i, j) = std::real((*this)(i, j));
            }
        }
        return real_mat;
    }
    // Extract imaginary part as real_matrix_owned
    [[nodiscard]] real_matrix_owned imag() const {
        real_matrix_owned imag_mat(this->rows_, this->cols_);
        for (size_t i = 0; i < this->rows_; ++i) {
            for (size_t j = 0; j < this->cols_; ++j) {
                imag_mat(i, j) = std::imag((*this)(i, j));
            }
        }
        return imag_mat;
    }
    // Convert to real_matrix_owned by taking magnitude
    [[nodiscard]] real_matrix_owned magnitude() const {
        real_matrix_owned mag_mat(this->rows_, this->cols_);
        for (size_t i = 0; i < this->rows_; ++i) {
            for (size_t j = 0; j < this->cols_; ++j) {
                mag_mat(i, j) = std::abs((*this)(i, j));
            }
        }
        return mag_mat;
    }
    // Convert to real_matrix_owned by taking phase (angle)
    [[nodiscard]] real_matrix_owned phase() const {
        real_matrix_owned phase_mat(this->rows_, this->cols_);
        for (size_t i = 0; i < this->rows_; ++i) {
            for (size_t j = 0; j < this->cols_; ++j) {
                phase_mat(i, j) = std::arg((*this)(i, j));
            }
        }
        return phase_mat;
    }
    // conjugate transpose
    [[nodiscard]] complex_matrix_owned conjugate_transpose() const {
        complex_matrix_owned adj_mat(this->cols_, this->rows_);
        for (size_t i = 0; i < this->rows_; ++i) {
            for (size_t j = 0; j < this->cols_; ++j) {
                adj_mat(j, i) = std::conj((*this)(i, j));
            }
        }
        return adj_mat;
    }

    // --- Sum of elements ---
    [[nodiscard]] inline complex_matrix_owned sum_columns() {
        complex_matrix_owned result(1, cols_, std::complex<double>{0.0, 0.0});
        for (size_t j = 0; j < cols_; ++j) {
            for (size_t i = 0; i < rows_; ++i) {
                result(0, j) += (*this)(i, j);
            }
        }
        return result;
    }

    [[nodiscard]] inline complex_matrix_owned sum_rows() {
        complex_matrix_owned result(rows_, 1, std::complex<double>{0.0, 0.0});
        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                result(i, 0) += (*this)(i, j);
            }
        }
        return result;
    }
};

typedef complex_matrix_owned matrixc;

// --- Free functions ---
inline void swap(complex_matrix_owned &a, complex_matrix_owned &b) noexcept {
    a.swap(b);
}

// Binary operators (return new matrix)
[[nodiscard]] inline complex_matrix_owned
operator+(const complex_matrix_owned &a, const complex_matrix_owned &b) {
    if (a.rows() != b.rows() || a.cols() != b.cols()) {
        throw std::invalid_argument(
            "Matrix: operands must have the same shape");
    }
    complex_matrix_owned result(a);
    result += b;
    return result;
}

[[nodiscard]] inline complex_matrix_owned
operator-(const complex_matrix_owned &a, const complex_matrix_owned &b) {
    if (a.rows() != b.rows() || a.cols() != b.cols()) {
        throw std::invalid_argument(
            "Matrix: operands must have the same shape");
    }
    complex_matrix_owned result(a);
    result -= b;
    return result;
}

[[nodiscard]] inline complex_matrix_owned
operator*(const complex_matrix_owned &a, const complex_matrix_owned &b) {
    if (a.cols() != b.rows()) {
        throw std::invalid_argument(
            "Matrix: inner dimensions must agree for multiplication");
    }
    complex_matrix_owned result(a.rows(), b.cols(), std::complex<double>{});
    for (size_t j = 0; j < b.cols(); ++j) {
        for (size_t i = 0; i < a.rows(); ++i) {
            std::complex<double> sum = std::complex<double>{};
            for (size_t k = 0; k < a.cols(); ++k) {
                sum += a(i, k) * b(k, j);
            }
            result(i, j) = sum;
        }
    }
    return result;
}

[[nodiscard]] inline complex_matrix_owned
operator-(const complex_matrix_owned &mat) {
    complex_matrix_owned result(mat);
    for (auto &val : result) {
        val = -val;
    }
    return result;
}

[[nodiscard]] inline complex_matrix_owned
operator*(const complex_matrix_owned &mat, const double &scalar) {
    complex_matrix_owned result(mat);
    result *= scalar;
    return result;
}

[[nodiscard]] inline complex_matrix_owned
operator*(const double &scalar, const complex_matrix_owned &mat) {
    return mat * scalar;
}

[[nodiscard]] inline complex_matrix_owned
operator/(const complex_matrix_owned &mat, const double &scalar) {
    complex_matrix_owned result(mat);
    result /= scalar;
    return result;
}

} // namespace msl::matrix

#endif // MSL_complex_MATRIX_OWNED_HPP