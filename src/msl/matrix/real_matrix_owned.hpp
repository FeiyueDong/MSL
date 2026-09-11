/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: real_matrix_owned.hpp
** -----
** File Created: Saturday, 11th October 2025 21:20:14
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Sunday, 14th December 2025 17:04:52
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_REAL_MATRIX_OWNED_HPP
#define MSL_REAL_MATRIX_OWNED_HPP

#include <algorithm>
#include <initializer_list>
#include <vector>

#include "real_matrix_base.hpp"

namespace msl::matrix {
class real_matrix_view; // Forward declaration

class real_matrix_owned : public real_matrix_base {
private:
    using Base = real_matrix_base;
    std::vector<double> storage_; // Actual data storage

    // Update base span when storage changes
    void update_span() { this->data_ = std::span<double>(storage_); }

public:
    // --- Constructors ---
    // Default constructor
    real_matrix_owned() : Base() {}

    // Size constructor (zero-initialized)
    real_matrix_owned(size_t rows, size_t cols)
        : Base(), storage_(rows * cols) {
        this->rows_ = rows;
        this->cols_ = cols;
        update_span();
    }

    // Size + value constructor
    real_matrix_owned(size_t rows, size_t cols, const double &value)
        : Base(), storage_(rows * cols, value) {
        this->rows_ = rows;
        this->cols_ = cols;
        update_span();
    }

    // Construct from span (deep copy)
    real_matrix_owned(size_t rows, size_t cols, std::span<const double> data)
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
    real_matrix_owned(size_t rows,
                      size_t cols,
                      std::initializer_list<double> values)
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
    real_matrix_owned(const real_matrix_owned &other)
        : Base(), storage_(other.storage_) {
        this->rows_ = other.rows_;
        this->cols_ = other.cols_;
        update_span();
    }
    explicit real_matrix_owned(const real_matrix_view &view);

    // --- Copy operations (deep copy) ---
    real_matrix_owned &operator=(const real_matrix_owned &other) {
        if (this != &other) {
            storage_ = other.storage_;
            this->rows_ = other.rows_;
            this->cols_ = other.cols_;
            update_span();
        }
        return *this;
    }
    real_matrix_owned &operator=(const real_matrix_view &view);

    // --- Move operations (no-throw guarantee) ---
    real_matrix_owned(real_matrix_owned &&other) noexcept
        : Base(), storage_(std::move(other.storage_)) {
        this->rows_ = other.rows_;
        this->cols_ = other.cols_;
        update_span();

        other.rows_ = 0;
        other.cols_ = 0;
        other.data_ = {};
    }
    real_matrix_owned &operator=(real_matrix_owned &&other) noexcept {
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
    ~real_matrix_owned() = default;

    // --- Size management ---
    void resize(size_t rows, size_t cols, const double &value = double()) {
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
    void fill(const double &value) {
        std::fill(storage_.begin(), storage_.end(), value);
    }

    void fill_zeros() { std::fill(storage_.begin(), storage_.end(), double{}); }

    // --- Row/Column extraction (returns owning copies) ---
    [[nodiscard]] real_matrix_owned row_copy(size_t i) const {
        if (i >= this->rows_) {
            throw std::out_of_range("Matrix: row index out of range");
        }
        real_matrix_owned row(1, this->cols_);
        for (size_t j = 0; j < this->cols_; ++j) {
            row(0, j) = (*this)(i, j);
        }
        return row;
    }
    [[nodiscard]] real_matrix_owned column_copy(size_t j) const {
        if (j >= this->cols_) {
            throw std::out_of_range("Matrix: column index out of range");
        }
        real_matrix_owned col(this->rows_, 1);
        const auto &col_span = this->column(j);
        std::copy(col_span.begin(), col_span.end(), col.storage_.begin());
        return col;
    }

    // --- Submatrix extraction (returns an owning copy) ---
    [[nodiscard]] real_matrix_owned submatrix_copy(size_t row_start,
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
        real_matrix_owned result(sub_rows, sub_cols);

        for (size_t j = 0; j < sub_cols; ++j) {
            for (size_t i = 0; i < sub_rows; ++i) {
                result(i, j) = (*this)(row_start + i, col_start + j);
            }
        }
        return result;
    }

    // --- In-place operations ---
    real_matrix_owned &operator+=(const real_matrix_base &other) {
        if (this->rows_ != other.rows() || this->cols_ != other.cols()) {
            throw std::invalid_argument(
                "Matrix: operands must have the same shape");
        }
        for (size_t i = 0; i < storage_.size(); ++i) {
            storage_[i] += other.data()[i];
        }
        return *this;
    }
    real_matrix_owned &operator-=(const real_matrix_base &other) {
        if (this->rows_ != other.rows() || this->cols_ != other.cols()) {
            throw std::invalid_argument(
                "Matrix: operands must have the same shape");
        }
        for (size_t i = 0; i < storage_.size(); ++i) {
            storage_[i] -= other.data()[i];
        }
        return *this;
    }
    real_matrix_owned &operator*=(const double &scalar) {
        for (auto &val : storage_) {
            val *= scalar;
        }
        return *this;
    }
    real_matrix_owned &operator/=(const double &scalar) {
        for (auto &val : storage_) {
            val /= scalar;
        }
        return *this;
    }

    // --- Swap ---
    void swap(real_matrix_owned &other) noexcept {
        using std::swap;
        swap(storage_, other.storage_);
        swap(this->rows_, other.rows_);
        swap(this->cols_, other.cols_);
        update_span();
        other.update_span();
    }

    // --- Factory methods ---
    [[nodiscard]] static real_matrix_owned zeros(size_t rows, size_t cols) {
        return real_matrix_owned(rows, cols, double{0});
    }
    [[nodiscard]] static real_matrix_owned ones(size_t rows, size_t cols) {
        return real_matrix_owned(rows, cols, double{1});
    }
    [[nodiscard]] static real_matrix_owned identity(size_t n) {
        real_matrix_owned result(n, n, double{});
        for (size_t i = 0; i < n; ++i) {
            result(i, i) = 1.0;
        }
        return result;
    }
    // Create diagonal matrix from vector
    [[nodiscard]] static real_matrix_owned
    diagonal(std::span<const double> diag) {
        size_t n = diag.size();
        real_matrix_owned result(n, n, 0.0);
        for (size_t i = 0; i < n; ++i) {
            result(i, i) = diag[i];
        }
        return result;
    }

    // --- Sum of elements ---
    [[nodiscard]] inline real_matrix_owned sum_columns() {
        real_matrix_owned result(1, cols_, 0.0);
        for (size_t j = 0; j < cols_; ++j) {
            for (size_t i = 0; i < rows_; ++i) {
                result(0, j) += (*this)(i, j);
            }
        }
        return result;
    }

    [[nodiscard]] inline real_matrix_owned sum_rows() {
        real_matrix_owned result(rows_, 1, 0.0);
        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                result(i, 0) += (*this)(i, j);
            }
        }
        return result;
    }
};

typedef real_matrix_owned matrixd;

// --- Free functions ---
inline void swap(real_matrix_owned &a, real_matrix_owned &b) noexcept {
    a.swap(b);
}

// Binary operators (return new matrix)
[[nodiscard]] inline real_matrix_owned operator+(const real_matrix_owned &a,
                                                 const real_matrix_owned &b) {
    if (a.rows() != b.rows() || a.cols() != b.cols()) {
        throw std::invalid_argument(
            "Matrix: operands must have the same shape");
    }
    real_matrix_owned result(a);
    result += b;
    return result;
}

[[nodiscard]] inline real_matrix_owned operator-(const real_matrix_owned &a,
                                                 const real_matrix_owned &b) {
    if (a.rows() != b.rows() || a.cols() != b.cols()) {
        throw std::invalid_argument(
            "Matrix: operands must have the same shape");
    }
    real_matrix_owned result(a);
    result -= b;
    return result;
}

[[nodiscard]] inline real_matrix_owned operator*(const real_matrix_owned &a,
                                                 const real_matrix_owned &b) {
    if (a.cols() != b.rows()) {
        throw std::invalid_argument(
            "Matrix: inner dimensions must agree for multiplication");
    }
    real_matrix_owned result(a.rows(), b.cols(), 0.0);
    for (size_t j = 0; j < b.cols(); ++j) {
        for (size_t i = 0; i < a.rows(); ++i) {
            double sum = 0.0;
            for (size_t k = 0; k < a.cols(); ++k) {
                sum += a(i, k) * b(k, j);
            }
            result(i, j) = sum;
        }
    }
    return result;
}

[[nodiscard]] inline real_matrix_owned operator-(const real_matrix_owned &mat) {
    real_matrix_owned result(mat);
    for (auto &val : result) {
        val = -val;
    }
    return result;
}

[[nodiscard]] inline real_matrix_owned operator*(const real_matrix_owned &mat,
                                                 const double &scalar) {
    real_matrix_owned result(mat);
    result *= scalar;
    return result;
}

[[nodiscard]] inline real_matrix_owned operator*(const double &scalar,
                                                 const real_matrix_owned &mat) {
    return mat * scalar;
}

[[nodiscard]] inline real_matrix_owned operator/(const real_matrix_owned &mat,
                                                 const double &scalar) {
    real_matrix_owned result(mat);
    result /= scalar;
    return result;
}
} // namespace msl::matrix

#endif // MSL_REAL_MATRIX_OWNED_HPP