/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: real_matrix_view.hpp
** -----
** File Created: Saturday, 11th October 2025 21:20:53
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Sunday, 14th December 2025 17:04:58
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_REAL_MATRIX_VIEW_HPP
#define MSL_REAL_MATRIX_VIEW_HPP

#include <cassert>

#include "real_matrix_base.hpp"

namespace msl::matrix {
class const_real_matrix_view; // Forward declaration
class real_matrix_owned;      // Forward declaration

class real_matrix_view : public real_matrix_base {
private:
    using Base = real_matrix_base;

public:
    // --- Constructors ---
    // Default constructor (empty view)
    real_matrix_view() : Base() {}

    // Construct from raw pointer
    real_matrix_view(double *data, size_t rows, size_t cols)
        : Base(rows, cols, std::span<double>(data, rows * cols)) {}

    // Construct from real_matrix_owned (most common use case)
    real_matrix_view(real_matrix_owned &owner);

    // Construct from span
    real_matrix_view(std::span<double> data, size_t rows, size_t cols)
        : Base(rows, cols, data) {
        if (data.size() != rows * cols) {
            throw std::invalid_argument(
                "Matrix view: data span size must equal rows * cols");
        }
    }

    // --- Copy semantics: Shallow copy (copies the view) ---
    real_matrix_view(const real_matrix_view &other) = default;
    real_matrix_view &operator=(const real_matrix_view &other) = default;

    // --- Move semantics ---
    real_matrix_view(real_matrix_view &&other) noexcept = default;
    real_matrix_view &operator=(real_matrix_view &&other) noexcept = default;

    // --- Cannot construct view from const owned (use const_matrix_view) ---
    real_matrix_view(const real_matrix_owned &) = delete;

    // --- Destructor ---
    ~real_matrix_view() = default;

    // --- Subview creation ---
    [[nodiscard]] real_matrix_view subview(size_t row_start,
                                           size_t row_end,
                                           size_t col_start,
                                           size_t col_end) {
        if (row_start >= row_end || row_end > this->rows_) {
            throw std::out_of_range("Matrix view: row range is out of bounds");
        }
        if (col_start >= col_end || col_end > this->cols_) {
            throw std::out_of_range(
                "Matrix view: column range is out of bounds");
        }
        if (row_start != 0 || row_end != this->rows_) {
            throw std::invalid_argument(
                "Matrix view: only full-column subviews are supported in "
                "column-major layout");
        }

        size_t sub_cols = col_end - col_start;
        double *sub_data = this->data() + col_start * this->rows_;
        return real_matrix_view(sub_data, this->rows_, sub_cols);
    }

    // --- Column subview (efficient) ---
    [[nodiscard]] real_matrix_view column_view(size_t j) {
        auto col_span = this->column(j);
        return real_matrix_view(col_span.data(), this->rows_, 1);
    }
    [[nodiscard]] const_real_matrix_view column_view(size_t j) const;

    // --- Fill operations ---
    void fill(const double &value) {
        std::fill(this->begin(), this->end(), value);
    }

    void fill_zeros() { std::fill(this->begin(), this->end(), 0.0); }

    // --- In-place operations ---
    real_matrix_view &operator+=(const real_matrix_base &other) {
        if (this->rows_ != other.rows() || this->cols_ != other.cols()) {
            throw std::invalid_argument(
                "Matrix: operands must have the same shape");
        }
        for (size_t i = 0; i < this->size(); ++i) {
            this->data_[i] += other.data()[i];
        }
        return *this;
    }
    real_matrix_view &operator-=(const real_matrix_base &other) {
        if (this->rows_ != other.rows() || this->cols_ != other.cols()) {
            throw std::invalid_argument(
                "Matrix: operands must have the same shape");
        }
        for (size_t i = 0; i < this->size(); ++i) {
            this->data_[i] -= other.data()[i];
        }
        return *this;
    }
    real_matrix_view &operator*=(const double &scalar) {
        for (size_t i = 0; i < this->size(); ++i) {
            this->data_[i] *= scalar;
        }
        return *this;
    }
    real_matrix_view &operator/=(const double &scalar) {
        for (size_t i = 0; i < this->size(); ++i) {
            this->data_[i] /= scalar;
        }
        return *this;
    }

    // --- Copy data from another matrix ---
    void copy_from(const real_matrix_base &src) {
        if (this->rows_ != src.rows() || this->cols_ != src.cols()) {
            throw std::invalid_argument(
                "Matrix: source must have the same shape");
        }
        std::copy(src.begin(), src.end(), this->begin());
    }

    // --- Convert to owned matrix (deep copy) ---
    [[nodiscard]] real_matrix_owned to_owned() const;

    // --- Swap (swaps the views, not the data) ---
    void swap(real_matrix_view &other) noexcept {
        using std::swap;
        swap(this->data_, other.data_);
        swap(this->rows_, other.rows_);
        swap(this->cols_, other.cols_);
    }
};

/**
 * @brief Const matrix view
 *
 * Similar to matrix_view but provides read-only access
 */
class const_real_matrix_view {
private:
    std::span<const double> data_;
    size_t rows_{0};
    size_t cols_{0};

public:
    // --- Constructors ---
    const_real_matrix_view() = default;

    const_real_matrix_view(const double *data, size_t rows, size_t cols)
        : data_(data, rows * cols), rows_(rows), cols_(cols) {}

    const_real_matrix_view(const real_matrix_owned &owner);

    const_real_matrix_view(const real_matrix_view &view)
        : data_(view.data_span()), rows_(view.rows()), cols_(view.cols()) {}

    const_real_matrix_view(std::span<const double> data,
                           size_t rows,
                           size_t cols)
        : data_(data), rows_(rows), cols_(cols) {
        if (data.size() != rows * cols) {
            throw std::invalid_argument(
                "Matrix view: data span size must equal rows * cols");
        }
    }

    // --- Copy and move ---
    const_real_matrix_view(const const_real_matrix_view &) = default;
    const_real_matrix_view &operator=(const const_real_matrix_view &) = default;
    const_real_matrix_view(const_real_matrix_view &&) noexcept = default;
    const_real_matrix_view &
    operator=(const_real_matrix_view &&) noexcept = default;

    // --- Size queries ---
    [[nodiscard]] size_t rows() const noexcept { return rows_; }
    [[nodiscard]] size_t cols() const noexcept { return cols_; }
    [[nodiscard]] size_t size() const noexcept { return rows_ * cols_; }
    [[nodiscard]] std::pair<size_t, size_t> shape() const noexcept {
        return {rows_, cols_};
    }

    // --- Element access (read-only) ---
    [[nodiscard]] const double &operator()(size_t i, size_t j) const noexcept {
        return data_[j * rows_ + i];
    }

    [[nodiscard]] const double &operator[](size_t idx) const noexcept {
        return data_[idx];
    }

    // --- Column access ---
    [[nodiscard]] std::span<const double> column(size_t j) const noexcept {
        return {data_.data() + j * rows_, rows_};
    }

    // --- Data access ---
    [[nodiscard]] const double *data() const noexcept { return data_.data(); }
    [[nodiscard]] std::span<const double> data_span() const noexcept {
        return data_;
    }

    // --- Iterators ---
    [[nodiscard]] auto begin() const noexcept { return data_.begin(); }
    [[nodiscard]] auto end() const noexcept { return data_.end(); }
    [[nodiscard]] auto cbegin() const noexcept { return data_.begin(); }
    [[nodiscard]] auto cend() const noexcept { return data_.end(); }

    // --- Convert to owned matrix ---
    [[nodiscard]] real_matrix_owned to_owned() const;
};

// --- Column subview (efficient) ---
inline const_real_matrix_view real_matrix_view::column_view(size_t j) const {
    auto col_span = this->column(j);
    return const_real_matrix_view(col_span.data(), this->rows_, 1);
}

typedef real_matrix_view matrixd_view;
typedef const_real_matrix_view const_matrixd_view;

// --- Helper functions ---
inline void swap(real_matrix_view &a, real_matrix_view &b) noexcept {
    a.swap(b);
}

// Create view from owned matrix
[[nodiscard]] inline real_matrix_view make_view(real_matrix_owned &mat) {
    return real_matrix_view(mat);
}

[[nodiscard]] inline const_real_matrix_view
make_view(const real_matrix_owned &mat) {
    return const_real_matrix_view(mat);
}

} // namespace msl::matrix

#endif // MSL_MATRIX_VIEW_HPP