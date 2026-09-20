/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: real_matrix_base.hpp
** -----
** File Created: Saturday, 11th October 2025 21:19:59
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Sunday, 14th December 2025 17:04:45
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_REAL_MATRIX_BASE_HPP
#define MSL_REAL_MATRIX_BASE_HPP

#include <cassert>
#include <span>
#include <stdexcept>
#include <utility>

namespace msl::matrix {
class real_matrix_base {
protected:
    std::span<double> data_; // Unified data access through span
    size_t rows_{0};
    size_t cols_{0};

public:
    real_matrix_base() = default;

    real_matrix_base(size_t rows, size_t cols, std::span<double> data)
        : data_(data), rows_(rows), cols_(cols) {
        if (data_.size() != rows * cols) {
            throw std::invalid_argument(
                "Matrix: data span size must equal rows * cols");
        }
    }

    real_matrix_base(const real_matrix_base &) = default;
    real_matrix_base(real_matrix_base &&) noexcept = default;
    real_matrix_base &operator=(const real_matrix_base &) = default;
    real_matrix_base &operator=(real_matrix_base &&) noexcept = default;

    // --- Size queries ---
    [[nodiscard]] size_t rows() const noexcept { return rows_; }
    [[nodiscard]] size_t cols() const noexcept { return cols_; }
    [[nodiscard]] size_t size() const noexcept { return rows_ * cols_; }
    [[nodiscard]] std::pair<size_t, size_t> shape() const noexcept {
        return {rows_, cols_};
    }
    [[nodiscard]] bool empty() const noexcept { return size() == 0; }

    // --- Element access (CRTP: no virtual function overhead) ---
    [[nodiscard]] double &operator()(size_t i, size_t j) noexcept {
        return data_[j * rows_ + i]; // Column-major order
    }

    [[nodiscard]] const double &operator()(size_t i, size_t j) const noexcept {
        return data_[j * rows_ + i];
    }

    // Linear indexing (useful for iteration)
    [[nodiscard]] double &operator[](size_t idx) noexcept { return data_[idx]; }

    [[nodiscard]] const double &operator[](size_t idx) const noexcept {
        return data_[idx];
    }

    // --- Column access (efficient in column-major layout) ---
    [[nodiscard]] std::span<double> column(size_t j) {
        if (j >= cols_) {
            throw std::out_of_range("Matrix: column index out of range");
        }
        return {data_.data() + j * rows_, rows_};
    }

    [[nodiscard]] std::span<const double> column(size_t j) const {
        if (j >= cols_) {
            throw std::out_of_range("Matrix: column index out of range");
        }
        return {data_.data() + j * rows_, rows_};
    }

    // --- Data access ---
    [[nodiscard]] double *data() noexcept { return data_.data(); }
    [[nodiscard]] const double *data() const noexcept { return data_.data(); }
    [[nodiscard]] std::span<double> data_span() noexcept { return data_; }
    [[nodiscard]] std::span<const double> data_span() const noexcept {
        return data_;
    }

    // --- Iterators ---
    [[nodiscard]] auto begin() noexcept { return data_.begin(); }
    [[nodiscard]] auto end() noexcept { return data_.end(); }
    [[nodiscard]] auto begin() const noexcept { return data_.begin(); }
    [[nodiscard]] auto end() const noexcept { return data_.end(); }
    //[[nodiscard]] auto cbegin() const noexcept { return data_.cbegin(); }
    //[[nodiscard]] auto cend() const noexcept { return data_.cend(); }

    // --- Apply function to each element ---
    template <typename Func>
    void apply(Func &&func) {
        for (auto &val : data_) {
            val = func(val);
        }
    }

    // --- Common operations ---
    [[nodiscard]] inline double sum() {
        double total = 0.0;
        for (const auto &val : data_) {
            total += val;
        }
        return total;
    }

    [[nodiscard]] inline double trace() {
        if (rows_ != cols_) {
            throw std::invalid_argument("Trace requires a square matrix");
        }
        double tr = 0.0;
        for (size_t i = 0; i < rows_; ++i) {
            tr += (*this)(i, i);
        }
        return tr;
    }

protected:
    // Helper for checking consistency
    [[nodiscard]] bool is_consistent() const noexcept {
        return data_.size() == rows_ * cols_;
    }
};

} // namespace msl::matrix

#endif // MSL_REAL_MATRIX_BASE_HPP