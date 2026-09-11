/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: polynominal.hpp
** -----
** File Created: Saturday, 13th December 2025 22:56:46
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Sunday, 14th December 2025 17:05:05
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_POLYNOMIAL_HPP
#define MSL_POLYNOMIAL_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <span>
#include <stdexcept>
#include <vector>

#include "matrix/martrix_decompose.hpp"
#include "matrix/matrix_operation.hpp"
#include "matrix/real_matrix_owned.hpp"

namespace msl::polynomial {
/**
 * @brief Polynomial model and least-squares fitting utilities.
 *
 * Coefficients are stored in ascending power order:
 * coeffs_[0] + coeffs_[1] * x + coeffs_[2] * x^2 + ...
 */
class Polynomial {
private:
    std::vector<double> coeffs_{};

public:
    Polynomial() = default;
    /**
     * @brief Construct polynomial directly from coefficients.
     *
     * @param coeffs Coefficients in ascending power order
     */
    Polynomial(const std::span<const double> &coeffs) {
        set_coefficients(coeffs);
    }

    /**
     * @brief Create polynomial by fitting to y values using implicit
     * x = [0, 1, 2, ...].
     *
     * @param y Sample values
     * @param n Polynomial degree
     * @return Fitted Polynomial
     */
    [[nodiscard]] static Polynomial from_fit(std::span<const double> y,
                                             std::size_t n) {
        std::vector<double> x(y.size());
        std::iota(x.begin(), x.end(), 0.0);
        Polynomial poly;
        poly.fit(x, y, n);
        return poly;
    }

    /**
     * @brief Create polynomial by fitting to data points.
     *
     * @param x Independent variable samples
     * @param y Dependent variable samples
     * @param n Polynomial degree
     * @return Fitted Polynomial
     */
    [[nodiscard]] static Polynomial from_fit(std::span<const double> x,
                                             std::span<const double> y,
                                             std::size_t n = 0) {
        Polynomial poly;
        poly.fit(x, y, n);
        return poly;
    }

    /**
     * @brief Evaluate polynomial at a single point.
     */
    double operator()(double x) const { return evaluate(x); }

    /**
     * @brief Evaluate polynomial at multiple points (zero-copy output).
     *
     * @param x_values Query points
     * @param result Output buffer (must have same size as @p x_values)
     */
    void operator()(std::span<const double> x_values,
                    std::span<double> result) const {
        if (x_values.size() != result.size()) {
            throw std::invalid_argument(
                "Polynomial: x_values and result spans must have same size");
        }

        for (std::size_t i = 0; i < x_values.size(); ++i) {
            result[i] = evaluate(x_values[i]);
        }
    }

    /**
     * @brief Evaluate polynomial at multiple points.
     *
     * @param x_values Query points
     * @return Evaluation results
     */
    std::vector<double> operator()(std::span<const double> x_values) const {
        std::vector<double> results(x_values.size());
        operator()(x_values, results);
        return results;
    }

    /**
     * @brief Replace polynomial coefficients.
     *
     * @param coeffs Coefficients in ascending power order
     */
    void set_coefficients(std::span<const double> coeffs) {
        if (coeffs.empty()) {
            throw std::invalid_argument(
                "Polynomial: coefficients cannot be empty");
        }

        coeffs_.assign(coeffs.begin(), coeffs.end());
    }

    /**
     * @brief Fit polynomial coefficients from data points.
     *
     * @param x Independent variable samples
     * @param y Dependent variable samples
     * @param n Polynomial degree
     */
    void fit(std::span<const double> x,
             std::span<const double> y,
             std::size_t n = 0) {
        calc_coefficients(x, y, n);
    }

    /**
     * @brief Get polynomial coefficients.
     */
    [[nodiscard]] std::vector<double> coefficients() const { return coeffs_; }

    /**
     * @brief Get polynomial degree.
     */
    [[nodiscard]] std::size_t degree() const {
        if (coeffs_.empty()) {
            throw std::runtime_error("Polynomial: coefficients are not set");
        }
        return coeffs_.size() - 1;
    }

    /**
     * @brief Evaluate first derivative at a single point.
     *
     * @param x Query point
     * @return Derivative value
     */
    [[nodiscard]] double derivative(double x) const {
        if (coeffs_.empty()) {
            throw std::runtime_error("Polynomial: coefficients are not set");
        }

        if (coeffs_.size() == 1) {
            return 0.0;
        }

        double result = 0.0;
        for (std::size_t i = coeffs_.size() - 1; i > 0; --i) {
            result = result * x + static_cast<double>(i) * coeffs_[i];
        }
        return result;
    }

private:
    /**
     * @brief Evaluate polynomial using Horner's method.
     */
    double evaluate(double x) const {
        if (coeffs_.empty()) {
            throw std::runtime_error("Polynomial: coefficients are not set");
        }

        double result = 0.0;
        for (std::size_t i = coeffs_.size(); i-- > 0;) {
            result = result * x + coeffs_[i];
        }

        return result;
    }

    /**
     * @brief Compute coefficients by least-squares fitting.
     *
     * Builds the Vandermonde system and solves it with an economy QR
     * factorization followed by triangular back-substitution.
     *
     * @throws std::invalid_argument if `x` and `y` differ in size or fewer
     * than `n + 1` samples are provided
     * @throws std::runtime_error if the system is rank deficient (for example
     * duplicate sample locations)
     */
    void calc_coefficients(std::span<const double> x,
                           std::span<const double> y,
                           std::size_t n = 0) {
        if (x.size() != y.size()) {
            throw std::invalid_argument(
                "Polynomial: x and y must have same size");
        }

        if (x.size() < n + 1) {
            throw std::invalid_argument(
                "Polynomial: not enough points to fit the polynomial");
        }

        if (n == 0) {
            coeffs_.resize(1);
            coeffs_[0] = std::accumulate(y.begin(), y.end(), 0.0)
                         / static_cast<double>(y.size());
            return;
        }

        // Least-squares fit via economy QR: solve R * c = Q^T * y by
        // back-substitution (no explicit inverse).
        matrix::matrixd A(x.size(), n + 1);
        for (std::size_t i = 0; i < x.size(); ++i) {
            double x_pow = 1.0;
            for (std::size_t j = 0; j <= n; ++j) {
                A(i, j) = x_pow;
                x_pow *= x[i];
            }
        }

        const auto qr_result = matrix::qr(A); // Q: m x k, R: k x k
        const auto &Q = qr_result[0];
        const auto &R = qr_result[1];
        const std::size_t k = n + 1;

        const auto y_column = matrix::matrixd(y.size(), 1, y);
        const auto Qt_y = matrix::transpose(Q) * y_column;

        // Rank check on the diagonal of R (column order is preserved)
        double max_pivot = 0.0;
        for (std::size_t i = 0; i < k; ++i) {
            max_pivot = std::max(max_pivot, std::abs(R(i, i)));
        }
        const double tolerance = std::numeric_limits<double>::epsilon()
                                 * static_cast<double>(std::max(x.size(), k))
                                 * max_pivot;

        std::vector<double> coefficients(k, 0.0);
        for (std::size_t i = k; i-- > 0;) {
            if (std::abs(R(i, i)) <= tolerance) {
                throw std::runtime_error(
                    "Polynomial: rank-deficient fit (duplicate or dependent "
                    "sample locations)");
            }
            double sum = Qt_y(i, 0);
            for (std::size_t j = i + 1; j < k; ++j) {
                sum -= R(i, j) * coefficients[j];
            }
            coefficients[i] = sum / R(i, i);
        }

        coeffs_ = std::move(coefficients);
    }
};

// --- Free functions ---

/**
 * @brief Polynomial fitting.
 *
 * @param x Independent variable samples
 * @param y Dependent variable samples
 * @param n Polynomial degree
 * @return Polynomial coefficients
 */
inline std::vector<double> polyfit(std::span<const double> x,
                                   std::span<const double> y,
                                   std::size_t n = 0) {
    return Polynomial::from_fit(x, y, n).coefficients();
}

/**
 * @brief Polynomial fitting (zero-copy output).
 *
 * @param x Independent variable samples
 * @param y Dependent variable samples
 * @param result Output buffer for coefficients (must have size n+1)
 * @param n Polynomial degree
 */
inline void polyfit(std::span<const double> x,
                    std::span<const double> y,
                    std::span<double> result,
                    std::size_t n = 0) {
    if (result.size() != n + 1) {
        throw std::invalid_argument(
            "polyfit: result span size must equal n + 1");
    }

    auto coeffs = Polynomial::from_fit(x, y, n).coefficients();
    std::copy(coeffs.begin(), coeffs.end(), result.begin());
}

/**
 * @brief Evaluate polynomial at one point.
 */
inline double polyval(const std::span<const double> &coeffs, double x) {
    return Polynomial(coeffs)(x);
}

/**
 * @brief Evaluate polynomial at multiple points.
 */
inline std::vector<double> polyval(const std::span<const double> &coeffs,
                                   const std::span<const double> &x_values) {
    return Polynomial(coeffs)(x_values);
}

/**
 * @brief Evaluate polynomial at multiple points (zero-copy output).
 *
 * @param coeffs Polynomial coefficients
 * @param x_values Query points
 * @param result Output buffer (must have same size as @p x_values)
 */
inline void polyval(std::span<const double> coeffs,
                    std::span<const double> x_values,
                    std::span<double> result) {
    Polynomial poly(coeffs);
    poly(x_values, result);
}

}; // namespace msl::polynomial

#endif // MSL_POLYNOMIAL_HPP
