/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2026, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: secant.hpp
** -----
** Author: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_SECANT_HPP
#define MSL_SECANT_HPP

#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include "root_options.hpp"
#include "root_result.hpp"

namespace msl::equation {

/**
 * @brief Find a root using the secant method.
 *
 * @param f Function whose root is sought
 * @param x0 First initial value
 * @param x1 Second initial value
 * @param tol Absolute tolerance for residual or step size
 * @param max_iter Maximum number of iterations
 * @return Root-finding result
 */
template <typename Func>
    requires std::is_invocable_r_v<double, Func, double>
inline root_result secant(Func &&f,
                          double x0,
                          double x1,
                          double tol = 1e-12,
                          size_t max_iter = 100) {
    if (x0 == x1) {
        throw std::invalid_argument("Secant: initial values must be distinct");
    }
    if (tol <= 0.0) {
        throw std::invalid_argument("Secant: tolerance must be positive");
    }
    if (max_iter == 0) {
        throw std::invalid_argument("Secant: max_iter must be positive");
    }

    auto func = std::forward<Func>(f);
    double f0 = std::invoke(func, x0);
    double f1 = std::invoke(func, x1);
    if (!std::isfinite(f0) || !std::isfinite(f1)) {
        return {x1, f1, 0, root_status::non_finite_value};
    }
    if (std::abs(f0) <= tol) {
        return {x0, f0, 0, root_status::converged};
    }
    if (std::abs(f1) <= tol) {
        return {x1, f1, 0, root_status::converged};
    }

    for (size_t iter = 1; iter <= max_iter; ++iter) {
        double denom = f1 - f0;
        if (std::abs(denom) <= std::numeric_limits<double>::epsilon()) {
            return {x1, f1, iter, root_status::zero_derivative};
        }

        double x2 = x1 - f1 * (x1 - x0) / denom;
        double f2 = std::invoke(func, x2);
        if (!std::isfinite(x2) || !std::isfinite(f2)) {
            return {x2, f2, iter, root_status::non_finite_value};
        }

        if (std::abs(f2) <= tol || std::abs(x2 - x1) <= tol) {
            return {x2, f2, iter, root_status::converged};
        }

        x0 = x1;
        f0 = f1;
        x1 = x2;
        f1 = f2;
    }

    return {x1, f1, max_iter, root_status::max_iterations};
}

/**
 * @brief Find a root with the secant method using root_options.
 *
 * @param f Function whose root is sought
 * @param x0 First initial value
 * @param x1 Second initial value
 * @param options Root-finding options
 * @return Root-finding result
 */
template <typename Func>
    requires std::is_invocable_r_v<double, Func, double>
inline root_result
secant(Func &&f, double x0, double x1, const root_options &options) {
    return secant(std::forward<Func>(f), x0, x1, options.tol, options.max_iter);
}

/**
 * @brief Convenience wrapper for the secant method that returns only the root.
 *
 * @throws std::runtime_error if the algorithm does not converge.
 */
template <typename Func>
    requires std::is_invocable_r_v<double, Func, double>
inline double secant_root(Func &&f,
                          double x0,
                          double x1,
                          double tol = 1e-12,
                          size_t max_iter = 100) {
    return secant(std::forward<Func>(f), x0, x1, tol, max_iter).value();
}

/**
 * @brief Convenience wrapper for the secant method with root_options.
 *
 * @throws std::runtime_error if the algorithm does not converge.
 */
template <typename Func>
    requires std::is_invocable_r_v<double, Func, double>
inline double
secant_root(Func &&f, double x0, double x1, const root_options &options) {
    return secant(std::forward<Func>(f), x0, x1, options).value();
}

} // namespace msl::equation

#endif // MSL_SECANT_HPP
