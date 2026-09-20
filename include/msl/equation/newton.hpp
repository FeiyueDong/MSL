/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2026, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: newton.hpp
** -----
** Author: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_NEWTON_HPP
#define MSL_NEWTON_HPP

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
 * @brief Find a root using Newton's method with an explicit derivative.
 *
 * @param f Function whose root is sought
 * @param df Derivative of f
 * @param x0 Initial guess
 * @param tol Absolute tolerance for residual or step size
 * @param max_iter Maximum number of iterations
 * @param derivative_tol Minimum acceptable derivative magnitude
 * @return Root-finding result
 */
template <typename Func, typename Deriv>
    requires std::is_invocable_r_v<double, Func, double>
             && std::is_invocable_r_v<double, Deriv, double>
inline root_result newton(Func &&f,
                          Deriv &&df,
                          double x0,
                          double tol = 1e-12,
                          size_t max_iter = 50,
                          double derivative_tol = 1e-14) {
    if (tol <= 0.0) {
        throw std::invalid_argument("Newton: tolerance must be positive");
    }
    if (max_iter == 0) {
        throw std::invalid_argument("Newton: max_iter must be positive");
    }
    if (derivative_tol <= 0.0) {
        throw std::invalid_argument(
            "Newton: derivative tolerance must be positive");
    }

    auto func = std::forward<Func>(f);
    auto deriv = std::forward<Deriv>(df);
    double x = x0;
    double fx = std::invoke(func, x);
    if (!std::isfinite(x) || !std::isfinite(fx)) {
        return {x, fx, 0, root_status::non_finite_value};
    }
    if (std::abs(fx) <= tol) {
        return {x, fx, 0, root_status::converged};
    }

    for (size_t iter = 1; iter <= max_iter; ++iter) {
        double dfx = std::invoke(deriv, x);
        if (!std::isfinite(dfx)) {
            return {x, fx, iter, root_status::non_finite_value};
        }
        if (std::abs(dfx) <= derivative_tol) {
            return {x, fx, iter, root_status::zero_derivative};
        }

        double x_next = x - fx / dfx;
        double f_next = std::invoke(func, x_next);
        if (!std::isfinite(x_next) || !std::isfinite(f_next)) {
            return {x_next, f_next, iter, root_status::non_finite_value};
        }

        if (std::abs(f_next) <= tol || std::abs(x_next - x) <= tol) {
            return {x_next, f_next, iter, root_status::converged};
        }

        x = x_next;
        fx = f_next;
    }

    return {x, fx, max_iter, root_status::max_iterations};
}

/**
 * @brief Find a root with Newton's method using root_options.
 *
 * @param f Function whose root is sought
 * @param df Derivative of f
 * @param x0 Initial guess
 * @param options Root-finding options
 * @return Root-finding result
 */
template <typename Func, typename Deriv>
    requires std::is_invocable_r_v<double, Func, double>
             && std::is_invocable_r_v<double, Deriv, double>
inline root_result
newton(Func &&f, Deriv &&df, double x0, const root_options &options) {
    return newton(std::forward<Func>(f),
                  std::forward<Deriv>(df),
                  x0,
                  options.tol,
                  options.max_iter,
                  options.derivative_tol);
}

/**
 * @brief Convenience wrapper for Newton's method that returns only the root.
 *
 * @throws std::runtime_error if the algorithm does not converge.
 */
template <typename Func, typename Deriv>
    requires std::is_invocable_r_v<double, Func, double>
             && std::is_invocable_r_v<double, Deriv, double>
inline double newton_root(Func &&f,
                          Deriv &&df,
                          double x0,
                          double tol = 1e-12,
                          size_t max_iter = 50,
                          double derivative_tol = 1e-14) {
    return newton(std::forward<Func>(f),
                  std::forward<Deriv>(df),
                  x0,
                  tol,
                  max_iter,
                  derivative_tol)
        .value();
}

/**
 * @brief Convenience wrapper for Newton's method with root_options.
 *
 * @throws std::runtime_error if the algorithm does not converge.
 */
template <typename Func, typename Deriv>
    requires std::is_invocable_r_v<double, Func, double>
             && std::is_invocable_r_v<double, Deriv, double>
inline double
newton_root(Func &&f, Deriv &&df, double x0, const root_options &options) {
    return newton(std::forward<Func>(f), std::forward<Deriv>(df), x0, options)
        .value();
}

} // namespace msl::equation

#endif // MSL_NEWTON_HPP
