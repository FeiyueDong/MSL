/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2026, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: regula_falsi.hpp
** -----
** Author: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_REGULA_FALSI_HPP
#define MSL_REGULA_FALSI_HPP

#include <algorithm>
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
 * @brief Find a root in a bracketing interval using the false-position method.
 *
 * @param f Function whose root is sought
 * @param a Lower interval endpoint
 * @param b Upper interval endpoint
 * @param tol Absolute tolerance for residual or interval width
 * @param max_iter Maximum number of iterations
 * @return Root-finding result
 */
template <typename Func>
    requires std::is_invocable_r_v<double, Func, double>
inline root_result regula_falsi(Func &&f,
                                double a,
                                double b,
                                double tol = 1e-12,
                                size_t max_iter = 100) {
    if (a >= b) {
        throw std::invalid_argument(
            "Regula falsi: lower endpoint must be less than upper endpoint");
    }
    if (tol <= 0.0) {
        throw std::invalid_argument("Regula falsi: tolerance must be positive");
    }
    if (max_iter == 0) {
        throw std::invalid_argument("Regula falsi: max_iter must be positive");
    }

    auto func = std::forward<Func>(f);
    double fa = std::invoke(func, a);
    double fb = std::invoke(func, b);
    if (!std::isfinite(fa) || !std::isfinite(fb)) {
        return {std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(),
                0,
                root_status::non_finite_value};
    }
    if (std::abs(fa) <= tol) {
        return {a, fa, 0, root_status::converged};
    }
    if (std::abs(fb) <= tol) {
        return {b, fb, 0, root_status::converged};
    }
    if (fa * fb > 0.0) {
        return {std::numeric_limits<double>::quiet_NaN(),
                std::min(std::abs(fa), std::abs(fb)),
                0,
                root_status::invalid_interval};
    }

    double x = a;
    double fx = fa;
    for (size_t iter = 1; iter <= max_iter; ++iter) {
        x = (a * fb - b * fa) / (fb - fa);
        fx = std::invoke(func, x);
        if (!std::isfinite(x) || !std::isfinite(fx)) {
            return {x, fx, iter, root_status::non_finite_value};
        }

        if (std::abs(fx) <= tol || std::abs(b - a) <= tol) {
            return {x, fx, iter, root_status::converged};
        }

        if (fa * fx <= 0.0) {
            b = x;
            fb = fx;
        } else {
            a = x;
            fa = fx;
        }
    }

    return {x, fx, max_iter, root_status::max_iterations};
}

/**
 * @brief Find a root with false-position using root_options.
 *
 * @param f Function whose root is sought
 * @param a Lower interval endpoint
 * @param b Upper interval endpoint
 * @param options Root-finding options
 * @return Root-finding result
 */
template <typename Func>
    requires std::is_invocable_r_v<double, Func, double>
inline root_result
regula_falsi(Func &&f, double a, double b, const root_options &options) {
    return regula_falsi(
        std::forward<Func>(f), a, b, options.tol, options.max_iter);
}

/**
 * @brief Convenience wrapper for false-position that returns only the root.
 *
 * @throws std::runtime_error if the algorithm does not converge.
 */
template <typename Func>
    requires std::is_invocable_r_v<double, Func, double>
inline double regula_falsi_root(Func &&f,
                                double a,
                                double b,
                                double tol = 1e-12,
                                size_t max_iter = 100) {
    return regula_falsi(std::forward<Func>(f), a, b, tol, max_iter).value();
}

/**
 * @brief Convenience wrapper for false-position with root_options.
 *
 * @throws std::runtime_error if the algorithm does not converge.
 */
template <typename Func>
    requires std::is_invocable_r_v<double, Func, double>
inline double
regula_falsi_root(Func &&f, double a, double b, const root_options &options) {
    return regula_falsi(std::forward<Func>(f), a, b, options).value();
}

} // namespace msl::equation

#endif // MSL_REGULA_FALSI_HPP
