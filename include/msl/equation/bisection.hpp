/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2026, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: bisection.hpp
** -----
** Author: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_BISECTION_HPP
#define MSL_BISECTION_HPP

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
 * @brief Find a root in a bracketing interval using the bisection method.
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
inline root_result bisection(Func &&f,
                             double a,
                             double b,
                             double tol = 1e-12,
                             size_t max_iter = 100) {
    if (a >= b) {
        throw std::invalid_argument(
            "Bisection: lower endpoint must be less than upper endpoint");
    }
    if (tol <= 0.0) {
        throw std::invalid_argument("Bisection: tolerance must be positive");
    }
    if (max_iter == 0) {
        throw std::invalid_argument("Bisection: max_iter must be positive");
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

    double mid = a;
    double fm = fa;
    for (size_t iter = 1; iter <= max_iter; ++iter) {
        mid = 0.5 * (a + b);
        fm = std::invoke(func, mid);
        if (!std::isfinite(fm)) {
            return {mid, fm, iter, root_status::non_finite_value};
        }

        if (std::abs(fm) <= tol || 0.5 * (b - a) <= tol) {
            return {mid, fm, iter, root_status::converged};
        }

        if (fa * fm <= 0.0) {
            b = mid;
            fb = fm;
        } else {
            a = mid;
            fa = fm;
        }
    }

    return {mid, fm, max_iter, root_status::max_iterations};
}

/**
 * @brief Find a root with bisection using root_options.
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
bisection(Func &&f, double a, double b, const root_options &options) {
    return bisection(
        std::forward<Func>(f), a, b, options.tol, options.max_iter);
}

/**
 * @brief Convenience wrapper for bisection that returns only the root.
 *
 * @throws std::runtime_error if the algorithm does not converge.
 */
template <typename Func>
    requires std::is_invocable_r_v<double, Func, double>
inline double bisection_root(Func &&f,
                             double a,
                             double b,
                             double tol = 1e-12,
                             size_t max_iter = 100) {
    return bisection(std::forward<Func>(f), a, b, tol, max_iter).value();
}

/**
 * @brief Convenience wrapper for bisection with root_options.
 *
 * @throws std::runtime_error if the algorithm does not converge.
 */
template <typename Func>
    requires std::is_invocable_r_v<double, Func, double>
inline double
bisection_root(Func &&f, double a, double b, const root_options &options) {
    return bisection(std::forward<Func>(f), a, b, options).value();
}

} // namespace msl::equation

#endif // MSL_BISECTION_HPP
