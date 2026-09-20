/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2026, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: brent.hpp
** -----
** Author: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_BRENT_HPP
#define MSL_BRENT_HPP

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
 * @brief Find a root in a bracketing interval using Brent's method.
 *
 * Brent's method combines bisection, secant, and inverse quadratic
 * interpolation. The interval endpoints must bracket a root.
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
inline root_result
brent(Func &&f, double a, double b, double tol = 1e-12, size_t max_iter = 100) {
    if (a >= b) {
        throw std::invalid_argument(
            "Brent: lower endpoint must be less than upper endpoint");
    }
    if (tol <= 0.0) {
        throw std::invalid_argument("Brent: tolerance must be positive");
    }
    if (max_iter == 0) {
        throw std::invalid_argument("Brent: max_iter must be positive");
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

    if (std::abs(fa) < std::abs(fb)) {
        std::swap(a, b);
        std::swap(fa, fb);
    }

    double c = a;
    double fc = fa;
    double d = b - a;
    double e = d;
    double root = b;
    double residual = fb;

    for (size_t iter = 1; iter <= max_iter; ++iter) {
        if (fb * fc > 0.0) {
            c = a;
            fc = fa;
            d = b - a;
            e = d;
        }

        if (std::abs(fc) < std::abs(fb)) {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }

        double tolerance =
            2.0 * std::numeric_limits<double>::epsilon() * std::abs(b) + tol;
        double midpoint = 0.5 * (c - b);
        root = b;
        residual = fb;

        if (std::abs(fb) <= tol || std::abs(midpoint) <= tolerance) {
            return {root, residual, iter - 1, root_status::converged};
        }

        if (std::abs(e) >= tolerance && std::abs(fa) > std::abs(fb)) {
            double s = fb / fa;
            double p = 0.0;
            double q = 0.0;

            if (a == c) {
                p = 2.0 * midpoint * s;
                q = 1.0 - s;
            } else {
                q = fa / fc;
                double r = fb / fc;
                p = s * (2.0 * midpoint * q * (q - r) - (b - a) * (r - 1.0));
                q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            }

            if (p > 0.0) {
                q = -q;
            } else {
                p = -p;
            }

            double min1 = 3.0 * midpoint * q - std::abs(tolerance * q);
            double min2 = std::abs(e * q);
            if (2.0 * p < std::min(min1, min2)) {
                e = d;
                d = p / q;
            } else {
                d = midpoint;
                e = d;
            }
        } else {
            d = midpoint;
            e = d;
        }

        a = b;
        fa = fb;
        if (std::abs(d) > tolerance) {
            b += d;
        } else {
            b += midpoint > 0.0 ? tolerance : -tolerance;
        }

        fb = std::invoke(func, b);
        if (!std::isfinite(fb)) {
            return {b, fb, iter, root_status::non_finite_value};
        }
    }

    return {root, residual, max_iter, root_status::max_iterations};
}

/**
 * @brief Find a root with Brent's method using root_options.
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
brent(Func &&f, double a, double b, const root_options &options) {
    return brent(std::forward<Func>(f), a, b, options.tol, options.max_iter);
}

/**
 * @brief Convenience wrapper for Brent's method that returns only the root.
 *
 * @throws std::runtime_error if the algorithm does not converge.
 */
template <typename Func>
    requires std::is_invocable_r_v<double, Func, double>
inline double brent_root(Func &&f,
                         double a,
                         double b,
                         double tol = 1e-12,
                         size_t max_iter = 100) {
    return brent(std::forward<Func>(f), a, b, tol, max_iter).value();
}

/**
 * @brief Convenience wrapper for Brent's method with root_options.
 *
 * @throws std::runtime_error if the algorithm does not converge.
 */
template <typename Func>
    requires std::is_invocable_r_v<double, Func, double>
inline double
brent_root(Func &&f, double a, double b, const root_options &options) {
    return brent(std::forward<Func>(f), a, b, options).value();
}

} // namespace msl::equation

#endif // MSL_BRENT_HPP
