/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2026, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: adaptive_simpson.hpp
** -----
** Author: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_ADAPTIVE_SIMPSON_HPP
#define MSL_ADAPTIVE_SIMPSON_HPP

#include <cmath>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace msl::integral {

// ============================================================================
// Function Quadrature (Adaptive Simpson's Rule)
// ============================================================================

namespace detail {
template <typename Func>
inline double adaptive_simpson_recursive(Func &f,
                                         double a,
                                         double b,
                                         double fa,
                                         double fb,
                                         double fc,
                                         double whole,
                                         double tol,
                                         size_t max_depth,
                                         size_t depth) {
    double c = 0.5 * (a + b);
    double left_mid = 0.5 * (a + c);
    double right_mid = 0.5 * (c + b);

    double fd = std::invoke(f, left_mid);
    double fe = std::invoke(f, right_mid);

    double left = (c - a) / 6.0 * (fa + 4.0 * fd + fc);
    double right = (b - c) / 6.0 * (fc + 4.0 * fe + fb);
    double refined = left + right;
    double error = (refined - whole) / 15.0;

    if (std::abs(error) <= tol || depth >= max_depth) {
        return refined + error;
    }

    return adaptive_simpson_recursive(
               f, a, c, fa, fc, fd, left, 0.5 * tol, max_depth, depth + 1)
           + adaptive_simpson_recursive(
               f, c, b, fc, fb, fe, right, 0.5 * tol, max_depth, depth + 1);
}
} // namespace detail

/**
 * @brief Adaptive Simpson's quadrature for function
 *
 * @param f Function to integrate
 * @param a Lower limit
 * @param b Upper limit
 * @param tol Absolute tolerance
 * @param max_depth Maximum recursion depth
 * @return Integral value
 */
template <typename Func>
    requires std::is_invocable_r_v<double, Func, double>
inline double adaptive_simpson(Func &&f,
                               double a,
                               double b,
                               double tol = 1e-8,
                               size_t max_depth = 50) {
    if (a >= b) {
        throw std::invalid_argument(
            "Adaptive Simpson: lower limit must be less than upper limit");
    }
    if (tol <= 0.0) {
        throw std::invalid_argument(
            "Adaptive Simpson: tolerance must be positive");
    }
    if (max_depth == 0) {
        throw std::invalid_argument(
            "Adaptive Simpson: max_depth must be positive");
    }

    auto func = std::forward<Func>(f);
    if constexpr (std::is_pointer_v<std::decay_t<Func>>) {
        if (func == nullptr) {
            throw std::invalid_argument(
                "Adaptive Simpson: function cannot be null");
        }
    }

    double c = 0.5 * (a + b);
    double fa = std::invoke(func, a);
    double fb = std::invoke(func, b);
    double fc = std::invoke(func, c);
    double whole = (b - a) / 6.0 * (fa + 4.0 * fc + fb);

    return detail::adaptive_simpson_recursive(
        func, a, b, fa, fb, fc, whole, tol, max_depth, 0);
}

/**
 * @brief Short alias for adaptive Simpson's quadrature
 *
 * @param f Function to integrate
 * @param a Lower limit
 * @param b Upper limit
 * @param tol Absolute tolerance
 * @param max_depth Maximum recursion depth
 * @return Integral value
 */
template <typename Func>
    requires std::is_invocable_r_v<double, Func, double>
inline double
quad(Func &&f, double a, double b, double tol = 1e-8, size_t max_depth = 50) {
    return adaptive_simpson(std::forward<Func>(f), a, b, tol, max_depth);
}

} // namespace msl::integral

#endif // MSL_ADAPTIVE_SIMPSON_HPP
