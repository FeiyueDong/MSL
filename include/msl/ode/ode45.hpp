/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2026, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: ode45.hpp
** -----
** Author: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_ODE45_HPP
#define MSL_ODE45_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include "detail/fixed_step.hpp"
#include "detail/state_ops.hpp"
#include "ode_options.hpp"
#include "ode_result.hpp"

namespace msl::ode {

namespace detail {
struct ode45_step_result {
    state y_next;
    double error_norm{};
};

inline double ode45_error_norm(const state &y,
                               const state &y_next,
                               const state &err,
                               const ode_options &options) {
    require_same_size(y, y_next);
    require_same_size(y, err);

    double sum = 0.0;
    for (size_t i = 0; i < y.size(); ++i) {
        double scale =
            options.atol
            + options.rtol * std::max(std::abs(y[i]), std::abs(y_next[i]));
        double ratio = err[i] / scale;
        sum += ratio * ratio;
    }

    return std::sqrt(sum / static_cast<double>(y.size()));
}

template <typename Func>
inline ode45_step_result ode45_step(Func &f,
                                    double t,
                                    const state &y,
                                    double h,
                                    const ode_options &options) {
    state k1 = eval_rhs(f, t, y);
    state y2 = add_scaled(y, k1, h * (1.0 / 5.0));
    state k2 = eval_rhs(f, t + h * (1.0 / 5.0), y2);

    state y3 = add_scaled(y, k1, h * (3.0 / 40.0), k2, h * (9.0 / 40.0));
    state k3 = eval_rhs(f, t + h * (3.0 / 10.0), y3);

    state y4(y.size());
    for (size_t i = 0; i < y.size(); ++i) {
        y4[i] = y[i]
                + h
                      * ((44.0 / 45.0) * k1[i] - (56.0 / 15.0) * k2[i]
                         + (32.0 / 9.0) * k3[i]);
    }
    state k4 = eval_rhs(f, t + h * (4.0 / 5.0), y4);

    state y5(y.size());
    for (size_t i = 0; i < y.size(); ++i) {
        y5[i] =
            y[i]
            + h
                  * ((19372.0 / 6561.0) * k1[i] - (25360.0 / 2187.0) * k2[i]
                     + (64448.0 / 6561.0) * k3[i] - (212.0 / 729.0) * k4[i]);
    }
    state k5 = eval_rhs(f, t + h * (8.0 / 9.0), y5);

    state y6(y.size());
    for (size_t i = 0; i < y.size(); ++i) {
        y6[i] = y[i]
                + h
                      * ((9017.0 / 3168.0) * k1[i] - (355.0 / 33.0) * k2[i]
                         + (46732.0 / 5247.0) * k3[i] + (49.0 / 176.0) * k4[i]
                         - (5103.0 / 18656.0) * k5[i]);
    }
    state k6 = eval_rhs(f, t + h, y6);

    state y_next(y.size());
    for (size_t i = 0; i < y.size(); ++i) {
        y_next[i] =
            y[i]
            + h
                  * ((35.0 / 384.0) * k1[i] + (500.0 / 1113.0) * k3[i]
                     + (125.0 / 192.0) * k4[i] - (2187.0 / 6784.0) * k5[i]
                     + (11.0 / 84.0) * k6[i]);
    }
    state k7 = eval_rhs(f, t + h, y_next);

    state y_fourth(y.size());
    for (size_t i = 0; i < y.size(); ++i) {
        y_fourth[i] =
            y[i]
            + h
                  * ((5179.0 / 57600.0) * k1[i] + (7571.0 / 16695.0) * k3[i]
                     + (393.0 / 640.0) * k4[i] - (92097.0 / 339200.0) * k5[i]
                     + (187.0 / 2100.0) * k6[i] + (1.0 / 40.0) * k7[i]);
    }

    state err(y.size());
    for (size_t i = 0; i < y.size(); ++i) {
        err[i] = y_next[i] - y_fourth[i];
    }

    double error_norm = ode45_error_norm(y, y_next, err, options);
    return {std::move(y_next), error_norm};
}

inline double
ode45_initial_step(double t0, double t1, const ode_options &options) {
    double span = t1 - t0;
    double h = options.initial_step > 0.0 ? options.initial_step : span / 100.0;
    if (options.max_step > 0.0) {
        h = std::min(h, options.max_step);
    }
    return std::max(h, options.min_step);
}

inline double
ode45_next_step(double h, double error_norm, const ode_options &options) {
    constexpr double safety = 0.9;
    constexpr double min_factor = 0.2;
    constexpr double max_factor = 5.0;

    double factor = max_factor;
    if (error_norm > 0.0) {
        factor = safety * std::pow(error_norm, -0.2);
        factor = std::clamp(factor, min_factor, max_factor);
    }

    double next_h = h * factor;
    if (options.max_step > 0.0) {
        next_h = std::min(next_h, options.max_step);
    }
    return next_h;
}

inline void validate_ode45_options(const ode_options &options) {
    if (options.rtol <= 0.0) {
        throw std::invalid_argument("ODE45: rtol must be positive");
    }
    if (options.atol <= 0.0) {
        throw std::invalid_argument("ODE45: atol must be positive");
    }
    if (options.min_step <= 0.0) {
        throw std::invalid_argument("ODE45: min_step must be positive");
    }
    if (options.max_step < 0.0) {
        throw std::invalid_argument("ODE45: max_step cannot be negative");
    }
    if (options.initial_step < 0.0) {
        throw std::invalid_argument("ODE45: initial_step cannot be negative");
    }
    if (options.max_steps == 0) {
        throw std::invalid_argument("ODE45: max_steps must be positive");
    }
}

template <typename Func>
inline ode_result ode45_solve_to(Func &&f,
                                 std::span<const double> y0,
                                 double t0,
                                 double t1,
                                 const ode_options &options,
                                 bool store_internal_steps) {
    validate_initial_value(y0);
    validate_ode45_options(options);
    if (!std::isfinite(t0) || !std::isfinite(t1)) {
        throw std::invalid_argument("ODE45: time bounds must be finite");
    }
    if (t0 >= t1) {
        throw std::invalid_argument(
            "ODE45: initial time must be less than final time");
    }

    auto func = std::forward<Func>(f);
    ode_result result;
    result.t.push_back(t0);
    result.y.push_back(to_state(y0));

    double t = t0;
    state y = result.y.back();
    double h = ode45_initial_step(t0, t1, options);

    while (t < t1) {
        if (result.steps + result.rejected_steps >= options.max_steps) {
            result.status = ode_status::max_steps;
            return result;
        }

        double remaining = t1 - t;
        h = std::min(h, remaining);
        if (h < options.min_step && remaining > options.min_step) {
            result.status = ode_status::step_size_underflow;
            return result;
        }

        auto step = ode45_step(func, t, y, h, options);
        if (!is_finite_state(step.y_next)) {
            result.t.push_back(t + h);
            result.y.push_back(std::move(step.y_next));
            ++result.steps;
            result.status = ode_status::non_finite_value;
            return result;
        }

        if (step.error_norm <= 1.0) {
            t += h;
            y = std::move(step.y_next);
            ++result.steps;
            if (store_internal_steps || t >= t1) {
                result.t.push_back(t);
                result.y.push_back(y);
            }
        } else {
            ++result.rejected_steps;
        }

        h = ode45_next_step(h, step.error_norm, options);
    }

    result.status = ode_status::success;
    return result;
}

} // namespace detail

/**
 * @brief Solve an initial value problem with ODE45.
 *
 * Uses the Dormand-Prince 5(4) embedded Runge-Kutta method with adaptive step
 * size control.
 *
 * @param f Right-hand side function f(t, y)
 * @param y0 Initial state
 * @param t0 Initial time
 * @param t1 Final time
 * @param options Solver options
 * @return Accepted internal time/state trajectory
 */
template <typename Func>
    requires std::is_invocable_r_v<state, Func, double, const state &>
inline ode_result ode45(Func &&f,
                        std::span<const double> y0,
                        double t0,
                        double t1,
                        const ode_options &options = {}) {
    return detail::ode45_solve_to(
        std::forward<Func>(f), y0, t0, t1, options, true);
}

/**
 * @brief Solve with ODE45 and return states at requested times.
 *
 * Internal adaptive steps are shortened to land exactly on each requested
 * output time. No dense-output interpolation is used in this first version.
 *
 * @param f Right-hand side function f(t, y)
 * @param y0 Initial state
 * @param t_eval Strictly increasing output times
 * @param options Solver options
 * @return Time/state trajectory at t_eval
 */
template <typename Func>
    requires std::is_invocable_r_v<state, Func, double, const state &>
inline ode_result ode45_eval(Func &&f,
                             std::span<const double> y0,
                             std::span<const double> t_eval,
                             const ode_options &options = {}) {
    detail::validate_eval_times(t_eval);

    auto func = std::forward<Func>(f);
    ode_result result;
    result.t.push_back(t_eval.front());
    result.y.push_back(detail::to_state(y0));

    state y = result.y.back();
    for (size_t i = 1; i < t_eval.size(); ++i) {
        ode_result segment = detail::ode45_solve_to(
            func, y, t_eval[i - 1], t_eval[i], options, false);
        result.steps += segment.steps;
        result.rejected_steps += segment.rejected_steps;
        if (!segment.success()) {
            result.status = segment.status;
            if (segment.t.size() > 1) {
                result.t.push_back(segment.t.back());
                result.y.push_back(segment.y.back());
            }
            return result;
        }

        y = segment.final_state();
        result.t.push_back(t_eval[i]);
        result.y.push_back(y);
    }

    result.status = ode_status::success;
    return result;
}

/**
 * @brief Solve with ODE45 and return only the final state.
 *
 * @throws std::runtime_error if the solver does not complete successfully.
 */
template <typename Func>
    requires std::is_invocable_r_v<state, Func, double, const state &>
inline state ode45_final(Func &&f,
                         std::span<const double> y0,
                         double t0,
                         double t1,
                         const ode_options &options = {}) {
    auto result = ode45(std::forward<Func>(f), y0, t0, t1, options);
    return result.value();
}

} // namespace msl::ode

#endif // MSL_ODE45_HPP
