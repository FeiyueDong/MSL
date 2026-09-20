/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2026, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: fixed_step.hpp
** -----
** Author: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_ODE_FIXED_STEP_HPP
#define MSL_ODE_FIXED_STEP_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <utility>

#include "../ode_options.hpp"
#include "../ode_result.hpp"
#include "state_ops.hpp"

namespace msl::ode::detail {

inline void validate_fixed_step_args(std::span<const double> y0,
                                     double t0,
                                     double t1,
                                     double dt,
                                     const ode_options &options) {
    validate_initial_value(y0);
    if (!std::isfinite(t0) || !std::isfinite(t1)) {
        throw std::invalid_argument("ODE: time bounds must be finite");
    }
    if (t0 >= t1) {
        throw std::invalid_argument(
            "ODE: initial time must be less than final time");
    }
    if (!std::isfinite(dt) || dt <= 0.0) {
        throw std::invalid_argument("ODE: time step must be positive");
    }
    if (options.max_steps == 0) {
        throw std::invalid_argument("ODE: max_steps must be positive");
    }
}

inline void validate_eval_times(std::span<const double> t_eval) {
    if (t_eval.size() < 2) {
        throw std::invalid_argument(
            "ODE: t_eval must contain at least two time points");
    }
    for (double t : t_eval) {
        if (!std::isfinite(t)) {
            throw std::invalid_argument("ODE: t_eval values must be finite");
        }
    }
    for (size_t i = 1; i < t_eval.size(); ++i) {
        if (t_eval[i] <= t_eval[i - 1]) {
            throw std::invalid_argument(
                "ODE: t_eval values must be strictly increasing");
        }
    }
}

template <typename Func, typename Stepper>
inline ode_result fixed_step_solve(Func &&f,
                                   std::span<const double> y0,
                                   double t0,
                                   double t1,
                                   double dt,
                                   const ode_options &options,
                                   Stepper &&stepper) {
    validate_fixed_step_args(y0, t0, t1, dt, options);

    auto func = std::forward<Func>(f);
    auto step = std::forward<Stepper>(stepper);

    ode_result result;
    result.t.push_back(t0);
    result.y.push_back(to_state(y0));

    double t = t0;
    state y = result.y.back();

    while (t < t1) {
        if (result.steps >= options.max_steps) {
            result.status = ode_status::max_steps;
            return result;
        }

        double h = std::min(dt, t1 - t);
        state y_next = step(func, t, y, h);
        double t_next = t + h;

        if (!std::isfinite(t_next) || !is_finite_state(y_next)) {
            result.t.push_back(t_next);
            result.y.push_back(std::move(y_next));
            ++result.steps;
            result.status = ode_status::non_finite_value;
            return result;
        }

        t = t_next;
        y = std::move(y_next);
        result.t.push_back(t);
        result.y.push_back(y);
        ++result.steps;
    }

    result.status = ode_status::success;
    return result;
}

template <typename Func, typename Stepper>
inline ode_result fixed_step_solve_eval(Func &&f,
                                        std::span<const double> y0,
                                        std::span<const double> t_eval,
                                        double dt,
                                        const ode_options &options,
                                        Stepper &&stepper) {
    validate_eval_times(t_eval);
    validate_fixed_step_args(y0, t_eval.front(), t_eval.back(), dt, options);

    auto func = std::forward<Func>(f);
    auto step = std::forward<Stepper>(stepper);

    ode_result result;
    result.t.push_back(t_eval.front());
    result.y.push_back(to_state(y0));

    double t = t_eval.front();
    state y = result.y.back();

    for (size_t i = 1; i < t_eval.size(); ++i) {
        double target = t_eval[i];
        while (t < target) {
            if (result.steps >= options.max_steps) {
                result.status = ode_status::max_steps;
                return result;
            }

            double h = std::min(dt, target - t);
            state y_next = step(func, t, y, h);
            double t_next = t + h;

            if (!std::isfinite(t_next) || !is_finite_state(y_next)) {
                result.t.push_back(t_next);
                result.y.push_back(std::move(y_next));
                ++result.steps;
                result.status = ode_status::non_finite_value;
                return result;
            }

            t = t_next;
            y = std::move(y_next);
            ++result.steps;
        }

        t = target;
        result.t.push_back(t);
        result.y.push_back(y);
    }

    result.status = ode_status::success;
    return result;
}

} // namespace msl::ode::detail

#endif // MSL_ODE_FIXED_STEP_HPP
