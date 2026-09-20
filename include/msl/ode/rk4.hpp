/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2026, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: rk4.hpp
** -----
** Author: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_ODE_RK4_HPP
#define MSL_ODE_RK4_HPP

#include <cstddef>
#include <functional>
#include <span>
#include <type_traits>
#include <utility>

#include "detail/fixed_step.hpp"
#include "detail/state_ops.hpp"
#include "ode_options.hpp"
#include "ode_result.hpp"

namespace msl::ode {

/**
 * @brief Solve an initial value problem with classical fourth-order RK.
 *
 * @param f Right-hand side function f(t, y)
 * @param y0 Initial state
 * @param t0 Initial time
 * @param t1 Final time
 * @param dt Time step
 * @param options Solver options
 * @return Full time/state trajectory
 */
template <typename Func>
    requires std::is_invocable_r_v<state, Func, double, const state &>
inline ode_result rk4(Func &&f,
                      std::span<const double> y0,
                      double t0,
                      double t1,
                      double dt,
                      const ode_options &options = {}) {
    auto stepper = [](auto &func, double t, const state &y, double h) {
        state k1 = detail::eval_rhs(func, t, y);
        state y2 = detail::add_scaled(y, k1, 0.5 * h);
        state k2 = detail::eval_rhs(func, t + 0.5 * h, y2);
        state y3 = detail::add_scaled(y, k2, 0.5 * h);
        state k3 = detail::eval_rhs(func, t + 0.5 * h, y3);
        state y4 = detail::add_scaled(y, k3, h);
        state k4 = detail::eval_rhs(func, t + h, y4);
        return detail::combine_rk4(y, k1, k2, k3, k4, h);
    };

    return detail::fixed_step_solve(
        std::forward<Func>(f), y0, t0, t1, dt, options, stepper);
}

/**
 * @brief Solve with RK4 and return states at requested times.
 *
 * The first entry of t_eval is the initial time for y0. Internal steps use dt;
 * the last step before each output time is shortened to land exactly on that
 * output time.
 *
 * @param f Right-hand side function f(t, y)
 * @param y0 Initial state
 * @param t_eval Strictly increasing output times
 * @param dt Internal time step
 * @param options Solver options
 * @return Time/state trajectory at t_eval
 */
template <typename Func>
    requires std::is_invocable_r_v<state, Func, double, const state &>
inline ode_result rk4_eval(Func &&f,
                           std::span<const double> y0,
                           std::span<const double> t_eval,
                           double dt,
                           const ode_options &options = {}) {
    auto stepper = [](auto &func, double t, const state &y, double h) {
        state k1 = detail::eval_rhs(func, t, y);
        state y2 = detail::add_scaled(y, k1, 0.5 * h);
        state k2 = detail::eval_rhs(func, t + 0.5 * h, y2);
        state y3 = detail::add_scaled(y, k2, 0.5 * h);
        state k3 = detail::eval_rhs(func, t + 0.5 * h, y3);
        state y4 = detail::add_scaled(y, k3, h);
        state k4 = detail::eval_rhs(func, t + h, y4);
        return detail::combine_rk4(y, k1, k2, k3, k4, h);
    };

    return detail::fixed_step_solve_eval(
        std::forward<Func>(f), y0, t_eval, dt, options, stepper);
}

/**
 * @brief Solve with RK4 and return only the final state.
 *
 * @throws std::runtime_error if the solver does not complete successfully.
 */
template <typename Func>
    requires std::is_invocable_r_v<state, Func, double, const state &>
inline state rk4_final(Func &&f,
                       std::span<const double> y0,
                       double t0,
                       double t1,
                       double dt,
                       const ode_options &options = {}) {
    auto result = rk4(std::forward<Func>(f), y0, t0, t1, dt, options);
    return result.value();
}

} // namespace msl::ode

#endif // MSL_ODE_RK4_HPP
