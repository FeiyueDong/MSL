/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2026, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: ode_options.hpp
** -----
** Author: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_ODE_OPTIONS_HPP
#define MSL_ODE_OPTIONS_HPP

#include <cstddef>

namespace msl::ode {

/**
 * @brief Options shared by ODE solvers.
 *
 * Fixed-step solvers currently use max_steps. Tolerances and step-size bounds
 * are reserved for adaptive solvers such as RK45.
 */
struct ode_options {
    double rtol = 1e-6;
    double atol = 1e-9;
    double initial_step = 0.0;
    double min_step = 1e-12;
    double max_step = 0.0;
    size_t max_steps = 100000;
};

} // namespace msl::ode

#endif // MSL_ODE_OPTIONS_HPP
