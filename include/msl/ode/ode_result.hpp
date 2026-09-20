/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2026, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: ode_result.hpp
** -----
** Author: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_ODE_RESULT_HPP
#define MSL_ODE_RESULT_HPP

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace msl::ode {

using state = std::vector<double>;

/**
 * @brief Status code for ODE solvers.
 */
enum class ode_status {
    success,
    max_steps,
    step_size_underflow,
    non_finite_value,
};

/**
 * @brief Result of an ODE solve.
 */
struct ode_result {
    std::vector<double> t;
    std::vector<state> y;
    size_t steps{};
    size_t rejected_steps{};
    ode_status status{ode_status::success};

    [[nodiscard]] constexpr bool success() const {
        return status == ode_status::success;
    }

    [[nodiscard]] constexpr explicit operator bool() const { return success(); }

    [[nodiscard]] const state &final_state() const {
        if (y.empty()) {
            throw std::runtime_error("ODE result contains no states");
        }
        return y.back();
    }

    [[nodiscard]] double final_time() const {
        if (t.empty()) {
            throw std::runtime_error("ODE result contains no time values");
        }
        return t.back();
    }

    [[nodiscard]] const state &value() const {
        if (!success()) {
            throw std::runtime_error(
                "ODE solver did not complete successfully");
        }
        return final_state();
    }
};

} // namespace msl::ode

#endif // MSL_ODE_RESULT_HPP
