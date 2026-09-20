/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2026, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: root_result.hpp
** -----
** Author: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_ROOT_RESULT_HPP
#define MSL_ROOT_RESULT_HPP

#include <cstddef>
#include <stdexcept>

namespace msl::equation {

/**
 * @brief Status code for one-dimensional root-finding algorithms.
 */
enum class root_status {
    converged,
    max_iterations,
    invalid_interval,
    zero_derivative,
    non_finite_value,
};

/**
 * @brief Result of a one-dimensional root-finding algorithm.
 */
struct root_result {
    double root{};
    double residual{};
    size_t iterations{};
    root_status status{root_status::max_iterations};

    [[nodiscard]] constexpr bool converged() const {
        return status == root_status::converged;
    }

    [[nodiscard]] constexpr explicit operator bool() const {
        return converged();
    }

    /**
     * @brief Return the root if the algorithm converged.
     *
     * @throws std::runtime_error if the algorithm did not converge.
     */
    [[nodiscard]] double value() const {
        if (!converged()) {
            throw std::runtime_error("Root result did not converge");
        }
        return root;
    }

    [[nodiscard]] explicit operator double() const { return value(); }
};
} // namespace msl::equation

#endif // MSL_ROOT_RESULT_HPP
