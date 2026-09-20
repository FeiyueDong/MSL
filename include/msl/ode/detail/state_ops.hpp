/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2026, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: state_ops.hpp
** -----
** Author: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_ODE_STATE_OPS_HPP
#define MSL_ODE_STATE_OPS_HPP

#include <cmath>
#include <cstddef>
#include <functional>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "../ode_result.hpp"

namespace msl::ode::detail {

inline bool is_finite(double value) { return std::isfinite(value); }

inline bool is_finite_state(std::span<const double> y) {
    for (double value : y) {
        if (!std::isfinite(value)) {
            return false;
        }
    }
    return true;
}

inline void validate_initial_value(std::span<const double> y0) {
    if (y0.empty()) {
        throw std::invalid_argument("ODE: initial state cannot be empty");
    }
    if (!is_finite_state(y0)) {
        throw std::invalid_argument("ODE: initial state must be finite");
    }
}

inline state to_state(std::span<const double> y) {
    return state(y.begin(), y.end());
}

inline void require_same_size(std::span<const double> lhs,
                              std::span<const double> rhs) {
    if (lhs.size() != rhs.size()) {
        throw std::invalid_argument(
            "ODE: state sizes must match (lhs=" + std::to_string(lhs.size())
            + ", rhs=" + std::to_string(rhs.size()) + ")");
    }
}

template <typename Func>
inline state eval_rhs(Func &f, double t, const state &y) {
    state dy = std::invoke(f, t, y);
    if (dy.size() != y.size()) {
        throw std::invalid_argument(
            "ODE: derivative size must match state size");
    }
    return dy;
}

inline state add_scaled(const state &y, const state &dy, double scale) {
    if (y.size() != dy.size()) {
        throw std::invalid_argument(
            "ODE: add_scaled size mismatch (y=" + std::to_string(y.size())
            + ", dy=" + std::to_string(dy.size()) + ")");
    }
    state result(y.size());
    for (size_t i = 0; i < y.size(); ++i) {
        result[i] = y[i] + scale * dy[i];
    }
    return result;
}

inline state add_scaled(const state &y,
                        const state &k1,
                        double a1,
                        const state &k2,
                        double a2) {
    if (y.size() != k1.size() || y.size() != k2.size()) {
        throw std::invalid_argument(
            "ODE: add_scaled size mismatch (y=" + std::to_string(y.size())
            + ", k1=" + std::to_string(k1.size())
            + ", k2=" + std::to_string(k2.size()) + ")");
    }
    state result(y.size());
    for (size_t i = 0; i < y.size(); ++i) {
        result[i] = y[i] + a1 * k1[i] + a2 * k2[i];
    }
    return result;
}

inline state combine_rk4(const state &y,
                         const state &k1,
                         const state &k2,
                         const state &k3,
                         const state &k4,
                         double h) {
    require_same_size(y, k1);
    require_same_size(y, k2);
    require_same_size(y, k3);
    require_same_size(y, k4);

    state result(y.size());
    for (size_t i = 0; i < y.size(); ++i) {
        result[i] =
            y[i] + h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
    return result;
}

} // namespace msl::ode::detail

#endif // MSL_ODE_STATE_OPS_HPP
