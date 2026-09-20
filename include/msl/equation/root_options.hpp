/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2026, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: root_options.hpp
** -----
** Author: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_ROOT_OPTIONS_HPP
#define MSL_ROOT_OPTIONS_HPP

#include <cstddef>

namespace msl::equation {

/**
 * @brief Options shared by one-dimensional root-finding algorithms.
 */
struct root_options {
    double tol = 1e-12;
    size_t max_iter = 100;
    double derivative_tol = 1e-14;
};

} // namespace msl::equation

#endif // MSL_ROOT_OPTIONS_HPP
