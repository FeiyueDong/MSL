/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: window.hpp
** -----
** File Created: Friday, 7th November 2025 17:53:39
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Sunday, 14th December 2025 17:05:44
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_WINDOW_HPP
#define MSL_WINDOW_HPP

#include <cmath>
#include <numbers>
#include <span>
#include <stdexcept>
#include <vector>

#include "matrix/real_matrix_owned.hpp"

namespace msl::signal {
// ========================================================================
// 1. Generate window functions
// ========================================================================
/**
 * @brief Generate Hamming window coefficients.
 *
 * @param window Output coefficient buffer.
 * @param n Window length.
 */
inline void hamming_window(std::span<double> window, size_t n) {
    if (window.size() != n) {
        throw std::invalid_argument(
            "Hamming window: output buffer size must be equal to n");
    }
    if (n <= 1) {
        if (n == 1) {
            window[0] = 1.0;
        }
        return;
    }

    for (size_t i = 0; i < n; ++i) {
        window[i] =
            0.54 - 0.46 * std::cos(2.0 * std::numbers::pi * i / (n - 1));
    }
}

/**
 * @brief Generate and return a Hamming window.
 *
 * @param n Window length.
 * @return Hamming window coefficients.
 */
inline std::vector<double> hamming_window(size_t n) {
    std::vector<double> window(n);
    hamming_window(window, n);
    return window;
}

/**
 * @brief Generate Hann window coefficients.
 *
 * @param window Output coefficient buffer.
 * @param n Window length.
 */
inline void hann_window(std::span<double> window, size_t n) {
    if (window.size() != n) {
        throw std::invalid_argument(
            "Hann window: output buffer size must be equal to n");
    }
    if (n <= 1) {
        if (n == 1) {
            window[0] = 1.0;
        }
        return;
    }

    for (size_t i = 0; i < n; ++i) {
        window[i] =
            0.5 * (1.0 - std::cos(2.0 * std::numbers::pi * i / (n - 1)));
    }
}

/**
 * @brief Generate and return a Hann window.
 *
 * @param n Window length.
 * @return Hann window coefficients.
 */
inline std::vector<double> hann_window(size_t n) {
    std::vector<double> window(n);
    hann_window(window, n);
    return window;
}

/**
 * @brief Generate Blackman window coefficients.
 *
 * @param window Output coefficient buffer.
 * @param n Window length.
 */
inline void blackman_window(std::span<double> window, size_t n) {
    if (window.size() != n) {
        throw std::invalid_argument(
            "Blackman window: output buffer size must be equal to n");
    }
    if (n <= 1) {
        if (n == 1) {
            window[0] = 1.0;
        }
        return;
    }

    for (size_t i = 0; i < n; ++i) {
        window[i] = 0.42 - 0.5 * std::cos(2.0 * std::numbers::pi * i / (n - 1))
                    + 0.08 * std::cos(4.0 * std::numbers::pi * i / (n - 1));
    }
}

/**
 * @brief Generate and return a Blackman window.
 *
 * @param n Window length.
 * @return Blackman window coefficients.
 */
inline std::vector<double> blackman_window(size_t n) {
    std::vector<double> window(n);
    blackman_window(window, n);
    return window;
}

/**
 * @brief Generate rectangular window coefficients.
 *
 * @param window Output coefficient buffer.
 * @param n Window length.
 */
inline void rectangular_window(std::span<double> window, size_t n) {
    if (window.size() != n) {
        throw std::invalid_argument(
            "Rectangular window: output buffer size must be equal to n");
    }
    for (size_t i = 0; i < n; ++i) {
        window[i] = 1.0;
    }
}

/**
 * @brief Generate and return a rectangular window.
 *
 * @param n Window length.
 * @return Rectangular window coefficients.
 */
inline std::vector<double> rectangular_window(size_t n) {
    std::vector<double> window(n);
    rectangular_window(window, n);
    return window;
}

// ========================================================================
// 2. Apply window functions
// ========================================================================
/**
 * @brief Apply a window to a signal in-place.
 *
 * @param signal Input/output signal.
 * @param window Window coefficients, same length as signal.
 */
inline void apply_window(std::span<double> signal,
                         std::span<const double> window) {
    if (signal.size() != window.size()) {
        throw std::invalid_argument(
            "Apply window: signal and window must have the same length");
    }

    for (size_t i = 0; i < signal.size(); ++i) {
        signal[i] *= window[i];
    }
}

/**
 * @brief Apply a window and write to a separate output buffer.
 *
 * @param input Input signal.
 * @param output Windowed output signal.
 * @param window Window coefficients.
 */
inline void apply_window(std::span<const double> input,
                         std::span<double> output,
                         std::span<const double> window) {
    if (input.size() != output.size() || input.size() != window.size()) {
        throw std::invalid_argument("Apply window: input, output, and window "
                                    "must have the same length");
    }

    for (size_t i = 0; i < input.size(); ++i) {
        output[i] = input[i] * window[i];
    }
}

/**
 * @brief Apply a window to each matrix column.
 *
 * Window length must equal matrix row count. Element `(i, j)` is multiplied
 * by `window[i]`.
 *
 * @param input Input matrix.
 * @param window Window coefficients.
 * @return Windowed matrix.
 */
inline matrix::matrixd
apply_window_columns(const matrix::real_matrix_base &input,
                     std::span<const double> window) {
    size_t n_rows = input.rows();
    size_t n_cols = input.cols();

    if (n_rows != window.size()) {
        throw std::invalid_argument(
            "Apply window: number of rows in input must match window length");
    }

    matrix::matrixd output(n_rows, n_cols);

    for (size_t j = 0; j < n_cols; ++j) {
        for (size_t i = 0; i < n_rows; ++i) {
            output(i, j) = input(i, j) * window[i];
        }
    }

    return output;
}

} // namespace msl::signal

#endif // MSL_WINDOW_HPP