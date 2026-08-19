/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2026, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: fft.hpp
** -----
** File Created: Friday, 9th January 2026 14:58:10
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Friday, 6th March 2026 10:49:32
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_FFT_HPP
#define MSL_FFT_HPP

#include <algorithm>
#include <complex>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <vector>

#include <eigen3/unsupported/Eigen/FFT>

#include "matrix/complex_matrix_base.hpp"
#include "matrix/complex_matrix_owned.hpp"
#include "matrix/real_matrix_base.hpp"
#include "matrix/real_matrix_owned.hpp"

namespace msl::signal {

// ============================================================================
// 1. 1D FFT for vectors
// ============================================================================

/**
 * @brief Compute forward FFT from real input into caller-provided buffer.
 *
 * @param input Real-valued input samples.
 * @param output Complex output buffer, size must be `nfft`.
 * @param nfft FFT length. `0` means `input.size()`.
 */
inline void fft(std::span<const double> input,
                std::span<std::complex<double>> output,
                std::size_t nfft = 0) {
    nfft = (nfft == 0) ? input.size() : nfft;
    if (output.size() != nfft) {
        throw std::invalid_argument(
            "FFT: output buffer size must be equal to nfft");
    }
    Eigen::FFT<double> fft_engine;
    if (input.size() < nfft) {
        std::vector<double> padded(nfft, 0.0);
        std::copy(input.begin(), input.end(), padded.begin());
        fft_engine.fwd(output.data(), padded.data(), nfft);
    } else {
        fft_engine.fwd(output.data(), input.data(), nfft);
    }
}

/**
 * @brief 1D Forward FFT: real vector -> complex vector
 *
 * @param input Real-valued input signal
 * @param nfft FFT length. `0` means `input.size()`.
 * @return Complex-valued frequency domain (size determined by nfft)
 */
inline std::vector<std::complex<double>> fft(std::span<const double> input,
                                             std::size_t nfft = 0) {
    const std::size_t fft_len = (nfft == 0) ? input.size() : nfft;
    std::vector<std::complex<double>> output(fft_len);
    fft(input, output, fft_len);
    return output;
}

/**
 * @brief Compute forward FFT from complex input into caller-provided buffer.
 *
 * Full complex-to-complex FFT (no symmetry assumptions)
 *
 * @param input Complex-valued input samples.
 * @param output Complex output buffer, size must be `nfft`.
 * @param nfft FFT length. `0` means `input.size()`.
 */
inline void fft(std::span<const std::complex<double>> input,
                std::span<std::complex<double>> output,
                std::size_t nfft = 0) {
    nfft = (nfft == 0) ? input.size() : nfft;
    if (output.size() != nfft) {
        throw std::invalid_argument(
            "FFT: output buffer size must be equal to nfft");
    }
    Eigen::FFT<double> fft_engine;
    if (input.size() < nfft) {
        std::vector<std::complex<double>> padded(nfft, {0.0, 0.0});
        std::copy(input.begin(), input.end(), padded.begin());
        fft_engine.fwd(output.data(), padded.data(), nfft);
    } else {
        fft_engine.fwd(output.data(), input.data(), nfft);
    }
}

/**
 * @brief 1D Forward FFT: complex vector -> complex vector
 *
 * Full complex-to-complex FFT (no symmetry assumptions)
 *
 * @param input Complex-valued input signal
 * @return Complex-valued frequency domain (same size as input)
 */
inline std::vector<std::complex<double>>
fft(std::span<const std::complex<double>> input, std::size_t nfft = 0) {
    const std::size_t fft_len = (nfft == 0) ? input.size() : nfft;
    std::vector<std::complex<double>> output(fft_len);
    fft(input, output, fft_len);
    return output;
}

// ============================================================================
// 2. 1D Inverse FFT for vectors
// ============================================================================

/**
 * @brief Compute inverse FFT into caller-provided complex buffer.
 *
 * @param input Complex frequency-domain input.
 * @param output Complex time-domain output buffer, size must be `nfft`.
 * @param nfft IFFT length. `0` means `input.size()`.
 */
inline void ifft(std::span<const std::complex<double>> input,
                 std::span<std::complex<double>> output,
                 std::size_t nfft = 0) {
    nfft = (nfft == 0) ? input.size() : nfft;
    if (output.size() != nfft) {
        throw std::invalid_argument(
            "IFFT: output buffer size must be equal to nfft");
    }
    if (input.size() < nfft) {
        throw std::invalid_argument(
            "IFFT: FFT length nfft cannot be greater than input size");
    }
    Eigen::FFT<double> fft_engine;

    fft_engine.inv(output.data(), input.data(), nfft);
}

/**
 * @brief 1D Inverse FFT: complex vector -> complex vector
 *
 * @param input Complex frequency domain data
 * @param nfft IFFT length. `0` means `input.size()`.
 * @return Complex time domain signal
 */
inline std::vector<std::complex<double>>
ifft(std::span<const std::complex<double>> input, std::size_t nfft = 0) {
    const std::size_t fft_len = (nfft == 0) ? input.size() : nfft;
    std::vector<std::complex<double>> output(fft_len);
    ifft(input, output, fft_len);
    return output;
}

/**
 * @brief Compute inverse FFT into caller-provided real buffer.
 *
 * Assumes the input spectrum corresponds to a real-valued time signal
 * (conjugate symmetry).
 *
 * @param input Complex frequency-domain input.
 * @param output Real time-domain output buffer, size must be `nfft`.
 * @param nfft IFFT length. `0` means `input.size()`.
 */
inline void ifft_real(std::span<const std::complex<double>> input,
                      std::span<double> output,
                      std::size_t nfft = 0) {
    nfft = (nfft == 0) ? input.size() : nfft;
    if (output.size() != nfft) {
        throw std::invalid_argument(
            "IFFT: output buffer size must be equal to nfft");
    }
    if (input.size() < nfft) {
        throw std::invalid_argument(
            "IFFT: FFT length nfft cannot be greater than input size");
    }
    Eigen::FFT<double> fft_engine;

    fft_engine.inv(output.data(), input.data(), nfft);
}
/**
 * @brief 1D Inverse FFT: complex vector -> real vector
 *
 * Assumes input has conjugate symmetry (from real FFT)
 *
 * @param input Complex frequency domain data
 * @param nfft IFFT length. `0` means `input.size()`.
 * @return Real-valued time domain signal
 */
inline std::vector<double>
ifft_real(std::span<const std::complex<double>> input, std::size_t nfft = 0) {
    const std::size_t fft_len = (nfft == 0) ? input.size() : nfft;
    std::vector<double> output(fft_len);
    ifft_real(input, output, fft_len);
    return output;
}

// ============================================================================
// 3. 2D FFT for matrices (column-wise FFT)
// ============================================================================

/**
 * @brief 2D FFT: real matrix -> complex matrix (column-wise)
 *
 * Applies 1D FFT to each column of the matrix independently.
 * Useful for processing multiple signals simultaneously.
 *
 * @param input Real-valued matrix (each column is a signal)
 * @param nfft FFT length per column; `0` means use `input.rows()`.
 * @return Complex matrix with FFT of each column
 */
inline matrix::matrixc fft_columns(const matrix::real_matrix_base &input,
                                   size_t nfft = 0) {
    size_t n_rows = input.rows();
    size_t n_cols = input.cols();

    nfft = (nfft == 0) ? n_rows : nfft;
    matrix::matrixc output(nfft, n_cols);

    // Process each column
    for (size_t j = 0; j < n_cols; ++j) {
        auto col_span = input.column(j);
        auto transformed = fft(col_span, nfft);
        auto col_fft = output.column(j);
        std::copy(transformed.begin(), transformed.end(), col_fft.begin());
    }

    return output;
}

/**
 * @brief 2D FFT: complex matrix -> complex matrix (column-wise)
 *
 * Full complex FFT for each column
 *
 * @param input Complex-valued matrix (each column is a signal)
 * @param nfft FFT length per column; `0` means use `input.rows()`.
 * @return Complex matrix with FFT of each column
 */
inline matrix::matrixc fft_columns(const matrix::complex_matrix_base &input,
                                   size_t nfft = 0) {
    size_t n_rows = input.rows();
    size_t n_cols = input.cols();

    nfft = (nfft == 0) ? n_rows : nfft;
    matrix::matrixc output(nfft, n_cols);

    for (size_t j = 0; j < n_cols; ++j) {
        auto col_span = input.column(j);
        auto transformed = fft(col_span, nfft);
        auto col_fft = output.column(j);
        std::copy(transformed.begin(), transformed.end(), col_fft.begin());
    }

    return output;
}

// ============================================================================
// 4. 2D Inverse FFT for matrices (column-wise)
// ============================================================================

/**
 * @brief 2D Inverse FFT: complex matrix -> complex matrix (column-wise)
 *
 * @param input Complex-valued matrix (each column is a signal in frequency
 * domain)
 * @param nfft IFFT length per column; `0` means use `input.rows()`.
 * @return Complex matrix with IFFT of each column (time domain)
 */
inline matrix::matrixc ifft_columns(const matrix::complex_matrix_base &input,
                                    size_t nfft = 0) {
    size_t n_rows = input.rows();
    size_t n_cols = input.cols();

    nfft = (nfft == 0) ? n_rows : nfft;
    if (n_rows < nfft) {
        throw std::invalid_argument(
            "IFFT columns: nfft cannot be greater than input rows");
    }

    matrix::matrixc output(nfft, n_cols);

    Eigen::FFT<double> fft_engine;

    for (size_t j = 0; j < n_cols; ++j) {
        auto col_span = input.column(j);
        auto col_ifft = output.column(j);
        fft_engine.inv(col_ifft.data(), col_span.data(), nfft);
    }

    return output;
}

/**
 * @brief 2D Inverse FFT: complex matrix -> real matrix (column-wise)
 *
 * Assumes input has conjugate symmetry (from real FFT)
 *
 * @param input Complex frequency domain matrix
 * @param nfft Inverse FFT length per column; `0` means use `input.rows()`.
 * @return Real-valued time domain matrix
 */
inline matrix::matrixd
ifft_columns_real(const matrix::complex_matrix_base &input, size_t nfft = 0) {
    size_t n_rows = input.rows();
    size_t n_cols = input.cols();

    nfft = (nfft == 0) ? n_rows : nfft;
    if (n_rows < nfft) {
        throw std::invalid_argument(
            "IFFT columns real: nfft cannot be greater than input rows");
    }

    matrix::matrixd output(nfft, n_cols);

    Eigen::FFT<double> fft_engine;

    for (size_t j = 0; j < n_cols; ++j) {
        auto col_span = input.column(j);
        auto col_ifft = output.column(j);
        fft_engine.inv(col_ifft.data(), col_span.data(), nfft);
    }

    return output;
}

// ============================================================================
// 5. 2D FFT for matrices (row-wise FFT)
// ============================================================================

/**
 * @brief 2D FFT: real matrix -> complex matrix (row-wise)
 *
 * Applies 1D FFT to each row of the matrix independently.
 *
 * @param input Real-valued matrix (each row is a signal)
 * @param nfft FFT length per row; `0` means use `input.cols()`.
 * @return Complex matrix with FFT of each row
 */
inline matrix::matrixc fft_rows(const matrix::real_matrix_base &input,
                                size_t nfft = 0) {
    size_t n_rows = input.rows();
    size_t n_cols = input.cols();

    nfft = (nfft == 0) ? n_cols : nfft;
    matrix::matrixc output(n_rows, nfft);

    for (size_t i = 0; i < n_rows; ++i) {
        std::vector<double> row_vec(n_cols);
        for (size_t j = 0; j < n_cols; ++j) {
            row_vec[j] = input(i, j);
        }
        auto transformed = fft(row_vec, nfft);
        for (size_t j = 0; j < nfft; ++j) {
            output(i, j) = transformed[j];
        }
    }

    return output;
}

/**
 * @brief 2D FFT: complex matrix -> complex matrix (row-wise)
 *
 * Applies full complex FFT to each row independently.
 *
 * @param input Complex-valued matrix (each row is a signal)
 * @param nfft FFT length per row; `0` means use `input.cols()
 * @return Complex matrix with FFT of each row
 */
inline matrix::matrixc fft_rows(const matrix::complex_matrix_base &input,
                                size_t nfft = 0) {
    size_t n_rows = input.rows();
    size_t n_cols = input.cols();

    nfft = (nfft == 0) ? n_cols : nfft;
    matrix::matrixc output(n_rows, nfft);

    for (size_t i = 0; i < n_rows; ++i) {
        std::vector<std::complex<double>> row_vec(n_cols);
        for (size_t j = 0; j < n_cols; ++j) {
            row_vec[j] = input(i, j);
        }
        auto transformed = fft(row_vec, nfft);
        for (size_t j = 0; j < nfft; ++j) {
            output(i, j) = transformed[j];
        }
    }

    return output;
}

/**
 * @brief 2D Inverse FFT: complex matrix -> complex matrix (row-wise)
 *
 * Full complex IFFT for each row independently.
 *
 * @param input Complex-valued matrix (each row is a signal in frequency
 * domain)
 * @param nfft IFFT length per row; `0` means use `input
 * @return Complex matrix with IFFT of each row (time domain)
 */
inline matrix::matrixc ifft_rows(const matrix::complex_matrix_base &input,
                                 size_t nfft = 0) {
    size_t n_rows = input.rows();
    size_t n_cols = input.cols();

    nfft = (nfft == 0) ? n_cols : nfft;
    if (n_cols < nfft) {
        throw std::invalid_argument(
            "IFFT rows: nfft cannot be greater than input columns");
    }

    matrix::matrixc output(n_rows, nfft);

    Eigen::FFT<double> fft_engine;

    for (size_t i = 0; i < n_rows; ++i) {
        std::vector<std::complex<double>> row_vec(input.cols());
        for (size_t j = 0; j < n_cols; ++j) {
            row_vec[j] = input(i, j);
        }
        auto row_ifft = output.get_row(i);
        fft_engine.inv(row_ifft.data(), row_vec.data(), nfft);
    }

    return output;
}

/**
 * @brief 2D Inverse FFT: complex matrix -> real matrix (row-wise)
 *
 * Assumes input has conjugate symmetry (from real FFT)
 *
 * @param input Complex frequency domain matrix
 * @param nfft Inverse FFT length per row; `0` means use `input.cols()`.
 * @return Real-valued time domain matrix
 */
inline matrix::matrixd ifft_rows_real(const matrix::complex_matrix_base &input,
                                      size_t nfft = 0) {
    size_t n_rows = input.rows();
    size_t n_cols = input.cols();

    nfft = (nfft == 0) ? n_cols : nfft;
    if (n_cols < nfft) {
        throw std::invalid_argument(
            "IFFT rows real: nfft cannot be greater than input columns");
    }

    matrix::matrixd output(n_rows, nfft);

    Eigen::FFT<double> fft_engine;

    for (size_t i = 0; i < n_rows; ++i) {
        std::vector<std::complex<double>> row_vec(input.cols());
        for (size_t j = 0; j < n_cols; ++j) {
            row_vec[j] = input(i, j);
        }
        auto row_ifft = output.get_row(i);
        fft_engine.inv(row_ifft.data(), row_vec.data(), nfft);
    }

    return output;
}

// ============================================================================
// 6. Utility functions
// ============================================================================

/**
 * @brief Fill frequency bins for an FFT grid.
 *
 * Frequencies are generated as `k * sample_rate / n` for `k = 0..n-1`.
 *
 * @param freqs Output frequency buffer.
 * @param n Number of FFT points.
 * @param sample_rate Sampling rate in Hz, must be greater than 0.
 */
inline void
fft_frequencies(std::span<double> freqs, size_t n, double sample_rate = 1.0) {
    if (freqs.size() != n) {
        throw std::invalid_argument(
            "FFT frequencies: output buffer size must be equal to n");
    }
    if (n == 0) {
        return;
    }
    if (sample_rate <= 0.0) {
        throw std::invalid_argument(
            "FFT frequencies: sample_rate must be greater than 0");
    }
    double df = sample_rate / n;

    for (size_t i = 0; i < freqs.size(); ++i) {
        freqs[i] = i * df;
    }
}

/**
 * @brief Compute FFT frequencies for a given sample rate
 *
 * @param n Number of samples
 * @param sample_rate Sampling rate (Hz)
 * @return Vector of frequency values
 */
inline std::vector<double> fft_frequencies(size_t n, double sample_rate = 1.0) {
    std::vector<double> freqs(n);
    fft_frequencies(freqs, n, sample_rate);
    return freqs;
}

/**
 * @brief Compute power spectrum into caller-provided output buffer.
 *
 * @param fft_result Complex FFT bins.
 * @param power Output buffer of power values `|X[k]|^2`.
 */
inline void power_spectrum(std::span<const std::complex<double>> fft_result,
                           std::span<double> power) {
    if (power.size() != fft_result.size()) {
        throw std::invalid_argument("Power spectrum: output buffer size must "
                                    "be equal to FFT result size");
    }

    for (size_t i = 0; i < fft_result.size(); ++i) {
        power[i] = std::norm(fft_result[i]); // |z|^2
    }
}

/**
 * @brief Compute power spectrum from FFT result
 *
 * @param fft_result Complex FFT output
 * @return Power spectrum (|FFT|^2)
 */
inline std::vector<double>
power_spectrum(std::span<const std::complex<double>> fft_result) {
    std::vector<double> power(fft_result.size());
    power_spectrum(fft_result, power);
    return power;
}

/**
 * @brief Compute magnitude spectrum into caller-provided output buffer.
 *
 * @param fft_result Complex FFT bins.
 * @param magnitude Output buffer of magnitudes `|X[k]|`.
 */
inline void magnitude_spectrum(std::span<const std::complex<double>> fft_result,
                               std::span<double> magnitude) {
    if (magnitude.size() != fft_result.size()) {
        throw std::invalid_argument("Magnitude spectrum: output buffer size "
                                    "must be equal to FFT result size");
    }

    for (size_t i = 0; i < fft_result.size(); ++i) {
        magnitude[i] = std::abs(fft_result[i]);
    }
}

/**
 * @brief Compute magnitude spectrum from FFT result
 *
 * @param fft_result Complex FFT output
 * @return Magnitude spectrum (|FFT|)
 */
inline std::vector<double>
magnitude_spectrum(std::span<const std::complex<double>> fft_result) {
    std::vector<double> magnitude(fft_result.size());
    magnitude_spectrum(fft_result, magnitude);
    return magnitude;
}

/**
 * @brief Compute phase spectrum into caller-provided output buffer.
 *
 * @param fft_result Complex FFT bins.
 * @param phase Output buffer of phase values `arg(X[k])` in radians.
 */
inline void phase_spectrum(std::span<const std::complex<double>> fft_result,
                           std::span<double> phase) {
    if (phase.size() != fft_result.size()) {
        throw std::invalid_argument("Phase spectrum: output buffer size must "
                                    "be equal to FFT result size");
    }

    for (size_t i = 0; i < fft_result.size(); ++i) {
        phase[i] = std::arg(fft_result[i]);
    }
}

/**
 * @brief Compute phase spectrum from FFT result
 *
 * @param fft_result Complex FFT output
 * @return Phase spectrum (arg(FFT))
 */
inline std::vector<double>
phase_spectrum(std::span<const std::complex<double>> fft_result) {
    std::vector<double> phase(fft_result.size());
    phase_spectrum(fft_result, phase);
    return phase;
}
} // namespace msl::signal

#endif // MSL_FFT_HPP
