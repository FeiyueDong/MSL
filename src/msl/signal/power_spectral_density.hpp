/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2026, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: power_spectral_density.hpp
** -----
** File Created: Friday, 9th January 2026 14:58:10
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Friday, 6th March 2026 09:45:11
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_POWER_SPECTRAL_DENSITY_HPP
#define MSL_POWER_SPECTRAL_DENSITY_HPP

#include <algorithm>
#include <complex>
#include <span>
#include <stdexcept>
#include <vector>

#include "fft.hpp"
#include "window.hpp"

namespace msl::signal {
// ========================================================================
// 1. Cross Power Spectral Density (CPSD) computation
// ========================================================================
/**
 * @brief Compute cross power spectral density using Welch's method.
 *
 * @param x First input signal.
 * @param y Second input signal.
 * Returns a full two-sided spectrum using `X * conj(Y)`. The density scale is
 * `1 / (sampling_rate * sum(window^2))` before segment averaging.
 *
 * @param output Output CPSD spectrum with size `nfft`.
 * @param window Segment window coefficients, size must equal `segment_length`.
 * @param noverlap Overlap samples between adjacent segments.
 * @param segment_length Number of input samples in each segment.
 * @param nfft FFT length; `0` means `segment_length`.
 * @param sampling_rate Sampling rate in Hz.
 */
inline void cpsd_welch(std::span<const double> x,
                       std::span<const double> y,
                       std::span<std::complex<double>> output,
                       const std::vector<double> &window = hann_window(1024),
                       size_t noverlap = 512,
                       size_t segment_length = 1024,
                       size_t nfft = 0,
                       double sampling_rate = 1.0) {
    if (segment_length == 0) {
        throw std::invalid_argument(
            "CPSD Welch: segment_length must be greater than 0");
    }
    nfft = nfft == 0 ? segment_length : nfft;
    if (nfft < segment_length) {
        throw std::invalid_argument(
            "CPSD Welch: nfft must be at least segment_length");
    }
    if (sampling_rate <= 0.0) {
        throw std::invalid_argument(
            "CPSD Welch: sampling_rate must be greater than 0");
    }
    if (noverlap >= segment_length) {
        throw std::invalid_argument(
            "CPSD Welch: noverlap must be less than segment_length");
    }
    if (x.size() != y.size()) {
        throw std::invalid_argument("Input signals must have the same length.");
    }
    if (window.size() != segment_length) {
        throw std::invalid_argument(
            "CPSD Welch: window size must be equal to segment_length");
    }
    if (output.size() != nfft) {
        throw std::invalid_argument(
            "CPSD Welch: output buffer size must be equal to nfft");
    }
    if (x.size() < segment_length || y.size() < segment_length) {
        throw std::invalid_argument(
            "CPSD Welch: input signals must be at least as long "
            "as segment_length.");
    }
    const size_t step = segment_length - noverlap;
    const size_t num_segments = 1 + (x.size() - segment_length) / step;

    std::vector<std::complex<double>> psd_accum(nfft,
                                                std::complex<double>(0.0, 0.0));
    double window_norm = 0.0;
    for (double w : window) {
        window_norm += w * w;
    }
    if (window_norm == 0.0) {
        throw std::invalid_argument(
            "CPSD Welch: window energy must be greater than 0");
    }

    for (size_t seg = 0; seg < num_segments; ++seg) {
        size_t start = seg * step;

        std::vector<double> x_segment(x.begin() + start,
                                      x.begin() + start + segment_length);
        std::vector<double> y_segment(y.begin() + start,
                                      y.begin() + start + segment_length);

        // Apply window
        for (size_t i = 0; i < segment_length; ++i) {
            x_segment[i] *= window[i];
            y_segment[i] *= window[i];
        }

        // Compute FFTs
        auto Xf = signal::fft(x_segment, nfft);
        auto Yf = signal::fft(y_segment, nfft);

        // Accumulate cross power
        for (size_t k = 0; k < nfft; ++k) {
            psd_accum[k] += Xf[k] * std::conj(Yf[k]);
        }
    }

    // Average and normalize
    const double scale =
        static_cast<double>(num_segments) * window_norm * sampling_rate;
    for (size_t k = 0; k < nfft; ++k) {
        output[k] = psd_accum[k] / scale;
    }
}

/**
 * @brief Return cross power spectral density using Welch's method.
 *
 * @param x First input signal.
 * @param y Second input signal.
 * @param window Segment window coefficients, size must equal `segment_length`.
 * @param noverlap Overlap samples between adjacent segments.
 * @param segment_length Number of input samples in each segment.
 * @param nfft FFT length; `0` means `segment_length`.
 * @param sampling_rate Sampling rate in Hz.
 * @return Full two-sided CPSD spectrum with size `nfft`.
 */
inline std::vector<std::complex<double>>
cpsd_welch(std::span<const double> x,
           std::span<const double> y,
           const std::vector<double> &window = hann_window(1024),
           size_t noverlap = 512,
           size_t segment_length = 1024,
           size_t nfft = 0,
           double sampling_rate = 1.0) {
    nfft = nfft == 0 ? segment_length : nfft;
    std::vector<std::complex<double>> output(nfft);
    cpsd_welch(
        x, y, output, window, noverlap, segment_length, nfft, sampling_rate);
    return output;
}

/**
 * @brief Compute cross power spectral density from one FFT frame.
 *
 * @param x First input signal.
 * @param y Second input signal.
 * @param output Output CPSD spectrum with size `nfft`.
 * @param nfft FFT length.
 */
inline void cpsd(std::span<const double> x,
                 std::span<const double> y,
                 std::span<std::complex<double>> output,
                 size_t nfft) {
    if (nfft == 0) {
        throw std::invalid_argument("CPSD: nfft must be greater than 0");
    }
    if (x.size() != y.size()) {
        throw std::invalid_argument("Input signals must have the same length.");
    }
    if (output.size() != nfft) {
        throw std::invalid_argument(
            "CPSD: output buffer size must be equal to nfft");
    }
    if (x.size() < nfft || y.size() < nfft) {
        throw std::invalid_argument(
            "CPSD: input signals must be at least as long as nfft.");
    }

    // Compute FFTs
    auto Xf = signal::fft(x, nfft);
    auto Yf = signal::fft(y, nfft);

    // Compute Cross Power Spectral Density
    for (size_t k = 0; k < nfft; ++k) {
        output[k] = Xf[k] * std::conj(Yf[k]);
    }
}

/**
 * @brief Return cross power spectral density from one FFT frame.
 *
 * @param x First input signal.
 * @param y Second input signal.
 * @param nfft FFT length.
 * @return CPSD spectrum with size `nfft`.
 */
inline std::vector<std::complex<double>>
cpsd(std::span<const double> x, std::span<const double> y, size_t nfft) {
    std::vector<std::complex<double>> output(nfft);
    cpsd(x, y, output, nfft);
    return output;
}

// ========================================================================
// 2. Power Spectral Density (PSD) computation
// ========================================================================
/**
 * @brief Compute power spectral density using Welch's method.
 *
 * @param x Input signal.
 * @param output Output PSD spectrum with size `nfft`.
 * @param window Segment window coefficients, size must equal `segment_length`.
 * @param noverlap Overlap samples between adjacent segments.
 * @param segment_length Number of input samples in each segment.
 * @param nfft FFT length; `0` means `segment_length`.
 * @param sampling_rate Sampling rate in Hz.
 */
inline void psd_welch(std::span<const double> x,
                      std::span<double> output,
                      const std::vector<double> &window = hann_window(1024),
                      size_t noverlap = 512,
                      size_t segment_length = 1024,
                      size_t nfft = 0,
                      double sampling_rate = 1.0) {
    nfft = nfft == 0 ? segment_length : nfft;
    if (output.size() != nfft) {
        throw std::invalid_argument(
            "PSD Welch: output buffer size must be equal to nfft");
    }

    std::vector<std::complex<double>> cpsd_result(nfft);
    cpsd_welch(x,
               x,
               cpsd_result,
               window,
               noverlap,
               segment_length,
               nfft,
               sampling_rate);

    for (size_t k = 0; k < nfft; ++k) {
        output[k] = std::real(cpsd_result[k]);
    }
}

/**
 * @brief Return power spectral density using Welch's method.
 *
 * @param x Input signal.
 * @param window Segment window coefficients, size must equal `segment_length`.
 * @param noverlap Overlap samples between adjacent segments.
 * @param segment_length Number of input samples in each segment.
 * @param nfft FFT length; `0` means `segment_length`.
 * @param sampling_rate Sampling rate in Hz.
 * @return Full two-sided PSD spectrum with size `nfft`.
 */
inline std::vector<double>
psd_welch(std::span<const double> x,
          const std::vector<double> &window = hann_window(1024),
          size_t noverlap = 512,
          size_t segment_length = 1024,
          size_t nfft = 0,
          double sampling_rate = 1.0) {
    nfft = nfft == 0 ? segment_length : nfft;
    std::vector<double> output(nfft);
    psd_welch(x, output, window, noverlap, segment_length, nfft, sampling_rate);
    return output;
}

/**
 * @brief Compute power spectral density from one FFT frame.
 *
 * @param x Input signal.
 * @param output Output PSD spectrum with size `nfft`.
 * @param nfft FFT length.
 */
inline void
psd(std::span<const double> x, std::span<double> output, size_t nfft) {
    if (nfft == 0) {
        throw std::invalid_argument("PSD: nfft must be greater than 0");
    }
    if (output.size() != nfft) {
        throw std::invalid_argument(
            "PSD: output buffer size must be equal to nfft");
    }

    std::vector<std::complex<double>> cpsd_result(nfft);
    cpsd(x, x, cpsd_result, nfft);

    for (size_t k = 0; k < nfft; ++k) {
        output[k] = std::real(cpsd_result[k]);
    }
}

/**
 * @brief Return power spectral density from one FFT frame.
 *
 * @param x Input signal.
 * @param nfft FFT length.
 * @return PSD spectrum with size `nfft`.
 */
inline std::vector<double> psd(std::span<const double> x, size_t nfft) {
    std::vector<double> output(nfft);
    psd(x, output, nfft);
    return output;
}

} // namespace msl::signal

#endif // MSL_POWER_SPECTRAL_DENSITY_HPP
