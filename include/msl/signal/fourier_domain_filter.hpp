/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2026, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: fourier_domain_filter.hpp
** -----
** File Created: Friday, 9th January 2026 14:58:10
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Friday, 6th March 2026 11:58:29
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_FOURIER_DOMAIN_FILTER_HPP
#define MSL_FOURIER_DOMAIN_FILTER_HPP

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <functional>
#include <numbers>
#include <span>
#include <stdexcept>
#include <vector>

#include "fft.hpp"
#include "filter_design.hpp"
#include "matrix.hpp"
#include "matrix/real_matrix_base.hpp"

namespace msl::signal {
/**
 * @brief Fourier domain filter designer
 *
 * Designs bandpass or bandstop filters in the Fourier domain.
 * The filter is designed by creating an ideal frequency response and
 * applying a window function to reduce ringing artifacts.
 */
class FourierDomainFilter {
public:
    /**
     * @brief Transition shaping used around cutoff edges.
     *
     * - rectangular: ideal sinc filter, sharp cutoff but high side lobes
     * - hamming: smoother transition, lower side lobes
     * - hanning: similar to Hamming, slightly wider main lobe
     * - blackman: even smoother, very low side lobes but wider transition
     */
    enum class WindowType { rectangular, hamming, hanning, blackman };

private:
    // Low cutoff frequency (0-1 where 1 = Nyquist)
    double fc_low_{0.0};
    // High cutoff frequency (0-1 where 1 = Nyquist)
    double fc_high_{1.0};

    // Filter type
    FilterType type_{FilterType::bandpass};

    // Window function type
    WindowType window_type_{WindowType::rectangular};
    // Transition band width (normalized)
    double transition_band_{0.0};

    // Window function generator based on current settings
    std::function<double(double)> window_function_;

public:
    FourierDomainFilter() = default;

    /**
     * @brief Create a lowpass Fourier-domain filter.
     *
     * @param fc_high Cutoff frequency (normalized: 0 < fc < 1)
     * @param window_type Window function type
     * @param transition_band Transition band width normalized of Nyquist
     * @return FourierDomainFilter with computed window function
     */
    [[nodiscard]] static FourierDomainFilter
    lowpass(double fc_high,
            WindowType window_type = WindowType::rectangular,
            double transition_band = 0.0) {
        FourierDomainFilter filter;
        filter.fc_low_ = 0.0;
        filter.fc_high_ = fc_high;
        filter.type_ = FilterType::lowpass;
        filter.window_type_ = window_type;
        filter.transition_band_ = transition_band;
        filter.design();
        return filter;
    }

    /**
     * @brief Create a highpass Fourier-domain filter.
     *
     * @param fc_low Cutoff frequency (normalized: 0 < fc < 1)
     * @param window_type Window function type
     * @param transition_band Transition band width normalized of Nyquist
     * @return FourierDomainFilter with computed window function
     */
    [[nodiscard]] static FourierDomainFilter
    highpass(double fc_low,
             WindowType window_type = WindowType::rectangular,
             double transition_band = 0.0) {
        FourierDomainFilter filter;
        filter.fc_low_ = fc_low;
        filter.fc_high_ = 1.0;
        filter.type_ = FilterType::highpass;
        filter.window_type_ = window_type;
        filter.transition_band_ = transition_band;
        filter.design();
        return filter;
    }

    /**
     * @brief Create a bandpass Fourier-domain filter.
     *
     * @param fc_low Low cutoff frequency (normalized: 0 < fc < 1)
     * @param fc_high High cutoff frequency (normalized: 0 < fc < 1)
     * @param window_type Window function type
     * @param transition_band Transition band width normalized of Nyquist
     * @return FourierDomainFilter with computed window function
     */
    [[nodiscard]] static FourierDomainFilter
    bandpass(double fc_low,
             double fc_high,
             WindowType window_type = WindowType::rectangular,
             double transition_band = 0.0) {
        FourierDomainFilter filter;
        filter.fc_low_ = fc_low;
        filter.fc_high_ = fc_high;
        filter.type_ = FilterType::bandpass;
        filter.window_type_ = window_type;
        filter.transition_band_ = transition_band;
        filter.design();
        return filter;
    }

    /**
     * @brief Create a bandstop Fourier-domain filter.
     *
     * @param fc_low Low cutoff frequency (normalized: 0 < fc < 1)
     * @param fc_high High cutoff frequency (normalized: 0 < fc < 1)
     * @param window_type Window function type
     * @param transition_band Transition band width normalized of Nyquist
     * @return FourierDomainFilter with computed window function
     */
    [[nodiscard]] static FourierDomainFilter
    bandstop(double fc_low,
             double fc_high,
             WindowType window_type = WindowType::rectangular,
             double transition_band = 0.0) {
        FourierDomainFilter filter;
        filter.fc_low_ = fc_low;
        filter.fc_high_ = fc_high;
        filter.type_ = FilterType::bandstop;
        filter.window_type_ = window_type;
        filter.transition_band_ = transition_band;
        filter.design();
        return filter;
    }

    /**
     * @brief Set window type
     * @param window_type Window function type
     */
    void set_window_type(WindowType window_type) {
        window_type_ = window_type;
        design();
    }

    /**
     * @brief Set transition band width
     * @param width Transition band width (normalized, e.g., 0.01 = 1% of
     * Nyquist)
     */
    void set_transition_band(double width) {

        transition_band_ = width;
        design();
    }

    /**
     * @brief Redesign the filter with new parameters, same as constructor
     *
     * @param fc_low Low cutoff frequency (normalized: 0 < fc < 1, where 1 =
     * Nyquist frequency)
     * @param fc_high High cutoff frequency (normalized: 0 < fc < 1, where 1 =
     * Nyquist frequency)
     * @param type Filter type (lowpass/highpass/bandpass/bandstop)
     * @param window_type Window function type
     * @param transition_band Transition band width normalized of Nyquist
     */
    void redesign(double fc_low,
                  double fc_high,
                  FilterType type = FilterType::bandpass,
                  WindowType window_type = WindowType::rectangular,
                  double transition_band = 0.0) {
        fc_low_ = fc_low;
        fc_high_ = fc_high;
        type_ = type;
        window_type_ = window_type;
        transition_band_ = transition_band;
        design();
    }

    /**
     * @brief Apply the filter to the input signal
     *
     * @param signal Input signal (time domain)
     * @param output Output buffer for filtered signal (time domain), must be
     * same size as input
     */
    void apply(std::span<const double> signal, std::span<double> output) const {
        if (signal.size() != output.size()) {
            throw std::invalid_argument("FourierDomainFilter: output buffer "
                                        "size must match input signal size");
        }
        if (signal.empty()) {
            return;
        }

        std::size_t fft_size = signal.size();

        std::vector<double> fft_input(fft_size, 0.0);
        std::copy(signal.begin(), signal.end(), fft_input.begin());

        auto fft_data = signal::fft(fft_input, fft_size);
        std::vector<std::complex<double>> filtered_fft(fft_data.size());

        for (std::size_t i = 0; i < fft_data.size(); ++i) {
            double freq = static_cast<double>(i) / fft_size;
            if (freq > 0.5)
                freq -= 1.0; // Map to [-0.5, 0.5]
            freq *= 2.0;     // Normalize to [-1, 1]
            double gain = window_function_(freq);
            filtered_fft[i] = fft_data[i] * gain;
        }

        signal::ifft_real(filtered_fft, output, fft_size);
    }

    /**
     * @brief Apply the filter to the input signal
     *
     * @param signal Input signal (time domain)
     * @return Filtered signal (time domain)
     */
    std::vector<double> apply(std::span<const double> signal) const {
        std::vector<double> output(signal.size());
        apply(signal, output);
        return output;
    }

    /**
     * @brief Apply the filter to the input signal matrix
     *
     * @param signal_matrix Input signal matrix (time domain)
     * @return Filtered signal matrix (time domain)
     */
    matrix::matrixd apply(const matrix::real_matrix_base &signal_matrix) const {
        matrix::matrixd output(signal_matrix.rows(), signal_matrix.cols());
        apply(signal_matrix, output);
        return output;
    }

    /**
     * @brief Apply the filter to all columns of an input matrix.
     *
     * @param signal_matrix Input matrix where each column is one signal.
     * @param output Output matrix. Must have the same shape as input.
     */
    void apply(const matrix::real_matrix_base &signal_matrix,
               matrix::matrixd &output) const {
        if (output.rows() != signal_matrix.rows()
            || output.cols() != signal_matrix.cols()) {
            throw std::invalid_argument("FourierDomainFilter: output matrix "
                                        "shape must match input matrix");
        }

        for (std::size_t j = 0; j < signal_matrix.cols(); ++j) {
            std::vector<double> column(signal_matrix.column(j).begin(),
                                       signal_matrix.column(j).end());
            auto filtered_column = apply(column);
            std::copy(filtered_column.begin(),
                      filtered_column.end(),
                      output.column(j).begin());
        }
    }


private:
    /**
     * @brief Validate filter parameters
     *
     * @throws std::invalid_argument if parameters are invalid
     */
    void validate_parameters() const {
        switch (type_) {
            case FilterType::lowpass:
                if (fc_high_ <= 0.0 || fc_high_ >= 1.0) {
                    throw std::invalid_argument(
                        "FourierDomainFilter: lowpass cutoff must satisfy 0 < "
                        "fc_high < 1");
                }
                break;

            case FilterType::highpass:
                if (fc_low_ <= 0.0 || fc_low_ >= 1.0) {
                    throw std::invalid_argument(
                        "FourierDomainFilter: highpass cutoff must satisfy 0 < "
                        "fc_low < 1");
                }
                break;

            case FilterType::bandpass:
            case FilterType::bandstop:
                if (fc_low_ <= 0.0 || fc_low_ >= 1.0 || fc_high_ <= 0.0
                    || fc_high_ >= 1.0) {
                    throw std::invalid_argument(
                        "FourierDomainFilter: band cutoffs must satisfy 0 < "
                        "fc_low < fc_high < 1");
                }
                if (fc_low_ >= fc_high_) {
                    throw std::invalid_argument(
                        "FourierDomainFilter: low frequency must be < high "
                        "frequency");
                }
                break;

            default:
                throw std::invalid_argument(
                    "FourierDomainFilter: unsupported filter type");
        }

        if (transition_band_ < 0.0 || transition_band_ >= 0.5) {
            throw std::invalid_argument("FourierDomainFilter: transition band "
                                        "must be in range [0, 0.5)");
        }
    }

    /**
     * @brief Design the filter
     */
    void design() {
        validate_parameters();

        window_function_ = [=, this](double f) -> double {
            f = std::abs(f); // 对称处理
            f = std::clamp(f, 0.0, 1.0);

            const double fl = std::clamp(fc_low_, 0.0, 1.0);
            const double fh = std::clamp(fc_high_, 0.0, 1.0);
            const double tw = std::clamp(transition_band_, 0.0, 0.5);

            if (tw == 0.0) {
                switch (type_) {
                    case FilterType::lowpass:
                        return (f <= fh) ? 1.0 : 0.0;
                    case FilterType::highpass:
                        return (f >= fl) ? 1.0 : 0.0;
                    case FilterType::bandpass:
                        return (f >= fl && f <= fh) ? 1.0 : 0.0;
                    case FilterType::bandstop:
                        return (f >= fl && f <= fh) ? 0.0 : 1.0;
                    default:
                        return 0.0;
                }
            }

            auto smooth = [&](double x) {
                double t = std::clamp((x + tw) / (2.0 * tw), 0.0, 1.0);
                switch (window_type_) {
                    case WindowType::rectangular:
                        return (t >= 0.5) ? 1.0 : 0.0;
                    case WindowType::hanning:
                        return 0.5 * (1 - std::cos(std::numbers::pi * t));
                    case WindowType::hamming:
                        return 0.54 - 0.46 * std::cos(std::numbers::pi * t);
                    case WindowType::blackman:
                        return 0.42 - 0.5 * std::cos(std::numbers::pi * t)
                               + 0.08 * std::cos(2 * std::numbers::pi * t);
                }
                return 0.0;
            };

            double gain = 0.0;

            switch (type_) {
                case FilterType::lowpass: {
                    const double pass_end = fh - tw;
                    const double stop_start = fh + tw;

                    if (f <= pass_end)
                        gain = 1.0;
                    else if (f >= stop_start)
                        gain = 0.0;
                    else
                        gain = smooth(stop_start - f);
                    break;
                }

                case FilterType::highpass: {
                    const double stop_end = fl - tw;
                    const double pass_start = fl + tw;

                    if (f <= stop_end)
                        gain = 0.0;
                    else if (f >= pass_start)
                        gain = 1.0;
                    else
                        gain = smooth(f - stop_end);
                    break;
                }

                case FilterType::bandpass: {
                    const double f_low_start = fl - tw;
                    const double f_low_end = fl + tw;
                    const double f_high_start = fh - tw;
                    const double f_high_end = fh + tw;

                    if (f < f_low_start || f > f_high_end)
                        gain = 0.0;
                    else if (f >= f_low_end && f <= f_high_start)
                        gain = 1.0;
                    else if (f >= f_low_start && f < f_low_end)
                        gain = smooth(f - f_low_start);
                    else if (f > f_high_start && f <= f_high_end)
                        gain = smooth(f_high_end - f);
                    break;
                }

                case FilterType::bandstop: {
                    const double f_low_start = fl - tw;
                    const double f_low_end = fl + tw;
                    const double f_high_start = fh - tw;
                    const double f_high_end = fh + tw;

                    if (f < f_low_start || f > f_high_end)
                        gain = 1.0;
                    else if (f >= f_low_end && f <= f_high_start)
                        gain = 0.0;
                    else if (f >= f_low_start && f < f_low_end)
                        gain = 1.0 - smooth(f - f_low_start);
                    else if (f > f_high_start && f <= f_high_end)
                        gain = 1.0 - smooth(f_high_end - f);
                    break;
                }
            }

            return std::clamp(gain, 0.0, 1.0);
        };
    }
};

// ============================================================================
// Convenience functions
// ============================================================================

/**
 * @brief Design and apply bandpass filter in Fourier domain
 *
 * @param signal Input signal
 * @param result Output buffer for filtered signal, must be same size as input
 * @param fc_low Low cutoff frequency (normalized, 0 < fc < 1)
 * @param fc_high High cutoff frequency (normalized, 0 < fc < 1)
 * @param window_type Window function type
 * @param transition_band Transition band width normalized of Nyquist
 */
inline void fourier_bandpass(std::span<const double> signal,
                             std::span<double> result,
                             double fc_low,
                             double fc_high,
                             FourierDomainFilter::WindowType window_type =
                                 FourierDomainFilter::WindowType::rectangular,
                             double transition_band = 0.0) {
    auto filter =
        FourierDomainFilter::bandpass(fc_low, fc_high, window_type, transition_band);
    filter.apply(signal, result);
}

/**
 * @brief Design and apply bandpass filter in Fourier domain
 *
 * @param signal Input signal
 * @param fc_low Low cutoff frequency (normalized, 0 < fc < 1)
 * @param fc_high High cutoff frequency (normalized, 0 < fc < 1)
 * @param window_type Window function type
 * @param transition_band Transition band width normalized of Nyquist
 * @return Filtered signal
 */
inline std::vector<double>
fourier_bandpass(std::span<const double> signal,
                 double fc_low,
                 double fc_high,
                 FourierDomainFilter::WindowType window_type =
                     FourierDomainFilter::WindowType::rectangular,
                 double transition_band = 0.0) {
    auto filter =
        FourierDomainFilter::bandpass(fc_low, fc_high, window_type, transition_band);
    return filter.apply(signal);
}

/**
 * @brief Design and apply bandstop filter in Fourier domain
 *
 * @param signal Input signal
 * @param result Output buffer for filtered signal, must be same size as input
 * @param fc_low Low cutoff frequency (normalized, 0 < fc < 1)
 * @param fc_high High cutoff frequency (normalized, 0 < fc < 1)
 * @param window_type Window function type
 * @param transition_band Transition band width normalized of Nyquist
 */
inline void fourier_bandstop(std::span<const double> signal,
                             std::span<double> result,
                             double fc_low,
                             double fc_high,
                             FourierDomainFilter::WindowType window_type =
                                 FourierDomainFilter::WindowType::rectangular,
                             double transition_band = 0.0) {
    auto filter =
        FourierDomainFilter::bandstop(fc_low, fc_high, window_type, transition_band);
    filter.apply(signal, result);
}

/**
 * @brief Design and apply bandstop filter in Fourier domain
 *
 * @param signal Input signal
 * @param fc_low Low cutoff frequency (normalized, 0 < fc < 1)
 * @param fc_high High cutoff frequency (normalized, 0 < fc < 1)
 * @param window_type Window function type
 * @param transition_band Transition band width normalized of Nyquist
 * @return Filtered signal
 */
inline std::vector<double>
fourier_bandstop(std::span<const double> signal,
                 double fc_low,
                 double fc_high,
                 FourierDomainFilter::WindowType window_type =
                     FourierDomainFilter::WindowType::rectangular,
                 double transition_band = 0.0) {
    auto filter =
        FourierDomainFilter::bandstop(fc_low, fc_high, window_type, transition_band);
    return filter.apply(signal);
}

/**
 * @brief Design and apply lowpass filter in Fourier domain
 *
 * @param signal Input signal
 * @param result Output buffer for filtered signal, must be same size as input
 * @param fc_high High cutoff frequency (normalized, 0 < fc < 1)
 * @param window_type Window function type
 * @param transition_band Transition band width normalized of Nyquist
 */
inline void fourier_lowpass(std::span<const double> signal,
                            std::span<double> result,
                            double fc_high,
                            FourierDomainFilter::WindowType window_type =
                                FourierDomainFilter::WindowType::rectangular,
                            double transition_band = 0.0) {
    auto filter =
        FourierDomainFilter::lowpass(fc_high, window_type, transition_band);
    filter.apply(signal, result);
}

/**
 * @brief Design and apply lowpass filter in Fourier domain
 *
 * @param signal Input signal
 * @param fc_high High cutoff frequency (normalized, 0 < fc < 1)
 * @param window_type Window function type
 * @param transition_band Transition band width normalized of Nyquist
 * @return Filtered signal
 */
inline std::vector<double>
fourier_lowpass(std::span<const double> signal,
                double fc_high,
                FourierDomainFilter::WindowType window_type =
                    FourierDomainFilter::WindowType::rectangular,
                double transition_band = 0.0) {
    auto filter =
        FourierDomainFilter::lowpass(fc_high, window_type, transition_band);
    return filter.apply(signal);
}

/**
 * @brief Design and apply highpass filter in Fourier domain
 *
 * @param signal Input signal
 * @param result Output buffer for filtered signal, must be same size as input
 * @param fc_low Low cutoff frequency (normalized, 0 < fc < 1)
 * @param window_type Window function type
 * @param transition_band Transition band width normalized of Nyquist
 */
inline void fourier_highpass(std::span<const double> signal,
                             std::span<double> result,
                             double fc_low,
                             FourierDomainFilter::WindowType window_type =
                                 FourierDomainFilter::WindowType::rectangular,
                             double transition_band = 0.0) {
    auto filter =
        FourierDomainFilter::highpass(fc_low, window_type, transition_band);
    filter.apply(signal, result);
}

/**
 * @brief Design and apply highpass filter in Fourier domain
 *
 * @param signal Input signal
 * @param fc_low Low cutoff frequency (normalized, 0 < fc < 1)
 * @param window_type Window function type
 * @param transition_band Transition band width normalized of Nyquist
 * @return Filtered signal
 */
inline std::vector<double>
fourier_highpass(std::span<const double> signal,
                 double fc_low,
                 FourierDomainFilter::WindowType window_type =
                     FourierDomainFilter::WindowType::rectangular,
                 double transition_band = 0.0) {
    auto filter =
        FourierDomainFilter::highpass(fc_low, window_type, transition_band);
    return filter.apply(signal);
}

} // namespace msl::signal

#endif // MSL_FOURIER_DOMAIN_FILTER_HPP