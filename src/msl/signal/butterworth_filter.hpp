/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2026, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: butterworth_filter.hpp
** -----
** File Created: Friday, 9th January 2026 14:58:10
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Friday, 6th March 2026 10:15:55
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_BUTTERWORTH_FILTER_HPP
#define MSL_BUTTERWORTH_FILTER_HPP

#include "filter.hpp"
#include "filter_design.hpp"
#include "filtfilt.hpp"

#include <cmath>
#include <complex>
#include <numbers>
#include <span>
#include <stdexcept>
#include <vector>

namespace msl::signal {
/**
 * @brief Butterworth filter designer
 *
 * Butterworth filters have maximally flat magnitude response in the passband.
 * No ripple in passband or stopband, monotonic response.
 *
 * The digital filter is designed using bilinear transformation from analog
 * prototype.
 */
class ButterworthFilter {
private:
    int order_ = 0;
    double fc_low_ =
        0.0; // Low cutoff frequency (normalized, 0-1 where 1 = Nyquist)
    double fc_high_ = 0.0; // High cutoff frequency
    FilterType type_ = FilterType::lowpass;

    FilterCoefficients coeffs_;

public:
    /**
     * @brief Default constructor
     */
    ButterworthFilter() = default;

    /**
     * @brief Create a lowpass Butterworth filter.
     *
     * @param order Filter order (must be positive)
     * @param fc Cutoff frequency (normalized: 0 < fc < 1)
     * @return ButterworthFilter with computed coefficients
     */
    [[nodiscard]] static ButterworthFilter lowpass(int order, double fc) {
        ButterworthFilter filter;
        filter.order_ = order;
        filter.type_ = FilterType::lowpass;
        filter.fc_low_ = fc;
        filter.fc_high_ = 0.0;
        filter.design();
        return filter;
    }

    /**
     * @brief Create a highpass Butterworth filter.
     *
     * @param order Filter order (must be positive)
     * @param fc Cutoff frequency (normalized: 0 < fc < 1)
     * @return ButterworthFilter with computed coefficients
     */
    [[nodiscard]] static ButterworthFilter highpass(int order, double fc) {
        ButterworthFilter filter;
        filter.order_ = order;
        filter.type_ = FilterType::highpass;
        filter.fc_low_ = 0.0;
        filter.fc_high_ = fc;
        filter.design();
        return filter;
    }

    /**
     * @brief Create a bandpass Butterworth filter.
     *
     * @param order Filter order (must be positive)
     * @param fc_low Low cutoff frequency (normalized)
     * @param fc_high High cutoff frequency (normalized)
     * @return ButterworthFilter with computed coefficients
     */
    [[nodiscard]] static ButterworthFilter
    bandpass(int order, double fc_low, double fc_high) {
        ButterworthFilter filter;
        filter.order_ = order;
        filter.fc_low_ = fc_low;
        filter.fc_high_ = fc_high;
        filter.type_ = FilterType::bandpass;
        filter.design();
        return filter;
    }

    /**
     * @brief Create a bandstop Butterworth filter.
     *
     * @param order Filter order (must be positive)
     * @param fc_low Low cutoff frequency (normalized)
     * @param fc_high High cutoff frequency (normalized)
     * @return ButterworthFilter with computed coefficients
     */
    [[nodiscard]] static ButterworthFilter
    bandstop(int order, double fc_low, double fc_high) {
        ButterworthFilter filter;
        filter.order_ = order;
        filter.fc_low_ = fc_low;
        filter.fc_high_ = fc_high;
        filter.type_ = FilterType::bandstop;
        filter.design();
        return filter;
    }

    /**
     * @brief Get filter coefficients
     *
     * @return FilterCoefficients struct containing numerator (b) and
     * denominator (a) coefficients
     */
    [[nodiscard]] const FilterCoefficients &coefficients() const {
        return coeffs_;
    }

    /**
     * @brief Get numerator coefficients
     *
     * @return Vector of numerator coefficients (b)
     */
    [[nodiscard]] const std::vector<double> &b() const { return coeffs_.b; }

    /**
     * @brief Get denominator coefficients
     *
     * @return Vector of denominator coefficients (a)
     */
    [[nodiscard]] const std::vector<double> &a() const { return coeffs_.a; }

    /** @brief Get filter order. */
    [[nodiscard]] int order() const { return order_; }
    /** @brief Get filter type. */
    [[nodiscard]] FilterType type() const { return type_; }
    /** @brief Get low cutoff frequency (0.0 for highpass). */
    [[nodiscard]] double fc_low() const { return fc_low_; }
    /** @brief Get high cutoff frequency (0.0 for lowpass). */
    [[nodiscard]] double fc_high() const { return fc_high_; }

    /**
     * @brief Redesign filter with new parameters, same with constructor
     *
     * @param order Filter order (must be positive)
     * @param fc Cutoff frequency (normalized: 0 < fc < 1, where 1 = Nyquist
     * frequency)
     * @param type Filter type (lowpass or highpass)
     */
    void redesign(int order, double fc, FilterType type = FilterType::lowpass) {
        if (type != FilterType::lowpass && type != FilterType::highpass) {
            throw std::invalid_argument(
                "Butterworth: single frequency redesign only for "
                "lowpass/highpass");
        }

        order_ = order;
        type_ = type;

        if (type == FilterType::lowpass) {
            fc_low_ = fc;
            fc_high_ = 0.0;
        } else // highpass
        {
            fc_low_ = 0.0;
            fc_high_ = fc;
        }

        design();
    }

    /**
     * @brief Redesign filter with new parameters, same with constructor
     *
     * @param order Filter order (must be positive)
     * @param fc_low Low cutoff frequency (normalized)
     * @param fc_high High cutoff frequency (normalized)
     * @param type Filter type (bandpass or bandstop)
     */
    void redesign(int order,
                  double fc_low,
                  double fc_high,
                  FilterType type = FilterType::bandpass) {
        if (type != FilterType::bandpass && type != FilterType::bandstop) {
            throw std::invalid_argument(
                "Butterworth: two-frequency redesign only for "
                "bandpass/bandstop");
        }

        order_ = order;
        fc_low_ = fc_low;
        fc_high_ = fc_high;
        type_ = type;

        design();
    }

private:
    // Validate parameters according to filter type.
    void validate_parameters() const {
        if (order_ <= 0) {
            throw std::invalid_argument(
                "Butterworth: filter order must be positive");
        }

        switch (type_) {
            case FilterType::lowpass:
                if (fc_low_ <= 0.0 || fc_low_ >= 1.0) {
                    throw std::invalid_argument(
                        "Butterworth: lowpass cutoff must satisfy 0 < fc < 1");
                }
                break;

            case FilterType::highpass:
                if (fc_high_ <= 0.0 || fc_high_ >= 1.0) {
                    throw std::invalid_argument(
                        "Butterworth: highpass cutoff must satisfy 0 < fc < 1");
                }
                break;

            case FilterType::bandpass:
            case FilterType::bandstop:
                if (fc_low_ <= 0.0 || fc_low_ >= 1.0) {
                    throw std::invalid_argument(
                        "Butterworth: low cutoff must satisfy 0 < fc_low < 1");
                }
                if (fc_high_ <= 0.0 || fc_high_ >= 1.0) {
                    throw std::invalid_argument("Butterworth: high cutoff must "
                                                "satisfy 0 < fc_high < 1");
                }
                if (fc_low_ >= fc_high_) {
                    throw std::invalid_argument(
                        "Butterworth: fc_low must be less than fc_high");
                }
                break;

            default:
                throw std::invalid_argument(
                    "Butterworth: unsupported filter type");
        }
    }

    /**
     * @brief Design the filter
     */
    void design() {
        validate_parameters();

        if (type_ == FilterType::bandpass || type_ == FilterType::bandstop) {
            design_bandpass_bandstop();
        } else {
            design_lowpass_highpass();
        }
    }

    /**
     * @brief Design lowpass or highpass filter
     */
    void design_lowpass_highpass() {
        const double fc = (type_ == FilterType::lowpass) ? fc_low_ : fc_high_;

        // Pre-warp cutoff frequency (T = 1 assumed)
        const double wc =
            2.0 * std::tan(std::numbers::pi * fc / 2.0); // <-- 修正

        // Analog prototype poles (normalized)
        auto analog_poles = get_analog_poles();

        // Map analog poles to digital poles using bilinear transform:
        // s = p * wc, z = (2 + s) / (2 - s)
        std::vector<std::complex<double>> digital_poles;
        digital_poles.reserve(order_);

        for (const auto &p : analog_poles) {
            std::complex<double> s = p * wc; // analog pole scaled
            std::complex<double> pz =
                (2.0 + s) / (2.0 - s); // <-- 修正: T=1 mapping
            digital_poles.push_back(pz);
        }

        // Denominator polynomial in z^{-1}: prod (1 - pz * z^{-1})
        coeffs_.a =
            poles_to_polynomial(digital_poles); // now returns a in z^{-k} order

        // Numerator: for lowpass the digital numerator corresponds to (1 +
        // z^{-1})^N (binomial)
        if (type_ == FilterType::lowpass) {
            coeffs_.b = compute_lowpass_numerator(); // returns coefficients in
                                                     // z^{-k} order (C(N,k))
        } else                                       // highpass
        {
            coeffs_.b = compute_highpass_numerator(); // also z^{-k} order
        }

        // Normalize: choose frequency point according to type
        if (type_ == FilterType::lowpass)
            normalize_gain(0.0); // unity at DC
        else
            normalize_gain(1.0); // unity at Nyquist (f=1 corresponds to ω=π)
    }

    /**
     * @brief Design bandpass or bandstop filter
     */
    void design_bandpass_bandstop() {
        // Pre-warp digital cutoff freqs (T = 1)
        const double wc1 = 2.0 * std::tan(std::numbers::pi * fc_low_ / 2.0);
        const double wc2 = 2.0 * std::tan(std::numbers::pi * fc_high_ / 2.0);

        const double w0 = std::sqrt(wc1 * wc2); // analog center freq
        const double bw = wc2 - wc1;            // analog bandwidth

        // Get analog lowpass prototype poles (order_ of them)
        auto lp_poles = get_analog_poles(); // returns poles of normalized LP

        // Transform each lowpass pole p into the two analog poles of the
        // requested band filter. Bandpass and bandstop use different
        // LP -> band transformations.
        std::vector<std::complex<double>> analog_poles;
        if (type_ == FilterType::bandpass) {
            analog_poles = transform_lowpass_to_bandpass(lp_poles, w0, bw);
        } else {
            analog_poles = transform_lowpass_to_bandstop(lp_poles, w0, bw);
        }

        // Bilinear transform each analog s pole to z-domain: z = (2 + s) / (2 -
        // s)
        std::vector<std::complex<double>> digital_poles;
        digital_poles.reserve(analog_poles.size());
        for (const auto &s : analog_poles) {
            std::complex<double> z = (2.0 + s) / (2.0 - s);
            digital_poles.push_back(z);
        }

        // Denominator polynomial from digital poles (z^{-1} representation)
        coeffs_.a = poles_to_polynomial(
            digital_poles); // returns [a0,a1,...,aM] for z^{-1} powers

        // Numerator: depends on filter type
        const int M = 2 * order_;
        coeffs_.b.assign(M + 1, 0.0);

        if (type_ == FilterType::bandpass) {
            // Bandpass: zeros at s=0 and s=∞ map to z=1 and z=-1
            // Digital numerator has factor (1 - z^{-2})^N
            // (1 - t^2)^N = sum_{m=0..N} C(N,m) * (-1)^m * t^{2m}
            auto binom = [](int n, int k) -> double {
                double r = 1.0;
                for (int i = 1; i <= k; ++i)
                    r *= double(n - (k - i)) / double(i);
                return r;
            };

            for (int m = 0; m <= order_; ++m) {
                double c = binom(order_, m)
                           * ((m % 2 == 0) ? 1.0 : -1.0); // (-1)^m * C(N,m)
                coeffs_.b[2 * m] = c;
            }

            // Normalize gain at the analog center frequency mapped back to
            // the digital domain through the inverse bilinear transform
            const double fc_center = analog_to_normalized_frequency(w0);
            normalize_gain(fc_center);
        } else // FilterType::bandstop
        {
            // Bandstop: the analog zeros at s = ±j*w0 map to digital zeros at
            // z = e^{±j*omega0}, where omega0 is the inverse-warped center.
            // For order N this gives N conjugate zero pairs and the digital
            // numerator ∏ (1 - 2*cos(omega0)*z^{-1} + z^{-2})^N.

            const double omega0 =
                2.0 * std::atan(w0 / 2.0); // digital center frequency
            const double two_cos_omega0 = 2.0 * std::cos(omega0);

            std::vector<double> poly(1, 1.0);
            for (int i = 0; i < order_; ++i) {
                std::vector<double> next(poly.size() + 2, 0.0);
                for (size_t k = 0; k < poly.size(); ++k) {
                    next[k] += poly[k]; // 1 * old_coeff
                    next[k + 1] -=
                        two_cos_omega0 * poly[k]; // -2cos(ω0) * old_coeff
                    next[k + 2] += poly[k];       // 1 * old_coeff
                }
                poly.swap(next);
            }

            coeffs_.b = poly;

            // Normalize the passband gain at DC. The symmetric analog
            // prototype also gives unity gain at Nyquist.
            normalize_gain(0.0);
        }
    }

    /**
     * @brief Map lowpass prototype poles to analog bandpass poles.
     *
     * Uses the LP -> BP transform `s_L = (s^2 + w0^2) / (bw * s)`, which for
     * each prototype pole `p` yields the roots of `s^2 - bw*p*s + w0^2 = 0`:
     *
     *   s = (bw/2)*p ± sqrt(((bw/2)*p)^2 - w0^2)
     *
     * @param lp_poles Normalized lowpass prototype poles.
     * @param w0 Analog center frequency.
     * @param bw Analog bandwidth.
     */
    static std::vector<std::complex<double>> transform_lowpass_to_bandpass(
        const std::vector<std::complex<double>> &lp_poles,
        double w0,
        double bw) {
        std::vector<std::complex<double>> poles;
        poles.reserve(2 * lp_poles.size());
        const std::complex<double> w0_squared(w0 * w0, 0.0);

        for (const auto &p : lp_poles) {
            const std::complex<double> alpha = (bw / 2.0) * p;
            const std::complex<double> beta =
                std::sqrt(alpha * alpha - w0_squared);
            poles.push_back(alpha + beta);
            poles.push_back(alpha - beta);
        }

        return poles;
    }

    /**
     * @brief Map lowpass prototype poles to analog bandstop poles.
     *
     * Uses the LP -> BS transform `s_L = (bw * s) / (s^2 + w0^2)`, which for
     * each prototype pole `p` yields the roots of `p*s^2 - bw*s + p*w0^2 = 0`:
     *
     *   s = bw/(2*p) ± sqrt((bw/(2*p))^2 - w0^2)
     *
     * @param lp_poles Normalized lowpass prototype poles.
     * @param w0 Analog center frequency.
     * @param bw Analog bandwidth.
     */
    static std::vector<std::complex<double>> transform_lowpass_to_bandstop(
        const std::vector<std::complex<double>> &lp_poles,
        double w0,
        double bw) {
        std::vector<std::complex<double>> poles;
        poles.reserve(2 * lp_poles.size());
        const std::complex<double> w0_squared(w0 * w0, 0.0);

        for (const auto &p : lp_poles) {
            const std::complex<double> alpha = bw / (2.0 * p);
            const std::complex<double> beta =
                std::sqrt(alpha * alpha - w0_squared);
            poles.push_back(alpha + beta);
            poles.push_back(alpha - beta);
        }

        return poles;
    }


    /**
     * @brief Get analog Butterworth poles
     *
     * @return poles of normalized lowpass Butterworth filter (cutoff = 1s^(-1))
     */
    std::vector<std::complex<double>> get_analog_poles() const {
        std::vector<std::complex<double>> poles;
        poles.reserve(order_);

        for (int k = 0; k < order_; ++k) {
            // Pole angle: (2k + 1)π / (2N) + π/2
            double angle = std::numbers::pi * (2.0 * k + 1.0) / (2.0 * order_)
                           + std::numbers::pi / 2.0;

            poles.emplace_back(std::cos(angle), std::sin(angle));
        }

        return poles;
    }

    /**
     * @brief Convert poles to polynomial coefficients
     *
     * Given poles p1, p2, ..., pN, compute coefficients of:
     * (z - p1)(z - p2)...(z - pN) = a[0] + a[1]*z + ... + a[N]*z^N
     *
     * @return Coefficients [a0, a1, ..., aN]
     */
    std::vector<double>
    poles_to_polynomial(const std::vector<std::complex<double>> &poles) const {
        // We compute polynomial in z^{-1} form:
        //   A(z) = ∏ (1 - p_i z^{-1}) = a0 + a1 z^{-1} + ... + aN z^{-N}
        // Start with poly = [1]
        std::vector<std::complex<double>> poly(1);
        poly[0] = std::complex<double>(1.0, 0.0);

        for (const auto &p : poles) {
            std::vector<std::complex<double>> next(
                poly.size() + 1, std::complex<double>(0.0, 0.0));
            for (size_t k = 0; k < poly.size(); ++k) {
                // multiply existing terms by 1 (coefficient for z^0) -> shift 0
                next[k] += poly[k];
                // multiply existing terms by (-p) and shift by one power
                // (z^{-1})
                next[k + 1] += -p * poly[k];
            }
            poly.swap(next);
        }

        // Convert to real (imag parts should be ~0)
        std::vector<double> result(poly.size());
        for (size_t i = 0; i < poly.size(); ++i)
            result[i] = poly[i].real();

        // Normalize so a[0] == 1.0 (should already be 1)
        double a0 = result[0];
        if (std::abs(a0) < 1e-300)
            throw std::runtime_error(
                "Butterworth: numerical instability in poles_to_polynomial");
        for (auto &c : result)
            c /= a0;

        return result; // coefficients [a0, a1, ..., aN] corresponding to
                       // z^{-0}, z^{-1}, ...
    }

    /**
     * @brief Compute numerator for lowpass filter
     *
     * @return Coefficients for (1 + z^{-1})^N expansion (binomial coefficients)
     */
    std::vector<double> compute_lowpass_numerator() const {
        std::vector<double> b(order_ + 1);

        // Binomial coefficients
        b[0] = 1.0;
        for (int i = 1; i <= order_; ++i) {
            b[i] = b[i - 1] * (order_ - i + 1) / i;
        }

        return b;
    }

    /**
     * @brief Compute numerator for highpass filter
     *
     * @return Coefficients for (1 - z^{-1})^N expansion (binomial coefficients
     * with alternating signs)
     */
    std::vector<double> compute_highpass_numerator() const {
        auto b = compute_lowpass_numerator();

        // Alternate signs for highpass
        for (int i = 0; i <= order_; ++i) {
            if (i % 2 == 1) {
                b[i] = -b[i];
            }
        }

        return b;
    }

    /**
     * @brief Map an analog frequency to a normalized digital frequency.
     *
     * Inverse of the pre-warped bilinear transform (T = 1):
     * f = (2 / pi) * atan(Omega / 2).
     *
     * @param omega Analog frequency (rad/s)
     * @return Normalized frequency in [0, 1], where 1 corresponds to Nyquist
     */
    static double analog_to_normalized_frequency(double omega) {
        return 2.0 / std::numbers::pi * std::atan(omega / 2.0);
    }

    /**
     * @brief Normalize filter gain at specified frequency
     *
     * @param f Normalized frequency (0 to 1, where 1 corresponds to Nyquist)
     */
    void normalize_gain(double f) {
        // Evaluate H(e^{jω}) at frequency f where f in [0,1], 1 => Nyquist => ω
        // = π
        const std::complex<double> j(0.0, 1.0);
        const double omega = std::numbers::pi * f; // ω = π * f
        const std::complex<double> inv_z =
            std::exp(-j * omega); // z^{-1} = e^{-j ω}

        // Evaluate numerator and denominator at z^{-k}
        std::complex<double> num(0.0, 0.0), den(0.0, 0.0);

        std::complex<double> zpow(1.0, 0.0); // z^0 = 1
        for (size_t k = 0; k < coeffs_.b.size(); ++k) {
            num += coeffs_.b[k] * zpow;
            zpow *= inv_z;
        }

        zpow = std::complex<double>(1.0, 0.0);
        for (size_t k = 0; k < coeffs_.a.size(); ++k) {
            den += coeffs_.a[k] * zpow;
            zpow *= inv_z;
        }

        double gain = std::abs(num / den);
        if (gain == 0.0)
            return; // avoid division by zero

        for (auto &c : coeffs_.b)
            c /= gain;
    }
};

// ============================================================================
// Convenience functions
// ============================================================================

/**
 * @brief Design lowpass Butterworth filter
 *
 * @param order Filter order
 * @param fc Cutoff frequency (normalized, 0 < fc < 1)
 * @return Filter coefficients
 */
inline FilterCoefficients butterworth_lowpass_design(int order, double fc) {
    return ButterworthFilter::lowpass(order, fc).coefficients();
}

/**
 * @brief Design highpass Butterworth filter
 *
 * @param order Filter order
 * @param fc Cutoff frequency (normalized, 0 < fc < 1)
 * @return Filter coefficients
 */
inline FilterCoefficients butterworth_highpass_design(int order, double fc) {
    return ButterworthFilter::highpass(order, fc).coefficients();
}

/**
 * @brief Design bandpass Butterworth filter
 *
 * @param order Filter order
 * @param fc_low Low cutoff frequency (normalized, 0 < fc_low < 1)
 * @param fc_high High cutoff frequency (normalized, 0 < fc_high < 1)
 * @return Filter coefficients
 */
inline FilterCoefficients
butterworth_bandpass_design(int order, double fc_low, double fc_high) {
    return ButterworthFilter::bandpass(order, fc_low, fc_high).coefficients();
}

/**
 * @brief Design bandstop Butterworth filter
 *
 * @param order Filter order
 * @param fc_low Low cutoff frequency (normalized, 0 < fc_low < 1)
 * @param fc_high High cutoff frequency (normalized, 0 < fc_high < 1)
 * @return Filter coefficients
 */
inline FilterCoefficients
butterworth_bandstop_design(int order, double fc_low, double fc_high) {
    return ButterworthFilter::bandstop(order, fc_low, fc_high).coefficients();
}

// ============================================================================
// Convenience functions - Apply filters
// ============================================================================

/**
 * @brief Design and apply lowpass filter in one step
 *
 * @param signal Input signal
 * @param order Filter order
 * @param fc Cutoff frequency (normalized)
 * @param zero_phase Use zero-phase filtering (default: true)
 * @return Filtered signal
 */
inline std::vector<double> butterworth_lowpass(std::span<const double> signal,
                                               int order,
                                               double fc,
                                               bool zero_phase = true) {
    auto coeffs = butterworth_lowpass_design(order, fc);
    return zero_phase ? filtfilt(signal, coeffs) : filter(signal, coeffs);
}

/**
 * @brief Design and apply lowpass filter in one step
 *
 * @param signal Input signal
 * @param result Output span for filtered signal (must be same size as input)
 * @param order Filter order
 * @param fc Cutoff frequency (normalized)
 * @param zero_phase Use zero-phase filtering (default: true)
 */
inline void butterworth_lowpass(std::span<const double> signal,
                                std::span<double> result,
                                int order,
                                double fc,
                                bool zero_phase = true) {
    if (result.size() != signal.size()) {
        throw std::invalid_argument(
            "Butterworth: result span must have the same size as input signal");
    }

    auto coeffs = butterworth_lowpass_design(order, fc);
    if (zero_phase)
        filtfilt(signal, result, coeffs);
    else
        filter(signal, result, coeffs);
}

/**
 * @brief Design and apply highpass filter
 */
inline std::vector<double> butterworth_highpass(std::span<const double> signal,
                                                int order,
                                                double fc,
                                                bool zero_phase = true) {
    auto coeffs = butterworth_highpass_design(order, fc);
    return zero_phase ? filtfilt(signal, coeffs) : filter(signal, coeffs);
}

/**
 * @brief Design and apply highpass filter in one step.
 *
 * @param signal Input signal.
 * @param result Output signal, must match input size.
 * @param order Filter order.
 * @param fc Normalized cutoff frequency in (0, 1).
 * @param zero_phase Use filtfilt when true, otherwise causal filter.
 */
inline void butterworth_highpass(std::span<const double> signal,
                                 std::span<double> result,
                                 int order,
                                 double fc,
                                 bool zero_phase = true) {
    if (result.size() != signal.size()) {
        throw std::invalid_argument(
            "Butterworth: result span must have the same size as input signal");
    }

    auto coeffs = butterworth_highpass_design(order, fc);
    if (zero_phase)
        filtfilt(signal, result, coeffs);
    else
        filter(signal, result, coeffs);
}

/**
 * @brief Design and apply bandpass filter
 */
inline std::vector<double> butterworth_bandpass(std::span<const double> signal,
                                                int order,
                                                double fc_low,
                                                double fc_high,
                                                bool zero_phase = true) {
    auto coeffs = butterworth_bandpass_design(order, fc_low, fc_high);
    return zero_phase ? filtfilt(signal, coeffs) : filter(signal, coeffs);
}

/**
 * @brief Design and apply bandpass filter in one step.
 *
 * @param signal Input signal.
 * @param result Output signal, must match input size.
 * @param order Filter order.
 * @param fc_low Normalized low cutoff in (0, 1).
 * @param fc_high Normalized high cutoff in (0, 1), must be greater than
 * `fc_low`.
 * @param zero_phase Use filtfilt when true, otherwise causal filter.
 */
inline void butterworth_bandpass(std::span<const double> signal,
                                 std::span<double> result,
                                 int order,
                                 double fc_low,
                                 double fc_high,
                                 bool zero_phase = true) {
    if (result.size() != signal.size()) {
        throw std::invalid_argument(
            "Butterworth: result span must have the same size as input signal");
    }

    auto coeffs = butterworth_bandpass_design(order, fc_low, fc_high);
    if (zero_phase)
        filtfilt(signal, result, coeffs);
    else
        filter(signal, result, coeffs);
}

/**
 * @brief Design and apply bandstop filter
 */
inline std::vector<double> butterworth_bandstop(std::span<const double> signal,
                                                int order,
                                                double fc_low,
                                                double fc_high,
                                                bool zero_phase = true) {
    auto coeffs = butterworth_bandstop_design(order, fc_low, fc_high);
    return zero_phase ? filtfilt(signal, coeffs) : filter(signal, coeffs);
}

/**
 * @brief Design and apply bandstop filter in one step.
 *
 * @param signal Input signal.
 * @param result Output signal, must match input size.
 * @param order Filter order.
 * @param fc_low Normalized low cutoff in (0, 1).
 * @param fc_high Normalized high cutoff in (0, 1), must be greater than
 * `fc_low`.
 * @param zero_phase Use filtfilt when true, otherwise causal filter.
 */
inline void butterworth_bandstop(std::span<const double> signal,
                                 std::span<double> result,
                                 int order,
                                 double fc_low,
                                 double fc_high,
                                 bool zero_phase = true) {
    if (result.size() != signal.size()) {
        throw std::invalid_argument(
            "Butterworth: result span must have the same size as input signal");
    }

    auto coeffs = butterworth_bandstop_design(order, fc_low, fc_high);
    if (zero_phase)
        filtfilt(signal, result, coeffs);
    else
        filter(signal, result, coeffs);
}
} // namespace msl::signal

#endif // MSL_BUTTERWORTH_FILTER_HPP