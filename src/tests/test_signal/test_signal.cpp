#include <algorithm>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <vector>

#include "../test_utils.hpp"
#include "matrix.hpp"
#include "signal.hpp"


using namespace msl;

static std::complex<double>
frequency_response(const signal::FilterCoefficients &coeffs, double f) {
    const std::complex<double> z_inv =
        std::exp(std::complex<double>(0.0, -std::numbers::pi * f));
    const auto evaluate = [&](const std::vector<double> &c) {
        std::complex<double> value = 0.0;
        std::complex<double> power = 1.0;
        for (const double coefficient : c) {
            value += coefficient * power;
            power *= z_inv;
        }
        return value;
    };
    return evaluate(coeffs.b) / evaluate(coeffs.a);
}

static double magnitude_response(const signal::FilterCoefficients &coeffs,
                                 double f) {
    return std::abs(frequency_response(coeffs, f));
}

static void
write_complex_vector(std::ofstream &file,
                     const std::vector<std::complex<double>> &values) {
    file << std::setprecision(17);
    for (const auto &value : values) {
        file << value.real() << " " << value.imag() << "\n";
    }
}

static void write_complex_matrix(std::ofstream &file,
                                 const matrix::matrixc &mat) {
    file << std::setprecision(17);
    for (size_t i = 0; i < mat.rows(); ++i) {
        for (size_t j = 0; j < mat.cols(); ++j) {
            file << mat(i, j).real() << " " << mat(i, j).imag()
                 << (j + 1 == mat.cols() ? '\n' : ' ');
        }
    }
}

static void write_matrix(std::ofstream &file, const matrix::matrixd &mat) {
    file << std::setprecision(17);
    for (size_t i = 0; i < mat.rows(); ++i) {
        for (size_t j = 0; j < mat.cols(); ++j) {
            file << mat(i, j) << (j + 1 == mat.cols() ? '\n' : ' ');
        }
    }
}

// Digital center corresponding to the analog geometric mean of the pre-warped
// cutoffs.
static double band_center_frequency(double low, double high) {
    const double w1 = 2.0 * std::tan(std::numbers::pi * low / 2.0);
    const double w2 = 2.0 * std::tan(std::numbers::pi * high / 2.0);
    return 2.0 / std::numbers::pi * std::atan(std::sqrt(w1 * w2) / 2.0);
}

// Check that all poles of the digital transfer function lie inside the unit
// circle, using the eigenvalues of the denominator companion matrix.
static bool filter_is_stable(const signal::FilterCoefficients &coeffs) {
    if (coeffs.a.empty()) {
        return false;
    }
    const size_t order = coeffs.a.size() - 1;
    if (order == 0) {
        return true;
    }

    matrix::matrixd companion(order, order, 0.0);
    for (size_t j = 0; j < order; ++j) {
        companion(0, j) = -coeffs.a[j + 1] / coeffs.a[0];
    }
    for (size_t i = 1; i < order; ++i) {
        companion(i, i - 1) = 1.0;
    }

    const auto eigenvalues = matrix::eig(companion)[1];
    for (size_t i = 0; i < order; ++i) {
        if (std::abs(eigenvalues(i, i)) >= 1.0) {
            return false;
        }
    }
    return true;
}

int test_fft() {
    std::vector<double> x{1.0, 2.0, 3.0, 4.0};
    auto X = signal::fft(x);
    EXPECT_EQ(X.size(), 4);
    EXPECT_CPLX_NEAR(X[0], std::complex<double>(10.0, 0.0), 1e-12);
    EXPECT_CPLX_NEAR(X[1], std::complex<double>(-2.0, 2.0), 1e-12);
    EXPECT_CPLX_NEAR(X[2], std::complex<double>(-2.0, 0.0), 1e-12);
    EXPECT_CPLX_NEAR(X[3], std::complex<double>(-2.0, -2.0), 1e-12);

    auto restored = signal::ifft_real(X);
    EXPECT_EQ(restored.size(), x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        EXPECT_NEAR(restored[i], x[i], 1e-12);
    }

    auto freqs = signal::fft_frequencies(4, 8.0);
    EXPECT_NEAR(freqs[0], 0.0, 1e-12);
    EXPECT_NEAR(freqs[1], 2.0, 1e-12);
    EXPECT_NEAR(freqs[2], 4.0, 1e-12);
    EXPECT_NEAR(freqs[3], 6.0, 1e-12);

    auto mag = signal::magnitude_spectrum(X);
    EXPECT_NEAR(mag[0], 10.0, 1e-12);
    EXPECT_NEAR(mag[1], std::sqrt(8.0), 1e-12);

    std::vector<double> short_input{1.0, -2.0, 0.5};
    std::vector<double> explicit_padding(8, 0.0);
    std::copy(short_input.begin(), short_input.end(), explicit_padding.begin());
    const auto automatic = signal::fft(short_input, 8);
    const auto explicit_result = signal::fft(explicit_padding);
    for (size_t i = 0; i < automatic.size(); ++i) {
        EXPECT_CPLX_NEAR(automatic[i], explicit_result[i], 1e-12);
    }
    const auto padded_round_trip = signal::ifft_real(automatic);
    for (size_t i = 0; i < short_input.size(); ++i) {
        EXPECT_NEAR(padded_round_trip[i], short_input[i], 1e-12);
    }
    for (size_t i = short_input.size(); i < padded_round_trip.size(); ++i) {
        EXPECT_NEAR(padded_round_trip[i], 0.0, 1e-12);
    }

    std::vector<std::complex<double>> complex_input{{1.0, 2.0}, {-0.5, 1.0}};
    std::vector<std::complex<double>> complex_padding(4, {0.0, 0.0});
    std::copy(
        complex_input.begin(), complex_input.end(), complex_padding.begin());
    const auto complex_automatic = signal::fft(complex_input, 4);
    const auto complex_explicit = signal::fft(complex_padding);
    for (size_t i = 0; i < complex_automatic.size(); ++i) {
        EXPECT_CPLX_NEAR(complex_automatic[i], complex_explicit[i], 1e-12);
    }

    matrix::matrixd samples(2, 2, {1.0, 2.0, 3.0, 4.0});
    const auto row_fft = signal::fft_rows(samples, 4);
    const auto expected_row = signal::fft(std::vector<double>{1.0, 2.0}, 4);
    for (size_t k = 0; k < 4; ++k) {
        EXPECT_CPLX_NEAR(row_fft(0, k), expected_row[k], 1e-12);
    }
    const auto column_fft = signal::fft_columns(samples, 4);
    const auto expected_column = signal::fft(std::vector<double>{1.0, 3.0}, 4);
    for (size_t k = 0; k < 4; ++k) {
        EXPECT_CPLX_NEAR(column_fft(k, 0), expected_column[k], 1e-12);
    }

    return 0;
}

int test_fft_matrix_roundtrip() {
    matrix::matrixd real_rows(3, 4);
    for (size_t i = 0; i < real_rows.rows(); ++i) {
        for (size_t j = 0; j < real_rows.cols(); ++j) {
            real_rows(i, j) = 1.0 + 0.5 * static_cast<double>(i)
                              + 0.25 * static_cast<double>(j);
        }
    }

    // Real row round-trip
    const auto real_spectrum = signal::fft_rows(real_rows);
    const auto real_restored = signal::ifft_rows_real(real_spectrum);
    EXPECT_EQ(real_restored.rows(), real_rows.rows());
    EXPECT_EQ(real_restored.cols(), real_rows.cols());
    for (size_t i = 0; i < real_rows.rows(); ++i) {
        for (size_t j = 0; j < real_rows.cols(); ++j) {
            EXPECT_NEAR(real_restored(i, j), real_rows(i, j), 1e-12);
        }
    }

    // Real row round-trip with zero padding
    const auto padded_spectrum = signal::fft_rows(real_rows, 8);
    const auto padded_restored = signal::ifft_rows_real(padded_spectrum, 8);
    EXPECT_EQ(padded_restored.rows(), real_rows.rows());
    EXPECT_EQ(padded_restored.cols(), 8);
    for (size_t i = 0; i < real_rows.rows(); ++i) {
        for (size_t j = 0; j < real_rows.cols(); ++j) {
            EXPECT_NEAR(padded_restored(i, j), real_rows(i, j), 1e-12);
        }
        for (size_t j = real_rows.cols(); j < 8; ++j) {
            EXPECT_NEAR(padded_restored(i, j), 0.0, 1e-12);
        }
    }

    // Complex row round-trip
    matrix::matrixc complex_rows(2, 3);
    for (size_t i = 0; i < complex_rows.rows(); ++i) {
        for (size_t j = 0; j < complex_rows.cols(); ++j) {
            complex_rows(i, j) =
                std::complex<double>(0.5 * static_cast<double>(i + 1),
                                     -0.25 * static_cast<double>(j + 1));
        }
    }
    const auto complex_spectrum = signal::fft_rows(complex_rows);
    const auto complex_restored = signal::ifft_rows(complex_spectrum);
    EXPECT_EQ(complex_restored.rows(), complex_rows.rows());
    EXPECT_EQ(complex_restored.cols(), complex_rows.cols());
    for (size_t i = 0; i < complex_rows.rows(); ++i) {
        for (size_t j = 0; j < complex_rows.cols(); ++j) {
            EXPECT_CPLX_NEAR(complex_restored(i, j), complex_rows(i, j), 1e-12);
        }
    }

    // Complex row round-trip with zero padding
    const auto complex_padded_spectrum = signal::fft_rows(complex_rows, 5);
    const auto complex_padded_restored =
        signal::ifft_rows(complex_padded_spectrum, 5);
    EXPECT_EQ(complex_padded_restored.rows(), complex_rows.rows());
    EXPECT_EQ(complex_padded_restored.cols(), 5);
    for (size_t i = 0; i < complex_rows.rows(); ++i) {
        for (size_t j = 0; j < complex_rows.cols(); ++j) {
            EXPECT_CPLX_NEAR(
                complex_padded_restored(i, j), complex_rows(i, j), 1e-12);
        }
        for (size_t j = complex_rows.cols(); j < 5; ++j) {
            EXPECT_CPLX_NEAR(complex_padded_restored(i, j),
                             std::complex<double>(0.0, 0.0),
                             1e-12);
        }
    }

    // Column round-trip for real and complex inputs
    const auto real_column_restored =
        signal::ifft_columns_real(signal::fft_columns(real_rows));
    for (size_t i = 0; i < real_rows.rows(); ++i) {
        for (size_t j = 0; j < real_rows.cols(); ++j) {
            EXPECT_NEAR(real_column_restored(i, j), real_rows(i, j), 1e-12);
        }
    }
    const auto complex_column_restored =
        signal::ifft_columns(signal::fft_columns(complex_rows));
    for (size_t i = 0; i < complex_rows.rows(); ++i) {
        for (size_t j = 0; j < complex_rows.cols(); ++j) {
            EXPECT_CPLX_NEAR(
                complex_column_restored(i, j), complex_rows(i, j), 1e-12);
        }
    }

    return 0;
}

int test_butterworth_bandpass() {
    struct Band {
        double low;
        double high;
    };
    const Band bands[] = {{0.05, 0.1}, {0.1, 0.3}, {0.004, 0.8}};
    const double cutoff_gain = 1.0 / std::sqrt(2.0);

    for (const auto &band : bands) {
        const auto coeffs =
            signal::butterworth_bandpass_design(2, band.low, band.high);
        const double center = band_center_frequency(band.low, band.high);

        EXPECT_NEAR(magnitude_response(coeffs, center), 1.0, 1e-9);
        EXPECT_NEAR(magnitude_response(coeffs, band.low), cutoff_gain, 1e-9);
        EXPECT_NEAR(magnitude_response(coeffs, band.high), cutoff_gain, 1e-9);
        EXPECT_TRUE(magnitude_response(coeffs, 0.0) < 1e-12);
        EXPECT_TRUE(magnitude_response(coeffs, 1.0) < 1e-12);
        EXPECT_TRUE(filter_is_stable(coeffs));
    }

    return 0;
}

int test_butterworth_bandstop() {
    struct Band {
        double low;
        double high;
    };
    const Band bands[] = {{0.05, 0.1}, {0.1, 0.3}, {0.004, 0.8}};
    const double cutoff_gain = 1.0 / std::sqrt(2.0);

    for (const auto &band : bands) {
        const auto coeffs =
            signal::butterworth_bandstop_design(2, band.low, band.high);
        const double notch = band_center_frequency(band.low, band.high);

        EXPECT_NEAR(magnitude_response(coeffs, notch), 0.0, 1e-9);
        EXPECT_NEAR(magnitude_response(coeffs, 0.0), 1.0, 1e-9);
        EXPECT_NEAR(magnitude_response(coeffs, 1.0), 1.0, 1e-9);
        EXPECT_NEAR(magnitude_response(coeffs, band.low), cutoff_gain, 1e-9);
        EXPECT_NEAR(magnitude_response(coeffs, band.high), cutoff_gain, 1e-9);
        EXPECT_TRUE(filter_is_stable(coeffs));
    }

    return 0;
}

int test_butterworth_lowpass_highpass() {
    const double cutoff_gain = 1.0 / std::sqrt(2.0);
    const double cutoffs[] = {0.004, 0.1, 0.3, 0.8};

    for (const double fc : cutoffs) {
        const auto lowpass = signal::butterworth_lowpass_design(2, fc);
        EXPECT_NEAR(magnitude_response(lowpass, 0.0), 1.0, 1e-9);
        EXPECT_NEAR(magnitude_response(lowpass, fc), cutoff_gain, 1e-9);
        EXPECT_TRUE(filter_is_stable(lowpass));

        const auto highpass = signal::butterworth_highpass_design(2, fc);
        EXPECT_NEAR(magnitude_response(highpass, 1.0), 1.0, 1e-9);
        EXPECT_NEAR(magnitude_response(highpass, fc), cutoff_gain, 1e-9);
        EXPECT_TRUE(filter_is_stable(highpass));
    }

    // Odd and higher orders
    const auto order3 = signal::butterworth_lowpass_design(3, 0.25);
    EXPECT_NEAR(magnitude_response(order3, 0.25), cutoff_gain, 1e-9);
    EXPECT_TRUE(filter_is_stable(order3));

    const auto order4 = signal::butterworth_highpass_design(4, 0.25);
    EXPECT_NEAR(magnitude_response(order4, 0.25), cutoff_gain, 1e-9);
    EXPECT_TRUE(filter_is_stable(order4));

    return 0;
}

int test_welch_and_covariance() {
    const std::vector<double> x{1.0, 2.0, 0.0, -1.0, 3.0, 2.0};
    const std::vector<double> y{0.0, 1.0, 2.0, 1.0, -1.0, 2.0};
    const std::vector<double> window(4, 1.0);
    constexpr size_t segment_length = 4;
    constexpr size_t overlap = 2;
    constexpr size_t nfft = 8;
    constexpr double sampling_rate = 8.0;

    const auto estimate = signal::cpsd_welch(
        x, y, window, overlap, segment_length, nfft, sampling_rate);
    const std::vector<std::complex<double>> expected_spectrum{
        {0.375, 0.0},
        {0.05981917382415922, 0.05445752147247765},
        {-0.125, 0.28125},
        {-0.02856917382415922, 0.21070752147247765},
        {-0.0625, 0.0},
        {-0.02856917382415922, -0.21070752147247765},
        {-0.125, -0.28125},
        {0.05981917382415922, -0.05445752147247765}};
    for (size_t k = 0; k < nfft; ++k) {
        EXPECT_CPLX_NEAR(estimate[k], expected_spectrum[k], 1e-12);
    }

    const auto reverse = signal::cpsd_welch(
        y, x, window, overlap, segment_length, nfft, sampling_rate);
    const auto auto_power = signal::psd_welch(
        x, window, overlap, segment_length, nfft, sampling_rate);
    for (size_t k = 0; k < nfft; ++k) {
        EXPECT_CPLX_NEAR(reverse[k], std::conj(estimate[k]), 1e-12);
        EXPECT_TRUE(auto_power[k] >= -1e-14);
    }
    const auto frequencies = signal::fft_frequencies(nfft, sampling_rate);
    EXPECT_NEAR(frequencies[1] - frequencies[0], 1.0, 1e-12);

    const std::vector<double> covariance_x{1.0, 2.0, 3.0};
    const std::vector<double> covariance_y{3.0, 1.0, 0.0};
    const auto covariance =
        signal::xcov_unbiased(covariance_x, covariance_y, 2);
    const std::vector<double> expected_covariance{
        4.0 / 3.0, 1.0 / 6.0, -1.0, -1.0 / 6.0, 5.0 / 3.0};
    EXPECT_EQ(covariance.size(), expected_covariance.size());
    for (size_t i = 0; i < covariance.size(); ++i) {
        EXPECT_NEAR(covariance[i], expected_covariance[i], 1e-12);
    }
    const auto reverse_covariance =
        signal::xcov_unbiased(covariance_y, covariance_x, 2);
    for (size_t i = 0; i < covariance.size(); ++i) {
        EXPECT_NEAR(covariance[i],
                    reverse_covariance[covariance.size() - 1 - i],
                    1e-12);
    }

    return 0;
}

int test_power_spectrum() {
    const std::vector<double> x{1.0, -2.0, 0.5, 3.0};
    const std::vector<double> y{0.5, 1.0, -1.0, 2.0};
    constexpr size_t nfft = 8;

    // Single-frame power spectrum equals |FFT|^2
    const auto power = signal::power_spectrum(x, nfft);
    EXPECT_EQ(power.size(), nfft);
    const auto spectrum = signal::fft(x, nfft);
    for (size_t k = 0; k < nfft; ++k) {
        EXPECT_NEAR(power[k], std::norm(spectrum[k]), 1e-12);
        EXPECT_TRUE(power[k] >= 0.0);
    }

    // Zero-padding: explicit padded input matches automatic padding
    std::vector<double> padded(nfft, 0.0);
    std::copy(x.begin(), x.end(), padded.begin());
    const auto explicit_power = signal::power_spectrum(padded, nfft);
    for (size_t k = 0; k < nfft; ++k) {
        EXPECT_NEAR(power[k], explicit_power[k], 1e-12);
    }

    // Cross power spectrum equals X * conj(Y) and is conjugate antisymmetric
    const auto cross = signal::cross_power_spectrum(x, y, nfft);
    const auto reverse_cross = signal::cross_power_spectrum(y, x, nfft);
    const auto y_spectrum = signal::fft(y, nfft);
    for (size_t k = 0; k < nfft; ++k) {
        EXPECT_CPLX_NEAR(
            cross[k], spectrum[k] * std::conj(y_spectrum[k]), 1e-12);
        EXPECT_CPLX_NEAR(reverse_cross[k], std::conj(cross[k]), 1e-12);
    }

    // Invalid nfft and mismatched input sizes are rejected
    bool threw_nfft = false;
    try {
        (void)signal::cross_power_spectrum(x, y, 0);
    } catch (const std::invalid_argument &) {
        threw_nfft = true;
    }
    EXPECT_TRUE(threw_nfft);

    bool threw_size = false;
    try {
        (void)signal::cross_power_spectrum(x, std::vector<double>{1.0}, nfft);
    } catch (const std::invalid_argument &) {
        threw_size = true;
    }
    EXPECT_TRUE(threw_size);

    return 0;
}

int test_windows_and_filter() {
    auto hann = signal::hann_window(5);
    std::vector<double> hann_expected{0.0, 0.5, 1.0, 0.5, 0.0};
    for (size_t i = 0; i < hann.size(); ++i) {
        EXPECT_NEAR(hann[i], hann_expected[i], 1e-12);
    }

    auto hamming = signal::hamming_window(5);
    std::vector<double> hamming_expected{0.08, 0.54, 1.0, 0.54, 0.08};
    for (size_t i = 0; i < hamming.size(); ++i) {
        EXPECT_NEAR(hamming[i], hamming_expected[i], 1e-12);
    }

    std::vector<double> input{1.0, 2.0, 4.0, 8.0};
    signal::FilterCoefficients coeffs({0.5, 0.5}, {1.0});
    auto filtered = signal::filter(input, coeffs);
    std::vector<double> expected{0.5, 1.5, 3.0, 6.0};
    for (size_t i = 0; i < filtered.size(); ++i) {
        EXPECT_NEAR(filtered[i], expected[i], 1e-12);
    }

    matrix::matrixd mat(4, 2);
    for (size_t i = 0; i < mat.rows(); ++i) {
        mat(i, 0) = input[i];
        mat(i, 1) = 2.0 * input[i];
    }
    auto filtered_cols = signal::filter_columns(mat, coeffs);
    EXPECT_EQ(filtered_cols.rows(), 4);
    EXPECT_EQ(filtered_cols.cols(), 2);
    EXPECT_NEAR(filtered_cols(3, 0), 6.0, 1e-12);
    EXPECT_NEAR(filtered_cols(3, 1), 12.0, 1e-12);

    bool threw = false;
    try {
        signal::FilterCoefficients bad({1.0}, {2.0});
        (void)signal::filter(input, bad);
    } catch (const std::invalid_argument &) {
        threw = true;
    }
    EXPECT_TRUE(threw);

    return 0;
}

int main() {
    int result = test_fft() + test_fft_matrix_roundtrip()
                 + test_butterworth_lowpass_highpass()
                 + test_butterworth_bandpass() + test_butterworth_bandstop()
                 + test_welch_and_covariance() + test_power_spectrum()
                 + test_windows_and_filter();

    auto compare_dir = msl::test::result_dir("signal");

    std::vector<double> x{1.0, 2.0, 3.0, 4.0};
    auto X = signal::fft(x);
    auto restored = signal::ifft_real(X);
    auto hann = signal::hann_window(5);
    auto hamming = signal::hamming_window(5);
    signal::FilterCoefficients coeffs({0.5, 0.5}, {1.0});
    auto filtered = signal::filter(x, coeffs);

    std::ofstream fft_file(compare_dir / "signal_fft.txt");
    fft_file << std::setprecision(17);
    for (size_t i = 0; i < X.size(); ++i) {
        fft_file << x[i] << " " << X[i].real() << " " << X[i].imag() << " "
                 << restored[i] << "\n";
    }

    std::ofstream window_file(compare_dir / "signal_windows.txt");
    window_file << std::setprecision(17);
    for (size_t i = 0; i < hann.size(); ++i) {
        window_file << hann[i] << " " << hamming[i] << "\n";
    }

    std::ofstream filter_file(compare_dir / "signal_filter.txt");
    filter_file << std::setprecision(17);
    for (size_t i = 0; i < filtered.size(); ++i) {
        filter_file << x[i] << " " << filtered[i] << "\n";
    }

    // ---- Butterworth comparison data for MATLAB validation ----
    struct ButterCase {
        int type; // 1 = lowpass, 2 = highpass, 3 = bandpass, 4 = bandstop
        int order;
        double low;
        double high;
    };
    const std::vector<ButterCase> butter_cases{{1, 2, 0.3, 0.0},
                                               {2, 2, 0.3, 0.0},
                                               {3, 2, 0.1, 0.3},
                                               {4, 2, 0.1, 0.3},
                                               {3, 2, 0.004, 0.8}};

    std::ofstream cases_file(compare_dir / "butterworth_cases.txt");
    cases_file << std::setprecision(17);
    for (const auto &c : butter_cases) {
        cases_file << c.type << " " << c.order << " " << c.low << " " << c.high
                   << "\n";
    }

    std::vector<signal::FilterCoefficients> butter_coeffs;
    butter_coeffs.reserve(butter_cases.size());
    for (const auto &c : butter_cases) {
        switch (c.type) {
            case 1:
                butter_coeffs.push_back(
                    signal::butterworth_lowpass_design(c.order, c.low));
                break;
            case 2:
                butter_coeffs.push_back(
                    signal::butterworth_highpass_design(c.order, c.low));
                break;
            case 3:
                butter_coeffs.push_back(signal::butterworth_bandpass_design(
                    c.order, c.low, c.high));
                break;
            default:
                butter_coeffs.push_back(signal::butterworth_bandstop_design(
                    c.order, c.low, c.high));
                break;
        }
    }

    constexpr size_t response_points = 101;
    std::ofstream response_file(compare_dir / "butterworth_response.txt");
    response_file << std::setprecision(17);
    for (size_t i = 0; i < response_points; ++i) {
        const double f =
            static_cast<double>(i) / static_cast<double>(response_points - 1);
        response_file << f;
        for (const auto &c : butter_coeffs) {
            response_file << " " << magnitude_response(c, f);
        }
        response_file << "\n";
    }

    constexpr size_t butter_signal_length = 64;
    std::vector<double> butter_input(butter_signal_length);
    for (size_t i = 0; i < butter_signal_length; ++i) {
        const double t = static_cast<double>(i);
        butter_input[i] = std::sin(2.0 * std::numbers::pi * 0.05 * t)
                          + 0.5 * std::sin(2.0 * std::numbers::pi * 0.35 * t)
                          + 0.01 * t;
    }

    std::vector<std::vector<double>> butter_filtered;
    std::vector<std::vector<double>> butter_filtfilt;
    for (const auto &c : butter_coeffs) {
        butter_filtered.push_back(signal::filter(butter_input, c));
        butter_filtfilt.push_back(signal::filtfilt(butter_input, c));
    }

    std::ofstream butter_filter_file(compare_dir / "butterworth_filter.txt");
    butter_filter_file << std::setprecision(17);
    std::ofstream butter_filtfilt_file(compare_dir
                                       / "butterworth_filtfilt.txt");
    butter_filtfilt_file << std::setprecision(17);
    for (size_t i = 0; i < butter_signal_length; ++i) {
        butter_filter_file << butter_input[i];
        butter_filtfilt_file << butter_input[i];
        for (size_t k = 0; k < butter_coeffs.size(); ++k) {
            butter_filter_file << " " << butter_filtered[k][i];
            butter_filtfilt_file << " " << butter_filtfilt[k][i];
        }
        butter_filter_file << "\n";
        butter_filtfilt_file << "\n";
    }

    // ---- Matrix FFT/IFFT comparison data ----
    matrix::matrixd fft_matrix(4, 3);
    for (size_t i = 0; i < fft_matrix.rows(); ++i) {
        for (size_t j = 0; j < fft_matrix.cols(); ++j) {
            fft_matrix(i, j) = std::sin(0.5 * static_cast<double>(i)
                                        + 0.7 * static_cast<double>(j))
                               + 0.1 * static_cast<double>(i * j);
        }
    }
    const auto fft_rows_result = signal::fft_rows(fft_matrix);
    const auto ifft_rows_result = signal::ifft_rows_real(fft_rows_result);
    const auto fft_columns_result = signal::fft_columns(fft_matrix);
    const auto ifft_columns_result =
        signal::ifft_columns_real(fft_columns_result);

    std::ofstream fft_matrix_file(compare_dir / "signal_fft_matrix_input.txt");
    write_matrix(fft_matrix_file, fft_matrix);
    std::ofstream fft_rows_file(compare_dir / "signal_fft_rows.txt");
    write_complex_matrix(fft_rows_file, fft_rows_result);
    std::ofstream ifft_rows_file(compare_dir / "signal_ifft_rows.txt");
    write_matrix(ifft_rows_file, ifft_rows_result);
    std::ofstream fft_columns_file(compare_dir / "signal_fft_columns.txt");
    write_complex_matrix(fft_columns_file, fft_columns_result);
    std::ofstream ifft_columns_file(compare_dir / "signal_ifft_columns.txt");
    write_matrix(ifft_columns_file, ifft_columns_result);

    // ---- Welch PSD/CPSD comparison data ----
    constexpr size_t welch_segment = 8;
    constexpr size_t welch_overlap = 4;
    constexpr size_t welch_nfft = 16;
    constexpr double welch_fs = 8.0;
    std::vector<double> welch_x(32);
    std::vector<double> welch_y(32);
    for (size_t n = 0; n < welch_x.size(); ++n) {
        const double t = static_cast<double>(n);
        welch_x[n] = std::sin(2.0 * std::numbers::pi * 0.1 * t)
                     + 0.5 * std::cos(2.0 * std::numbers::pi * 0.25 * t) + 0.1;
        welch_y[n] = std::cos(2.0 * std::numbers::pi * 0.1 * t)
                     - 0.25 * std::sin(2.0 * std::numbers::pi * 0.25 * t);
    }
    const auto welch_window = signal::hann_window(welch_segment);
    const auto welch_psd = signal::psd_welch(welch_x,
                                             welch_window,
                                             welch_overlap,
                                             welch_segment,
                                             welch_nfft,
                                             welch_fs);
    const auto welch_cpsd = signal::cpsd_welch(welch_x,
                                               welch_y,
                                               welch_window,
                                               welch_overlap,
                                               welch_segment,
                                               welch_nfft,
                                               welch_fs);
    std::ofstream welch_x_file(compare_dir / "signal_welch_x.txt");
    welch_x_file << std::setprecision(17);
    std::ofstream welch_y_file(compare_dir / "signal_welch_y.txt");
    welch_y_file << std::setprecision(17);
    for (size_t n = 0; n < welch_x.size(); ++n) {
        welch_x_file << welch_x[n] << "\n";
        welch_y_file << welch_y[n] << "\n";
    }
    std::ofstream welch_window_file(compare_dir / "signal_welch_window.txt");
    welch_window_file << std::setprecision(17);
    for (const double value : welch_window) {
        welch_window_file << value << "\n";
    }
    std::ofstream welch_psd_file(compare_dir / "signal_welch_psd.txt");
    welch_psd_file << std::setprecision(17);
    for (const double value : welch_psd) {
        welch_psd_file << value << "\n";
    }
    std::ofstream welch_cpsd_file(compare_dir / "signal_welch_cpsd.txt");
    write_complex_vector(welch_cpsd_file, welch_cpsd);

    // ---- Unbiased cross-covariance comparison data ----
    std::vector<double> xcov_x(8);
    std::vector<double> xcov_y(8);
    for (size_t n = 0; n < xcov_x.size(); ++n) {
        const double t = static_cast<double>(n);
        xcov_x[n] = std::sin(0.4 * t) + 0.05 * t;
        xcov_y[n] = std::cos(0.3 * t) - 0.02 * t * t;
    }
    const auto xcov_result = signal::xcov_unbiased(xcov_x, xcov_y, 3);
    std::ofstream xcov_input_file(compare_dir / "signal_xcov_inputs.txt");
    xcov_input_file << std::setprecision(17);
    for (size_t n = 0; n < xcov_x.size(); ++n) {
        xcov_input_file << xcov_x[n] << " " << xcov_y[n] << "\n";
    }
    std::ofstream xcov_file(compare_dir / "signal_xcov.txt");
    xcov_file << std::setprecision(17);
    for (const double value : xcov_result) {
        xcov_file << value << "\n";
    }

    // ---- Detrend comparison data ----
    std::vector<double> detrend_input(16);
    for (size_t n = 0; n < detrend_input.size(); ++n) {
        const double t = static_cast<double>(n);
        detrend_input[n] =
            0.5 + 0.05 * t + std::sin(2.0 * std::numbers::pi * 0.2 * t);
    }
    const auto detrend_mean_removed = signal::detrend(detrend_input, 0);
    const auto detrend_linear_removed = signal::detrend(detrend_input, 1);
    std::ofstream detrend_file(compare_dir / "signal_detrend.txt");
    detrend_file << std::setprecision(17);
    for (size_t n = 0; n < detrend_input.size(); ++n) {
        detrend_file << detrend_input[n] << " " << detrend_mean_removed[n]
                     << " " << detrend_linear_removed[n] << "\n";
    }

    return result + msl::test::summary();
}
