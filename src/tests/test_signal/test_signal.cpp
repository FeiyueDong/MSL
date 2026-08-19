#include <algorithm>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "matrix.hpp"
#include "signal.hpp"
#include "signal/filter.hpp"

using namespace msl;

static int g_failures = 0;
static int g_total = 0;

#define EXPECT_TRUE(cond)                                                      \
    do {                                                                       \
        ++g_total;                                                             \
        if (!(cond)) {                                                         \
            ++g_failures;                                                      \
            std::cerr << "[FAIL] " << __FILE__ << ":" << __LINE__ << " - "     \
                      << #cond << "\n";                                        \
        }                                                                      \
    } while (0)

#define EXPECT_EQ(a, b) EXPECT_TRUE((a) == (b))

#define EXPECT_NEAR(a, b, eps)                                                 \
    do {                                                                       \
        ++g_total;                                                             \
        if (std::fabs((a) - (b)) > (eps)) {                                    \
            ++g_failures;                                                      \
            std::cerr << "[FAIL] " << __FILE__ << ":" << __LINE__ << " - "     \
                      << #a << " ~= " << #b << " (got: " << (a) << " vs "      \
                      << (b) << ")\n";                                         \
        }                                                                      \
    } while (0)

#define EXPECT_CPLX_NEAR(a, b, eps)                                            \
    do {                                                                       \
        const auto actual = (a);                                               \
        const auto expected = (b);                                             \
        EXPECT_NEAR(actual.real(), expected.real(), eps);                      \
        EXPECT_NEAR(actual.imag(), expected.imag(), eps);                      \
    } while (0)

std::filesystem::path project_root() {
    auto path = std::filesystem::current_path();
    while (!path.empty()) {
        if (std::filesystem::exists(path / "msl")
            && std::filesystem::exists(path / "tests")) {
            return path;
        }
        auto parent = path.parent_path();
        if (parent == path) {
            break;
        }
        path = parent;
    }
    return std::filesystem::current_path();
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
    int result =
        test_fft() + test_welch_and_covariance() + test_windows_and_filter();

    auto compare_dir =
        project_root() / "test_result" / "signal" / "matlab_compare";
    std::filesystem::create_directories(compare_dir);

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

    std::cout << "Total checks: " << g_total << ", failures: " << g_failures
              << "\n";
    return result + (g_failures == 0 ? 0 : 1);
}
