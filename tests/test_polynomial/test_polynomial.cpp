#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "../test_utils.hpp"
#include "polynomial.hpp"
#include "polynomial/polynomial.hpp"


using namespace msl;

int test_polynomial() {
    polynomial::Polynomial direct(std::vector<double>{1.0, 2.0, 3.0});
    EXPECT_EQ(direct.degree(), 2);
    EXPECT_NEAR(direct(2.0), 17.0, 1e-12);
    EXPECT_NEAR(direct.derivative(2.0), 14.0, 1e-12);

    std::vector<double> query{0.0, 1.0, 2.0};
    auto values = direct(query);
    EXPECT_EQ(values.size(), query.size());
    EXPECT_NEAR(values[0], 1.0, 1e-12);
    EXPECT_NEAR(values[1], 6.0, 1e-12);
    EXPECT_NEAR(values[2], 17.0, 1e-12);

    std::vector<double> out(query.size());
    polynomial::polyval(std::vector<double>{1.0, 2.0, 3.0}, query, out);
    EXPECT_NEAR(out[2], 17.0, 1e-12);

    std::vector<double> x{-2.0, -1.0, 0.0, 1.0, 2.0};
    std::vector<double> y;
    for (double xi : x) {
        y.push_back(1.0 - 2.0 * xi + 0.5 * xi * xi);
    }

    auto coeffs = polynomial::polyfit(x, y, 2);
    EXPECT_EQ(coeffs.size(), 3);
    EXPECT_NEAR(coeffs[0], 1.0, 1e-10);
    EXPECT_NEAR(coeffs[1], -2.0, 1e-10);
    EXPECT_NEAR(coeffs[2], 0.5, 1e-10);

    std::vector<double> coeffs_out(3);
    polynomial::polyfit(x, y, coeffs_out, 2);
    EXPECT_NEAR(coeffs_out[0], 1.0, 1e-10);
    EXPECT_NEAR(coeffs_out[1], -2.0, 1e-10);
    EXPECT_NEAR(coeffs_out[2], 0.5, 1e-10);

    auto constant_fit = polynomial::Polynomial::from_fit(y, 0);
    EXPECT_EQ(constant_fit.degree(), 0);

    bool threw = false;
    try {
        polynomial::Polynomial bad;
        (void)bad.degree();
    } catch (const std::runtime_error &) {
        threw = true;
    }
    EXPECT_TRUE(threw);

    // Noisy least-squares fit stays close to the underlying quadratic
    std::vector<double> noisy_x(20);
    std::vector<double> noisy_y(20);
    for (size_t i = 0; i < noisy_x.size(); ++i) {
        noisy_x[i] = -2.0 + 0.2 * static_cast<double>(i);
        const double xi = noisy_x[i];
        noisy_y[i] =
            1.0 - 2.0 * xi + 0.5 * xi * xi + 0.01 * std::sin(10.0 * xi);
    }
    const auto noisy_coeffs = polynomial::polyfit(noisy_x, noisy_y, 2);
    EXPECT_EQ(noisy_coeffs.size(), 3);
    EXPECT_NEAR(noisy_coeffs[0], 1.0, 0.02);
    EXPECT_NEAR(noisy_coeffs[1], -2.0, 0.02);
    EXPECT_NEAR(noisy_coeffs[2], 0.5, 0.02);

    // Degree too large for the sample count is rejected
    bool threw_degree = false;
    try {
        (void)polynomial::polyfit(std::vector<double>{0.0, 1.0, 2.0},
                                  std::vector<double>{1.0, 2.0, 3.0},
                                  3);
    } catch (const std::invalid_argument &) {
        threw_degree = true;
    }
    EXPECT_TRUE(threw_degree);

    // Duplicate sample locations make the Vandermonde system rank deficient
    bool threw_rank = false;
    try {
        (void)polynomial::polyfit(std::vector<double>{0.0, 0.0, 0.0},
                                  std::vector<double>{1.0, 2.0, 3.0},
                                  2);
    } catch (const std::runtime_error &) {
        threw_rank = true;
    }
    EXPECT_TRUE(threw_rank);

    auto compare_dir = msl::test::result_dir("polynomial");
    std::ofstream coeff_file(compare_dir / "polynomial_coefficients.txt");
    coeff_file << std::setprecision(17);
    for (double c : coeffs) {
        coeff_file << c << "\n";
    }

    std::ofstream values_file(compare_dir / "polynomial_values.txt");
    values_file << std::setprecision(17);
    for (size_t i = 0; i < query.size(); ++i) {
        values_file << query[i] << " " << values[i] << " " << out[i] << "\n";
    }

    return msl::test::summary();
}

int main() { return test_polynomial(); }
