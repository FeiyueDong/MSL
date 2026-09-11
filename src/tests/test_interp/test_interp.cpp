#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "interp.hpp"
#include "test_utils.hpp"

using namespace msl;

static void expect_vector_near(const std::vector<double> &actual,
                               const std::vector<double> &expected,
                               double eps) {
    EXPECT_EQ(actual.size(), expected.size());
    for (size_t i = 0; i < actual.size() && i < expected.size(); ++i) {
        EXPECT_NEAR(actual[i], expected[i], eps);
    }
}

int test_interp() {
    std::vector<double> x{0.0, 1.0, 2.0, 3.0};
    std::vector<double> y{1.0, 3.0, 5.0, 7.0};
    std::vector<double> xq{-1.0, 0.5, 1.5, 3.0, 4.0};
    std::vector<double> linear_expected{-1.0, 2.0, 4.0, 7.0, 9.0};

    expect_vector_near(
        interp::interp1_linear(x, y, xq), linear_expected, 1e-12);
    expect_vector_near(interp::interp1_cubic(x, y, xq), linear_expected, 1e-10);
    expect_vector_near(interp::interp1_pchip(x, y, xq), linear_expected, 1e-12);
    expect_vector_near(interp::interp1_akima(x, y, xq), linear_expected, 1e-12);

    std::vector<double> near_xq{0.25, 1.4, 2.6};
    std::vector<double> near_expected{1.0, 3.0, 7.0};
    expect_vector_near(
        interp::interp1_near(x, y, near_xq), near_expected, 1e-12);

    std::vector<double> out(xq.size());
    interp::interp1_linear(x, y, xq, out);
    expect_vector_near(out, linear_expected, 1e-12);

    bool threw = false;
    try {
        (void)interp::interp1_linear(std::vector<double>{0.0, 0.0},
                                     std::vector<double>{1.0, 2.0},
                                     std::vector<double>{0.5});
    } catch (const std::invalid_argument &) {
        threw = true;
    }
    EXPECT_TRUE(threw);

    auto compare_dir = msl::test::result_dir("interp");

    // Nonlinear, non-uniform dataset that distinguishes the methods
    const std::vector<double> reference_x{0.0, 0.7, 1.8, 3.0, 5.0};
    const std::vector<double> reference_y{0.0, 1.2, 0.8, 2.5, 1.0};
    for (size_t i = 0; i < reference_x.size(); ++i) {
        const std::vector<double> point{reference_x[i]};
        EXPECT_NEAR(interp::interp1_linear(reference_x, reference_y, point)[0],
                    reference_y[i],
                    1e-12);
        EXPECT_NEAR(interp::interp1_cubic(reference_x, reference_y, point)[0],
                    reference_y[i],
                    1e-12);
        EXPECT_NEAR(interp::interp1_pchip(reference_x, reference_y, point)[0],
                    reference_y[i],
                    1e-12);
        EXPECT_NEAR(interp::interp1_akima(reference_x, reference_y, point)[0],
                    reference_y[i],
                    1e-12);
        EXPECT_NEAR(
            interp::interp1_polynomial(reference_x, reference_y, point)[0],
            reference_y[i],
            1e-10);
    }

    std::vector<double> dense_xq(101);
    for (size_t i = 0; i < dense_xq.size(); ++i) {
        dense_xq[i] = 5.0 * static_cast<double>(i)
                      / static_cast<double>(dense_xq.size() - 1);
    }
    const auto dense_linear =
        interp::interp1_linear(reference_x, reference_y, dense_xq);
    const auto dense_cubic =
        interp::interp1_cubic(reference_x, reference_y, dense_xq);
    const auto dense_pchip =
        interp::interp1_pchip(reference_x, reference_y, dense_xq);
    const auto dense_akima =
        interp::interp1_akima(reference_x, reference_y, dense_xq);
    const auto dense_polynomial =
        interp::interp1_polynomial(reference_x, reference_y, dense_xq);
    const auto dense_nearest =
        interp::interp1_near(reference_x, reference_y, dense_xq);

    std::ofstream data_file(compare_dir / "interp_data.txt");
    data_file << std::setprecision(17);
    for (size_t i = 0; i < reference_x.size(); ++i) {
        data_file << reference_x[i] << " " << reference_y[i] << "\n";
    }

    std::ofstream interp_file(compare_dir / "interp_results.txt");
    interp_file << std::setprecision(17);
    for (size_t i = 0; i < dense_xq.size(); ++i) {
        interp_file << dense_xq[i] << " " << dense_linear[i] << " "
                    << dense_cubic[i] << " " << dense_pchip[i] << " "
                    << dense_akima[i] << " " << dense_polynomial[i] << "\n";
    }

    std::ofstream nearest_file(compare_dir / "interp_nearest.txt");
    nearest_file << std::setprecision(17);
    for (size_t i = 0; i < dense_xq.size(); ++i) {
        nearest_file << dense_xq[i] << " " << dense_nearest[i] << "\n";
    }

    return 0;
}

int test_interp_state() {
    // Default-constructed interpolators reject evaluation before set_data
    interp::CubicSpline cubic;
    EXPECT_EQ(cubic.size(), 0);
    bool threw_eval = false;
    try {
        (void)cubic(0.5);
    } catch (const std::runtime_error &) {
        threw_eval = true;
    }
    EXPECT_TRUE(threw_eval);

    bool threw_range = false;
    try {
        (void)cubic.range();
    } catch (const std::runtime_error &) {
        threw_range = true;
    }
    EXPECT_TRUE(threw_range);

    bool threw_in_range = false;
    try {
        (void)cubic.in_range(0.5);
    } catch (const std::runtime_error &) {
        threw_in_range = true;
    }
    EXPECT_TRUE(threw_in_range);

    interp::AkimaSpline akima;
    bool threw_akima = false;
    try {
        (void)akima(0.5);
    } catch (const std::runtime_error &) {
        threw_akima = true;
    }
    EXPECT_TRUE(threw_akima);

    // Invalid input data is rejected
    const auto throws_invalid = [](auto &&func) {
        try {
            func();
        } catch (const std::invalid_argument &) {
            return true;
        } catch (...) {
            return false;
        }
        return false;
    };

    interp::CubicSpline spline;
    EXPECT_TRUE(throws_invalid([&] {
        spline.set_data(std::vector<double>{0.0, 1.0, 2.0},
                        std::vector<double>{1.0, 2.0});
    }));
    EXPECT_TRUE(throws_invalid([&] {
        spline.set_data(std::vector<double>{0.0}, std::vector<double>{1.0});
    }));
    EXPECT_TRUE(throws_invalid([&] {
        spline.set_data(std::vector<double>{0.0, 1.0, 1.0},
                        std::vector<double>{1.0, 2.0, 3.0});
    }));

    // Not-a-knot is the default; Natural differs on nonlinear data
    const std::vector<double> x{0.0, 1.0, 2.0, 3.0, 4.0};
    const std::vector<double> y{0.0, 1.0, 0.0, 1.0, 0.0};
    auto default_spline = interp::CubicSpline::from_data(x, y);
    auto natural_spline = interp::CubicSpline::from_data(
        x, y, interp::CubicSpline::BoundaryCondition::Natural);
    const double default_value = default_spline(0.5);
    const double natural_value = natural_spline(0.5);
    EXPECT_TRUE(std::abs(default_value - natural_value) > 1e-6);

    // Both boundary conditions interpolate the data points exactly
    for (size_t i = 0; i < x.size(); ++i) {
        EXPECT_NEAR(default_spline(x[i]), y[i], 1e-12);
        EXPECT_NEAR(natural_spline(x[i]), y[i], 1e-12);
    }

    return 0;
}

int main() {
    test_interp();
    test_interp_state();
    return msl::test::summary();
}
