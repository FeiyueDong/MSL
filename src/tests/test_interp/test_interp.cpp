#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "interp.hpp"

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

static void expect_vector_near(const std::vector<double> &actual,
                               const std::vector<double> &expected,
                               double eps) {
    EXPECT_EQ(actual.size(), expected.size());
    for (size_t i = 0; i < actual.size() && i < expected.size(); ++i) {
        EXPECT_NEAR(actual[i], expected[i], eps);
    }
}

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

    auto compare_dir =
        project_root() / "test_result" / "interp" / "matlab_compare";
    std::filesystem::create_directories(compare_dir);

    auto linear = interp::interp1_linear(x, y, xq);
    auto cubic = interp::interp1_cubic(x, y, xq);
    auto pchip = interp::interp1_pchip(x, y, xq);
    auto akima = interp::interp1_akima(x, y, xq);
    std::ofstream interp_file(compare_dir / "interp_results.txt");
    interp_file << std::setprecision(17);
    for (size_t i = 0; i < xq.size(); ++i) {
        interp_file << xq[i] << " " << linear[i] << " " << cubic[i] << " "
                    << pchip[i] << " " << akima[i] << "\n";
    }

    auto nearest = interp::interp1_near(x, y, near_xq);
    std::ofstream nearest_file(compare_dir / "interp_nearest.txt");
    nearest_file << std::setprecision(17);
    for (size_t i = 0; i < near_xq.size(); ++i) {
        nearest_file << near_xq[i] << " " << nearest[i] << "\n";
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
    std::cout << "Total checks: " << g_total << ", failures: " << g_failures
              << "\n";
    return g_failures == 0 ? 0 : 1;
}
