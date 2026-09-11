#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "equation.hpp"
#include "test_utils.hpp"

using namespace msl;

int test_equation() {
    auto square_minus_two = [](double x) { return x * x - 2.0; };
    auto square_minus_two_deriv = [](double x) { return 2.0 * x; };
    double sqrt2 = std::sqrt(2.0);

    auto bisection_result =
        equation::bisection(square_minus_two, 0.0, 2.0, 1e-12);
    EXPECT_TRUE(bisection_result);
    EXPECT_EQ(bisection_result.status, equation::root_status::converged);
    EXPECT_NEAR(bisection_result.root, sqrt2, 1e-10);
    EXPECT_NEAR(bisection_result.residual, 0.0, 1e-10);
    EXPECT_NEAR(bisection_result.value(), sqrt2, 1e-10);
    EXPECT_NEAR(static_cast<double>(bisection_result), sqrt2, 1e-10);

    double bisection_root =
        equation::bisection_root(square_minus_two, 0.0, 2.0, 1e-12);
    EXPECT_NEAR(bisection_root, sqrt2, 1e-10);

    equation::root_options options;
    options.tol = 1e-12;
    options.max_iter = 100;
    options.derivative_tol = 1e-14;

    auto bisection_options =
        equation::bisection(square_minus_two, 0.0, 2.0, options);
    EXPECT_TRUE(bisection_options.converged());
    EXPECT_NEAR(bisection_options.root, sqrt2, 1e-10);
    EXPECT_NEAR(equation::bisection_root(square_minus_two, 0.0, 2.0, options),
                sqrt2,
                1e-10);

    auto brent_result = equation::brent(square_minus_two, 0.0, 2.0, 1e-12);
    EXPECT_TRUE(brent_result.converged());
    EXPECT_NEAR(brent_result.root, sqrt2, 1e-12);
    EXPECT_NEAR(
        equation::brent_root(square_minus_two, 0.0, 2.0, 1e-12), sqrt2, 1e-12);
    EXPECT_NEAR(equation::brent(square_minus_two, 0.0, 2.0, options).root,
                sqrt2,
                1e-12);
    EXPECT_NEAR(equation::brent_root(square_minus_two, 0.0, 2.0, options),
                sqrt2,
                1e-12);

    auto newton_result =
        equation::newton(square_minus_two, square_minus_two_deriv, 1.0, 1e-12);
    EXPECT_TRUE(newton_result.converged());
    EXPECT_NEAR(newton_result.root, sqrt2, 1e-12);
    EXPECT_NEAR(equation::newton_root(
                    square_minus_two, square_minus_two_deriv, 1.0, 1e-12),
                sqrt2,
                1e-12);
    EXPECT_NEAR(
        equation::newton(square_minus_two, square_minus_two_deriv, 1.0, options)
            .root,
        sqrt2,
        1e-12);
    EXPECT_NEAR(equation::newton_root(
                    square_minus_two, square_minus_two_deriv, 1.0, options),
                sqrt2,
                1e-12);

    auto secant_result = equation::secant(square_minus_two, 1.0, 2.0, 1e-12);
    EXPECT_TRUE(secant_result.converged());
    EXPECT_NEAR(secant_result.root, sqrt2, 1e-10);
    EXPECT_NEAR(
        equation::secant_root(square_minus_two, 1.0, 2.0, 1e-12), sqrt2, 1e-10);
    EXPECT_NEAR(equation::secant(square_minus_two, 1.0, 2.0, options).root,
                sqrt2,
                1e-10);
    EXPECT_NEAR(equation::secant_root(square_minus_two, 1.0, 2.0, options),
                sqrt2,
                1e-10);

    auto falsi_result =
        equation::regula_falsi(square_minus_two, 0.0, 2.0, 1e-12);
    EXPECT_TRUE(falsi_result.converged());
    EXPECT_NEAR(falsi_result.root, sqrt2, 1e-10);
    EXPECT_NEAR(equation::regula_falsi_root(square_minus_two, 0.0, 2.0, 1e-12),
                sqrt2,
                1e-10);
    EXPECT_NEAR(
        equation::regula_falsi(square_minus_two, 0.0, 2.0, options).root,
        sqrt2,
        1e-10);
    EXPECT_NEAR(
        equation::regula_falsi_root(square_minus_two, 0.0, 2.0, options),
        sqrt2,
        1e-10);

    auto fixed_point = equation::bisection(
        [](double x) { return std::cos(x) - x; }, 0.0, 1.0, 1e-12);
    EXPECT_TRUE(fixed_point.converged());
    EXPECT_NEAR(fixed_point.root, 0.7390851332151607, 1e-10);

    auto invalid_interval =
        equation::bisection([](double x) { return x * x + 1.0; }, -1.0, 1.0);
    EXPECT_EQ(invalid_interval.status, equation::root_status::invalid_interval);
    auto brent_invalid_interval =
        equation::brent([](double x) { return x * x + 1.0; }, -1.0, 1.0);
    EXPECT_EQ(brent_invalid_interval.status,
              equation::root_status::invalid_interval);
    bool value_threw = false;
    try {
        (void)invalid_interval.value();
    } catch (const std::runtime_error &) {
        value_threw = true;
    }
    EXPECT_TRUE(value_threw);

    auto zero_derivative = equation::newton(
        [](double x) { return x * x + 1.0; }, [](double) { return 0.0; }, 1.0);
    EXPECT_EQ(zero_derivative.status, equation::root_status::zero_derivative);

    auto max_iterations =
        equation::secant(square_minus_two, 1.0, 2.0, 1e-16, 1);
    EXPECT_EQ(max_iterations.status, equation::root_status::max_iterations);

    equation::root_options limited;
    limited.max_iter = 1;
    auto max_iterations_options =
        equation::secant(square_minus_two, 1.0, 2.0, limited);
    EXPECT_EQ(max_iterations_options.status,
              equation::root_status::max_iterations);
    auto brent_max_iterations =
        equation::brent(square_minus_two, 0.0, 2.0, limited);
    EXPECT_EQ(brent_max_iterations.status,
              equation::root_status::max_iterations);

    equation::root_options derivative_limited;
    derivative_limited.derivative_tol = 10.0;
    auto zero_derivative_options = equation::newton(
        square_minus_two, square_minus_two_deriv, 1.0, derivative_limited);
    EXPECT_EQ(zero_derivative_options.status,
              equation::root_status::zero_derivative);

    bool threw = false;
    try {
        (void)equation::bisection(square_minus_two, 2.0, 0.0);
    } catch (const std::invalid_argument &) {
        threw = true;
    }
    EXPECT_TRUE(threw);

    auto compare_dir = msl::test::result_dir("equation");
    std::ofstream roots_file(compare_dir / "equation_roots.txt");
    roots_file << std::setprecision(17);
    roots_file << bisection_result.root << " " << bisection_result.residual
               << " " << bisection_result.iterations << "\n";
    roots_file << brent_result.root << " " << brent_result.residual << " "
               << brent_result.iterations << "\n";
    roots_file << newton_result.root << " " << newton_result.residual << " "
               << newton_result.iterations << "\n";
    roots_file << secant_result.root << " " << secant_result.residual << " "
               << secant_result.iterations << "\n";
    roots_file << falsi_result.root << " " << falsi_result.residual << " "
               << falsi_result.iterations << "\n";
    roots_file << fixed_point.root << " " << fixed_point.residual << " "
               << fixed_point.iterations << "\n";

    return msl::test::summary();
}

int main() { return test_equation(); }
