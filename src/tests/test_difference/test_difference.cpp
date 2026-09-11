#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "../test_utils.hpp"
#include "difference.hpp"
#include "matrix.hpp"


using namespace msl;

static void write_matrix(std::ofstream &file, const matrix::matrixd &mat) {
    file << std::setprecision(17);
    for (size_t i = 0; i < mat.rows(); ++i) {
        for (size_t j = 0; j < mat.cols(); ++j) {
            file << mat(i, j) << (j + 1 == mat.cols() ? '\n' : ' ');
        }
    }
}

int test_difference() {
    std::vector<double> y{1.0, 4.0, 9.0, 16.0};
    auto d = difference::diff(y);
    EXPECT_EQ(d.size(), 3);
    EXPECT_NEAR(d[0], 3.0, 1e-12);
    EXPECT_NEAR(d[1], 5.0, 1e-12);
    EXPECT_NEAR(d[2], 7.0, 1e-12);

    auto fg = difference::forward_gradient(y, 0.5);
    EXPECT_EQ(fg.size(), 3);
    EXPECT_NEAR(fg[0], 6.0, 1e-12);
    EXPECT_NEAR(fg[1], 10.0, 1e-12);
    EXPECT_NEAR(fg[2], 14.0, 1e-12);

    auto cg = difference::central_gradient(y, 1.0);
    EXPECT_EQ(cg.size(), y.size());
    EXPECT_NEAR(cg[0], 3.0, 1e-12);
    EXPECT_NEAR(cg[1], 4.0, 1e-12);
    EXPECT_NEAR(cg[2], 6.0, 1e-12);
    EXPECT_NEAR(cg[3], 7.0, 1e-12);

    matrix::matrixd mat(3, 3);
    for (size_t i = 0; i < mat.rows(); ++i) {
        for (size_t j = 0; j < mat.cols(); ++j) {
            mat(i, j) = 10.0 * static_cast<double>(i) + static_cast<double>(j);
        }
    }

    auto row_diff = difference::diff(mat, 0);
    EXPECT_EQ(row_diff.rows(), 2);
    EXPECT_EQ(row_diff.cols(), 3);
    for (size_t i = 0; i < row_diff.rows(); ++i) {
        for (size_t j = 0; j < row_diff.cols(); ++j) {
            EXPECT_NEAR(row_diff(i, j), 10.0, 1e-12);
        }
    }

    auto col_diff = difference::diff(mat, 1);
    EXPECT_EQ(col_diff.rows(), 3);
    EXPECT_EQ(col_diff.cols(), 2);
    for (size_t i = 0; i < col_diff.rows(); ++i) {
        for (size_t j = 0; j < col_diff.cols(); ++j) {
            EXPECT_NEAR(col_diff(i, j), 1.0, 1e-12);
        }
    }

    auto lap = difference::laplacian(mat);
    EXPECT_EQ(lap.rows(), 3);
    EXPECT_EQ(lap.cols(), 3);
    EXPECT_NEAR(lap(1, 1), 0.0, 1e-12);

    bool threw = false;
    try {
        (void)difference::diff(std::vector<double>{1.0});
    } catch (const std::invalid_argument &) {
        threw = true;
    }
    EXPECT_TRUE(threw);

    // Non-uniform forward gradient is exact for f(x) = x^2
    const std::vector<double> nu_x{0.0, 1.0, 3.0, 4.0};
    const std::vector<double> nu_y{0.0, 1.0, 9.0, 16.0};
    const auto nu_forward = difference::forward_gradient(nu_y, nu_x);
    EXPECT_NEAR(nu_forward[0], 1.0, 1e-12);
    EXPECT_NEAR(nu_forward[1], 4.0, 1e-12);
    EXPECT_NEAR(nu_forward[2], 7.0, 1e-12);

    const auto nu_central = difference::central_gradient(nu_x, nu_y);
    EXPECT_NEAR(nu_central[1], 2.0 * nu_x[1], 1e-12);
    EXPECT_NEAR(nu_central[2], 2.0 * nu_x[2], 1e-12);

    // Documented example: x = [0, 1, 3], f = x^2, interior derivative = 2
    const std::vector<double> nu_small_x{0.0, 1.0, 3.0};
    const std::vector<double> nu_small_y{0.0, 1.0, 9.0};
    const auto nu_small_central =
        difference::central_gradient(nu_small_x, nu_small_y);
    EXPECT_NEAR(nu_small_central[1], 2.0, 1e-12);

    // Non-uniform matrix central gradient matches the vector formulation
    matrix::matrixd nu_rows(4, 3);
    for (size_t i = 0; i < nu_rows.rows(); ++i) {
        for (size_t j = 0; j < nu_rows.cols(); ++j) {
            nu_rows(i, j) = nu_x[i] * nu_x[i]; // f = x^2, gradient along rows
        }
    }
    const auto nu_rows_grad = difference::central_gradient(nu_rows, nu_x, 0);
    for (size_t i = 1; i + 1 < nu_rows.rows(); ++i) {
        for (size_t j = 0; j < nu_rows.cols(); ++j) {
            EXPECT_NEAR(nu_rows_grad(i, j), 2.0 * nu_x[i], 1e-12);
        }
    }

    matrix::matrixd nu_cols(3, 4);
    for (size_t i = 0; i < nu_cols.rows(); ++i) {
        for (size_t j = 0; j < nu_cols.cols(); ++j) {
            nu_cols(i, j) = nu_x[j] * nu_x[j]; // f = x^2, along columns
        }
    }
    const auto nu_cols_grad = difference::central_gradient(nu_cols, nu_x, 1);
    for (size_t i = 0; i < nu_cols.rows(); ++i) {
        for (size_t j = 1; j + 1 < nu_cols.cols(); ++j) {
            EXPECT_NEAR(nu_cols_grad(i, j), 2.0 * nu_x[j], 1e-12);
        }
    }

    // Second derivative is exact for a quadratic at every point
    const std::vector<double> quad_y{0.0, 1.0, 4.0, 9.0, 16.0};
    const auto quad_second = difference::central_gradient2(quad_y, 1.0);
    for (const double value : quad_second) {
        EXPECT_NEAR(value, 2.0, 1e-12);
    }

    // 2D analytic fields on f(i, j) = i^2 + j^2
    matrix::matrixd field(5, 5);
    for (size_t i = 0; i < field.rows(); ++i) {
        for (size_t j = 0; j < field.cols(); ++j) {
            field(i, j) = static_cast<double>(i * i + j * j);
        }
    }
    const auto gradient_2d = difference::central_gradient2d(field);
    for (size_t i = 1; i + 1 < field.rows(); ++i) {
        for (size_t j = 1; j + 1 < field.cols(); ++j) {
            EXPECT_NEAR(
                gradient_2d.first(i, j), 2.0 * static_cast<double>(j), 1e-12);
            EXPECT_NEAR(
                gradient_2d.second(i, j), 2.0 * static_cast<double>(i), 1e-12);
        }
    }

    const auto laplacian_field = difference::laplacian(field);
    for (size_t i = 1; i + 1 < field.rows(); ++i) {
        for (size_t j = 1; j + 1 < field.cols(); ++j) {
            EXPECT_NEAR(laplacian_field(i, j), 4.0, 1e-12);
        }
    }

    matrix::matrixd vector_x(4, 4);
    matrix::matrixd vector_y(4, 4);
    for (size_t i = 0; i < 4; ++i) {
        for (size_t j = 0; j < 4; ++j) {
            vector_x(i, j) = static_cast<double>(j); // x
            vector_y(i, j) = static_cast<double>(i); // y
        }
    }
    const auto divergence_field = difference::divergence(vector_x, vector_y);
    for (size_t i = 0; i < 4; ++i) {
        for (size_t j = 0; j < 4; ++j) {
            EXPECT_NEAR(divergence_field(i, j), 2.0, 1e-12);
        }
    }

    matrix::matrixd curl_x(4, 4);
    matrix::matrixd curl_y(4, 4);
    for (size_t i = 0; i < 4; ++i) {
        for (size_t j = 0; j < 4; ++j) {
            curl_x(i, j) = -static_cast<double>(i); // -y
            curl_y(i, j) = static_cast<double>(j);  // x
        }
    }
    const auto curl_field = difference::curl(curl_x, curl_y);
    for (size_t i = 0; i < 4; ++i) {
        for (size_t j = 0; j < 4; ++j) {
            EXPECT_NEAR(curl_field(i, j), 2.0, 1e-12);
        }
    }

    // Savitzky-Golay derivative is exact for a quadratic (order 2 fit)
    std::vector<double> sg_quadratic(11);
    for (size_t i = 0; i < sg_quadratic.size(); ++i) {
        const double xi = static_cast<double>(i);
        sg_quadratic[i] = xi * xi;
    }
    const auto sg_derivative =
        difference::savgol_gradient(sg_quadratic, 5, 2, 1.0);
    for (size_t i = 0; i < sg_quadratic.size(); ++i) {
        EXPECT_NEAR(sg_derivative[i], 2.0 * static_cast<double>(i), 1e-10);
    }

    // A cubic fit reproduces a cubic exactly
    std::vector<double> sg_cubic(11);
    for (size_t i = 0; i < sg_cubic.size(); ++i) {
        const double xi = static_cast<double>(i);
        sg_cubic[i] = xi * xi * xi;
    }
    const auto sg_cubic_derivative =
        difference::savgol_gradient(sg_cubic, 7, 3, 1.0);
    for (size_t i = 0; i < sg_cubic.size(); ++i) {
        const double xi = static_cast<double>(i);
        EXPECT_NEAR(sg_cubic_derivative[i], 3.0 * xi * xi, 1e-10);
    }

    bool threw_savgol = false;
    try {
        (void)difference::savgol_gradient(
            std::vector<double>{1.0, 2.0, 3.0, 4.0, 5.0}, 4, 2, 1.0);
    } catch (const std::invalid_argument &) {
        threw_savgol = true;
    }
    EXPECT_TRUE(threw_savgol);

    // Invalid inputs must throw instead of underflowing or dividing by zero
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

    EXPECT_TRUE(
        throws_invalid([] { (void)difference::diff(std::vector<double>{}); }));
    EXPECT_TRUE(throws_invalid(
        [] { (void)difference::forward_gradient(std::vector<double>{}); }));
    EXPECT_TRUE(throws_invalid([] {
        (void)difference::forward_gradient(std::vector<double>{1.0, 2.0}, 0.0);
    }));
    EXPECT_TRUE(throws_invalid([] {
        (void)difference::central_gradient(std::vector<double>{1.0, 2.0}, 0.0);
    }));
    EXPECT_TRUE(throws_invalid([] {
        (void)difference::central_gradient2(std::vector<double>{0.0, 1.0, 4.0},
                                            0.0);
    }));

    // Non-uniform coordinates must be strictly increasing
    EXPECT_TRUE(throws_invalid([] {
        (void)difference::forward_gradient(std::vector<double>{1.0, 2.0, 4.0},
                                           std::vector<double>{0.0, 1.0, 1.0});
    }));
    EXPECT_TRUE(throws_invalid([] {
        (void)difference::forward_gradient(std::vector<double>{1.0, 2.0, 4.0},
                                           std::vector<double>{0.0, 3.0, 1.0});
    }));
    EXPECT_TRUE(throws_invalid([] {
        (void)difference::central_gradient(std::vector<double>{0.0, 1.0, 1.0},
                                           std::vector<double>{0.0, 1.0, 4.0});
    }));
    EXPECT_TRUE(throws_invalid([] {
        (void)difference::central_gradient(std::vector<double>{0.0, 1.0, 0.5},
                                           std::vector<double>{0.0, 1.0, 6.0});
    }));

    const matrix::matrixd validation_matrix(3, 3, 1.0);
    EXPECT_TRUE(throws_invalid([&] {
        (void)difference::forward_gradient(
            validation_matrix, std::vector<double>{0.0, 1.0, 1.0}, 0);
    }));
    EXPECT_TRUE(throws_invalid([&] {
        (void)difference::central_gradient(
            validation_matrix, std::vector<double>{0.0, 1.0, 0.5}, 1);
    }));
    EXPECT_TRUE(throws_invalid([&] {
        (void)difference::central_gradient(validation_matrix, 0.0, 0);
    }));
    EXPECT_TRUE(throws_invalid(
        [&] { (void)difference::laplacian(validation_matrix, 1.0, 0.0); }));

    auto compare_dir = msl::test::result_dir("difference");
    std::ofstream vec_file(compare_dir / "difference_vector.txt");
    vec_file << std::setprecision(17);
    for (size_t i = 0; i < y.size(); ++i) {
        vec_file << y[i] << " " << (i < d.size() ? d[i] : 0.0) << " " << cg[i]
                 << "\n";
    }

    std::ofstream mat_file(compare_dir / "difference_matrix.txt");
    mat_file << std::setprecision(17);
    for (size_t i = 0; i < mat.rows(); ++i) {
        for (size_t j = 0; j < mat.cols(); ++j) {
            mat_file << mat(i, j) << (j + 1 == mat.cols() ? '\n' : ' ');
        }
    }

    std::ofstream row_file(compare_dir / "difference_row_diff.txt");
    row_file << std::setprecision(17);
    for (size_t i = 0; i < row_diff.rows(); ++i) {
        for (size_t j = 0; j < row_diff.cols(); ++j) {
            row_file << row_diff(i, j)
                     << (j + 1 == row_diff.cols() ? '\n' : ' ');
        }
    }

    std::ofstream col_file(compare_dir / "difference_col_diff.txt");
    col_file << std::setprecision(17);
    for (size_t i = 0; i < col_diff.rows(); ++i) {
        for (size_t j = 0; j < col_diff.cols(); ++j) {
            col_file << col_diff(i, j)
                     << (j + 1 == col_diff.cols() ? '\n' : ' ');
        }
    }

    // ---- Non-uniform gradient comparison data ----
    std::ofstream nonuniform_input_file(compare_dir
                                        / "difference_nonuniform_input.txt");
    nonuniform_input_file << std::setprecision(17);
    for (size_t i = 0; i < nu_x.size(); ++i) {
        nonuniform_input_file << nu_x[i] << " " << nu_y[i] << "\n";
    }
    std::ofstream forward_nonuniform_file(
        compare_dir / "difference_forward_nonuniform.txt");
    forward_nonuniform_file << std::setprecision(17);
    for (const double value : nu_forward) {
        forward_nonuniform_file << value << "\n";
    }
    std::ofstream central_nonuniform_file(
        compare_dir / "difference_central_nonuniform.txt");
    central_nonuniform_file << std::setprecision(17);
    for (const double value : nu_central) {
        central_nonuniform_file << value << "\n";
    }

    // ---- Second derivative comparison data ----
    std::ofstream second_derivative_file(compare_dir
                                         / "difference_second_derivative.txt");
    second_derivative_file << std::setprecision(17);
    for (size_t i = 0; i < quad_y.size(); ++i) {
        second_derivative_file << static_cast<double>(i) << " " << quad_y[i]
                               << " " << quad_second[i] << "\n";
    }

    // ---- 2D gradient / Laplacian comparison data ----
    std::ofstream field_file(compare_dir / "difference_field.txt");
    write_matrix(field_file, field);
    std::ofstream grad_x_file(compare_dir / "difference_grad_x.txt");
    write_matrix(grad_x_file, gradient_2d.first);
    std::ofstream grad_y_file(compare_dir / "difference_grad_y.txt");
    write_matrix(grad_y_file, gradient_2d.second);
    std::ofstream laplacian_file(compare_dir / "difference_laplacian.txt");
    write_matrix(laplacian_file, laplacian_field);

    // ---- Divergence / curl comparison data ----
    std::ofstream vector_x_file(compare_dir / "difference_vector_x.txt");
    write_matrix(vector_x_file, vector_x);
    std::ofstream vector_y_file(compare_dir / "difference_vector_y.txt");
    write_matrix(vector_y_file, vector_y);
    std::ofstream divergence_file(compare_dir / "difference_divergence.txt");
    write_matrix(divergence_file, divergence_field);
    std::ofstream curl_file(compare_dir / "difference_curl.txt");
    write_matrix(curl_file, curl_field);

    // ---- Savitzky-Golay comparison data ----
    std::vector<double> savgol_input(15);
    for (size_t i = 0; i < savgol_input.size(); ++i) {
        const double xi = static_cast<double>(i);
        savgol_input[i] = 0.5 * xi * xi * xi - 2.0 * xi * xi + 0.3 * xi;
    }
    const auto savgol_result =
        difference::savgol_gradient(savgol_input, 5, 2, 1.0);
    std::ofstream savgol_file(compare_dir / "difference_savgol.txt");
    savgol_file << std::setprecision(17);
    for (size_t i = 0; i < savgol_input.size(); ++i) {
        savgol_file << static_cast<double>(i) << " " << savgol_input[i] << " "
                    << savgol_result[i] << "\n";
    }

    return msl::test::summary();
}

int main() { return test_difference(); }
