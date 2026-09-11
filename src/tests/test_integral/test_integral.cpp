#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "integral.hpp"
#include "matrix.hpp"

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

int test_integral() {
    const double pi = std::acos(-1.0);

    std::vector<double> y{0.0, 1.0, 4.0, 9.0};
    auto ct = integral::cumtrapz(y, 1.0);
    EXPECT_EQ(ct.size(), y.size());
    EXPECT_NEAR(ct[0], 0.0, 1e-12);
    EXPECT_NEAR(ct[1], 0.5, 1e-12);
    EXPECT_NEAR(ct[2], 3.0, 1e-12);
    EXPECT_NEAR(ct[3], 9.5, 1e-12);
    EXPECT_NEAR(integral::trapz(y, 1.0), 9.5, 1e-12);

    std::vector<double> x{0.0, 0.5, 2.0, 3.0};
    auto ct_nonuniform = integral::cumtrapz(x, y);
    EXPECT_NEAR(ct_nonuniform[1], 0.25, 1e-12);
    EXPECT_NEAR(ct_nonuniform[2], 4.0, 1e-12);
    EXPECT_NEAR(ct_nonuniform[3], 10.5, 1e-12);
    EXPECT_NEAR(integral::trapz(x, y), 10.5, 1e-12);

    std::vector<double> parabola{0.0, 1.0, 4.0};
    EXPECT_NEAR(integral::simpson(parabola, 1.0), 8.0 / 3.0, 1e-12);
    auto cs = integral::cumsimpson(parabola, 1.0);
    EXPECT_EQ(cs.size(), parabola.size());
    EXPECT_NEAR(cs[0], 0.0, 1e-12);
    EXPECT_NEAR(cs[1], 1.0 / 3.0, 1e-12);
    EXPECT_NEAR(cs[2], 8.0 / 3.0, 1e-12);

    std::vector<double> x_parabola{0.0, 0.5, 2.0};
    std::vector<double> y_parabola{0.0, 0.25, 4.0};
    EXPECT_NEAR(integral::simpson(x_parabola, y_parabola), 8.0 / 3.0, 1e-12);
    auto cs_nonuniform = integral::cumsimpson(x_parabola, y_parabola);
    EXPECT_EQ(cs_nonuniform.size(), y_parabola.size());
    EXPECT_NEAR(cs_nonuniform[1], 1.0 / 24.0, 1e-12);
    EXPECT_NEAR(cs_nonuniform[2], 8.0 / 3.0, 1e-12);

    // Cumulative Simpson is exact at every point for a quadratic
    std::vector<double> quadratic_y{0.0, 1.0, 4.0, 9.0, 16.0};
    auto cumulative = integral::cumsimpson(quadratic_y, 1.0);
    for (size_t i = 0; i < quadratic_y.size(); ++i) {
        const double exact = static_cast<double>(i) * static_cast<double>(i)
                             * static_cast<double>(i) / 3.0;
        EXPECT_NEAR(cumulative[i], exact, 1e-12);
    }

    std::vector<double> quadratic_x{0.0, 0.5, 1.5, 3.0, 4.0};
    std::vector<double> quadratic_y_nonuniform{0.0, 0.25, 2.25, 9.0, 16.0};
    auto cumulative_nonuniform =
        integral::cumsimpson(quadratic_x, quadratic_y_nonuniform);
    for (size_t i = 0; i < quadratic_x.size(); ++i) {
        const double exact =
            quadratic_x[i] * quadratic_x[i] * quadratic_x[i] / 3.0;
        EXPECT_NEAR(cumulative_nonuniform[i], exact, 1e-12);
    }

    // Uniform and non-uniform formulations agree on a uniform grid
    std::vector<double> uniform_axis{0.0, 1.0, 2.0, 3.0, 4.0};
    auto cumulative_axis = integral::cumsimpson(uniform_axis, quadratic_y);
    for (size_t i = 0; i < cumulative.size(); ++i) {
        EXPECT_NEAR(cumulative_axis[i], cumulative[i], 1e-12);
    }

    // Even number of points: the final interval falls back to trapezoidal
    std::vector<double> even_y{0.0, 1.0, 4.0, 9.0, 16.0, 25.0};
    auto even_cumulative = integral::cumsimpson(even_y, 1.0);
    for (size_t i = 0; i + 1 < even_y.size(); ++i) {
        const double exact = static_cast<double>(i) * static_cast<double>(i)
                             * static_cast<double>(i) / 3.0;
        EXPECT_NEAR(even_cumulative[i], exact, 1e-12);
    }
    EXPECT_NEAR(
        even_cumulative[5], even_cumulative[4] + 0.5 * (16.0 + 25.0), 1e-12);

    EXPECT_NEAR(
        integral::trapz([](double x) { return x; }, 0.0, 1.0, 16), 0.5, 1e-12);
    EXPECT_NEAR(integral::simpson([](double x) { return x * x; }, 0.0, 2.0, 16),
                8.0 / 3.0,
                1e-12);

    std::vector<double> romberg_x{0.0, 1.0, 2.0, 3.0, 4.0};
    std::vector<double> romberg_y{0.0, 1.0, 4.0, 9.0, 16.0};
    EXPECT_NEAR(integral::romberg(romberg_y, 1.0), 64.0 / 3.0, 1e-12);
    EXPECT_NEAR(integral::romberg(romberg_x, romberg_y), 64.0 / 3.0, 1e-12);
    EXPECT_NEAR(
        integral::romberg([](double x) { return std::sin(x); }, 0.0, pi, 1e-12),
        2.0,
        1e-10);
    EXPECT_NEAR(integral::adaptive_simpson(
                    [](double x) { return std::sin(x); }, 0.0, pi, 1e-12),
                2.0,
                1e-10);
    EXPECT_NEAR(integral::quad([](double x) { return x * x; }, 0.0, 2.0),
                8.0 / 3.0,
                1e-8);

    matrix::matrixd mat(4, 2);
    for (size_t i = 0; i < mat.rows(); ++i) {
        mat(i, 0) = static_cast<double>(i);
        mat(i, 1) = 2.0 * static_cast<double>(i);
    }
    auto col_trapz = integral::trapz(mat, 1.0);
    EXPECT_EQ(col_trapz.size(), 2);
    EXPECT_NEAR(col_trapz[0], 4.5, 1e-12);
    EXPECT_NEAR(col_trapz[1], 9.0, 1e-12);

    auto col_trapz_nonuniform = integral::trapz(x, mat);
    EXPECT_EQ(col_trapz_nonuniform.size(), 2);
    EXPECT_NEAR(col_trapz_nonuniform[0], 5.0, 1e-12);
    EXPECT_NEAR(col_trapz_nonuniform[1], 10.0, 1e-12);

    auto mat_cum = integral::cumtrapz(mat, 1.0);
    EXPECT_EQ(mat_cum.rows(), 4);
    EXPECT_EQ(mat_cum.cols(), 2);
    EXPECT_NEAR(mat_cum(3, 0), 4.5, 1e-12);
    EXPECT_NEAR(mat_cum(3, 1), 9.0, 1e-12);

    bool threw = false;
    try {
        (void)integral::trapz(std::vector<double>{1.0}, 1.0);
    } catch (const std::invalid_argument &) {
        threw = true;
    }
    EXPECT_TRUE(threw);

    bool threw_cumsimpson = false;
    try {
        (void)integral::cumsimpson(std::vector<double>{0.0, 1.0, 1.0},
                                   std::vector<double>{0.0, 1.0, 4.0});
    } catch (const std::invalid_argument &) {
        threw_cumsimpson = true;
    }
    EXPECT_TRUE(threw_cumsimpson);

    auto compare_dir =
        project_root() / "test_result" / "integral" / "matlab_compare";
    std::filesystem::create_directories(compare_dir);

    std::ofstream summary_file(compare_dir / "integral_summary.txt");
    summary_file << std::setprecision(17);
    summary_file << integral::trapz(y, 1.0) << "\n";
    summary_file << integral::trapz(x, y) << "\n";
    summary_file << integral::simpson(parabola, 1.0) << "\n";
    summary_file << integral::simpson(x_parabola, y_parabola) << "\n";
    summary_file << integral::romberg(romberg_y, 1.0) << "\n";
    summary_file << integral::adaptive_simpson(
        [](double x) { return std::sin(x); }, 0.0, pi, 1e-12)
                 << "\n";

    std::ofstream cumulative_file(compare_dir / "integral_cumulative.txt");
    cumulative_file << std::setprecision(17);
    for (size_t i = 0; i < y.size(); ++i) {
        cumulative_file << x[i] << " " << y[i] << " " << ct[i] << " "
                        << ct_nonuniform[i] << "\n";
    }

    // ---- Cumulative Simpson against analytic primitives ----
    const std::vector<double> cs_x{0.0, 1.0, 2.0, 3.0, 4.0};
    const std::vector<double> cs_y{0.0, 1.0, 4.0, 9.0, 16.0};
    const auto cs_uniform = integral::cumsimpson(cs_y, 1.0);
    std::ofstream cs_uniform_file(compare_dir
                                  / "integral_cumsimpson_uniform.txt");
    cs_uniform_file << std::setprecision(17);
    for (size_t i = 0; i < cs_x.size(); ++i) {
        cs_uniform_file << cs_x[i] << " " << cs_y[i] << " " << cs_uniform[i]
                        << "\n";
    }

    const std::vector<double> cs_x_nonuniform{0.0, 0.5, 1.5, 3.0, 4.0};
    const std::vector<double> cs_y_nonuniform{0.0, 0.25, 2.25, 9.0, 16.0};
    const auto cs_nonuniform_values =
        integral::cumsimpson(cs_x_nonuniform, cs_y_nonuniform);
    std::ofstream cs_nonuniform_file(compare_dir
                                     / "integral_cumsimpson_nonuniform.txt");
    cs_nonuniform_file << std::setprecision(17);
    for (size_t i = 0; i < cs_x_nonuniform.size(); ++i) {
        cs_nonuniform_file << cs_x_nonuniform[i] << " " << cs_y_nonuniform[i]
                           << " " << cs_nonuniform_values[i] << "\n";
    }

    const std::vector<double> cs_even_y{0.0, 1.0, 4.0, 9.0, 16.0, 25.0};
    const auto cs_even = integral::cumsimpson(cs_even_y, 1.0);
    std::ofstream cs_even_file(compare_dir / "integral_cumsimpson_even.txt");
    cs_even_file << std::setprecision(17);
    for (size_t i = 0; i < cs_even_y.size(); ++i) {
        cs_even_file << static_cast<double>(i) << " " << cs_even_y[i] << " "
                     << cs_even[i] << "\n";
    }

    // ---- Matrix integration against analytic columns ----
    matrix::matrixd integrand(5, 2);
    for (size_t i = 0; i < integrand.rows(); ++i) {
        const double xi = static_cast<double>(i);
        integrand(i, 0) = xi * xi;
        integrand(i, 1) = xi * xi * xi;
    }
    const auto simpson_columns = integral::simpson(integrand, 1.0);
    const auto trapz_columns = integral::trapz(integrand, 1.0);
    std::ofstream integrand_file(compare_dir / "integral_matrix.txt");
    integrand_file << std::setprecision(17);
    for (size_t i = 0; i < integrand.rows(); ++i) {
        integrand_file << integrand(i, 0) << " " << integrand(i, 1) << "\n";
    }
    std::ofstream matrix_results_file(compare_dir
                                      / "integral_matrix_results.txt");
    matrix_results_file << std::setprecision(17);
    matrix_results_file << simpson_columns[0] << " " << simpson_columns[1]
                        << "\n";
    matrix_results_file << trapz_columns[0] << " " << trapz_columns[1] << "\n";

    std::cout << "Total checks: " << g_total << ", failures: " << g_failures
              << "\n";
    return g_failures == 0 ? 0 : 1;
}

int main() { return test_integral(); }
