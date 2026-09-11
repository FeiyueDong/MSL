/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2026, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: test_utils.hpp
** -----
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Shared helpers and assertion macros for the MSL test executables.
*/

#ifndef MSL_TESTS_TEST_UTILS_HPP
#define MSL_TESTS_TEST_UTILS_HPP

#include <cmath>
#include <complex>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

namespace msl::test {

// Global check counters shared by the assertion macros.
inline int g_failures = 0;
inline int g_total = 0;

/**
 * @brief Locate the repository root from the current working directory.
 *
 * Works both when tests are launched from the repository root and when they
 * are launched from a build directory by `xmake run`.
 */
inline std::filesystem::path project_root() {
    auto path = std::filesystem::current_path();
    while (!path.empty()) {
        if (std::filesystem::exists(path / "src" / "msl")
            && std::filesystem::exists(path / "src" / "tests")) {
            return path;
        }
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

/**
 * @brief Create and return the MATLAB comparison directory for a module.
 *
 * The directory is `<repo_root>/test_result/<module>/matlab_compare`.
 */
inline std::filesystem::path result_dir(const std::string &module) {
    auto directory = project_root() / "test_result" / module / "matlab_compare";
    std::filesystem::create_directories(directory);
    return directory;
}

inline void expect_complex_near(const std::complex<double> &actual,
                                const std::complex<double> &expected,
                                double eps,
                                const char *file,
                                int line) {
    ++g_total;
    if (std::abs(actual.real() - expected.real()) > eps
        || std::abs(actual.imag() - expected.imag()) > eps) {
        ++g_failures;
        std::cerr << "[FAIL] " << file << ":" << line
                  << " - complex not near (got: " << actual << " vs "
                  << expected << ", eps=" << eps << ")\n";
    }
}

/**
 * @brief Print the final check summary and return a process exit code.
 */
inline int summary() {
    std::cout << "Total checks: " << g_total << ", failures: " << g_failures
              << "\n";
    return g_failures == 0 ? 0 : 1;
}

} // namespace msl::test

#define TEST_CASE(name)                                                        \
    std::cout << "=== Test: " << (name) << " ===\n";                           \
    do

#define EXPECT_TRUE(cond)                                                      \
    do {                                                                       \
        ++msl::test::g_total;                                                  \
        if (!(cond)) {                                                         \
            ++msl::test::g_failures;                                           \
            std::cerr << "[FAIL] " << __FILE__ << ":" << __LINE__ << " - "     \
                      << #cond << "\n";                                        \
        }                                                                      \
    } while (0)

#define EXPECT_EQ(a, b) EXPECT_TRUE((a) == (b))

#define EXPECT_NEAR(a, b, eps)                                                 \
    do {                                                                       \
        ++msl::test::g_total;                                                  \
        const double expect_near_actual = static_cast<double>(a);              \
        const double expect_near_expected = static_cast<double>(b);            \
        if (!(std::fabs(expect_near_actual - expect_near_expected)             \
              <= (eps))) {                                                     \
            ++msl::test::g_failures;                                           \
            std::cerr << "[FAIL] " << __FILE__ << ":" << __LINE__ << " - "     \
                      << #a << " ~= " << #b << " (got: " << expect_near_actual \
                      << " vs " << expect_near_expected << ", eps=" << (eps)   \
                      << ")\n";                                                \
        }                                                                      \
    } while (0)

#define EXPECT_CPLX_NEAR(a, b, eps)                                            \
    msl::test::expect_complex_near((a), (b), (eps), __FILE__, __LINE__)

#endif // MSL_TESTS_TEST_UTILS_HPP
