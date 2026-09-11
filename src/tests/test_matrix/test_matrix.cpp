// tests/test_matrix_all.cpp
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

#include "matrix.hpp"

using namespace msl;

// ---------- 简易断言框架 ----------
static int g_failures = 0;
static int g_total = 0;

#define TEST_CASE(name)                                                        \
    std::cout << "=== Test: " << name << " ===\n";                             \
    do

#define EXPECT_TRUE(cond)                                                      \
    do {                                                                       \
        ++g_total;                                                             \
        if (!(cond)) {                                                         \
            ++g_failures;                                                      \
            std::cerr << "[FAIL] " << __FILE__ << ":" << __LINE__ << " - "     \
                      << #cond << "\n";                                        \
        }                                                                      \
    } while (0)

#define EXPECT_EQ(a, b)                                                        \
    do {                                                                       \
        ++g_total;                                                             \
        if (!((a) == (b))) {                                                   \
            ++g_failures;                                                      \
            std::cerr << "[FAIL] " << __FILE__ << ":" << __LINE__ << " - "     \
                      << #a << " == " << #b << " (got: " << (a) << " vs "      \
                      << (b) << ")\n";                                         \
        }                                                                      \
    } while (0)

#define EXPECT_NEAR(a, b, eps)                                                 \
    do {                                                                       \
        ++g_total;                                                             \
        if (!(std::fabs((double)(a) - (double)(b)) <= (eps))) {                \
            ++g_failures;                                                      \
            std::cerr << "[FAIL] " << __FILE__ << ":" << __LINE__ << " - "     \
                      << #a << " ~= " << #b << " (got: " << (a) << " vs "      \
                      << (b) << ", eps=" << (eps) << ")\n";                    \
        }                                                                      \
    } while (0)

// For complex near
static void expect_complex_near(const std::complex<double> &a,
                                const std::complex<double> &b,
                                double eps,
                                const char *file,
                                int line) {
    ++g_total;
    if (std::abs(a.real() - b.real()) > eps
        || std::abs(a.imag() - b.imag()) > eps) {
        ++g_failures;
        std::cerr << "[FAIL] " << file << ":" << line
                  << " - complex not near (got: " << a << " vs " << b
                  << ", eps=" << eps << ")\n";
    }
}
#define EXPECT_CPLX_NEAR(a, b, eps)                                            \
    expect_complex_near((a), (b), (eps), __FILE__, __LINE__)

// ---------- 辅助函数 ----------
static bool approx_equal_double_vector(const std::vector<double> &a,
                                       const std::vector<double> &b,
                                       double eps = 1e-12) {
    if (a.size() != b.size())
        return false;
    for (size_t i = 0; i < a.size(); ++i)
        if (std::fabs(a[i] - b[i]) > eps)
            return false;
    return true;
}

static bool matrix_near(const matrix::real_matrix_base &actual,
                        const matrix::real_matrix_base &expected,
                        double tolerance) {
    if (actual.rows() != expected.rows() || actual.cols() != expected.cols()) {
        return false;
    }
    for (size_t j = 0; j < actual.cols(); ++j) {
        for (size_t i = 0; i < actual.rows(); ++i) {
            if (std::abs(actual(i, j) - expected(i, j)) > tolerance) {
                return false;
            }
        }
    }
    return true;
}

static bool columns_are_orthonormal(const matrix::real_matrix_base &vectors,
                                    double tolerance) {
    for (size_t i = 0; i < vectors.cols(); ++i) {
        for (size_t j = 0; j < vectors.cols(); ++j) {
            double dot = 0.0;
            for (size_t row = 0; row < vectors.rows(); ++row) {
                dot += vectors(row, i) * vectors(row, j);
            }
            const double expected = i == j ? 1.0 : 0.0;
            if (std::abs(dot - expected) > tolerance) {
                return false;
            }
        }
    }
    return true;
}

static bool
complex_columns_are_orthonormal(const matrix::complex_matrix_base &vectors,
                                double tolerance) {
    for (size_t i = 0; i < vectors.cols(); ++i) {
        for (size_t j = 0; j < vectors.cols(); ++j) {
            std::complex<double> dot = 0.0;
            for (size_t row = 0; row < vectors.rows(); ++row) {
                dot += std::conj(vectors(row, i)) * vectors(row, j);
            }
            const double expected = i == j ? 1.0 : 0.0;
            if (std::abs(dot - expected) > tolerance) {
                return false;
            }
        }
    }
    return true;
}

static bool lower_unit_triangular(const matrix::real_matrix_base &A,
                                  double tolerance) {
    if (A.rows() != A.cols()) {
        return false;
    }
    for (size_t j = 0; j < A.cols(); ++j) {
        for (size_t i = 0; i < A.rows(); ++i) {
            if (i == j) {
                if (std::abs(A(i, j) - 1.0) > tolerance) {
                    return false;
                }
            } else if (i < j && std::abs(A(i, j)) > tolerance) {
                return false;
            }
        }
    }
    return true;
}

static bool upper_triangular(const matrix::real_matrix_base &A,
                             double tolerance) {
    for (size_t j = 0; j < A.cols(); ++j) {
        for (size_t i = j + 1; i < A.rows(); ++i) {
            if (std::abs(A(i, j)) > tolerance) {
                return false;
            }
        }
    }
    return true;
}

static double singular_triplet_relative_residual(
    const matrix::real_matrix_base &A,
    const matrix::truncated_svd_result &decomposition,
    size_t index) {
    double residual_squared = 0.0;
    double matrix_norm_squared = 0.0;
    for (const double value : A) {
        matrix_norm_squared += value * value;
    }
    for (size_t row = 0; row < A.rows(); ++row) {
        double av = 0.0;
        for (size_t column = 0; column < A.cols(); ++column) {
            av += A(row, column) * decomposition.V(column, index);
        }
        const double residual = av
                                - decomposition.singular_values[index]
                                      * decomposition.U(row, index);
        residual_squared += residual * residual;
    }
    return std::sqrt(residual_squared)
           / std::max(1.0, std::sqrt(matrix_norm_squared));
}

static matrix::matrixd
diagonal_test_matrix(size_t rows, size_t cols, size_t numerical_rank) {
    matrix::matrixd result(rows, cols);
    const size_t diagonal_size = std::min(rows, cols);
    for (size_t i = 0; i < std::min(diagonal_size, numerical_rank); ++i) {
        result(i, i) = i == 0 ? 12.0 : (i == 1 ? 8.0 : 5.0);
    }
    for (size_t i = numerical_rank; i < diagonal_size; ++i) {
        result(i, i) = 0.1 / static_cast<double>(i + 1);
    }
    return result;
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

void write_matrix(std::ofstream &file, const matrix::matrixd &mat) {
    file << std::setprecision(17);
    for (size_t i = 0; i < mat.rows(); ++i) {
        for (size_t j = 0; j < mat.cols(); ++j) {
            file << mat(i, j) << (j + 1 == mat.cols() ? '\n' : ' ');
        }
    }
}

int test_base_op() {
    std::cout << "MSL matrix unit tests\n";

    // 1) 基本构造与 factory
    TEST_CASE("identity / zeros / ones / diagonal") {
        auto I = matrix::matrixd::identity(3); // matrixd = matrix_real alias
        EXPECT_EQ(I.rows(), 3);
        EXPECT_EQ(I.cols(), 3);
        for (size_t i = 0; i < 3; ++i)
            for (size_t j = 0; j < 3; ++j)
                EXPECT_EQ(I(i, j), (i == j ? 1.0 : 0.0));

        auto Z = matrix::matrixd::zeros(2, 3);
        EXPECT_EQ(Z.rows(), 2);
        EXPECT_EQ(Z.cols(), 3);
        for (auto v : Z)
            EXPECT_EQ(v, 0.0);

        auto O = matrix::matrixd::ones(2, 2);
        for (auto v : O)
            EXPECT_EQ(v, 1.0);

        std::vector<double> diag = {1.5, 2.5, 3.5};
        auto temp = matrix::matrixd::ones(3, 3);
        temp(0, 0) = diag[0];
        temp(1, 1) = diag[1];
        temp(2, 2) = diag[2];
        auto D =
            temp.diagonal(std::span<const double>(diag.data(), diag.size()));
        EXPECT_EQ(D.rows(), 3);
        EXPECT_EQ(D.cols(), 3);
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                if (i == j)
                    EXPECT_EQ(D(i, j), diag[i]);
                else
                    EXPECT_EQ(D(i, j), 0.0);
            }
        }
    }
    while (0)
        ;

    // 2) element access, resize, fill
    TEST_CASE("element access / resize / fill") {
        matrix::matrixd A(3, 2);
        EXPECT_EQ(A.rows(), 3);
        EXPECT_EQ(A.cols(), 2);
        A(0, 0) = 1.0;
        A(1, 0) = 2.0;
        A(2, 0) = 3.0;
        A(0, 1) = 4.0;
        A(1, 1) = 5.0;
        A(2, 1) = 6.0;
        EXPECT_EQ(A(2, 1), 6.0);
        A.fill(7.0);
        for (auto v : A)
            EXPECT_EQ(v, 7.0);

        A.resize(2, 2);
        EXPECT_EQ(A.rows(), 2);
        EXPECT_EQ(A.cols(), 2);
        A.fill_zeros();
        for (auto v : A)
            EXPECT_EQ(v, 0.0);
    }
    while (0)
        ;

    // 3) copy / move / swap
    TEST_CASE("copy / move / swap") {
        matrix::matrixd A(2, 2, 1.0);
        matrix::matrixd B = A; // copy
        EXPECT_EQ(B(0, 0), 1.0);
        matrix::matrixd C = std::move(B); // move
        EXPECT_EQ(C(1, 1), 1.0);
        // after move B may be empty (implementation-dependent) but should be
        // safe to resize
        C(0, 0) = 9.0;
        A.swap(C);
        EXPECT_EQ(A(0, 0), 9.0);
    }
    while (0)
        ;

    // 4) arithmetic operators
    TEST_CASE("operators + - * /") {
        matrix::matrixd A(2, 2);
        A.fill(2.0);
        matrix::matrixd B(2, 2);
        B.fill(1.0);
        auto S = A + B; // should be 3's
        for (auto v : S)
            EXPECT_EQ(v, 3.0);
        auto D = A - B; // should be 1's
        for (auto v : D)
            EXPECT_EQ(v, 1.0);
        auto scaled = A * 0.5;
        for (auto v : scaled)
            EXPECT_EQ(v, 1.0);
        auto scaled2 = 0.5 * A;
        for (auto v : scaled2)
            EXPECT_EQ(v, 1.0);
        auto divv = A / 2.0;
        for (auto v : divv)
            EXPECT_EQ(v, 1.0);

        // in-place ops
        A += B;
        for (auto v : A)
            EXPECT_EQ(v, 3.0);
        A -= B;
        for (auto v : A)
            EXPECT_EQ(v, 2.0);
        A *= 2.0;
        for (auto v : A)
            EXPECT_EQ(v, 4.0);
        A /= 2.0;
        for (auto v : A)
            EXPECT_EQ(v, 2.0);
    }
    while (0)
        ;

    // 5) row/column/submatrix/get_row/get_column
    TEST_CASE("row/col/submatrix") {
        matrix::matrixd M(3, 4);
        // fill with M(i,j) = i + 10*j  (col-major aware)
        for (size_t j = 0; j < M.cols(); ++j)
            for (size_t i = 0; i < M.rows(); ++i)
                M(i, j) = double(i + 10 * j);

        auto row1 = M.get_row(1);
        EXPECT_EQ(row1.rows(), 1);
        EXPECT_EQ(row1.cols(), M.cols());
        for (size_t j = 0; j < M.cols(); ++j)
            EXPECT_EQ(row1(0, j), M(1, j));

        auto col2 = M.get_column(2);
        EXPECT_EQ(col2.rows(), M.rows());
        EXPECT_EQ(col2.cols(), 1);
        for (size_t i = 0; i < M.rows(); ++i)
            EXPECT_EQ(col2(i, 0), M(i, 2));

        auto sub = M.submatrix(1, 3, 1, 4); // rows 1..2, cols 1..3
        EXPECT_EQ(sub.rows(), 2);
        EXPECT_EQ(sub.cols(), 3);
        for (size_t j = 0; j < sub.cols(); ++j)
            for (size_t i = 0; i < sub.rows(); ++i)
                EXPECT_EQ(sub(i, j), M(1 + i, 1 + j));
    }
    while (0)
        ;

    // 6) iterators and range-for / operator[] if provided
    TEST_CASE("iterators and operator[]") {
        matrix::matrixd A(3, 2);
        double val = 0.0;
        for (size_t j = 0; j < A.cols(); ++j)
            for (size_t i = 0; i < A.rows(); ++i)
                A(i, j) = val++;
        // sum via iterator
        double s1 = 0.0;
        for (auto &x : A)
            s1 += x;
        double s2 = 0.0;
        for (size_t k = 0; k < A.size(); ++k)
            s2 += A[k]; // operator[] assumed
        EXPECT_EQ(s1, s2);
    }
    while (0)
        ;

    // 7) matrix_view zero-copy and construction from external data
    TEST_CASE("matrix_view zero-copy") {
        const size_t R = 3, C = 2;
        double raw[R * C] = {
            1, 2, 3, 4, 5, 6}; // row-major-like literal BUT we'll treat as
                               // column-major when filling
        // Fill raw in column-major to match matrix_view expectation (col-major
        // index j*R + i)
        for (size_t j = 0; j < C; ++j)
            for (size_t i = 0; i < R; ++i)
                raw[j * R + i] = double(j * 10 + i);

        matrix::matrixd_view V(raw, R, C);
        EXPECT_EQ(V.rows(), R);
        EXPECT_EQ(V.cols(), C);
        // modify through view and check raw changed
        V(1, 1) = 99.0;
        EXPECT_EQ(raw[1 + 1 * R], 99.0);

        // construct owning matrix from view
        matrix::matrixd M = matrix::matrixd(V); // should copy data
        EXPECT_EQ(M.rows(), R);
        EXPECT_EQ(M.cols(), C);
        for (size_t j = 0; j < C; ++j)
            for (size_t i = 0; i < R; ++i)
                EXPECT_EQ(M(i, j), V(i, j));
    }
    while (0)
        ;

    // 8) multiplication (matrix * matrix) and (matrix * vector)
    TEST_CASE("matrix multiplication") {
        matrix::matrixd A(2, 3); // 2x3
        matrix::matrixd B(3, 2); // 3x2
        // A = [1 2 3; 4 5 6] in column-major -> careful filling
        // We'll set A(i,j) = i + j*10 for clarity
        for (size_t j = 0; j < A.cols(); ++j)
            for (size_t i = 0; i < A.rows(); ++i)
                A(i, j) = double(i + j * 10);
        for (size_t j = 0; j < B.cols(); ++j)
            for (size_t i = 0; i < B.rows(); ++i)
                B(i, j) = double(1 + i + j * 10);
        auto A_eigen = msl::matrix::eigen_interface::as_eigen(A);
        auto B_eigen = msl::matrix::eigen_interface::as_eigen(B);

        auto C_eigen = A_eigen * B_eigen; // should be 2x2
        auto C = msl::matrix::eigen_interface::from_eigen(C_eigen);
        EXPECT_EQ(C.rows(), 2);
        EXPECT_EQ(C.cols(), 2);

        // naive recompute
        std::vector<double> expected(4, 0.0);
        for (size_t i = 0; i < C.rows(); ++i)
            for (size_t j = 0; j < C.cols(); ++j)
                for (size_t k = 0; k < A.cols(); ++k)
                    expected[j * C.rows() + i] += A(i, k) * B(k, j);

        for (size_t j = 0; j < C.cols(); ++j)
            for (size_t i = 0; i < C.rows(); ++i)
                EXPECT_EQ(C(i, j), expected[j * C.rows() + i]);

        // matrix * std::vector
        std::vector<double> x(B.cols(), 1.0); // length 2
        Eigen::VectorXd y_eigen =
            B_eigen * Eigen::Map<Eigen::VectorXd>(x.data(), x.size()).eval();
        auto y = msl::matrix::eigen_interface::from_eigen(y_eigen);
        EXPECT_EQ(y.size(), B.rows());
        // y_i = sum_j B(i,j) * x[j]
        for (size_t i = 0; i < B.rows(); ++i) {
            double s = 0.0;
            for (size_t j = 0; j < B.cols(); ++j)
                s += B(i, j) * x[j];
            EXPECT_EQ(y[i], s);
        }
    }
    while (0)
        ;

    // 9) complex matrix basic use
    TEST_CASE("complex matrix basics") {
        matrix::matrixc C(2, 2);
        C(0, 0) = std::complex<double>(1.0, 2.0);
        C(1, 0) = std::complex<double>(3.0, -1.0);
        EXPECT_CPLX_NEAR(C(0, 0), std::complex<double>(1.0, 2.0), 1e-12);
        // conjugate if method exists (we will call if present)
        // use SFINAE would be robust, but for test we simply attempt and skip
        // compile-time if not present
    }
    while (0)
        ;

    // 总结
    std::cout << "Total checks: " << g_total << ", failures: " << g_failures
              << "\n";
    if (g_failures == 0) {
        std::cout << "All tests passed\n";
        return 0;
    } else {
        std::cerr << "Some tests failed\n";
        return 1;
    }
}

// ---------- decomposition tests (SVD, Eigen, etc.) ----------
int test_decompositions() {
    // Real matrix SVD test
    TEST_CASE("Real matrix SVD") {
        matrix::matrixd A(3, 2);
        A(0, 0) = 3.0;
        A(1, 0) = 2.0;
        A(2, 0) = 1.0;
        A(0, 1) = 4.0;
        A(1, 1) = 5.0;
        A(2, 1) = 6.0;

        auto svd = matrix::svd(A);
        auto U = svd[0];
        auto S = svd[1];
        auto V = svd[2];

        // Reconstruct A from U * S * V^T
        auto US = U * S;
        auto A_reconstructed = US * matrix::transpose(V);

        // Check reconstruction
        for (size_t j = 0; j < A.cols(); ++j)
            for (size_t i = 0; i < A.rows(); ++i)
                EXPECT_NEAR(A(i, j), A_reconstructed(i, j), 1e-10);

        // Singular values are non-negative and in descending order
        for (size_t i = 0; i < std::min(A.rows(), A.cols()); ++i) {
            EXPECT_TRUE(S(i, i) >= 0.0);
            if (i > 0) {
                EXPECT_TRUE(S(i - 1, i - 1) >= S(i, i));
            }
        }

        // Left and right singular vectors are orthonormal
        EXPECT_TRUE(columns_are_orthonormal(U, 1e-10));
        EXPECT_TRUE(columns_are_orthonormal(V, 1e-10));
    }
    while (0)
        ;

    // Complex matrix SVD test
    TEST_CASE("Complex matrix SVD") {
        matrix::matrixc A(3, 2);
        A(0, 0) = std::complex<double>(3.0, 1.0);
        A(1, 0) = std::complex<double>(2.0, -1.0);
        A(2, 0) = std::complex<double>(1.0, 1.0);
        A(0, 1) = std::complex<double>(4.0, -1.0);
        A(1, 1) = std::complex<double>(5.0, 1.0);
        A(2, 1) = std::complex<double>(6.0, -1.0);

        auto svd = matrix::svd(A);
        auto U = svd[0];
        auto S = svd[1];
        auto V = svd[2];

        // Reconstruct A from U * S * V^H
        auto US = U * S;
        auto A_reconstructed = US * matrix::conjugate_transpose(V);

        // Check reconstruction
        for (size_t j = 0; j < A.cols(); ++j)
            for (size_t i = 0; i < A.rows(); ++i)
                EXPECT_CPLX_NEAR(A(i, j), A_reconstructed(i, j), 1e-10);

        // Singular values are real, non-negative and descending
        for (size_t i = 0; i < std::min(A.rows(), A.cols()); ++i) {
            EXPECT_NEAR(S(i, i).imag(), 0.0, 1e-12);
            EXPECT_TRUE(S(i, i).real() >= 0.0);
            if (i > 0) {
                EXPECT_TRUE(S(i - 1, i - 1).real() >= S(i, i).real());
            }
        }

        // Left and right singular vectors are orthonormal
        EXPECT_TRUE(complex_columns_are_orthonormal(U, 1e-10));
        EXPECT_TRUE(complex_columns_are_orthonormal(V, 1e-10));
    }
    while (0)
        ;

    TEST_CASE("LU decomposition preserves partial pivoting") {
        const auto check_real = [](const matrix::matrixd &A, double tol) {
            const auto result = matrix::lu(A);
            const size_t n = A.rows();
            EXPECT_EQ(result.permutation.size(), n);
            EXPECT_EQ(result.L.rows(), n);
            EXPECT_EQ(result.L.cols(), n);
            EXPECT_EQ(result.U.rows(), n);
            EXPECT_EQ(result.U.cols(), n);

            matrix::matrixd PA(n, n);
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    PA(i, j) = A(result.permutation[i], j);
                }
            }
            EXPECT_TRUE(matrix_near(PA, result.L * result.U, tol));
            EXPECT_TRUE(lower_unit_triangular(result.L, tol));
            EXPECT_TRUE(upper_triangular(result.U, tol));
        };

        // No pivoting required
        check_real(matrix::matrixd(2, 2, {2.0, 1.0, 1.0, 3.0}), 1e-12);
        // Pivoting required
        check_real(matrix::matrixd(2, 2, {0.0, 1.0, 1.0, 1.0}), 1e-12);
        check_real(matrix::matrixd(
                       3, 3, {0.0, 2.0, 1.0, 1.0, 1.0, 1.0, 2.0, 0.0, 1.0}),
                   1e-12);

        // Complex pivoting case
        matrix::matrixc A_complex(
            2, 2, {{0.0, 0.0}, {0.0, 1.0}, {1.0, 0.0}, {0.0, 0.0}});
        const auto complex_result = matrix::lu(A_complex);
        EXPECT_EQ(complex_result.permutation.size(), 2);
        EXPECT_EQ(complex_result.L.rows(), 2);
        EXPECT_EQ(complex_result.U.rows(), 2);
        matrix::matrixc complex_PA(2, 2);
        for (size_t i = 0; i < 2; ++i) {
            for (size_t j = 0; j < 2; ++j) {
                complex_PA(i, j) = A_complex(complex_result.permutation[i], j);
            }
        }
        const auto complex_LU = complex_result.L * complex_result.U;
        for (size_t i = 0; i < 2; ++i) {
            for (size_t j = 0; j < 2; ++j) {
                EXPECT_CPLX_NEAR(complex_PA(i, j), complex_LU(i, j), 1e-12);
            }
        }

        // Non-square matrices are rejected
        bool threw_real = false;
        try {
            (void)matrix::lu(matrix::matrixd(2, 3));
        } catch (const std::invalid_argument &) {
            threw_real = true;
        }
        EXPECT_TRUE(threw_real);
        bool threw_complex = false;
        try {
            (void)matrix::lu(matrix::matrixc(2, 3));
        } catch (const std::invalid_argument &) {
            threw_complex = true;
        }
        EXPECT_TRUE(threw_complex);
    }
    while (0)
        ;

    TEST_CASE("QR decomposition supports tall, square and wide matrices") {
        const auto check_economy = [](const matrix::matrixd &A, double tol) {
            const size_t m = A.rows();
            const size_t n = A.cols();
            const size_t k = std::min(m, n);
            const auto result = matrix::qr(A);
            const auto &Q = result[0];
            const auto &R = result[1];

            EXPECT_EQ(Q.rows(), m);
            EXPECT_EQ(Q.cols(), k);
            EXPECT_EQ(R.rows(), k);
            EXPECT_EQ(R.cols(), n);
            EXPECT_TRUE(matrix_near(Q * R, A, tol));
            EXPECT_TRUE(columns_are_orthonormal(Q, tol));
            EXPECT_TRUE(upper_triangular(R, tol));
        };

        const auto check_full = [](const matrix::matrixd &A, double tol) {
            const size_t m = A.rows();
            const size_t n = A.cols();
            const auto result = matrix::qr(A, true);

            EXPECT_EQ(result[0].rows(), m);
            EXPECT_EQ(result[0].cols(), m);
            EXPECT_EQ(result[1].rows(), m);
            EXPECT_EQ(result[1].cols(), n);
            EXPECT_TRUE(matrix_near(result[0] * result[1], A, tol));
            EXPECT_TRUE(columns_are_orthonormal(result[0], tol));
            EXPECT_TRUE(upper_triangular(result[1], tol));
        };

        const matrix::matrixd tall(
            4, 2, {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0});
        const matrix::matrixd square(
            3, 3, {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0});
        const matrix::matrixd wide(
            2, 4, {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0});
        const matrix::matrixd rank_deficient(
            2, 3, {1.0, 2.0, 3.0, 2.0, 4.0, 6.0});

        check_economy(tall, 1e-12);
        check_economy(square, 1e-12);
        check_economy(wide, 1e-12);
        check_economy(rank_deficient, 1e-12);

        check_full(tall, 1e-12);
        check_full(wide, 1e-12);
    }
    while (0)
        ;

    TEST_CASE("Real nonsymmetric eigenvalues retain conjugate pairs") {
        constexpr double a = 0.8;
        constexpr double b = 0.3;
        matrix::matrixd A(2, 2, {a, -b, b, a});
        const auto decomposition = matrix::eig(A);
        const auto &V = decomposition[0];
        const auto &D = decomposition[1];

        bool found_positive = false;
        bool found_negative = false;
        const double matrix_norm = std::sqrt(2.0 * (a * a + b * b));
        for (size_t column = 0; column < 2; ++column) {
            const auto lambda = D(column, column);
            EXPECT_NEAR(lambda.real(), a, 1e-12);
            found_positive = found_positive || lambda.imag() > 0.0;
            found_negative = found_negative || lambda.imag() < 0.0;

            double residual_squared = 0.0;
            double vector_norm_squared = 0.0;
            for (size_t i = 0; i < 2; ++i) {
                std::complex<double> av = 0.0;
                for (size_t j = 0; j < 2; ++j) {
                    av += A(i, j) * V(j, column);
                }
                residual_squared += std::norm(av - lambda * V(i, column));
                vector_norm_squared += std::norm(V(i, column));
            }
            const double relative_residual =
                std::sqrt(residual_squared)
                / (matrix_norm * std::sqrt(vector_norm_squared));
            EXPECT_TRUE(relative_residual < 1e-12);
        }
        EXPECT_TRUE(found_positive);
        EXPECT_TRUE(found_negative);

        auto B = matrix::matrixd::identity(2) * 2.0;
        const auto generalized = matrix::eig(A, B);
        EXPECT_NEAR(std::abs(generalized[1](0, 0).imag()), b / 2.0, 1e-12);
        EXPECT_NEAR(std::abs(generalized[1](1, 1).imag()), b / 2.0, 1e-12);
    }
    while (0)
        ;

    TEST_CASE("Moore-Penrose pseudoinverse") {
        matrix::matrixd tall(3, 2, {1.0, 2.0, 2.0, 4.0, 3.0, 7.0});
        const auto tall_pinv = matrix::pinv(tall);
        EXPECT_EQ(tall_pinv.rows(), 2);
        EXPECT_EQ(tall_pinv.cols(), 3);
        EXPECT_TRUE(matrix_near((tall * tall_pinv) * tall, tall, 1e-11));
        EXPECT_TRUE(
            matrix_near((tall_pinv * tall) * tall_pinv, tall_pinv, 1e-11));

        const auto wide = matrix::transpose(tall);
        const auto wide_pinv = matrix::pinv(wide);
        EXPECT_EQ(wide_pinv.rows(), 3);
        EXPECT_EQ(wide_pinv.cols(), 2);
        EXPECT_TRUE(matrix_near((wide * wide_pinv) * wide, wide, 1e-11));

        matrix::matrixd deficient(3, 2, {1.0, 2.0, 2.0, 4.0, 3.0, 6.0});
        const auto deficient_pinv = matrix::pinv(deficient);
        EXPECT_TRUE(matrix_near(
            (deficient * deficient_pinv) * deficient, deficient, 1e-11));

        const auto zero = matrix::matrixd::zeros(2, 3);
        const auto zero_pinv = matrix::pinv(zero);
        EXPECT_EQ(zero_pinv.rows(), 3);
        EXPECT_EQ(zero_pinv.cols(), 2);
        for (const double value : zero_pinv) {
            EXPECT_EQ(value, 0.0);
        }

        matrix::matrixd ill_conditioned(2, 2, {1.0, 0.0, 0.0, 1e-10});
        const auto reduced_rank = matrix::pinv(ill_conditioned, 1e-8);
        const auto full_rank = matrix::pinv(ill_conditioned, 1e-12);
        EXPECT_NEAR(reduced_rank(0, 0), 1.0, 1e-12);
        EXPECT_EQ(reduced_rank(1, 1), 0.0);
        EXPECT_NEAR(full_rank(1, 1), 1e10, 1e-3);
    }
    while (0)
        ;

    TEST_CASE("Leading singular triplets") {
        matrix::matrixd A(3, 2, {3.0, 0.0, 0.0, 2.0, 0.0, 0.0});
        const auto leading = matrix::truncated_svd(A, 1);
        EXPECT_EQ(leading.U.rows(), 3);
        EXPECT_EQ(leading.U.cols(), 1);
        EXPECT_EQ(leading.V.rows(), 2);
        EXPECT_EQ(leading.V.cols(), 1);
        EXPECT_EQ(leading.singular_values.size(), 1);
        EXPECT_NEAR(leading.singular_values[0], 3.0, 1e-12);

        matrix::matrixd reconstructed(3, 2);
        for (size_t i = 0; i < reconstructed.rows(); ++i) {
            for (size_t j = 0; j < reconstructed.cols(); ++j) {
                reconstructed(i, j) = leading.U(i, 0)
                                      * leading.singular_values[0]
                                      * leading.V(j, 0);
            }
        }
        EXPECT_TRUE(
            matrix_near(reconstructed,
                        matrix::matrixd(3, 2, {3.0, 0.0, 0.0, 0.0, 0.0, 0.0}),
                        1e-12));

        const auto left_only = matrix::truncated_svd(A, 1, false);
        EXPECT_EQ(left_only.V.rows(), 0);
        EXPECT_EQ(left_only.V.cols(), 0);
    }
    while (0)
        ;

    TEST_CASE("Deterministic randomized truncated SVD") {
        const auto check_shape = [](size_t rows, size_t cols) {
            constexpr size_t requested_rank = 3;
            auto A = diagonal_test_matrix(rows, cols, requested_rank);
            const auto full = matrix::svd(A);
            const auto leading = matrix::truncated_svd(A, requested_rank);
            const auto repeated = matrix::truncated_svd(A, requested_rank);

            EXPECT_EQ(leading.U.rows(), rows);
            EXPECT_EQ(leading.U.cols(), requested_rank);
            EXPECT_EQ(leading.V.rows(), cols);
            EXPECT_EQ(leading.V.cols(), requested_rank);
            EXPECT_EQ(leading.singular_values.size(), requested_rank);
            EXPECT_TRUE(columns_are_orthonormal(leading.U, 1e-10));
            EXPECT_TRUE(columns_are_orthonormal(leading.V, 1e-10));
            EXPECT_TRUE(matrix_near(leading.U, repeated.U, 1e-13));
            EXPECT_TRUE(matrix_near(leading.V, repeated.V, 1e-13));

            for (size_t i = 0; i < requested_rank; ++i) {
                const double relative_error =
                    std::abs(leading.singular_values[i] - full[1](i, i))
                    / full[1](i, i);
                EXPECT_TRUE(relative_error < 1e-6);
                EXPECT_TRUE(i == 0
                            || leading.singular_values[i - 1]
                                   >= leading.singular_values[i]);
                EXPECT_TRUE(singular_triplet_relative_residual(A, leading, i)
                            < 1e-6);
                EXPECT_NEAR(leading.singular_values[i],
                            repeated.singular_values[i],
                            1e-13);
            }
        };

        check_shape(24, 24);
        check_shape(32, 20);
        check_shape(20, 32);

        matrix::matrixd repeated_values(24, 24);
        repeated_values(0, 0) = 10.0;
        repeated_values(1, 1) = 10.0;
        for (size_t i = 2; i < 24; ++i) {
            repeated_values(i, i) = 0.1 / static_cast<double>(i + 1);
        }
        const auto repeated_subspace =
            matrix::truncated_svd(repeated_values, 2);
        EXPECT_NEAR(repeated_subspace.singular_values[0], 10.0, 1e-8);
        EXPECT_NEAR(repeated_subspace.singular_values[1], 10.0, 1e-8);
        EXPECT_TRUE(columns_are_orthonormal(repeated_subspace.U, 1e-10));
        EXPECT_TRUE(columns_are_orthonormal(repeated_subspace.V, 1e-10));
        for (size_t i = 0; i < 2; ++i) {
            EXPECT_TRUE(singular_triplet_relative_residual(
                            repeated_values, repeated_subspace, i)
                        < 1e-8);
        }
    }
    while (0)
        ;

    TEST_CASE("Randomized truncated SVD rank-deficient and zero matrices") {
        matrix::truncated_svd_options options;
        options.oversampling = 2;
        options.power_iterations = 2;

        matrix::matrixd deficient(12, 10);
        deficient(0, 0) = 4.0;
        deficient(1, 1) = 2.0;
        const auto deficient_result =
            matrix::truncated_svd(deficient, 4, options);
        EXPECT_NEAR(deficient_result.singular_values[0], 4.0, 1e-12);
        EXPECT_NEAR(deficient_result.singular_values[1], 2.0, 1e-12);
        EXPECT_NEAR(deficient_result.singular_values[2], 0.0, 1e-12);
        EXPECT_NEAR(deficient_result.singular_values[3], 0.0, 1e-12);
        EXPECT_TRUE(columns_are_orthonormal(deficient_result.U, 1e-10));
        EXPECT_TRUE(columns_are_orthonormal(deficient_result.V, 1e-10));
        for (size_t i = 0; i < 4; ++i) {
            EXPECT_TRUE(singular_triplet_relative_residual(
                            deficient, deficient_result, i)
                        < 1e-10);
        }

        const auto zero = matrix::matrixd::zeros(12, 10);
        const auto zero_result = matrix::truncated_svd(zero, 4, options);
        EXPECT_TRUE(columns_are_orthonormal(zero_result.U, 1e-10));
        EXPECT_TRUE(columns_are_orthonormal(zero_result.V, 1e-10));
        for (const double value : zero_result.singular_values) {
            EXPECT_EQ(value, 0.0);
        }

        options.compute_right_vectors = false;
        const auto left_only = matrix::truncated_svd(deficient, 4, options);
        EXPECT_EQ(left_only.U.rows(), deficient.rows());
        EXPECT_EQ(left_only.U.cols(), 4);
        EXPECT_EQ(left_only.V.rows(), 0);
        EXPECT_EQ(left_only.V.cols(), 0);
        EXPECT_TRUE(columns_are_orthonormal(left_only.U, 1e-10));
    }
    while (0)
        ;

    TEST_CASE("Randomized truncated SVD validates rank and finite input") {
        matrix::matrixd A(4, 3);
        bool threw_zero_rank = false;
        bool threw_large_rank = false;
        bool threw_nan = false;
        bool threw_infinity = false;
        try {
            (void)matrix::truncated_svd(A, 0);
        } catch (const std::invalid_argument &) {
            threw_zero_rank = true;
        }
        try {
            (void)matrix::truncated_svd(A, 4);
        } catch (const std::invalid_argument &) {
            threw_large_rank = true;
        }
        A(0, 0) = std::numeric_limits<double>::quiet_NaN();
        try {
            (void)matrix::truncated_svd(A, 1);
        } catch (const std::invalid_argument &) {
            threw_nan = true;
        }
        A(0, 0) = std::numeric_limits<double>::infinity();
        try {
            (void)matrix::truncated_svd(A, 1);
        } catch (const std::invalid_argument &) {
            threw_infinity = true;
        }
        EXPECT_TRUE(threw_zero_rank);
        EXPECT_TRUE(threw_large_rank);
        EXPECT_TRUE(threw_nan);
        EXPECT_TRUE(threw_infinity);
    }
    while (0)
        ;

    // 总结
    std::cout << "Total decomposition checks: " << g_total
              << ", failures: " << g_failures << "\n";
    if (g_failures == 0) {
        std::cout << "All decomposition tests passed\n";
        return 0;
    } else {
        std::cerr << "Some decomposition tests failed\n";
        return 1;
    }
}

// ---------- 主程序 ----------
int main() {
    int res1 = test_base_op();
    int res2 = test_decompositions();

    auto compare_dir =
        project_root() / "test_result" / "matrix" / "matlab_compare";
    std::filesystem::create_directories(compare_dir);

    matrix::matrixd A(2, 3);
    matrix::matrixd B(3, 2);
    for (size_t i = 0; i < A.rows(); ++i) {
        for (size_t j = 0; j < A.cols(); ++j) {
            A(i, j) = static_cast<double>(1 + i + 2 * j);
        }
    }
    for (size_t i = 0; i < B.rows(); ++i) {
        for (size_t j = 0; j < B.cols(); ++j) {
            B(i, j) = static_cast<double>(1 + 3 * i + j);
        }
    }
    auto product = A * B;

    std::ofstream a_file(compare_dir / "matrix_A.txt");
    write_matrix(a_file, A);
    std::ofstream b_file(compare_dir / "matrix_B.txt");
    write_matrix(b_file, B);
    std::ofstream product_file(compare_dir / "matrix_product.txt");
    write_matrix(product_file, product);

    matrix::matrixd svd_input(3, 2);
    svd_input(0, 0) = 3.0;
    svd_input(1, 0) = 2.0;
    svd_input(2, 0) = 1.0;
    svd_input(0, 1) = 4.0;
    svd_input(1, 1) = 5.0;
    svd_input(2, 1) = 6.0;
    auto svd = matrix::svd(svd_input);
    auto reconstructed = (svd[0] * svd[1]) * matrix::transpose(svd[2]);
    std::ofstream svd_input_file(compare_dir / "matrix_svd_input.txt");
    write_matrix(svd_input_file, svd_input);
    std::ofstream svd_file(compare_dir / "matrix_svd_reconstruction.txt");
    write_matrix(svd_file, reconstructed);

    std::cout << "Test results: base operations = " << res1
              << ", decompositions = " << res2 << "\n";

    return res1 + res2;
}
