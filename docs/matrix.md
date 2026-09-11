# Matrix Module

**Namespace**: `msl::matrix`  
**Aggregator**: `msl/matrix.hpp`  
**Files**: 11 headers in `msl/matrix/`

## Overview

The matrix module provides real and complex matrix types with column-major storage, arithmetic operations, linear algebra (transpose, inverse, determinant), matrix decompositions (QR, LU, SVD, eigenvalue), and an Eigen3 bridge for zero-copy interop.

## Class Hierarchy

```
real_matrix_base          (std::span<double>, no virtual methods)
  |-- real_matrix_owned   (owns std::vector<double> = matrixd)
  |-- real_matrix_view    (non-owning mutable view = matrixd_view)

complex_matrix_base       (std::span<std::complex<double>>)
  |-- complex_matrix_owned (owns std::vector<complex<double>> = matrixc)
  |-- complex_matrix_view  (non-owning mutable view)

const_real_matrix_view    (read-only view = const_matrixd_view)
```

### Type Aliases

| Alias | Full Type |
|-------|-----------|
| `matrixd` | `real_matrix_owned` |
| `matrixc` | `complex_matrix_owned` |
| `matrixd_view` | `real_matrix_view` |
| `const_matrixd_view` | `const_real_matrix_view` |

## Column-Major Storage

All matrices use column-major storage: element `(i, j)` is at `data_[j * rows_ + i]`. This is compatible with Fortran, MATLAB, and Eigen's default layout.

```
Matrix A = [a b]     Storage: [a, c, b, d]
           [c d]     (2 rows, 2 cols)
```

## real_matrix_base

The base class for real matrices. Wraps a `std::span<double>` with row/column dimensions.

### Construction

```cpp
// From span (no ownership)
double buf[6] = {1,2,3,4,5,6};
real_matrix_base base(2, 3, buf);
```

### Element Access

```cpp
double &val = mat(i, j);           // (row, col) indexing
double &val = mat[idx];            // linear indexing (column-major)
std::span<double> col = mat.column(j); // column view
double *ptr = mat.data();          // raw pointer
```

### Size Queries

```cpp
size_t r = mat.rows();
size_t c = mat.cols();
size_t n = mat.size();            // rows * cols
auto [r, c] = mat.shape();
bool is_empty = mat.empty();
```

### Iteration

```cpp
for (auto &val : mat) { ... }     // column-major traversal
mat.apply([](double x) { return x * 2; }); // in-place transform
double s = mat.sum();             // sum of all elements
double t = mat.trace();           // trace (must be square)
```

## real_matrix_owned (matrixd)

Owns its storage via `std::vector<double>`. Full value semantics (deep copy, move).

### Construction

```cpp
matrixd A;                             // default (0x0)
matrixd A(3, 4);                       // 3x4, zero-initialized
matrixd A(3, 4, 1.0);                  // 3x4, filled with 1.0
matrixd A(2, 3, {1,2,3,4,5,6});       // from initializer list (row-major input → column-major storage)
matrixd A(2, 3, data_span);            // deep copy from span

matrixd A = matrixd::zeros(3, 3);      // factory: zeros
matrixd B = matrixd::ones(3, 3);       // factory: ones
matrixd I = matrixd::identity(4);      // factory: identity
matrixd D = matrixd::diagonal({1,2,3});// factory: diagonal
```

### Arithmetic

```cpp
matrixd C = A + B;
matrixd D = A - B;
matrixd E = A * B;    // matrix multiplication
matrixd F = A * 2.0;  // scalar multiplication
matrixd G = A / 2.0;  // scalar division
A += B;               // in-place add
A -= B;               // in-place subtract
A *= 2.0;             // in-place scalar multiply
```

### Resizing and Submatrix

```cpp
A.resize(5, 5, 0.0);       // resize, optional fill value
A.clear();                  // clear to 0x0
auto row = A.get_row(i);    // extract row (as new matrixd)
auto col = A.get_column(j); // extract column (as new matrixd)
auto sub = A.submatrix(r0, r1, c0, c1); // submatrix view
```

## real_matrix_view (matrixd_view)

Non-owning mutable view into an existing matrix. No copy — modifies the original.

```cpp
matrixd A(3, 3);
matrixd_view view(A);       // view into A
view(0, 0) = 42.0;          // modifies A(0,0)
matrixd copy = view.to_owned(); // deep copy
```

## const_real_matrix_view (const_matrixd_view)

Read-only view. Can be constructed from `const matrixd&`.

```cpp
void inspect(const const_matrixd_view &mat) {
    double val = mat(0, 0); // OK
    // mat(0, 0) = 42;   // ERROR: const view
}
```

## Complex Matrix Types

The same hierarchy exists for `std::complex<double>`:

```cpp
complex_matrix_base
  |-- complex_matrix_owned (matrixc)
  |-- complex_matrix_view
```

Identical API to real matrices, with `std::complex<double>` element type.

## Matrix Operations

Free functions in `msl::matrix` namespace:

| Function | Signature | Description |
|----------|-----------|-------------|
| `transpose` | `matrixd transpose(const real_matrix_base &A)` | Matrix transpose |
| `conjugate_transpose` | `matrixc conjugate_transpose(const complex_matrix_base &A)` | Conjugate transpose (Hermitian) |
| `trace` | `double trace(const real_matrix_base &A)` | Trace of square matrix |
| `inverse` | `matrixd inverse(const real_matrix_base &A)` | Matrix inverse (via Eigen) |
| `adjoint` | `matrixd adjoint(const matrixd &A)` | Classical adjoint (adjugate) |
| `determinant` | `double determinant(const real_matrix_base &A)` | Determinant (via Eigen) |

## Matrix Decompositions

All in `msl/matrix/martrix_decompose.hpp` (note the typo in the filename). All delegate to Eigen3.

### QR Decomposition

```cpp
auto [Q, R] = msl::matrix::qr(A);           // thin/economy QR (default)
auto [Q, R] = msl::matrix::qr(A, true);      // full QR
```
Returns `std::array<matrixd, 2>` — can use structured bindings.

For economy QR, with `k = min(rows, cols)`, `Q` is `m×k` and `R` is `k×n`,
satisfying `A ≈ Q*R` and `Q^T*Q ≈ I_k`. Full QR returns `Q` as `m×m` and `R`
as `m×n`.

### LU Decomposition

```cpp
auto result = msl::matrix::lu(A);
// result.permutation: row permutation
// result.L: unit lower triangular
// result.U: upper triangular
```
`A` must be square. Partial pivoting is used, so the decomposition satisfies
`P * A = L * U`, where row `i` of `P * A` is row `result.permutation[i]` of
`A`.

### SVD Decomposition

```cpp
auto [U, S, V] = msl::matrix::svd(A);
// A = U * S * V^T
```
`U` is `m×m`, `S` is `m×n` diagonal, `V` is `n×n`.

For callers that only consume the leading singular triplets:

```cpp
auto leading = msl::matrix::truncated_svd(A, 30, false);
// leading.U: m×30, leading.singular_values: 30 values
// leading.V is omitted because compute_right_vectors=false
```

`truncated_svd` uses deterministic randomized SVD. The default sample size is
`min(rank + 10, min(m,n))`, with two reorthogonalized power iterations and a
fixed seed. Only the projected small matrix is decomposed with JacobiSVD, and
the returned singular values are explicitly ordered from largest to smallest.

Advanced callers can override the defaults without changing the result type:

```cpp
msl::matrix::truncated_svd_options options;
options.oversampling = 10;
options.power_iterations = 2;
options.seed = 0x5eedULL;
options.compute_right_vectors = false;
auto leading = msl::matrix::truncated_svd(A, 30, options);
```

When `rank + oversampling` covers `min(m,n)`, the implementation falls back to
an exact thin SVD. Non-finite input and invalid ranks are rejected.

### Moore-Penrose Pseudoinverse

```cpp
auto A_plus = msl::matrix::pinv(A);        // MATLAB-like default tolerance
auto A_rank = msl::matrix::pinv(A, 1e-8); // explicit singular-value tolerance
```

The output is `n×m` for an `m×n` input. The default tolerance is
`max(m,n) * epsilon * largest_singular_value`; singular values less than or
equal to the selected tolerance are discarded.

### Eigenvalue Decomposition

```cpp
auto [V, D] = msl::matrix::eig(A);              // standard eigenvalue
auto [V, D] = msl::matrix::eig(A, B);           // generalized eigenvalue
// A * V = V * D  or  A * V = B * V * D
```
Returns complex matrices (`matrixc`). `V` contains eigenvectors, `D` is diagonal of eigenvalues.
For a real nonsymmetric input, conjugate eigenvalues and eigenvectors retain
their complex components.

## Eigen3 Interface

`msl/matrix/eigen_interface.hpp` provides zero-copy interop with Eigen3.

### Zero-Copy Views

```cpp
// MSL → Eigen (zero-copy)
auto eig_map = msl::matrix::eigen_interface::as_eigen(A);
// eig_map is an Eigen::Map<Eigen::MatrixXd>

// Eigen → MSL (deep copy)
auto msl_mat = msl::matrix::eigen_interface::from_eigen(eig_matrix);
```

### Convenience Functions

```cpp
auto C = msl::matrix::matmul(A, B);     // matrix multiply
auto y = msl::matrix::matvec(A, x);     // matrix-vector multiply
auto x = msl::matrix::solve(A, b);      // solve Ax = b (via Eigen)
```

The `view_from_eigen()` function creates MSL views from Eigen data. Rvalue overloads are deleted to prevent dangling references.
