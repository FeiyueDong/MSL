# MSL — Agent Guide

## Build

MSL is a C++20 header-only library. Xmake 3.0.3 or later is the build entrypoint;
Eigen3 is the single external C++ dependency, declared through
`add_requires("eigen")` in `xmake/products.lua`.

```sh
xmake f -p linux -m release -y
xmake
```

Choose `-p windows` for MSVC, `-p mingw` for MinGW, or `-p macosx` for macOS.
MinGW SDK discovery can be overridden with `MSYS2_ROOT`. For a separate global
package cache, set `XMAKE_GLOBALDIR` consistently for all Xmake commands.

The root `xmake.lua` builds `msl` and the development tests. External consumers
include only `xmake/products.lua` and depend on `msl`; this publishes both
`include/` and `include/msl/` and the Eigen package, without test targets.

## Tests

```sh
xmake run test_matrix
xmake run test_signal
xmake run test_difference
xmake run test_integral
xmake run test_interp
xmake run test_polynomial
xmake run test_equation
xmake run test_ode
```

- Shared assertions are in `tests/test_utils.hpp`: `TEST_CASE`, `EXPECT_TRUE`,
  `EXPECT_EQ`, `EXPECT_NEAR`, `EXPECT_CPLX_NEAR`. Finish with
  `msl::test::summary()` and a failing exit code when checks fail.
- Prefer deterministic in-memory fixtures and mathematical invariants. Check
  errors and boundaries as well as expected values.
- Tests generate optional comparison data under
  `test_result/<module>/matlab_compare/` (gitignored).
- After refreshing C++ outputs, run optional MATLAB comparisons with
  `matlab -batch "run('matlab/validation/run_all.m')"`. Toolbox-dependent
  sections can be skipped; report that explicitly. C++ tests require neither
  MATLAB nor Python.
- CI defines Linux GCC, Linux Clang, and Windows MSVC jobs. Distinguish local
  validation from remote CI or MATLAB evidence.

## Structure

| Path | Purpose |
|------|---------|
| `include/msl/<module>.hpp` | Public aggregator headers |
| `include/msl/<module>/` | Header-only implementations |
| `xmake/products.lua` | Product target and public dependency/includes |
| `tests/` | Eight test targets and shared helpers |
| `matlab/validation/` | Independent numerical comparisons |
| `docs/` | Current API, compatibility, and development guides |
| `resource/` | Reference data and legacy MATLAB utilities |

Namespaces mirror modules: `msl::matrix`, `msl::signal`, `msl::difference`,
`msl::integral`, `msl::interp`, `msl::polynomial`, `msl::equation`,
`msl::ode`, `msl::utils`.

## Numerical and API conventions

- **Column-major storage**: `data()[j * rows() + i]`.
- **Ownership**: `row_copy`, `column_copy`, `submatrix_copy` own their results;
  `column()` and views borrow data. Views require their owners to stay alive.
- **Runtime validation**: public matrix operations throw
  `std::invalid_argument` or `std::out_of_range`; element indexing is unchecked.
- **Decompositions**: LU returns permutation plus factors with `P*A = L*U`;
  real nonsymmetric eigendecomposition preserves complex eigenpairs.
- **Spectra**: single-frame power spectra are unnormalized; Welch PSD/CPSD
  are density-normalized, full two-sided spectra using `X * conj(Y)`.
- **Filtering**: cutoffs are normalized to Nyquist. `filtfilt` uses odd
  reflection with `3 * filter_order` padding and a direct initial-state solve.
- **Polynomial coefficients**: ascending powers, unlike MATLAB `polyfit`.
- **Typo filename**: `include/msl/matrix/martrix_decompose.hpp`; prefer the
  public aggregator instead of renaming it during unrelated work.
- **Eigen bridge**: `include/msl/matrix/eigen_interface.hpp` provides maps and
  `matmul`, `matvec`, `solve` in `msl::matrix::eigen_interface`.
- **Language**: builds use C++20; clangd selects C++23 for indexing.
- **Format/lint**: `.clang-format` (LLVM, 80 columns, indent 4), `.clang-tidy`.
- **Header author**: Dong Feiyue `<FeiyueDong@outlook.com>`.

Keep algorithm contracts in the formal module documentation and update the
[compatibility guide](docs/compatibility.md) for breaking public API changes.
See the [developer guide](docs/developer-guide.md) for consumer examples and
validation procedures, and [third-party notices](NOTICE.md) for dependency
licensing.
