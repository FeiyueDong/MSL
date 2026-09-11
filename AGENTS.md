# MSL — Agent Guide

## Build

Two build systems coexist:

- **xmake** (cross-platform, default platform: mingw):
  `xmake f -p mingw && xmake` — builds the header-only `msl` target and the test
  executables. MSYS2/GCC on `PATH` is auto-detected; set `MSYS2_ROOT` only when
  the SDK lives elsewhere. Eigen is found through the compiler default include
  paths (override with `EIGEN_INCLUDE_DIR` if needed).
- **Visual Studio 2022** (Windows only): open `MSL.sln`. Uses v143 toolchain, C++20 (`stdcpp20`).

Single external dependency: Eigen (header-only). The VS solution uses the
vcpkg manifest (`vcpkg_installed/`); xmake uses system Eigen on MSYS2/Linux.

## Run Tests

```sh
xmake run test_matrix        # or test_signal, test_difference, test_integral, test_interp, test_polynomial, test_equation, test_ode
```

- No test framework — shared macros in `src/tests/test_utils.hpp`: `TEST_CASE`, `EXPECT_TRUE`, `EXPECT_EQ`, `EXPECT_NEAR`, `EXPECT_CPLX_NEAR`; call `msl::test::summary()` at the end.
- Tests write output to `test_result/<module>/` (gitignored).
- Prefer deterministic unit assertions over output-file-only checks. Legacy signal/interp/integral/difference tests may read sample files, but new coverage should use small in-memory fixtures where practical.
- MATLAB cross-validation lives in `matlab/validation/`; run everything with `matlab -batch "run('matlab/validation/run_all.m')"` after the C++ tests have refreshed `test_result/`.

## Structure

| Path | Purpose |
|---|---|
| `src/msl/<module>.hpp` | Aggregator headers (public API entrypoints) |
| `src/msl/<module>/` | Header-only implementations (C++20, `inline`) |
| `project/` | Visual Studio `.vcxproj` static lib targets |
| `src/tests/` | xmake test executables and `test_utils.hpp` |
| `matlab/` | MATLAB cross-validation scripts |

Namespaces mirror directories: `msl::matrix`, `msl::signal`, `msl::difference`, `msl::integral`, `msl::interp`, `msl::polynomial`, `msl::equation`, `msl::ode`, `msl::utils`.

**Column-major** storage (`data_[j * rows_ + i]`).

## Gotchas

- **Typo filename**: `src/msl/matrix/martrix_decompose.hpp` (not `matrix`).
- **C++20 for builds, C++23 for indexing**: xmake uses `c++20`; clangd and VS Code IntelliSense use `c++23`/`gnu++23`.
- **Eigen bridge**: `src/msl/matrix/eigen_interface.hpp` provides zero-copy `Eigen::Map` views and `matmul()`, `matvec()`, `solve()` wrappers. FFT uses `unsupported/Eigen/FFT`.
- **Runtime validation**: public matrix operations throw `std::invalid_argument`/`std::out_of_range`; only element access `operator()` stays unchecked.
- **Format**: `.clang-format` (LLVM-based, column 80, indent 4). **Lint**: `.clang-tidy` (bugprone, cppcoreguidelines, modernize, performance, readability).
- **Header template**: psi-header with author "Dong Feiyue" `<FeiyueDong@outlook.com>`.
- **License**: MIT (MSL) + MPL 2.0 (Eigen3 dependency). See `licenses/LICENSE`.
