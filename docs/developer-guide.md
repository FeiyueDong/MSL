# Developer Guide

## Build requirements

MSL uses C++20 and Xmake 3.0.3 or later. Its only external C++ dependency is
Eigen3, declared with `add_requires("eigen")` in
[`xmake/products.lua`](../xmake/products.lua). Xmake resolves that package for
the selected toolchain. Linux installations commonly provide Eigen through
`libeigen3-dev`.

## Standalone build

Run from the repository root and select the platform explicitly:

```sh
xmake f -p linux -m release -y    # Linux, GCC
xmake                          # msl plus all eight test executables
xmake run test_matrix
xmake run test_signal
```

Other supported configurations:

```sh
xmake f -p linux -m release --toolchain=clang -y
xmake f -p windows -a x64 -m release -y
xmake f -p mingw -m release -y
xmake f -p macosx -m release -y
```

Windows uses MSVC through Xmake. For MinGW, make the GCC toolchain available
to Xmake; set `MSYS2_ROOT` to its SDK directory when automatic discovery is
insufficient. The repository's build entry is `xmake.lua`.

Xmake's global package cache can be redirected with `XMAKE_GLOBALDIR` when a
separate writable cache is needed. Keep that value consistent between
configuration, build, and test commands; a new cache may need dependency
network access.

## Embed the header-only target

The standalone root imports both products and development tests. An embedding
project should import only `xmake/products.lua`, which defines the `msl`
header-only target and its Eigen dependency. It resolves paths from its own
location, so it can be included from another project directory.

For this consumer layout:

```text
consumer/
├── xmake.lua
├── main.cpp
└── vendor/
    └── msl/
        ├── include/
        └── xmake/products.lua
```

Use:

```lua
set_project("msl-consumer")
set_languages("c++20")
add_rules("mode.debug", "mode.release")

includes("vendor/msl/xmake/products.lua")

target("example")
    set_kind("binary")
    add_files("main.cpp")
    add_deps("msl")
```

A minimal `main.cpp` is:

```cpp
#include <vector>
#include <msl/matrix.hpp>
#include <msl/signal.hpp>

int main() {
    const auto identity = msl::matrix::matrixd::identity(2);
    const std::vector<double> input{1.0, 2.0};
    const auto values = msl::matrix::eigen_interface::matvec(identity, input);
    const auto spectrum = msl::signal::fft(values, 4);
    return spectrum.size() == 4 ? 0 : 1;
}
```

Configure and build from `consumer/`:

```sh
xmake f -P . -p linux -m release -y
xmake show -l targets
xmake build example
xmake run example
```

The target list contains `msl` and `example`; MSL's development tests are not
imported. The dependency publishes both `include/` and `include/msl/` along
with Eigen's include configuration. When building without Xmake, supply both
MSL include roots, Eigen's include path, and C++20 explicitly.

## Testing

The eight C++ test targets are:

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

[`tests/test_utils.hpp`](../tests/test_utils.hpp) provides `TEST_CASE`,
`EXPECT_TRUE`, `EXPECT_EQ`, `EXPECT_NEAR`, and `EXPECT_CPLX_NEAR`. Tests finish
with `msl::test::summary()` and return a nonzero exit code on failure.

Prefer small deterministic in-memory fixtures. Numerical tests should check
mathematical invariants as well as fixed reference values:

| Module | Representative checks |
|--------|-----------------------|
| Matrix | `P*A = L*U`, QR dimensions and orthogonality, SVD reconstruction, complex eigenpairs, pseudoinverse identities, truncated SVD accuracy |
| Signal | Row/column FFT round trips, zero padding, Butterworth gains and poles, Welch normalization, covariance lag order, `filtfilt` startup residual and frozen SciPy output |
| Difference | Non-uniform matrix/vector gradients, analytic vector fields, Savitzky-Golay derivatives, invalid spacing and coordinates |
| Integral | Cumulative Simpson values against polynomial antiderivatives, quadrature accuracy |
| Interpolation | Nonlinear data, spline endpoint derivatives, extrapolation and uninitialized use |
| Polynomial | Coefficient order, QR fit accuracy, insufficient samples and rank deficiency |
| Equation / ODE | Analytic solutions, convergence/status reporting, failure paths |

C++ tests also write comparison data to
`test_result/<module>/matlab_compare/` (gitignored). Run tests from the project
root through `xmake run` before using the
[MATLAB validation scripts](../matlab/validation/README.md). MATLAB and Python
are not required to run the C++ suite; frozen SciPy values are embedded in
`test_signal`.

The [CI workflow](../.github/workflows/ci.yml) defines Linux GCC, Linux Clang,
and Windows MSVC jobs that build and run all eight C++ targets. MATLAB checks
are a separate workflow run locally; configured CI jobs are not evidence of a
particular revision having passed remotely.

## Code conventions

- Public aggregators live in `include/msl/<module>.hpp`; implementation headers
  live in `include/msl/<module>/`. Functions defined in headers are `inline`
  where required by the one-definition rule.
- Use column-major storage and keep ownership explicit. Views and spans borrow
  data; allocating return values and `*_copy` accessors own it.
- Public operations use exceptions for invalid sizes and ranges. Preserve each
  module's documented numerical and diagnostic contracts.
- Prefer a direct linear solve to an explicit inverse when solving `A*x = b`.
- Document mathematical conventions, normalization, output shapes, boundaries,
  and exceptions next to the API.
- Follow `.clang-format` (LLVM, four-space indentation, 80 columns) and the
  checks in `.clang-tidy`. Format only files within the intended change.
- Preserve the psi-header author convention: Dong Feiyue
  `<FeiyueDong@outlook.com>`.

The decomposition header retains the filename `martrix_decompose.hpp`.
Use `msl/matrix.hpp` to include it through the public aggregator. Builds use
C++20; the checked-in clangd configuration selects C++23 for indexing.

## Directory structure

```text
include/msl/       Public aggregators and implementation headers
xmake/products.lua Header-only product and dependency declaration
xmake.lua          Standalone build including development tests
tests/             Eight test targets and shared assertion helpers
matlab/validation/ MATLAB comparison scripts and run_all.m
resource/          Reference data and legacy MATLAB utilities
docs/              Current API and developer documentation
licenses/          Bundled dependency license text
```

## Add a module

1. Add `include/msl/<module>.hpp` and implementation headers under
   `include/msl/<module>/`, in namespace `msl::<module>`.
2. Add a deterministic test in `tests/test_<module>/` with its `xmake.lua`.
3. Include that directory from `tests/xmake.lua` and add the test to CI.
4. Document the public API and link it from the module table.
5. If independent numerical validation is useful, add a script under
   `matlab/validation/<module>/` and register it in `run_all.m`.
