# MSL — Modern Scientific Library

MSL is a C++20 header-only library of general numerical algorithms. Its only
external C++ dependency is Eigen3, used for linear algebra and FFT. MATLAB is
used for optional cross-validation; selected C++ regressions contain frozen
SciPy reference values and require no Python installation.

## API conventions

- Matrices store columns contiguously: element `(i, j)` is at
  `data()[j * rows() + i]`.
- Owning containers and non-owning views are distinct. Matrix `*_copy`
  accessors allocate independent storage; `column()` and matrix views borrow
  existing storage and require the owner to remain alive.
- Many vector algorithms provide an allocating overload and an overload that
  writes to a caller-provided `std::span`. Supplying an output buffer does not
  imply that the algorithm uses no temporary allocations.
- Public high-level operations validate shapes and required inputs. Element
  indexing remains unchecked; individual modules document their error rules.
- Polynomial coefficients use ascending powers. Filter coefficients describe
  ascending powers of `z^-1`. Filter cutoffs are normalized to Nyquist.

## Modules

| Module | Namespace | Public header | Purpose |
|--------|-----------|---------------|---------|
| [Matrix](matrix.md) | `msl::matrix` | `msl/matrix.hpp` | Real/complex matrices, views, decompositions, Eigen interop |
| [Signal](signal.md) | `msl::signal` | `msl/signal.hpp` | FFT/IFFT, Butterworth and Fourier filters, spectra, covariance, detrending |
| [Difference](difference.md) | `msl::difference` | `msl/difference.hpp` | Differences, gradients, vector calculus, Savitzky-Golay derivatives |
| [Integral](integral.md) | `msl::integral` | `msl/integral.hpp` | Trapezoidal, Simpson, Romberg, adaptive quadrature |
| [Interpolation](interp.md) | `msl::interp` | `msl/interp.hpp` | Linear, cubic spline, Akima, PCHIP, nearest, polynomial interpolation |
| [Polynomial](polynomial.md) | `msl::polynomial` | `msl/polynomial.hpp` | Least-squares fitting, evaluation, derivatives |
| [Equation](equation.md) | `msl::equation` | `msl/equation.hpp` | Scalar nonlinear root finding with diagnostics |
| [ODE](ode.md) | `msl::ode` | `msl/ode.hpp` | Fixed-step and adaptive initial value solvers |
| [Utilities](utils.md) | `msl::utils`, `msl` | Individual headers | Text data I/O and assertion handler |

## Quick start

Save this example as `main.cpp`:

```cpp
#include <iostream>
#include <vector>

#include <msl/difference.hpp>
#include <msl/equation.hpp>
#include <msl/integral.hpp>
#include <msl/interp.hpp>
#include <msl/matrix.hpp>
#include <msl/ode.hpp>
#include <msl/polynomial.hpp>
#include <msl/signal.hpp>

int main() {
    msl::matrix::matrixd A(2, 2, {3.0, 1.0, 1.0, 2.0});
    const std::vector<double> b{9.0, 8.0};
    const auto solution = msl::matrix::eigen_interface::solve(A, b);
    std::cout << solution[0] << ", " << solution[1] << '\n'; // 2, 3

    const std::vector<double> x{0.0, 1.0, 2.0, 3.0};
    const std::vector<double> y{0.0, 1.0, 4.0, 9.0};
    const auto gradient = msl::difference::central_gradient(x, y);
    const double area = msl::integral::trapz(x, y);
    std::cout << gradient[1] << ", " << area << '\n'; // 2, 9.5

    const auto interpolator = msl::interp::Linear::from_data(x, y);
    const auto coefficients = msl::polynomial::polyfit(x, y, 2);
    std::cout << interpolator(1.5) << ", "
              << msl::polynomial::polyval(coefficients, 1.5) << '\n';
    // 2.5, 2.25

    const auto spectrum = msl::signal::fft(y, 8); // zero-padded FFT
    std::cout << spectrum.size() << '\n'; // 8

    const double root = msl::equation::brent_root(
        [](double value) { return value * value - 2.0; }, 0.0, 2.0);
    const auto final_state = msl::ode::rk4_final(
        [](double, const msl::ode::state &state) {
            return msl::ode::state{state[0]};
        },
        std::vector<double>{1.0}, 0.0, 1.0, 0.01);
    std::cout << root << ", " << final_state[0] << '\n'; // sqrt(2), exp(1)
}
```

For a direct GCC build from the repository root with system Eigen on Linux:

```sh
c++ -std=c++20 -O2 -Iinclude -Iinclude/msl -I/usr/include/eigen3 main.cpp -o example
./example
```

Both MSL include roots are needed by the current headers. Xmake propagates
these and the Eigen dependency through the `msl` target; see the
[developer guide](developer-guide.md) for standalone and embedded builds.

## Further reading

- [API compatibility](compatibility.md): replaced names and numerical behavior
  to consider when updating existing callers.
- [Testing](developer-guide.md#testing): unit assertions and numerical invariants.
- [MATLAB validation](../matlab/validation/README.md): comparison data, runner,
  and optional toolbox checks.
- [Third-party notices](../NOTICE.md) and [Eigen license text](../licenses/LICENSE).
