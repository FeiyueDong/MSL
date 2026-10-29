# MATLAB Validation

These scripts independently compare MSL outputs with analytic solutions and
MATLAB numerical routines. MATLAB is optional for library consumers and for
the C++ unit tests.

## Generate comparison data

From the repository root, configure and build MSL, then run the C++ tests:

```sh
xmake f -p linux -m release -y
xmake
xmake run test_matrix
xmake run test_signal
xmake run test_integral
xmake run test_difference
xmake run test_interp
xmake run test_polynomial
xmake run test_equation
xmake run test_ode
```

Use the appropriate platform option on other systems. The C++ tests generate
`test_result/<module>/matlab_compare/`; generated files are gitignored. Refresh
them after changing code or switching revisions so MATLAB does not compare
against stale results.

## Run MATLAB

Run all module scripts from the repository root:

```sh
matlab -batch "run('matlab/validation/run_all.m')"
```

Or run one module after its C++ test:

```sh
xmake run test_signal
matlab -batch "run('matlab/validation/signal/validate_signal.m')"
```

`run_all.m` reports per-module results and raises an error if a module fails.
A module can still report PASS while a toolbox-specific section is skipped:
Butterworth, Welch, and Savitzky-Golay comparisons require Signal Processing
Toolbox. Report those skips separately from checks that actually ran.

## Coverage

| Module script | Main comparisons |
|---------------|------------------|
| [matrix](matrix/validate_matrix.m) | Matrix product, LU permutation, tall/wide QR, SVD reconstruction, complex eigenpairs, pseudoinverse, inverse/determinant, truncated SVD |
| [signal](signal/validate_signal.m) | Vector and row/column FFT/IFFT, windows, Butterworth response, causal and zero-phase filtering, Welch PSD/CPSD, covariance, detrending |
| [difference](difference/validate_difference.m) | Gradients, second derivatives, analytic vector fields, Savitzky-Golay derivatives |
| [integral](integral/validate_integral.m) | Total and cumulative integration, polynomial antiderivatives |
| [interp](interp/validate_interp.m) | Nonlinear/non-uniform interpolation, spline, PCHIP, modified Akima |
| [polynomial](polynomial/validate_polynomial.m) | Polynomial fitting, coefficient ordering, evaluation |
| [equation](equation/validate_equation.m) | Roots against analytic values and `fzero` |
| [ode](ode/validate_ode.m) | Initial value solutions against analytic values and `ode45` |

Numerical conventions belong to the [module guides](../../docs/index.md).
In particular, MSL polynomial coefficients ascend in power, LU includes a row
permutation, and CPSD uses `X * conj(Y)`; the scripts account for these
conventions rather than comparing differently defined quantities.

The sensitive fourth-order `filtfilt` startup regression is additionally
checked inside `test_signal`: a residual test and a frozen 64-sample SciPy
reference are compiled into the executable. Those checks do not require a
Python or MATLAB runtime. They supplement, rather than replace, this MATLAB
workflow.
