# API Compatibility

This guide helps update callers written against earlier MSL interfaces. The
module documentation describes the current API; the replacements below use
that API directly, without compatibility aliases.

## Replaced names and return types

| Earlier call | Current API | Contract |
|--------------|-------------|----------|
| `A.get_row(i)` | `A.row_copy(i)` | Owning matrix copy of a row |
| `A.get_column(j)` | `A.column_copy(j)` | Owning matrix copy of a column |
| `A.submatrix(r0, r1, c0, c1)` | `A.submatrix_copy(r0, r1, c0, c1)` | Owning copy using half-open ranges |
| `matrix::adjoint(A)` for real matrices | `matrix::transpose(A)` | Ordinary transpose |
| `matrix::adjoint(A)` for complex matrices | `matrix::conjugate_transpose(A)` | Hermitian transpose, including rectangular matrices |
| `signal::psd(x, nfft)` | `signal::power_spectrum(x, nfft)` | Unnormalized single-frame power |
| `signal::cpsd(x, y, nfft)` | `signal::cross_power_spectrum(x, y, nfft)` | Unnormalized `X * conj(Y)` |
| `auto [L, U] = matrix::lu(A)` | `auto factors = matrix::lu(A)` | Permutation plus `L` and `U`, satisfying `P*A = L*U` |

The output-buffer overloads of the single-frame spectrum functions use the
same new names. `psd_welch` and `cpsd_welch` retain their names and density
normalization. The overload `power_spectrum(complex_fft_bins)` computes
`|X[k]|^2` directly from existing FFT bins.

`A.column(j)` still returns a borrowed `std::span`. Changing an owning copy
does not update `A`; changing a mutable column span does. See
[matrix ownership](matrix.md#real_matrix_owned-matrixd).

## LU and QR dimensions

LU returns `lu_result` for real matrices and `complex_lu_result` for complex
matrices. For both, `(P*A)(i, j) = A(factors.permutation[i], j)`. A manual
triangular solve must apply the same row permutation to the right-hand side.

For real economy QR, `k = min(m, n)`: `Q` is `m x k` and `R` is `k x n`.
In particular, a wide input no longer produces oversized economy factors.
Real full QR (`qr(A, true)`) and complex QR return `m x m` and `m x n` factors.
SVD reconstruction is `U*S*transpose(V)` for real matrices and
`U*S*conjugate_transpose(V)` for complex matrices.

## Numerical behavior to review

- **Butterworth band filters** use the inverse-bilinear image of the analog
  geometric center. Bandpass and bandstop have distinct pole transforms.
  Recompute output baselines produced by older band-filter coefficients.
- **Zero-phase filtering** solves the initial-state linear system directly
  with LU. The public overloads, odd reflection, and `3 * filter_order`
  padding remain the contract. Numerical outputs can change for poorly
  conditioned startup systems; the C++ suite includes a fixed SciPy oracle.
- **Cumulative Simpson integration** integrates a three-point parabola over
  each half interval. Interior cumulative points are not averages of adjacent
  accumulated integrals. With an even number of samples, the final interval
  uses a trapezoid.
- **Savitzky-Golay gradients** use local polynomial least-squares derivatives,
  including fitted boundary windows, rather than averaged finite differences.
- **Polynomial fitting** uses QR and a triangular solve. Too few samples throw
  `std::invalid_argument`; rank-deficient fits throw `std::runtime_error`.
- **Eigendecomposition** of real nonsymmetric matrices retains complex
  eigenpairs. A caller must not assume real-valued output.
- **FFT** zero-pads when a requested transform length exceeds the input.
  Welch segment length, FFT length, and sampling rate are separate parameters.

See [signal](signal.md), [integral](integral.md),
[difference](difference.md), [polynomial](polynomial.md), and
[matrix](matrix.md) for formulas, sizes, and examples.

## Errors and boundary handling

Matrix shape mismatches and invalid construction sizes throw
`std::invalid_argument`; copy-accessor and view ranges throw
`std::out_of_range`. `column()` performs bounds checking and is not
`noexcept`. Element indexing remains unchecked. Contiguous subviews support
whole columns; use an owning submatrix copy for an arbitrary rectangle.

Difference routines reject insufficient samples before allocating shortened
outputs. Spacing must be positive and non-uniform coordinates strictly
increasing. The non-uniform matrix central gradient uses the same weighted
three-point formula as the vector overload along either axis.

An interpolator must have data before evaluation, `in_range()`, or `range()`;
otherwise it throws `std::runtime_error`. Cubic spline derivatives at the
right endpoint use the last valid interval. Default cubic and Akima behavior
is described in the [interpolation guide](interp.md).

## Validate an update

Build and run the affected C++ tests, including numerical invariants and
invalid inputs. Optional MATLAB comparison scripts use freshly generated
test outputs; see the [developer guide](developer-guide.md#testing) and
[MATLAB validation guide](../matlab/validation/README.md). Keep runtime,
toolbox, platform, and remote-CI results separate when reporting validation.
