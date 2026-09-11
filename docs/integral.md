# Integral Module

**Namespace**: `msl::integral`  
**Aggregator**: `msl/integral.hpp`  
**Files**: algorithm headers in `msl/integral/`

## Overview

The integral module provides numerical integration for 1D vectors and 2D matrices: trapezoidal rule (total and cumulative), Simpson's rule, Romberg integration, and adaptive Simpson quadrature.

## Trapz (Trapezoidal Rule)

Total integral using the trapezoidal rule.

```cpp
// Uniform spacing
double I = msl::integral::trapz(y, dx);

// Non-uniform spacing
double I = msl::integral::trapz(x, y);

// Matrix: integrate each column (uniform spacing)
auto integrals = msl::integral::trapz(mat, dx);  // returns vector<double>
msl::integral::trapz(mat, result, dx);           // zero-copy output

// Matrix: integrate each column (non-uniform spacing)
auto integrals = msl::integral::trapz(x, mat);

// Function
double I = msl::integral::trapz(f, a, b, intervals);
```

## Cumulative Trapz

Cumulative trapezoidal integration: `result[i] = ∫[0 to i] y dx`, with `result[0] = 0`.

```cpp
// Uniform spacing
auto cum = msl::integral::cumtrapz(y, dx);
msl::integral::cumtrapz(y, cum, dx);  // zero-copy

// Non-uniform spacing
auto cum = msl::integral::cumtrapz(x, y);

// Matrix: integrate each column
auto cum_mat = msl::integral::cumtrapz(mat, dx);
auto cum_mat = msl::integral::cumtrapz(x, mat);  // non-uniform
```

## Simpson's Rule

Higher accuracy than trapezoidal for smooth functions.

```cpp
// Total integral (requires at least 3 points)
// If even number of points, uses trapezoidal for last segment
double I = msl::integral::simpson(y, dx);

// Non-uniform spacing
double I = msl::integral::simpson(x, y);

// Function
double I = msl::integral::simpson(f, a, b, intervals);

// Cumulative Simpson
auto cum = msl::integral::cumsimpson(y, dx);
msl::integral::cumsimpson(y, cum, dx);  // zero-copy
auto cum = msl::integral::cumsimpson(x, y);

// Matrix: integrate each column
auto integrals = msl::integral::simpson(mat, dx);
auto integrals = msl::integral::simpson(x, mat);
```

Cumulative Simpson integrates the three-point interpolating parabola over the
two-interval span at even indices and over the preceding single interval at odd
indices, so every cumulative value is exact for quadratics. If the number of
points is even, the last interval falls back to the trapezoidal rule.

## Romberg Integration

Richardson extrapolation on the trapezoidal rule for high-precision integration.

```cpp
// Vector: y must have 2^k + 1 points
double I = msl::integral::romberg(y, dx, tol);  // tol default: 1e-10

// Uniformly spaced x/y samples
double I = msl::integral::romberg(x, y, tol);

// Function
double I = msl::integral::romberg(f, a, b, tol, max_iter);

// Matrix: integrate each column
auto integrals = msl::integral::romberg(mat, dx, tol);
auto integrals = msl::integral::romberg(x, mat, tol);
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `tol` | `1e-10` | Convergence tolerance |

## Adaptive Simpson Quadrature

Adaptive recursive Simpson's rule for arbitrary functions.

```cpp
double I = msl::integral::adaptive_simpson(f, a, b, tol, max_depth);
double I = msl::integral::quad(f, a, b, tol, max_depth);
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `f` | — | Callable compatible with `double(double)` |
| `a` | — | Lower limit |
| `b` | — | Upper limit (must be > a) |
| `tol` | `1e-8` | Tolerance |
| `max_depth` | `50` | Maximum recursion depth |

The function recursively subdivides the interval until the error estimate `|S_acb - S_ab| / 15 < tol`. Tolerance is halved at each recursion level.

`adaptive_simpson_quadrature.hpp` and `cumtrapz.hpp` are compatibility headers.
New code should include `adaptive_simpson.hpp` and `trapz.hpp` directly when
including algorithm-specific headers.
