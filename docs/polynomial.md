# Polynomial Module

**Namespace**: `msl::polynomial`  
**Aggregator**: `msl/polynomial.hpp`  
**File**: `msl/polynomial/polynomial.hpp`

## Overview

The polynomial module provides a `Polynomial` class with least-squares fitting (via QR decomposition), Horner evaluation, derivatives, and convenience `polyfit`/`polyval` free functions.

Coefficients are stored in **ascending power order**:
```
coeffs_[0] + coeffs_[1] * x + coeffs_[2] * x^2 + ... + coeffs_[n] * x^n
```

## Polynomial Class

### Construction

```cpp
// From coefficients
msl::polynomial::Polynomial p({1.0, 2.0, 3.0});  // 1 + 2x + 3x^2

// Fit to data (implicit x = 0, 1, 2, ...)
msl::polynomial::Polynomial p(y, degree);

// Fit to data (explicit x values)
msl::polynomial::Polynomial p(x, y, degree);
```

### Evaluation

```cpp
// Single point (Horner's method)
double val = p(2.5);

// Multiple points (allocating)
auto vals = p(x_values);

// Multiple points (zero-copy)
p(x_values, result);
```

### Fitting

```cpp
// Fit to new data
p.fit(x, y, degree);

// Set coefficients directly
p.set_coefficients({1.0, 2.0, 3.0});
```

### Accessors

```cpp
auto coeffs = p.coefficients();     // vector of coefficients
size_t deg = p.degree();            // polynomial degree
auto deriv = p.derivative(2.5);     // first derivative at x
```

## Free Functions

### polyfit

Least-squares polynomial fitting.

```cpp
// Allocating
auto coeffs = msl::polynomial::polyfit(x, y, n);

// Zero-copy (result must have size n+1)
msl::polynomial::polyfit(x, y, result, n);
```

### polyval

Polynomial evaluation.

```cpp
// Single point
double val = msl::polynomial::polyval(coeffs, 2.5);

// Multiple points (allocating)
auto vals = msl::polynomial::polyval(coeffs, x_values);

// Multiple points (zero-copy)
msl::polynomial::polyval(coeffs, x_values, result);
```

## Algorithm Details

### Least-Squares Fitting

1. Build Vandermonde matrix `A` where `A[i][j] = x[i]^j` (size m × (n+1))
2. Compute economy QR decomposition of A
3. Solve `R * c = Q^T * y` via back-substitution (no explicit inverse)
4. Returns coefficients in ascending power order

The fit rejects rank-deficient systems (for example duplicate sample
locations) with `std::runtime_error`, and throws `std::invalid_argument` when
fewer than `n + 1` samples are provided.

### Evaluation

Horner's method (nested multiplication):
```
result = coeffs[n]
for i = n-1 down to 0:
    result = result * x + coeffs[i]
```

### Derivative

Horner-based derivative evaluation:
```
result = n * coeffs[n]
for i = n-1 down to 1:
    result = result * x + i * coeffs[i]
```
