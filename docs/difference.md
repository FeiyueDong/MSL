# Difference Module

**Namespace**: `msl::difference`  
**Aggregator**: `msl/difference.hpp`  
**File**: `msl/difference/difference.hpp`

## Overview

The difference module provides numerical differentiation for 1D vectors and 2D matrices: first-order differences, forward/central gradients, second derivatives, Laplacian, divergence, curl, and Savitzky-Golay smoothed derivatives.

Vector routines provide output-span and allocating overloads; matrix routines
return owning matrices. Output sizes depend on the selected operation.

## Input validation

`diff`, forward gradients, and central gradients require at least two samples
along the differentiated axis. Second derivatives require at least three.
Uniform spacing must be positive; non-uniform coordinates must be strictly
increasing and match the selected axis length. Invalid sizes, coordinates,
spacing, or axes throw `std::invalid_argument`. Shortened-output overloads
validate the input length before allocating.

## Difference (First-Order)

```cpp
// Vector: diff[i] = y[i+1] - y[i], output size = n-1
auto dy = msl::difference::diff(y);
msl::difference::diff(y, result);  // zero-copy

// Matrix: diff along axis
auto dy = msl::difference::diff(mat, 0);  // axis=0: vertical diff (rows)
auto dy = msl::difference::diff(mat, 1);  // axis=1: horizontal diff (cols)
```

## Gradient (Forward Difference)

Forward difference approximation: `grad[i] = (y[i+1] - y[i]) / dx`. Output size is `n-1`.

```cpp
// Uniform spacing
auto grad = msl::difference::forward_gradient(y, dx);
msl::difference::forward_gradient(y, grad, dx);  // zero-copy

// Non-uniform spacing
auto grad = msl::difference::forward_gradient(y, x);

// Matrix
auto grad = msl::difference::forward_gradient(mat, dx, axis);
auto grad = msl::difference::forward_gradient(mat, x, axis);
```

## Gradient (Central Difference)

Central difference with forward/backward at boundaries. Output size equals input size.

```cpp
// Uniform: forward at start, central at interior, backward at end
auto grad = msl::difference::central_gradient(y, dx);

// Non-uniform: weighted central differences
auto grad = msl::difference::central_gradient(x, y);

// Matrix
auto grad = msl::difference::central_gradient(mat, dx, axis);
auto grad = msl::difference::central_gradient(mat, x, axis);
```

**Boundary treatment**:
- First point: forward difference `(y[1] - y[0]) / dx`
- Interior points: central difference `(y[i+1] - y[i-1]) / (2*dx)`
- Last point: backward difference `(y[n-1] - y[n-2]) / dx`

For non-uniform spacing, interior points use weighted central differences:
```
grad[i] = (dx_left * (y[i+1] - y[i]) / dx_right
         + dx_right * (y[i] - y[i-1]) / dx_left)
         / (dx_left + dx_right)
```

The matrix overload uses the same formula for both axes: `axis=0`
differentiates down rows using one coordinate per row; `axis=1` differentiates
across columns using one coordinate per column. Interior derivatives are exact
for a quadratic sampled on a non-uniform grid.

For vector inputs, the argument order is `forward_gradient(y, x)` but
`central_gradient(x, y)`; keep this distinction when changing methods.

## 2D Gradient

Returns both horizontal and vertical gradients as a pair.

```cpp
auto [grad_x, grad_y] = msl::difference::central_gradient2d(mat, dx, dy);
// grad_x: horizontal derivative ∂f/∂x
// grad_y: vertical derivative ∂f/∂y
```

## Second Derivative

Central difference approximation: `d²y/dx² ≈ (y[i+1] - 2*y[i] + y[i-1]) / dx²`. Output size equals input size.

```cpp
auto d2y = msl::difference::central_gradient2(y, dx);
msl::difference::central_gradient2(y, d2y, dx);  // zero-copy
```
Requires at least 3 points.

## Laplacian

2D Laplacian operator: `∇²f = ∂²f/∂x² + ∂²f/∂y²`. Boundaries are set to zero.

```cpp
auto lap = msl::difference::laplacian(mat, dx, dy);
```
Requires at least 3×3 matrix.

## Divergence

2D divergence of a vector field: `div(F) = ∂Fx/∂x + ∂Fy/∂y`.

```cpp
auto div = msl::difference::divergence(Fx, Fy, dx, dy);
```
`Fx` and `Fy` must have the same dimensions.

## Curl

2D curl (z-component): `curl(F) = ∂Fy/∂x - ∂Fx/∂y`.

```cpp
auto curl_z = msl::difference::curl(Fx, Fy, dx, dy);
```

## Savitzky-Golay Gradient

Smoothed first derivative using a local polynomial least-squares fit. Interior
points use centered derivative weights; endpoints differentiate a polynomial
fit to the first or last complete window. This is a polynomial derivative,
rather than an average of adjacent finite differences.

```cpp
auto grad = msl::difference::savgol_gradient(y, window_size, poly_order, dx);
msl::difference::savgol_gradient(y, grad, window_size, poly_order, dx);  // zero-copy
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `window_size` | `5` | Odd window size, at least 5 and no larger than the input |
| `poly_order` | `2` | Non-negative polynomial order, less than `window_size` |
| `dx` | `1.0` | Spacing |
