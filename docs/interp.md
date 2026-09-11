# Interpolation Module

**Namespace**: `msl::interp`  
**Aggregator**: `msl/interp.hpp`  
**Files**: 7 headers in `msl/interp/`

## Overview

The interpolation module provides 1D interpolation with six methods: nearest neighbor, linear, cubic C2 spline, Akima, PCHIP (shape-preserving), and polynomial interpolation. All interpolators inherit from a shared base class with extensible extrapolation modes.

## Interpolator Hierarchy

```
InterpolatorBase (abstract)
  |-- Linear       : piecewise linear interpolation
  |-- CubicSpline  : cubic C2 spline (Not-a-knot boundary by default)
  |-- AkimaSpline  : modified Akima / MAKIMA by default
  |-- PchipSpline  : shape-preserving piecewise cubic Hermite
  |-- Near         : nearest-neighbor (stepwise constant)
  |-- Polynomial   : polynomial interpolation
```

## Base Class (InterpolatorBase)

All interpolators share a common API:

```cpp
// Set data points
interp.set_data(x, y);

// Single point evaluation
double val = interp(3.14);

// Multiple point evaluation (allocating)
auto vals = interp(x_new);

// Multiple point evaluation (zero-copy)
interp(x_new, result);

// Extrapolation control
interp.set_extrapolation(ExtrapolationMode::Linear);
auto mode = interp.extrapolation_mode();

// Data access
auto x_data = interp.x_data();
auto y_data = interp.y_data();
size_t n = interp.size();
bool in = interp.in_range(x);
auto [xmin, xmax] = interp.range();
```

Evaluating an interpolator before `set_data()` (or `from_data()`) has been
called throws `std::runtime_error`; `in_range()` and `range()` do the same
instead of dereferencing empty storage.

## Extrapolation Modes

| Mode | Behavior |
|------|----------|
| `None` | Throw `std::out_of_range` |
| `Constant` | Return boundary value |
| `Linear` | Extrapolate using boundary slope |
| `Nearest` | Return nearest boundary point |
| `Periodic` | Periodic extension of the domain |
| `Polynomial` | Default; delegates to derived class for boundary polynomial |

## Interpolation Methods

### Linear

Piecewise linear interpolation between consecutive data points.

```cpp
msl::interp::Linear interp;
interp.set_data({0.0, 1.0, 2.0}, {0.0, 1.0, 4.0});
double val = interp(1.5);  // ≈ 2.5
```

### Spline

Cubic C2-continuous spline. The default boundary condition is Not-a-knot,
matching MATLAB `interp1(..., 'spline')`; Natural and Clamped boundaries are
available through `CubicSpline::BoundaryCondition`. The tridiagonal system is
solved in O(n).

```cpp
msl::interp::CubicSpline interp;
interp.set_data(x, y);
double val = interp(3.14);
```

### Akima

Akima interpolation — a piecewise cubic Hermite method with less oscillation
than standard splines for irregular data. The default uses modified Akima
weights, matching MATLAB `makima`; pass `modified_akima = false` to
`AkimaSpline::from_data` for the original Akima scheme.

```cpp
msl::interp::AkimaSpline interp;
interp.set_data(x, y);
```

### PCHIP

Piecewise Cubic Hermite Interpolating Polynomial — shape-preserving. Preserves monotonicity and local shape of the data.

```cpp
msl::interp::PchipSpline interp;
interp.set_data(x, y);
```

### Nearest

Nearest-neighbor (stepwise constant) interpolation. Returns the value of the closest data point.

```cpp
msl::interp::Near interp;
interp.set_data(x, y);
```

### Polynomial Interpolation

Interpolates all data points with a single polynomial. Uses Neville's algorithm for evaluation.

```cpp
msl::interp::Polynomial interp;
interp.set_data(x, y);  // n points → degree n-1 polynomial
```

**Note**: High-degree polynomial interpolation can oscillate wildly for non-uniformly spaced data (Runge's phenomenon). Consider spline or PCHIP for smooth interpolation of many points.

## Batch Evaluation

All interpolators support efficient batch evaluation:

```cpp
// Allocating
auto y_new = interp(x_new);

// Zero-copy
std::vector<double> result(x_new.size());
interp(x_new, result);
```

The batch evaluation iterates and calls `operator()` for each point. The base class handles extrapolation dispatch before calling the derived `interpolate()` method.
