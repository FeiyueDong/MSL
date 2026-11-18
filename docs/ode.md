# ODE Module

**Namespace**: `msl::ode`  
**Aggregator**: `msl/ode.hpp`

## Overview

The ODE module solves initial value problems

```text
y' = f(t, y),  y(t0) = y0
```

with `std::vector<double>` states:

```cpp
using msl::ode::state; // std::vector<double>
```

The right-hand side function has the form:

```cpp
auto f = [](double t, const msl::ode::state& y) {
    return msl::ode::state{y[0]};
};
```

## Result Type

```cpp
enum class ode_status {
    success,
    max_steps,
    step_size_underflow,
    non_finite_value
};

struct ode_result {
    std::vector<double> t;
    std::vector<state> y;
    size_t steps;
    size_t rejected_steps;
    ode_status status;

    bool success() const;
    explicit operator bool() const;
    const state& final_state() const;
    double final_time() const;
    const state& value() const; // throws if not successful
};
```

`step_size_underflow` also covers a step too small to advance the current
time in double precision. The solver returns before evaluating the
right-hand side or counting that step.

## Options

```cpp
struct ode_options {
    double rtol = 1e-6;
    double atol = 1e-9;
    double initial_step = 0.0;
    double min_step = 1e-12;
    double max_step = 0.0;
    size_t max_steps = 100000;
};
```

Fixed-step solvers mainly use `max_steps`. `ode45` uses all tolerance and
step-size options.

## Fixed-Step Solvers

Each fixed-step method has three forms:

```cpp
auto r = msl::ode::rk4(f, y0, t0, t1, dt);
auto y_end = msl::ode::rk4_final(f, y0, t0, t1, dt);
auto r_eval = msl::ode::rk4_eval(f, y0, t_eval, dt);
```

The `*_eval` form treats `t_eval.front()` as the initial time and returns states
exactly at the requested times. The last internal step before each requested
time is shortened; no interpolation is used.

Available fixed-step methods:

```cpp
msl::ode::euler(...)
msl::ode::heun(...)
msl::ode::rk4(...)
```

## ODE45

`ode45` uses the Dormand-Prince 5(4) embedded Runge-Kutta method with adaptive
step-size control.

```cpp
msl::ode::ode_options opt;
opt.rtol = 1e-8;
opt.atol = 1e-10;

auto r = msl::ode::ode45(f, y0, 0.0, 10.0, opt);
auto y_end = msl::ode::ode45_final(f, y0, 0.0, 10.0, opt);
auto r_eval = msl::ode::ode45_eval(f, y0, t_eval, opt);
```

The current `ode45_eval` implementation lands exactly on each requested output
time by shortening adaptive steps. Dense-output interpolation can be added later
if intermediate output should avoid restarting at each interval.

## Example

```cpp
auto f = [](double, const msl::ode::state& y) {
    return msl::ode::state{y[0]};
};

std::vector<double> y0{1.0};
auto r = msl::ode::rk4(f, y0, 0.0, 1.0, 0.01);

if (r) {
    double y1 = r.final_state()[0];
}
```

## MATLAB Validation

Run the C++ test first to generate comparison data:

```sh
xmake run test_ode
```

Then run:

```sh
matlab -batch "run('matlab/validation/ode/validate_ode.m')"
```

The script compares MSL results against analytic solutions and MATLAB `ode45`.
