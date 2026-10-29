# Utilities Module

**Namespaces**: `msl::utils` (data I/O), `msl` (assert handler)
**Files**: 2 headers in `msl/utils/`

## Overview

The utils module provides small utility functions used across the library: file I/O for numerical data and a custom assert handler.

## Data I/O (`data_io.hpp`)

Text files contain rows of whitespace-separated values; the in-memory vector
is column-major. Include `msl/utils/data_io.hpp` directly.

### Reading

```cpp
// Read data file (space/newline separated values)
auto data = msl::utils::ReadData("input.txt", num_cols, num_rows);
// Returns vector<double> in column-major order
```

### Writing

```cpp
// Write matrix data (column-major)
msl::utils::WriteData("output.txt", data, num_cols, num_rows);

// Write complex matrix data
msl::utils::WriteComplexData("output.txt", complex_data, num_cols, num_rows);
```

### JSON

```cpp
// Read entire JSON file as string
std::string json = msl::utils::ReadJson("config.json");
```

`ReadJson` reads text; it does not parse or validate JSON. If a file cannot be
opened, the readers log to `stderr` and return an empty vector or string;
writers log the failure and return. These legacy helpers do not throw an
exception on an open failure.

## Assert Handler (`msl_assert.hpp`)

Custom assertion handler that logs to stderr and a file, then calls `std::abort()`.

```cpp
msl::handle_assert(expr_str, message, file, line, function);
```

This is used internally by the MSL assertion mechanism. It writes to both `stderr` and `msl_error.log`, providing a record of assertion failures for debugging.
