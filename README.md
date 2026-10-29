# MSL — Modern Scientific Library

MSL is a C++20 header-only numerical library for real and complex matrices,
signal processing, differentiation, integration, interpolation, polynomial
fitting, scalar equations, and ordinary differential equations. Eigen3 is its
only external C++ dependency.

- [Overview and quick start](docs/index.md)
- [Build, embedding, and testing](docs/developer-guide.md)
- [API compatibility](docs/compatibility.md)
- [MATLAB cross-validation](matlab/validation/README.md)

## Build and test

With Xmake 3.0.3 or later and Eigen3 available, run from the repository root:

```sh
xmake f -p linux -m release -y
xmake
xmake run test_matrix
xmake run test_signal
```

Use `-p windows` for MSVC, `-p mingw` for MinGW, or `-p macosx` for macOS.
The [developer guide](docs/developer-guide.md) lists all eight test targets
and explains how to embed only the header-only `msl` target in another project.

## Dependencies and licensing

MSL is distributed under the MIT license. Eigen3 supplies linear algebra and
FFT support. Its MPL 2.0 license text is included in
[licenses/LICENSE](licenses/LICENSE); see also the
[third-party notices](NOTICE.md).
