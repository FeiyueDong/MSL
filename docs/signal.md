# Signal Module

**Namespace**: `msl::signal`  
**Aggregator**: `msl/signal.hpp`
**Files**: 10 headers in `include/msl/signal/`

## Overview

The signal module provides vector and matrix FFT/IFFT, Butterworth IIR filter
design, zero-phase forward-backward filtering, Fourier-domain filtering,
spectra and spectral densities, unbiased cross-covariance, windows, and
detrending.

`#include <msl/signal.hpp>` includes every header listed below, including
`filter.hpp`, `filtfilt.hpp`, and `cross_covariance.hpp`. Individual headers can
also be included directly.

## File Listing

| File | Contents |
|------|----------|
| `fft.hpp` | 1D/2D FFT, IFFT, frequency bins, power/magnitude/phase spectra |
| `butterworth_filter.hpp` | Butterworth IIR filter design and application |
| `filter.hpp` | Direct Form II Transposed IIR filter implementation |
| `filtfilt.hpp` | Zero-phase forward-backward filtering |
| `filter_design.hpp` | Filter type enums and coefficient struct |
| `fourier_domain_filter.hpp` | Frequency-domain filtering with window shaping |
| `power_spectral_density.hpp` | Single-frame power spectra and Welch PSD/CPSD |
| `cross_covariance.hpp` | Mean-removed unbiased real cross-covariance |
| `window.hpp` | Hamming, Hann, Blackman, rectangular window functions |
| `detrend.hpp` | Polynomial trend removal |

## FFT (`fft.hpp`)

### 1D FFT

```cpp
// Real → Complex (convenience)
auto X = msl::signal::fft(input);            // default length = input.size()
auto X = msl::signal::fft(input, 1024);      // zero-padded to 1024

// Real → Complex (zero-copy)
std::vector<std::complex<double>> Xbuf(1024);
msl::signal::fft(input, Xbuf, 1024);

// Complex → Complex
auto X = msl::signal::fft(complex_input);
```

For real and complex input, `nfft=0` selects the input length. A larger
transform length appends zeros; a smaller positive length transforms the first
`nfft` input samples. The caller-provided output must contain exactly `nfft`
elements (or the input length when `nfft=0`).

### 1D IFFT

```cpp
auto x = msl::signal::ifft(X);               // complex → complex
auto x = msl::signal::ifft_real(X);          // complex → real (conjugate symmetric)
```

### Matrix FFT (Column-wise / Row-wise)

```cpp
// Column-wise: each column is a signal
auto X = msl::signal::fft_columns(mat);      // matrix → matrixc

// Row-wise: each row is a signal
auto X = msl::signal::fft_rows(mat);         // matrix → matrixc

auto x = msl::signal::ifft_columns(X);       // inverse
auto x = msl::signal::ifft_rows_real(X);     // inverse to real
```

These apply independent one-dimensional transforms along the selected axis;
they are not a full two-dimensional transform. Row-wise inverse transforms
write the reconstructed values into the output matrix. Optional `nfft`
arguments control the transformed row or column length.

### Utility

```cpp
auto freqs = msl::signal::fft_frequencies(n, sample_rate); // frequency bins
auto power = msl::signal::power_spectrum(fft_result);      // |X[k]|^2
auto magn  = msl::signal::magnitude_spectrum(fft_result);  // |X[k]|
auto phase = msl::signal::phase_spectrum(fft_result);      // arg(X[k])
```

## Butterworth Filter (`butterworth_filter.hpp`)

Butterworth filters have maximally flat passband response. Design uses a
bilinear transformation from an analog prototype. Cutoffs are normalized to
Nyquist: `fc = frequency_hz / (sampling_rate_hz / 2)`, with `0 < fc < 1` and
`fc_low < fc_high` for a band filter. The requested order is the lowpass
prototype order; bandpass and bandstop transfer functions have twice that
order.

### Class API

```cpp
// Lowpass/highpass
auto filt = msl::signal::ButterworthFilter::lowpass(order, fc);
auto highpass = msl::signal::ButterworthFilter::highpass(order, fc);

// Bandpass/bandstop
auto bandpass =
    msl::signal::ButterworthFilter::bandpass(order, fc_low, fc_high);
auto bandstop =
    msl::signal::ButterworthFilter::bandstop(order, fc_low, fc_high);

// Get coefficients
auto coeffs = filt.coefficients();  // FilterCoefficients {b, a}
auto b = filt.b();                  // numerator
auto a = filt.a();                  // denominator

// Redesign with new parameters
filt.redesign(order, fc, msl::signal::FilterType::highpass);
```

### Convenience: Design-Only Functions

```cpp
auto coeffs = msl::signal::butterworth_lowpass_design(order, fc);
auto coeffs = msl::signal::butterworth_highpass_design(order, fc);
auto coeffs = msl::signal::butterworth_bandpass_design(order, fc_low, fc_high);
auto coeffs = msl::signal::butterworth_bandstop_design(order, fc_low, fc_high);
```

### Convenience: Design-and-Apply Functions

```cpp
// Apply with zero-phase filtering (default, recommended)
auto filtered = msl::signal::butterworth_lowpass(signal, order, fc);

// Causal filtering
auto filtered = msl::signal::butterworth_lowpass(signal, order, fc, false);

// Zero-copy output
msl::signal::butterworth_lowpass(signal, result, order, fc);

// Available: butterworth_lowpass, _highpass, _bandpass, _bandstop
```

### Band-filter normalization

With normalized cutoffs `f1` and `f2`, the prewarped analog cutoffs are
`Omega_i = 2*tan(pi*f_i/2)`. Their geometric center
`Omega_0 = sqrt(Omega_1*Omega_2)` maps back to the digital frequency
`f0 = (2/pi)*atan(Omega_0/2)`.

Bandpass gain is normalized at `f0`. Bandstop uses its own lowpass-to-bandstop
pole transform and places conjugate zeros at the mapped notch frequency; its
passband is normalized at DC. The center is not the arithmetic mean of the
digital cutoffs. Design checks the realness of the coefficient polynomial and
rejects an unusable normalization gain rather than silently returning
unnormalized coefficients.

## Filter Application (`filter.hpp`)

Causal Direct Form II Transposed IIR filtering with zero initial state.
The filter normalizes numerator and denominator by `a[0]` when needed.

```cpp
auto filtered = msl::signal::filter(signal, coeffs);  // convenience

// Zero-copy
msl::signal::filter(signal, result, coeffs);

// Apply to matrix columns
auto filtered_mat = msl::signal::filter_columns(signal_matrix, coeffs);
```

## Zero-Phase Filtering (`filtfilt.hpp`)

Applies the filter forward and backward to obtain a zero-phase result with
the squared magnitude response. The transfer-function order is doubled by
the two passes. The examples use normalized coefficients (`a[0] == 1`) from
the Butterworth designer.

```cpp
auto filtered = msl::signal::filtfilt(signal, coeffs);  // convenience

// Zero-copy
msl::signal::filtfilt(signal, result, coeffs);

// Apply to matrix columns
auto filtered_mat = msl::signal::filtfilt_columns(signal_matrix, coeffs);
```

Let `p = max(a.size(), b.size()) - 1`. The input length must exceed `3*p`;
otherwise `filtfilt` throws `std::invalid_argument`. Both ends use odd
reflection with `3*p` samples. The output span must match the input length.

The internal initial state solves the steady-state equation `(I-A)*zi = B`,
where `A` is the transposed direct-form state matrix and
`B[i] = b[i+1] - b[0]*a[i+1]`. It uses
`matrix::eigen_interface::solve` (LU) directly, avoiding an explicit inverse
that can amplify startup rounding error in a poorly conditioned system.
Each pass scales this state by its first padded input value.

For a fourth-order Butterworth prototype at 50 Hz with a 0.1-20 Hz passband:

```cpp
const double sampling_rate = 50.0;
const auto coeffs = msl::signal::butterworth_bandpass_design(
    4, 0.1 / (sampling_rate / 2.0), 20.0 / (sampling_rate / 2.0));
const auto filtered = msl::signal::filtfilt(input, coeffs);
```

The bandpass transfer-function order is 8, so its padding length is 24.
When comparing the same input and coefficients with SciPy, specify
`scipy.signal.filtfilt(b, a, input, padtype="odd", padlen=24)` explicitly;
the default padding length differs. The C++ regressions check a startup
relative residual below `1e-12` and all 64 samples of a frozen SciPy reference
within `1e-9`. These are test tolerances for that fixture, rather than a
universal accuracy guarantee for arbitrary IIR coefficients.

## Filter Design Types (`filter_design.hpp`)

```cpp
enum class FilterType { lowpass, highpass, bandpass, bandstop };
enum class FilterMethod { butterworth /*, chebyshev1, chebyshev2, elliptic */ };

struct FilterCoefficients {
    std::vector<double> b;  // numerator
    std::vector<double> a;  // denominator (a[0] = 1.0)
    size_t order() const;   // max(a.size(), b.size()) - 1
};
```

## Fourier Domain Filter (`fourier_domain_filter.hpp`)

Filters signals by multiplying with a designed frequency response in the Fourier domain. Supports smooth transition bands using window functions.

```cpp
// Class API
auto filt = msl::signal::FourierDomainFilter::bandpass(
    fc_low, fc_high,
    msl::signal::FourierDomainFilter::WindowType::hamming, transition_band);

auto filtered = filt.apply(signal);

// Convenience functions
auto bandpassed = msl::signal::fourier_bandpass(signal, fc_low, fc_high);
auto lowpassed = msl::signal::fourier_lowpass(signal, fc_high);
auto highpassed = msl::signal::fourier_highpass(signal, fc_low);
auto bandstopped = msl::signal::fourier_bandstop(signal, fc_low, fc_high);
```

Window types: `rectangular`, `hamming`, `hanning`, `blackman`. Transition band creates smooth rolloff between passband and stopband.

## Power Spectral Density (`power_spectral_density.hpp`)

### Welch's Method

```cpp
// PSD (auto-power)
auto psd = msl::signal::psd_welch(x);     // default: Hann window, 1024 pts, 50% overlap
auto psd = msl::signal::psd_welch(x, window, noverlap, segment_length);

// Cross-PSD
auto cpsd = msl::signal::cpsd_welch(
    x, y, window, noverlap, segment_length);

// 256-sample segments, 2048-point FFT, density in power/Hz at 100 Hz
auto cpsd = msl::signal::cpsd_welch(
    x, y, window, noverlap, 256, 2048, 100.0);
auto frequencies = msl::signal::fft_frequencies(2048, 100.0);

// Single-frame spectra (unnormalized)
auto power = msl::signal::power_spectrum(x, nfft);          // |X|^2
auto cross = msl::signal::cross_power_spectrum(x, y, nfft); // X * conj(Y)
```

When `nfft` exceeds the input length, single-frame functions zero-pad the
input before the FFT.

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `window` | `hann_window(1024)` | Window coefficients |
| `noverlap` | `512` | Overlap samples (default: 50%) |
| `segment_length` | `1024` | Number of windowed samples per segment |
| `nfft` | `0` | FFT length; `0` selects `segment_length`, otherwise it must be at least the segment length |
| `sampling_rate` | `1.0` | Samples per second used for power/Hz density scaling |

Welch functions return the full two-sided `nfft` spectrum using
`X * conj(Y)`. The normalization is the segment average divided by
`sampling_rate * sum(window^2)`. A MATLAB-style one-sided spectrum for real
signals can be formed from bins `0..nfft/2`, doubling interior positive
frequency bins but not DC or the Nyquist bin.

That Nyquist exception applies to even `nfft`; for odd lengths every retained
positive-frequency bin is doubled. Only complete segments contribute:
`1 + floor((N - segment_length) / (segment_length - noverlap))`. The input
must contain at least one full segment, `window.size()` must equal the
segment length, window energy must be positive, and
`0 <= noverlap < segment_length`. Segments are not automatically detrended.

Cross-spectrum conventions vary between packages. A reference using
`conj(X) * Y` needs conjugation or swapped inputs to match MSL's
`X * conj(Y)` convention; compare complex values as well as magnitudes.

`psd_welch` and `cpsd_welch` are density-normalized (power/Hz). The
single-frame `power_spectrum` and `cross_power_spectrum` are unnormalized
(`|X|^2` and `X * conj(Y)`), so they are deliberately not named PSD.

## Unbiased Cross-Covariance (`cross_covariance.hpp`)

```cpp
auto covariance = msl::signal::xcov_unbiased(x, y, max_lag);
```

The result is ordered from `-max_lag` through `+max_lag`. Both means are
removed, and each lag is divided by `N - abs(lag)`, matching MATLAB
`xcov(x, y, max_lag, "unbiased")` for equal-length real inputs.

For lag `k`, MSL sums `(x[n+k] - mean(x)) * (y[n] - mean(y))` over valid
overlapping indices. Inputs must be nonempty and equally sized, with
`max_lag < N`. Omitting `max_lag` returns all lags `-(N-1)..N-1`.

## Window Functions (`window.hpp`)

```cpp
auto w = msl::signal::hamming_window(n);
auto w = msl::signal::hann_window(n);
auto w = msl::signal::blackman_window(n);
auto w = msl::signal::rectangular_window(n);

// Apply window
msl::signal::apply_window(signal, window);                // in-place
msl::signal::apply_window(input, output, window);         // to output
auto mat = msl::signal::apply_window_columns(mat, window); // to matrix columns
```

## Detrend (`detrend.hpp`)

Removes polynomial trends from signals.

```cpp
// Remove linear trend (default)
auto detrended = msl::signal::detrend(data);

// Remove mean (n=0)
auto detrended = msl::signal::detrend(data, 0);

// With explicit x-values
auto detrended = msl::signal::detrend(x, data, n);

// Matrix: detrend each column
auto detrended = msl::signal::detrend(data_matrix, n);
```

Internally uses `msl::polynomial::Polynomial` for least-squares fitting via QR decomposition.
