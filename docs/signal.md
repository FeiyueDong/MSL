# Signal Module

**Namespace**: `msl::signal`  
**Aggregator**: `msl/signal.hpp` (incomplete — see note below)  
**Files**: 9 headers in `msl/signal/`

## Overview

The signal module provides digital signal processing capabilities: FFT/IFFT (1D/2D), Butterworth IIR filter design, zero-phase forward-backward filtering, Fourier-domain filtering, power spectral density estimation, window functions, and detrending.

> **Important**: `msl/signal.hpp` does NOT include `filtfilt.hpp` or `filter.hpp`. Include them explicitly:
> ```cpp
> #include "msl/signal/filter.hpp"
> #include "msl/signal/filtfilt.hpp"
> ```

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

### 1D IFFT

```cpp
auto x = msl::signal::ifft(X);               // complex → complex
auto x = msl::signal::ifft_real(X);          // complex → real (conjugate symmetric)
```

### 2D FFT (Column-wise / Row-wise)

```cpp
// Column-wise: each column is a signal
auto X = msl::signal::fft_columns(mat);      // matrix → matrixc

// Row-wise: each row is a signal
auto X = msl::signal::fft_rows(mat);         // matrix → matrixc

auto x = msl::signal::ifft_columns(X);       // inverse
auto x = msl::signal::ifft_rows_real(X);     // inverse to real
```

### Utility

```cpp
auto freqs = msl::signal::fft_frequencies(n, sample_rate); // frequency bins
auto power = msl::signal::power_spectrum(fft_result);      // |X[k]|^2
auto magn  = msl::signal::magnitude_spectrum(fft_result);  // |X[k]|
auto phase = msl::signal::phase_spectrum(fft_result);      // arg(X[k])
```

## Butterworth Filter (`butterworth_filter.hpp`)

Butterworth filters have maximally flat passband response. Design uses bilinear transformation from analog prototype.

### Class API

```cpp
// Lowpass/highpass
msl::signal::ButterworthFilter filt(order, fc, FilterType::lowpass);

// Bandpass/bandstop
msl::signal::ButterworthFilter filt(order, fc_low, fc_high, FilterType::bandpass);

// Get coefficients
auto coeffs = filt.coefficients();  // FilterCoefficients {b, a}
auto b = filt.b();                  // numerator
auto a = filt.a();                  // denominator

// Redesign with new parameters
filt.redesign(order, fc, FilterType::highpass);
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

## Filter Application (`filter.hpp`)

Direct Form II Transposed IIR filter — the most numerically stable form.

```cpp
auto filtered = msl::signal::filter(signal, coeffs);  // convenience

// Zero-copy
msl::signal::filter(signal, result, coeffs);

// Apply to matrix columns
auto filtered_mat = msl::signal::filter_columns(signal_matrix, coeffs);
```

## Zero-Phase Filtering (`filtfilt.hpp`)

Applies the filter forward and backward, eliminating phase distortion. Effective filter order is doubled. Equivalent to MATLAB's `filtfilt()`.

```cpp
auto filtered = msl::signal::filtfilt(signal, coeffs);  // convenience

// Zero-copy
msl::signal::filtfilt(signal, result, coeffs);

// Apply to matrix columns
auto filtered_mat = msl::signal::filtfilt_columns(signal_matrix, coeffs);
```

The implementation includes proper initial condition computation (`compute_filtfilt_zi`) and reflection padding (3x filter order on each side) to minimize startup transients.

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
msl::signal::FourierDomainFilter filt(
    fc_low, fc_high, FilterType::bandpass,
    FourierDomainFilter::WindowType::hamming, transition_band);

auto filtered = filt.apply(signal);

// Convenience functions
auto y = msl::signal::fourier_bandpass(signal, fc_low, fc_high);
auto y = msl::signal::fourier_lowpass(signal, fc_high);
auto y = msl::signal::fourier_highpass(signal, fc_low);
auto y = msl::signal::fourier_bandstop(signal, fc_low, fc_high);
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
