# fsk-baud — Blind symbol-rate estimation for FSK

`fsk-baud` estimates the **symbol (baud) rate** of a frequency-shift-keyed (FSK)
signal directly from its complex-baseband (IQ) samples, with **no prior
knowledge** of the modulation parameters — number of tones, tone spacing,
modulation index, carrier offset or pulse shaping are all unknown.

It works for 2-FSK and M-FSK, is robust to a residual carrier offset and to a
few dB of wideband noise (for typical oversampled captures), and is accurate to
well under 1 % of the true symbol rate.

```text
  true Rs  order    SNR    estimate   error %  method
------------------------------------------------------------
     2400      2  clean     2400.01    +0.000  auto(agree)
     4800      2   10dB     4799.88    -0.002  auto(agree)
     9600      4   12dB     9599.22    -0.008  auto(agree)
     1200      2    8dB     1200.01    +0.001  auto(agree)
```

## How it works

The instantaneous frequency (IF) of an FSK signal is ideally *piecewise
constant*: it only changes at symbol boundaries. The raw IF behaves like a
random NRZ sequence, which has a **continuous** spectrum with no line at the
symbol rate. The estimator therefore:

1. **Band-limits** the signal to its occupied band. FSK is usually heavily
   oversampled, so most of the captured bandwidth is noise; detecting the
   occupied band from the averaged power spectrum and band-pass filtering to it
   dramatically improves the SNR seen by the (noise-sensitive) frequency
   discriminator.
2. Forms a **transition signal**: the squared first difference of the
   mean-removed instantaneous frequency. This produces a non-negative pulse at
   every symbol transition, whose periodic mean creates discrete
   **cyclostationary lines** at integer multiples of the symbol rate.
3. **Locates the fundamental** of those lines:
   - the **spectral** estimator uses a long periodogram and a harmonic-comb
     search that resolves the harmonic / sub-harmonic ambiguity, then refines
     the peak with parabolic interpolation;
   - the **autocorrelation** estimator finds the smallest lag whose integer
     multiples are all autocorrelation peaks, then refines it against a far
     harmonic peak (dividing the lag-quantisation error by the harmonic index).

The `auto` method runs both and cross-checks them.

## Installation

```bash
pip install -e .            # from a clone of this repository
# or just install the single runtime dependency:
pip install numpy
```

Optional extras: `matplotlib` (for `examples/plot_diagnostics.py`) and `pytest`
(for the test suite).

## Library usage

```python
import numpy as np
from fsk_baud import estimate_symbol_rate, generate_fsk

# Synthesise a noisy 2-FSK signal (or load your own IQ samples).
sig = generate_fsk(symbol_rate=9600.0, fs=192_000.0, num_symbols=6000,
                   order=2, snr_db=10.0, carrier_offset=2_000.0, seed=0)

est = estimate_symbol_rate(sig.iq, fs=192_000.0, method="auto")
print(est.symbol_rate)   # ~9600 Hz
print(est.confidence)    # peak-to-floor ratio of the detected feature
```

`estimate_symbol_rate` returns a `SymbolRateEstimate` with the fields
`symbol_rate`, `method`, `confidence` and a `diagnostics` dictionary (spectra,
autocorrelation, detected occupied band, ...). It also casts to `float`:

```python
rate_hz = float(est)
```

### Choosing a method

| `method`     | description                                                |
|--------------|------------------------------------------------------------|
| `"spectral"` | frequency-domain cyclostationary line (default, accurate)  |
| `"autocorr"` | time-domain periodicity of the transition signal           |
| `"auto"`     | runs both and combines / cross-checks them                 |

Useful keyword arguments (forwarded to the estimators):

- `rs_min`, `rs_max` — restrict the symbol-rate search band (Hz).
- `bandlimit=False` — disable the occupied-band pre-filter.
- `transform="square"|"abs"` — nonlinearity used to build the transition signal.

## Command-line interface

```bash
# Estimate from a recorded capture (NumPy .npy of complex samples):
python -m fsk_baud estimate --input capture.npy --fs 1e6

# Interleaved float32 I/Q (I0 Q0 I1 Q1 ...):
python -m fsk_baud estimate --input capture.cf32 --format cf32 --fs 1e6

# Self-contained demo (synthesise + estimate):
python -m fsk_baud demo --fs 1e6 --symbol-rate 9600 --order 2 --snr 10
```

Supported input formats: `npy`, `cf32`, `cf64`, `ci16`, `ci8`.

## Examples

```bash
python examples/demo.py              # accuracy table across rates / SNRs
python examples/plot_diagnostics.py  # plot the transition spectrum (needs matplotlib)
```

## Testing

```bash
pip install pytest
pytest -q
```

## Module layout

| Module                  | Contents                                              |
|-------------------------|-------------------------------------------------------|
| `fsk_baud.estimation`   | estimators, occupied-band detection, band-limiting    |
| `fsk_baud.signal_gen`   | continuous-phase M-FSK generator and AWGN helper      |
| `fsk_baud.dsp`          | instantaneous frequency, autocorrelation, peak interp |
| `fsk_baud.cli`          | command-line interface                                |

## Limitations / notes

- Designed for complex-baseband (IQ) input. The signal should contain at least a
  few tens of symbols.
- Robustness depends on the in-band SNR, which improves with oversampling and
  longer captures; very low SNR or near-critically-sampled signals are harder.
- The default search band excludes symbol rates above `fs/4`; pass `rs_max` to
  override.
