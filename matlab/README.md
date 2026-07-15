# Blind FSK symbol-rate estimation — MATLAB

A MATLAB/Octave port of the `fsk_baud` algorithm: it estimates the **symbol
(baud) rate** of an FSK signal blindly from complex-baseband (IQ) samples, with
no prior knowledge of the tone count, tone spacing, modulation index or carrier
offset.

Only base MATLAB is required — no Signal Processing or Statistics toolboxes (the
windows, percentile and autocorrelation are implemented directly). The code also
runs unmodified under GNU Octave.

## Quick start

```matlab
addpath('matlab');                       % make the functions visible

% Synthesise a noisy 2-FSK signal (or load your own IQ vector instead).
sig = generateFSK(9600, 192000, 6000, 'order', 2, 'snrDb', 10, ...
                  'carrierOffset', 2000, 'seed', 0);

est = estimateSymbolRate(sig.iq, 192000, 'method', 'auto');
fprintf('Rs = %.1f Hz  (confidence %.1f, %s)\n', ...
        est.symbolRate, est.confidence, est.method);
```

Run the bundled demo and self-test:

```matlab
demoFskBaud      % accuracy table across rates / orders / SNRs
runTests         % assertion-based test suite
```

## How it works

1. **Band-limit** to the occupied band (detected from the averaged power
   spectrum) to raise the effective SNR before the noise-sensitive frequency
   discriminator.
2. Build a **transition signal** — the squared first difference of the
   mean-removed instantaneous frequency — whose periodic mean produces discrete
   cyclostationary lines at integer multiples of the symbol rate.
3. **Find the fundamental** either in the frequency domain (`spectral`,
   harmonic-comb search + parabolic refinement) or in the lag domain
   (`autocorr`, smallest periodic lag + far-harmonic refinement). `auto` runs
   both and cross-checks.

## HF 2-FSK blind symbol-rate estimator (`blindSymbolRateFSK.m`)

`blindSymbolRateFSK.m` is a self-contained estimator tailored to the HF 2-FSK /
Watterson scenario (e.g. `Fs=9600`, tones 600/900 Hz). It is a clean, toolbox-
free packaging of the IF-hard-limit → transition-clock-FFT approach, and is a
drop-in replacement for a hand-rolled estimation block:

```matlab
Rs_hat = blindSymbolRateFSK(rxChanNoisy, Fs, Rs_min, Rs_max);
% Rs == Delta f case (e.g. Rs=300 with tones 600/900): disable the notch
Rs_hat = blindSymbolRateFSK(rx, Fs, Rs_min, Rs_max, 'notchToneSpacing', false);
```

`blind2FSK_rx_test.m` is the full receiver simulation (2FSK + Watterson channel +
BER + Monte-Carlo + plots) wired to call this estimator, and `testBlindRsHF.m`
is a fast Monte-Carlo accuracy harness.

**Practical notes / limitations (HF channel).** Blind symbol-rate estimation on
a fast-fading multipath HF channel is genuinely hard:

- On multipath the instantaneous frequency **beats at the tone spacing `Δf`**
  (when direct/echo paths carry different tones), producing a strong spurious
  line at `Δf`; it is notched using a blind `Δf` estimate.
- `Rs = Δf` (e.g. 300 baud with 300 Hz spacing) is an **inherent ambiguity**
  (a symbol clock vs the tone beat) — disable the notch for that case.
- The modulator in the reference script is **discontinuous-phase** FSK, which
  injects extra tone-related lines into the transition signal; a
  continuous-phase (CPFSK) source is markedly easier to estimate blindly (the
  `generateFSK.m` / `estimateSymbolRate.m` tools above handle CPFSK very well).
- Accuracy improves with record length, higher `Eb/N0`, milder channels, and
  Monte-Carlo averaging; tune via `'notchHz'`, `'smoothDiv'`, `'subharmThresh'`.

## Functions

| File                              | Purpose                                        |
|-----------------------------------|------------------------------------------------|
| `estimateSymbolRate.m`            | High-level entry point (`method` dispatch)     |
| `estimateSymbolRateSpectral.m`    | Frequency-domain estimator                     |
| `estimateSymbolRateAutocorr.m`    | Time-domain (autocorrelation) estimator        |
| `generateFSK.m`                   | Continuous-phase M-FSK generator               |
| `addAWGN.m`                       | Complex AWGN at a target SNR                    |
| `detectOccupiedBand.m`            | Occupied-band detection                        |
| `bandlimitSignal.m`               | Brick-wall FFT band-pass filter                |
| `instantaneousFrequency.m`        | Phase-difference instantaneous frequency       |
| `transitionSignal.m`              | Squared-difference transition signal           |
| `selectFundamental.m`             | Harmonic-comb fundamental selection            |
| `autocorrelation.m`, `localMaxima.m`, `quadraticPeakInterp.m`, `analysisWindow.m`, `percentileValue.m`, `resolveSearchBand.m`, `parseOptions.m`, `preprocessSignal.m` | helpers |
| `demoFskBaud.m`, `runTests.m`     | demo and tests                                 |

## Options

`estimateSymbolRate(iq, fs, 'Name', Value, ...)`:

- `'method'`    — `'spectral'` (default), `'autocorr'`, `'auto'`.
- `'rsMin'`, `'rsMax'` — symbol-rate search band (Hz); empty for auto defaults.
- `'transform'` — `'square'` (default) or `'abs'`.
- `'window'`    — `'hann'` (default), `'hamming'`, `'blackman'`, `'none'`.
- `'bandlimit'` — `true` (default) to band-pass to the occupied band first.

The returned struct has fields `symbolRate`, `method`, `confidence` and
`diagnostics` (spectra, autocorrelation, detected occupied band, ...).
