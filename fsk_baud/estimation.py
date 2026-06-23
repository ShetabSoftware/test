"""Blind symbol-rate (baud-rate) estimation for FSK signals.

Overview
--------
The instantaneous frequency (IF) of an FSK signal is, ideally, piecewise
constant: its level changes only at symbol boundaries. The raw IF waveform
behaves like a random NRZ sequence and therefore has a *continuous* power
spectral density with no discrete component at the symbol rate. To expose the
symbol rate we apply a memoryless nonlinearity to a difference of the IF. The
squared first difference produces a train of non-negative pulses located exactly
at the symbol transitions; the mean of that pulse train is periodic with the
symbol period, which creates discrete cyclostationary lines at integer multiples
of the symbol rate.

Pipeline
--------
1. **Band-limiting.** FSK captures are usually heavily oversampled, so most of
   the bandwidth is noise. The occupied band is detected from the averaged power
   spectrum and the signal is band-pass filtered to it, dramatically improving
   the effective SNR before the (noise-sensitive) frequency discriminator.
2. **Transition signal.** The squared first difference of the mean-removed
   instantaneous frequency.
3. **Fundamental detection.**
   - ``spectral``: a single long periodogram of the transition signal; the
     fundamental line is chosen from the harmonic comb (resolving the
     harmonic/sub-harmonic ambiguity) and refined by parabolic interpolation.
   - ``autocorr``: the autocorrelation of the transition signal; the symbol
     period is the smallest lag whose integer multiples are all peaks, refined
     against a far harmonic peak for fine resolution.

Both estimators are robust down to a few dB of (wideband) SNR for typical
oversampled FSK and accurate to well under 1% of the true symbol rate.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from .dsp import (
    autocorrelation,
    instantaneous_frequency,
    next_fast_len,
    quadratic_peak_interpolation,
)


@dataclass
class SymbolRateEstimate:
    """Result of a blind symbol-rate estimation.

    Attributes
    ----------
    symbol_rate:
        Estimated symbol rate in symbols/second (hertz).
    method:
        Name of the estimator that produced the result.
    confidence:
        A heuristic confidence score (peak-to-floor ratio of the detected
        feature). Larger values indicate a more reliable estimate.
    diagnostics:
        Optional dictionary with intermediate quantities (spectra, frequency
        axes, candidate peaks, the detected occupied band, ...) useful for
        plotting and debugging.
    """

    symbol_rate: float
    method: str
    confidence: float
    diagnostics: dict = field(default_factory=dict)

    def __float__(self) -> float:
        return float(self.symbol_rate)

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return (
            f"SymbolRateEstimate(symbol_rate={self.symbol_rate:.6g} Hz, "
            f"method={self.method!r}, confidence={self.confidence:.3g})"
        )


# ---------------------------------------------------------------------------
# Front-end: occupied-band detection and band-pass filtering
# ---------------------------------------------------------------------------
def _make_window(window: str, n: int) -> np.ndarray:
    """Construct an analysis window of length ``n``."""
    window = (window or "none").lower()
    if window in ("none", "rect", "boxcar"):
        return np.ones(n)
    if window == "hann":
        return np.hanning(n)
    if window == "hamming":
        return np.hamming(n)
    if window == "blackman":
        return np.blackman(n)
    raise ValueError(f"unknown window: {window!r}")


def detect_occupied_band(
    iq: np.ndarray,
    fs: float,
    threshold_mult: float = 4.0,
) -> tuple[float, float]:
    """Estimate the centre frequency and width of the occupied band.

    A Welch power spectrum is thresholded relative to a robust noise-floor
    estimate (the 25th percentile of the PSD). The lowest and highest
    frequencies whose smoothed power exceeds ``threshold_mult`` times the floor
    delimit the occupied band. Even at a few dB of wideband SNR the in-band PSD
    sits many times above the floor for oversampled FSK, so this is reliable.

    Returns
    -------
    (center_hz, bandwidth_hz)
    """
    x = np.asarray(iq, dtype=np.complex128)
    n = x.size
    nperseg = int(min(n, max(256, n // 32)))
    nperseg = max(nperseg, 16)
    noverlap = nperseg // 2
    step = max(nperseg - noverlap, 1)
    win = np.hanning(nperseg)
    win_norm = float(np.sum(win ** 2)) or 1.0

    accum = np.zeros(nperseg)
    count = 0
    for s in range(0, n - nperseg + 1, step):
        seg = x[s : s + nperseg]
        accum += np.abs(np.fft.fft(seg * win)) ** 2 / (fs * win_norm)
        count += 1
    if count == 0:
        accum = np.abs(np.fft.fft(x, nperseg)) ** 2 / (fs * win_norm)
    psd = np.fft.fftshift(accum / max(count, 1))
    freqs = np.fft.fftshift(np.fft.fftfreq(nperseg, d=1.0 / fs))

    floor = float(np.percentile(psd, 25)) or 1e-30
    mask = psd > threshold_mult * floor
    if not np.any(mask):
        return 0.0, fs
    idx = np.flatnonzero(mask)
    f_lo, f_hi = freqs[idx[0]], freqs[idx[-1]]
    return 0.5 * (f_lo + f_hi), float(f_hi - f_lo)


def bandlimit(
    iq: np.ndarray,
    fs: float,
    center: float,
    bandwidth: float,
    margin: float = 1.5,
) -> tuple[np.ndarray, float]:
    """Band-pass the signal to the occupied band via a brick-wall FFT filter.

    Returns the filtered signal and the (two-sided) pass bandwidth actually
    applied, which the autocorrelation estimator uses to skip the filter's
    correlation main lobe.
    """
    x = np.asarray(iq, dtype=np.complex128)
    n = x.size
    nfft = next_fast_len(n)
    X = np.fft.fftshift(np.fft.fft(x, nfft))
    freqs = np.fft.fftshift(np.fft.fftfreq(nfft, d=1.0 / fs))
    half = max(bandwidth / 2.0 * margin, fs / nfft * 8)
    half = min(half, fs / 2.0)
    out_of_band = (freqs < center - half) | (freqs > center + half)
    X[out_of_band] = 0.0
    y = np.fft.ifft(np.fft.ifftshift(X))[:n]
    return y, 2.0 * half


def _preprocess(
    iq: np.ndarray, fs: float, do_bandlimit: bool
) -> tuple[np.ndarray, float, float, float]:
    """Optionally band-limit the signal; return (y, band_center, band_bw, eff_bw)."""
    if not do_bandlimit:
        return np.asarray(iq, dtype=np.complex128), 0.0, fs, fs
    center, bw = detect_occupied_band(iq, fs)
    if bw <= 0 or bw >= fs * 0.98:
        # Already wideband / band detection inconclusive: skip filtering.
        return np.asarray(iq, dtype=np.complex128), center, bw, fs
    y, eff_bw = bandlimit(iq, fs, center, bw)
    return y, center, bw, eff_bw


# ---------------------------------------------------------------------------
# Transition signal and search-band helpers
# ---------------------------------------------------------------------------
def _transition_signal(
    iq: np.ndarray, fs: float, transform: str = "square"
) -> np.ndarray:
    """Build the symbol-transition pulse train from the instantaneous frequency.

    The instantaneous frequency is differenced once to emphasise transitions and
    to remove any constant carrier offset, then passed through a memoryless
    nonlinearity (``square`` or ``abs``) so that transitions of either sign add
    constructively.
    """
    inst_freq = instantaneous_frequency(iq, fs)
    inst_freq = inst_freq - np.mean(inst_freq)
    diff = np.diff(inst_freq)
    if transform == "square":
        trans = diff * diff
    elif transform == "abs":
        trans = np.abs(diff)
    else:
        raise ValueError("transform must be 'square' or 'abs'")
    # Remove the DC component so it does not dominate the spectral search.
    return trans - np.mean(trans)


def _resolve_band(
    fs: float,
    n: int,
    rs_min: Optional[float],
    rs_max: Optional[float],
    min_symbols: int = 16,
) -> tuple[float, float]:
    """Determine the symbol-rate search band, applying sensible defaults."""
    if rs_max is None:
        # The transition pulse train is real; its informative band is below
        # fs/2. FSK is normally generously oversampled, so cap the default well
        # below Nyquist to avoid spurious high-frequency peaks.
        rs_max = fs / 4.0
    if rs_min is None:
        # Require at least ``min_symbols`` symbols within the record so that a
        # spectral line can actually form.
        rs_min = max(min_symbols * fs / max(n, 1), fs / 1e6)
    if rs_min <= 0:
        raise ValueError("rs_min must be positive")
    if rs_max <= rs_min:
        raise ValueError("rs_max must be greater than rs_min")
    if rs_max >= fs / 2.0:
        rs_max = fs / 2.0 * 0.999
    return float(rs_min), float(rs_max)


def _local_maxima(y: np.ndarray) -> np.ndarray:
    """Indices of (non-strict) local maxima of a 1-D array."""
    if y.size < 3:
        return np.arange(y.size)
    greater_left = y[1:-1] >= y[:-2]
    greater_right = y[1:-1] >= y[2:]
    return np.flatnonzero(greater_left & greater_right) + 1


# ---------------------------------------------------------------------------
# Spectral estimator
# ---------------------------------------------------------------------------
def _select_fundamental(
    spectrum: np.ndarray,
    fs: float,
    nfft: int,
    rs_min: float,
    rs_max: float,
    peak_frac: float = 0.2,
    max_harmonics: int = 8,
    family_tol: float = 0.9,
) -> int:
    """Find the FFT bin of the fundamental symbol-rate line.

    The transition pulse train produces lines at every integer multiple of the
    symbol rate, frequently with comparable amplitude, so the global maximum is
    often a harmonic. Each candidate peak is scored by the mean energy of its
    harmonic comb (noise floor removed): the true fundamental and every harmonic
    that happens to be a clean line all attain a similar score, whereas spurious
    peaks and sub-harmonics score poorly because their combs fall on the noise
    floor. The fundamental is then the *lowest-frequency* member of the family
    of candidates whose scores are tied with the best.
    """
    df = fs / nfft
    k_min = max(int(np.ceil(rs_min / df)), 1)
    k_max = min(int(np.floor(rs_max / df)), spectrum.size - 1)
    if k_max <= k_min:
        raise ValueError("search band contains no FFT bins; widen the band")

    band = spectrum[k_min : k_max + 1]
    noise_floor = float(np.median(band))
    positive = np.clip(spectrum - noise_floor, 0.0, None)

    band_pos = positive[k_min : k_max + 1]
    band_max = float(np.max(band_pos))
    if band_max <= 0:
        return int(k_min + np.argmax(band))

    rel_maxima = _local_maxima(band_pos)
    candidates = [
        k_min + int(i)
        for i in rel_maxima
        if band_pos[i] >= peak_frac * band_max
    ]
    if not candidates:
        candidates = [int(k_min + np.argmax(band_pos))]

    def comb_score(bin0: int) -> float:
        f0_frac, _ = quadratic_peak_interpolation(positive, bin0)
        n_harm = min(max_harmonics, int(np.floor(k_max / f0_frac)))
        if n_harm < 2:
            return 0.0
        total = 0.0
        for k in range(1, n_harm + 1):
            center = int(round(k * f0_frac))
            w = max(2, int(round(0.004 * center)))
            lo = max(center - w, 0)
            hi = min(center + w + 1, positive.size)
            total += float(np.max(positive[lo:hi]))
        return total / n_harm

    scores = {b: comb_score(b) for b in candidates}
    best = max(scores.values())
    if best <= 0:
        return int(max(candidates, key=lambda b: positive[b]))
    family = [b for b, s in scores.items() if s >= family_tol * best]
    return int(min(family))


def estimate_symbol_rate_spectral(
    iq: np.ndarray,
    fs: float,
    rs_min: Optional[float] = None,
    rs_max: Optional[float] = None,
    nfft: Optional[int] = None,
    transform: str = "square",
    window: str = "hann",
    bandlimit: bool = True,
) -> SymbolRateEstimate:
    """Estimate the FSK symbol rate from the cyclostationary spectral line.

    Parameters
    ----------
    iq:
        Complex baseband samples of the FSK signal.
    fs:
        Sample rate in hertz.
    rs_min, rs_max:
        Lower/upper bounds of the symbol-rate search band (hertz). Sensible
        defaults are chosen from ``fs`` and the record length when omitted.
    nfft:
        FFT length. Defaults to the next power of two >= signal length.
    transform:
        Nonlinearity applied to the differenced instantaneous frequency,
        ``"square"`` (default) or ``"abs"``.
    window:
        Window applied before the FFT (``"hann"``, ``"hamming"``,
        ``"blackman"`` or ``"none"``).
    bandlimit:
        If true (default) the signal is band-pass filtered to its occupied band
        before processing, which greatly improves robustness at low SNR.

    Returns
    -------
    SymbolRateEstimate
    """
    if fs <= 0:
        raise ValueError("fs must be positive")
    y, band_center, band_bw, _ = _preprocess(iq, fs, bandlimit)

    trans = _transition_signal(y, fs, transform=transform)
    n = trans.size
    if n < 8:
        raise ValueError("signal is too short for spectral estimation")

    rs_min, rs_max = _resolve_band(fs, n, rs_min, rs_max)

    win = _make_window(window, n)
    if nfft is None:
        nfft = next_fast_len(n)
        nfft = min(nfft, 1 << 22)
    spectrum = np.abs(np.fft.rfft(trans * win, nfft))
    freqs = np.fft.rfftfreq(nfft, d=1.0 / fs)

    peak_idx = _select_fundamental(spectrum, fs, nfft, rs_min, rs_max)
    interp_idx, _ = quadratic_peak_interpolation(spectrum, peak_idx)
    symbol_rate = interp_idx * fs / nfft

    band = (freqs >= rs_min) & (freqs <= rs_max)
    med = float(np.median(spectrum[band])) or 1e-30
    confidence = float(spectrum[peak_idx] / med)

    return SymbolRateEstimate(
        symbol_rate=float(symbol_rate),
        method="spectral",
        confidence=confidence,
        diagnostics={
            "freqs": freqs,
            "spectrum": spectrum,
            "search_band": (rs_min, rs_max),
            "peak_index": peak_idx,
            "nfft": nfft,
            "occupied_band": (band_center, band_bw),
        },
    )


# ---------------------------------------------------------------------------
# Autocorrelation estimator
# ---------------------------------------------------------------------------
def estimate_symbol_rate_autocorr(
    iq: np.ndarray,
    fs: float,
    rs_min: Optional[float] = None,
    rs_max: Optional[float] = None,
    transform: str = "square",
    bandlimit: bool = True,
    max_harmonics: int = 8,
    family_tol: float = 0.9,
) -> SymbolRateEstimate:
    """Estimate the FSK symbol rate from the transition-signal autocorrelation.

    The squared-difference transition signal has autocorrelation peaks at every
    integer multiple of the symbol period. The symbol period is the smallest lag
    whose harmonic comb (the lag and its multiples) is as strong as the best
    candidate -- this resolves both the period/multiple ambiguity and any
    correlation main lobe introduced by band-limiting. The coarse period is then
    refined against a far harmonic peak, which divides the lag-quantisation error
    by the harmonic number.
    """
    if fs <= 0:
        raise ValueError("fs must be positive")
    y, band_center, band_bw, eff_bw = _preprocess(iq, fs, bandlimit)

    trans = _transition_signal(y, fs, transform=transform)
    n = trans.size
    if n < 8:
        raise ValueError("signal is too short for autocorrelation estimation")

    rs_min, rs_max = _resolve_band(fs, n, rs_min, rs_max)

    # Skip the correlation main lobe introduced by the band-pass filter: its
    # width is roughly the reciprocal of the pass bandwidth.
    lag_hump = int(np.ceil(1.3 * fs / eff_bw)) if eff_bw < fs else 1
    lag_min = max(int(np.floor(fs / rs_max)), lag_hump, 2)
    lag_max = min(int(np.ceil(fs / rs_min)), n // 3)
    if lag_max <= lag_min + 2:
        raise ValueError("search band yields an empty autocorrelation lag range")

    max_lag = min(n - 1, lag_max * max_harmonics)
    acf = autocorrelation(trans, max_lag=max_lag)
    zero_lag = acf[0] if acf[0] != 0 else 1e-30
    a = np.clip(acf / zero_lag, 0.0, None)

    scores = np.zeros(lag_max + 1)
    for tau in range(lag_min, lag_max + 1):
        m = min(max_harmonics, max_lag // tau)
        if m < 2:
            continue
        scores[tau] = float(np.mean([a[k * tau] for k in range(1, m + 1)]))
    best = float(np.max(scores))
    if best <= 0:
        raise ValueError("no periodic structure found in autocorrelation")
    family = np.flatnonzero(scores >= family_tol * best)
    tau0 = int(family[family >= lag_min][0])

    # Fine refinement using the highest reliable harmonic peak.
    kmax = min(max_harmonics, max_lag // tau0)
    refined_period = float(tau0)
    for k in range(kmax, 0, -1):
        center = k * tau0
        w = max(2, tau0 // 3)
        lo = max(center - w, 1)
        hi = min(center + w, len(acf) - 2)
        if hi <= lo:
            continue
        local = lo + int(np.argmax(acf[lo : hi + 1]))
        interp_lag, _ = quadratic_peak_interpolation(acf, local)
        refined_period = interp_lag / k
        break

    symbol_rate = fs / refined_period
    confidence = float(a[tau0])

    return SymbolRateEstimate(
        symbol_rate=float(symbol_rate),
        method="autocorr",
        confidence=confidence,
        diagnostics={
            "acf": acf,
            "scores": scores,
            "lag_range": (lag_min, lag_max),
            "coarse_period": tau0,
            "refined_period": refined_period,
            "occupied_band": (band_center, band_bw),
        },
    )


# ---------------------------------------------------------------------------
# High-level dispatcher
# ---------------------------------------------------------------------------
def estimate_symbol_rate(
    iq: np.ndarray,
    fs: float,
    method: str = "spectral",
    rs_min: Optional[float] = None,
    rs_max: Optional[float] = None,
    **kwargs,
) -> SymbolRateEstimate:
    """High-level blind symbol-rate estimator.

    Parameters
    ----------
    iq:
        Complex baseband samples.
    fs:
        Sample rate in hertz.
    method:
        ``"spectral"`` (default), ``"autocorr"``, or ``"auto"``. ``"auto"`` runs
        both estimators; when they agree the two are averaged and the confidence
        is boosted, otherwise the more confident estimate is returned.
    rs_min, rs_max:
        Symbol-rate search band bounds (hertz).
    **kwargs:
        Forwarded to the underlying estimator(s) (e.g. ``transform``,
        ``bandlimit``, ``window``).

    Returns
    -------
    SymbolRateEstimate
    """
    if method == "spectral":
        return estimate_symbol_rate_spectral(
            iq, fs, rs_min=rs_min, rs_max=rs_max, **kwargs
        )
    if method == "autocorr":
        ac_kwargs = {k: v for k, v in kwargs.items() if k != "window"}
        return estimate_symbol_rate_autocorr(
            iq, fs, rs_min=rs_min, rs_max=rs_max, **ac_kwargs
        )
    if method == "auto":
        spec = estimate_symbol_rate_spectral(
            iq, fs, rs_min=rs_min, rs_max=rs_max, **kwargs
        )
        ac_kwargs = {k: v for k, v in kwargs.items() if k != "window"}
        ac = estimate_symbol_rate_autocorr(
            iq, fs, rs_min=rs_min, rs_max=rs_max, **ac_kwargs
        )
        rel_diff = abs(spec.symbol_rate - ac.symbol_rate) / max(
            spec.symbol_rate, 1e-30
        )
        if rel_diff < 0.02:
            combined = 0.5 * (spec.symbol_rate + ac.symbol_rate)
            return SymbolRateEstimate(
                symbol_rate=combined,
                method="auto(agree)",
                confidence=spec.confidence + ac.confidence,
                diagnostics={"spectral": spec, "autocorr": ac},
            )
        chosen = max((spec, ac), key=lambda e: e.confidence)
        return SymbolRateEstimate(
            symbol_rate=chosen.symbol_rate,
            method=f"auto({chosen.method})",
            confidence=chosen.confidence,
            diagnostics={"spectral": spec, "autocorr": ac},
        )
    raise ValueError(f"unknown method: {method!r}")
