"""Low-level DSP helpers used by the symbol-rate estimators."""

from __future__ import annotations

from typing import Tuple

import numpy as np


def _as_complex(iq: np.ndarray) -> np.ndarray:
    """Return ``iq`` as a 1-D complex128 array.

    Real-valued input is accepted and treated as the in-phase component with a
    zero quadrature component (i.e. an analytic signal is *not* formed here).
    """
    arr = np.asarray(iq)
    if arr.ndim != 1:
        arr = arr.reshape(-1)
    if not np.iscomplexobj(arr):
        arr = arr.astype(np.float64)
    return arr.astype(np.complex128, copy=False)


def instantaneous_phase(iq: np.ndarray) -> np.ndarray:
    """Unwrapped instantaneous phase (radians) of a complex signal."""
    x = _as_complex(iq)
    return np.unwrap(np.angle(x))


def instantaneous_frequency(iq: np.ndarray, fs: float) -> np.ndarray:
    """Estimate the instantaneous frequency of a complex baseband signal.

    The instantaneous frequency is the time derivative of the unwrapped phase,
    approximated here with a first difference. For an input of length ``N`` the
    returned array has length ``N - 1`` and is expressed in hertz.

    Parameters
    ----------
    iq:
        Complex baseband samples.
    fs:
        Sample rate in hertz.

    Returns
    -------
    numpy.ndarray
        Instantaneous frequency in hertz.
    """
    if fs <= 0:
        raise ValueError("fs must be positive")
    x = _as_complex(iq)
    if x.size < 2:
        raise ValueError("at least two samples are required")
    # Differentiate the phase directly from the product x[n] * conj(x[n-1]).
    # This is numerically equivalent to diff(unwrap(angle(x))) but avoids the
    # need for an explicit unwrap and is robust to large hops.
    product = x[1:] * np.conj(x[:-1])
    dphi = np.angle(product)
    return dphi * (fs / (2.0 * np.pi))


def quadratic_peak_interpolation(
    y: np.ndarray, k: int
) -> Tuple[float, float]:
    """Refine the location of a discrete peak with parabolic interpolation.

    Fits a parabola through the three samples centred on index ``k`` and returns
    the interpolated (fractional) peak location and the interpolated peak value.

    Parameters
    ----------
    y:
        Sample values (e.g. a magnitude spectrum).
    k:
        Integer index of the local maximum.

    Returns
    -------
    (float, float)
        ``(interpolated_index, interpolated_value)``.
    """
    y = np.asarray(y, dtype=np.float64)
    n = y.size
    if k <= 0 or k >= n - 1:
        # Cannot interpolate at the edges; return the sample as-is.
        return float(k), float(y[k])
    ym1, y0, yp1 = y[k - 1], y[k], y[k + 1]
    denom = ym1 - 2.0 * y0 + yp1
    if denom == 0.0:
        return float(k), float(y0)
    delta = 0.5 * (ym1 - yp1) / denom
    # Clamp to the valid interpolation range for numerical safety.
    delta = float(np.clip(delta, -1.0, 1.0))
    value = y0 - 0.25 * (ym1 - yp1) * delta
    return float(k) + delta, float(value)


def next_fast_len(n: int) -> int:
    """Smallest power of two that is >= ``n`` (fast and dependency-free)."""
    if n <= 1:
        return 1
    return 1 << (int(n - 1).bit_length())


def autocorrelation(x: np.ndarray, max_lag: int | None = None) -> np.ndarray:
    """Biased autocorrelation of a real sequence computed via the FFT.

    Returns the non-negative-lag portion of the autocorrelation, starting at
    lag 0.
    """
    x = np.asarray(x, dtype=np.float64)
    x = x - np.mean(x)
    n = x.size
    if n == 0:
        return np.zeros(0)
    nfft = next_fast_len(2 * n)
    spec = np.fft.rfft(x, nfft)
    acf_full = np.fft.irfft(spec * np.conj(spec), nfft)[:n]
    if max_lag is not None:
        acf_full = acf_full[: max_lag + 1]
    return acf_full
