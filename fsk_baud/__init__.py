"""Blind symbol (baud) rate estimation for FSK signals.

This package provides tools to estimate the symbol rate of a frequency-shift
keyed (FSK) signal directly from its complex baseband (IQ) samples, without any
prior knowledge of the modulation parameters (number of tones, tone spacing,
carrier offset, pulse shaping, etc.).

The core idea is that the instantaneous frequency of an FSK signal is a
piecewise-constant waveform whose level changes only at symbol boundaries.
A random NRZ-like waveform has no discrete spectral line at the symbol rate, so
a nonlinear transform (the squared first difference of the instantaneous
frequency) is applied to generate a strong cyclostationary spectral line at the
symbol rate. The location of that line is the estimated symbol rate.

Public API
----------
- :func:`fsk_baud.estimation.estimate_symbol_rate`
- :func:`fsk_baud.estimation.estimate_symbol_rate_spectral`
- :func:`fsk_baud.estimation.estimate_symbol_rate_autocorr`
- :class:`fsk_baud.estimation.SymbolRateEstimate`
- :func:`fsk_baud.signal_gen.generate_fsk`
- :func:`fsk_baud.dsp.instantaneous_frequency`
"""

from .dsp import (
    instantaneous_frequency,
    instantaneous_phase,
    quadratic_peak_interpolation,
)
from .estimation import (
    SymbolRateEstimate,
    bandlimit,
    detect_occupied_band,
    estimate_symbol_rate,
    estimate_symbol_rate_autocorr,
    estimate_symbol_rate_spectral,
)
from .signal_gen import add_awgn, generate_fsk

__version__ = "0.1.0"

__all__ = [
    "instantaneous_frequency",
    "instantaneous_phase",
    "quadratic_peak_interpolation",
    "SymbolRateEstimate",
    "estimate_symbol_rate",
    "estimate_symbol_rate_spectral",
    "estimate_symbol_rate_autocorr",
    "detect_occupied_band",
    "bandlimit",
    "generate_fsk",
    "add_awgn",
]
