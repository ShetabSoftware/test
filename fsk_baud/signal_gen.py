"""Continuous-phase FSK signal generation utilities.

These helpers are primarily intended for testing, demonstrations and
benchmarking of the blind symbol-rate estimators, but they are also a perfectly
usable standalone CPFSK modulator.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Sequence

import numpy as np


@dataclass
class FSKSignal:
    """Container for a generated FSK signal and its ground-truth parameters."""

    iq: np.ndarray
    fs: float
    symbol_rate: float
    symbols: np.ndarray
    tone_frequencies: np.ndarray
    carrier_offset: float = 0.0
    snr_db: Optional[float] = None
    meta: dict = field(default_factory=dict)

    @property
    def samples_per_symbol(self) -> float:
        return self.fs / self.symbol_rate

    @property
    def duration(self) -> float:
        return self.iq.size / self.fs


def _tone_frequencies(
    order: int,
    symbol_rate: float,
    modulation_index: float,
    freq_separation: Optional[float],
) -> np.ndarray:
    """Compute the set of tone frequencies for an M-ary FSK constellation.

    The tones are centred about zero (DC). If ``freq_separation`` is given it
    overrides ``modulation_index`` (which otherwise defines the spacing as
    ``modulation_index * symbol_rate``).
    """
    if order < 2:
        raise ValueError("FSK order must be >= 2")
    if (order & (order - 1)) != 0:
        # Not strictly required, but FSK orders are conventionally powers of two.
        pass
    if freq_separation is None:
        freq_separation = modulation_index * symbol_rate
    if freq_separation <= 0:
        raise ValueError("tone separation must be positive")
    # Symmetric tones: e.g. order=2 -> [-0.5, +0.5] * sep, order=4 -> [-1.5..1.5]
    indices = np.arange(order) - (order - 1) / 2.0
    return indices * freq_separation


def generate_fsk(
    symbol_rate: float,
    fs: float,
    num_symbols: int,
    order: int = 2,
    modulation_index: float = 1.0,
    freq_separation: Optional[float] = None,
    carrier_offset: float = 0.0,
    snr_db: Optional[float] = None,
    symbols: Optional[Sequence[int]] = None,
    tone_frequencies: Optional[Sequence[float]] = None,
    amplitude: float = 1.0,
    seed: Optional[int] = None,
) -> FSKSignal:
    """Generate a continuous-phase M-FSK signal at complex baseband.

    Parameters
    ----------
    symbol_rate:
        Symbol (baud) rate in symbols/second.
    fs:
        Sample rate in hertz. Must be greater than the occupied bandwidth.
    num_symbols:
        Number of symbols to transmit.
    order:
        Modulation order ``M`` (number of tones). Defaults to binary FSK.
    modulation_index:
        CPFSK modulation index ``h``; tone spacing is ``h * symbol_rate`` unless
        ``freq_separation`` is supplied.
    freq_separation:
        Explicit spacing in hertz between adjacent tones. Overrides
        ``modulation_index`` when provided.
    carrier_offset:
        Residual carrier frequency offset in hertz applied to the whole signal.
    snr_db:
        If given, complex AWGN is added to achieve this signal-to-noise ratio
        (in dB), measured over the full sampled bandwidth ``fs``.
    symbols:
        Optional explicit symbol sequence (values in ``range(order)``). If not
        given a random sequence is drawn.
    tone_frequencies:
        Optional explicit tone frequencies (hertz). Overrides the computed set.
    amplitude:
        Linear amplitude of the (noise-free) signal.
    seed:
        Seed for the random number generator (symbols and noise).

    Returns
    -------
    FSKSignal
        The generated signal together with its ground-truth parameters.
    """
    if symbol_rate <= 0:
        raise ValueError("symbol_rate must be positive")
    if fs <= 0:
        raise ValueError("fs must be positive")
    if num_symbols < 1:
        raise ValueError("num_symbols must be >= 1")

    rng = np.random.default_rng(seed)

    if tone_frequencies is not None:
        tones = np.asarray(tone_frequencies, dtype=np.float64)
        order = tones.size
    else:
        tones = _tone_frequencies(
            order, symbol_rate, modulation_index, freq_separation
        )

    if symbols is not None:
        sym = np.asarray(symbols, dtype=int)
        if sym.size != num_symbols:
            num_symbols = sym.size
        if sym.min() < 0 or sym.max() >= order:
            raise ValueError("symbol values must lie in range(order)")
    else:
        sym = rng.integers(0, order, size=num_symbols)

    sps = fs / symbol_rate
    total_samples = int(round(num_symbols * sps))
    if total_samples < 2:
        raise ValueError("signal is too short; increase num_symbols or fs")

    # Map each output sample to its underlying symbol index. Using floor on the
    # fractional samples-per-symbol grid correctly supports non-integer sps.
    sample_index = np.arange(total_samples)
    symbol_of_sample = np.minimum(
        (sample_index / sps).astype(int), num_symbols - 1
    )
    inst_freq = tones[sym[symbol_of_sample]] + carrier_offset

    # Continuous-phase integration of the instantaneous frequency.
    phase = 2.0 * np.pi * np.cumsum(inst_freq) / fs
    iq = amplitude * np.exp(1j * phase)

    if snr_db is not None:
        iq = add_awgn(iq, snr_db, rng=rng)

    return FSKSignal(
        iq=iq.astype(np.complex128),
        fs=float(fs),
        symbol_rate=float(symbol_rate),
        symbols=sym,
        tone_frequencies=tones,
        carrier_offset=float(carrier_offset),
        snr_db=snr_db,
        meta={
            "order": int(order),
            "modulation_index": float(modulation_index),
            "samples_per_symbol": float(sps),
        },
    )


def add_awgn(
    iq: np.ndarray,
    snr_db: float,
    rng: Optional[np.random.Generator] = None,
) -> np.ndarray:
    """Add complex additive white Gaussian noise at the requested SNR.

    SNR is defined as the ratio of measured signal power to noise power over the
    full sampled bandwidth.
    """
    if rng is None:
        rng = np.random.default_rng()
    x = np.asarray(iq, dtype=np.complex128)
    signal_power = float(np.mean(np.abs(x) ** 2))
    if signal_power <= 0:
        return x.copy()
    snr_linear = 10.0 ** (snr_db / 10.0)
    noise_power = signal_power / snr_linear
    noise = np.sqrt(noise_power / 2.0) * (
        rng.standard_normal(x.shape) + 1j * rng.standard_normal(x.shape)
    )
    return x + noise
