import numpy as np
import pytest

from fsk_baud.dsp import (
    autocorrelation,
    instantaneous_frequency,
    instantaneous_phase,
    next_fast_len,
    quadratic_peak_interpolation,
)


def test_instantaneous_frequency_pure_tone():
    fs = 1000.0
    f0 = 123.0
    n = np.arange(4096)
    x = np.exp(1j * 2 * np.pi * f0 * n / fs)
    inst = instantaneous_frequency(x, fs)
    assert np.allclose(np.mean(inst), f0, atol=1e-6)
    assert inst.size == x.size - 1


def test_instantaneous_frequency_requires_two_samples():
    with pytest.raises(ValueError):
        instantaneous_frequency(np.array([1 + 0j]), 1.0)


def test_instantaneous_frequency_rejects_bad_fs():
    with pytest.raises(ValueError):
        instantaneous_frequency(np.ones(4, dtype=complex), 0.0)


def test_instantaneous_phase_monotonic_for_positive_tone():
    fs = 1000.0
    n = np.arange(100)
    x = np.exp(1j * 2 * np.pi * 50 * n / fs)
    phase = instantaneous_phase(x)
    assert np.all(np.diff(phase) > 0)


def test_quadratic_peak_interpolation_centered():
    # Symmetric parabola peaks exactly at the centre sample.
    y = np.array([0.0, 1.0, 4.0, 1.0, 0.0])
    idx, val = quadratic_peak_interpolation(y, 2)
    assert idx == pytest.approx(2.0, abs=1e-9)
    assert val == pytest.approx(4.0, abs=1e-9)


def test_quadratic_peak_interpolation_offset():
    # Parabola y = -(x - 2.3)^2 + 10 sampled on the integer grid.
    x = np.arange(6)
    y = -((x - 2.3) ** 2) + 10
    k = int(np.argmax(y))
    idx, val = quadratic_peak_interpolation(y, k)
    assert idx == pytest.approx(2.3, abs=1e-6)
    assert val == pytest.approx(10.0, abs=1e-6)


def test_quadratic_peak_interpolation_edge():
    y = np.array([5.0, 1.0, 0.0])
    idx, val = quadratic_peak_interpolation(y, 0)
    assert idx == 0.0
    assert val == 5.0


def test_next_fast_len_powers_of_two():
    assert next_fast_len(1) == 1
    assert next_fast_len(5) == 8
    assert next_fast_len(1024) == 1024
    assert next_fast_len(1025) == 2048


def test_autocorrelation_periodic_signal():
    fs = 1.0
    period = 20
    n = np.arange(2000)
    x = np.sin(2 * np.pi * n / period)
    acf = autocorrelation(x, max_lag=60)
    # The autocorrelation of a sinusoid peaks at multiples of the period.
    assert acf[0] == pytest.approx(np.max(acf))
    peak_near_period = period + int(np.argmax(acf[period - 2 : period + 3])) - 2
    assert peak_near_period == period
