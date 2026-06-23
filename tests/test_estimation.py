import numpy as np
import pytest

from fsk_baud import (
    SymbolRateEstimate,
    estimate_symbol_rate,
    estimate_symbol_rate_autocorr,
    estimate_symbol_rate_spectral,
    generate_fsk,
)
from fsk_baud.estimation import detect_occupied_band


def _rel_err(est, truth):
    return abs(float(est) - truth) / truth


def test_generate_fsk_basic_shape():
    sig = generate_fsk(symbol_rate=1000.0, fs=20000.0, num_symbols=200, seed=1)
    assert sig.iq.dtype == np.complex128
    assert sig.samples_per_symbol == pytest.approx(20.0)
    # Constant-envelope (no noise) signal.
    assert np.allclose(np.abs(sig.iq), 1.0, atol=1e-9)


@pytest.mark.parametrize("symbol_rate", [800.0, 2400.0, 9600.0])
def test_spectral_noise_free_binary(symbol_rate):
    fs = 192000.0
    sig = generate_fsk(
        symbol_rate=symbol_rate, fs=fs, num_symbols=3000, order=2, seed=7
    )
    est = estimate_symbol_rate_spectral(sig.iq, fs)
    assert _rel_err(est, symbol_rate) < 0.02


@pytest.mark.parametrize("symbol_rate", [800.0, 2400.0, 9600.0])
def test_autocorr_noise_free_binary(symbol_rate):
    fs = 192000.0
    sig = generate_fsk(
        symbol_rate=symbol_rate, fs=fs, num_symbols=3000, order=2, seed=7
    )
    est = estimate_symbol_rate_autocorr(sig.iq, fs)
    assert _rel_err(est, symbol_rate) < 0.02


def test_spectral_4fsk():
    fs = 100000.0
    rs = 5000.0
    sig = generate_fsk(
        symbol_rate=rs, fs=fs, num_symbols=4000, order=4, seed=3
    )
    est = estimate_symbol_rate_spectral(sig.iq, fs)
    assert _rel_err(est, rs) < 0.02


def test_robust_to_carrier_offset():
    fs = 100000.0
    rs = 4000.0
    sig = generate_fsk(
        symbol_rate=rs,
        fs=fs,
        num_symbols=4000,
        order=2,
        carrier_offset=7777.0,
        seed=11,
    )
    est = estimate_symbol_rate(sig.iq, fs, method="auto")
    assert _rel_err(est, rs) < 0.02


@pytest.mark.parametrize("snr_db", [20.0, 10.0, 7.0])
def test_spectral_with_noise(snr_db):
    fs = 192000.0
    rs = 4800.0
    sig = generate_fsk(
        symbol_rate=rs,
        fs=fs,
        num_symbols=6000,
        order=2,
        snr_db=snr_db,
        seed=42,
    )
    est = estimate_symbol_rate_spectral(sig.iq, fs)
    assert _rel_err(est, rs) < 0.03


@pytest.mark.parametrize("snr_db", [20.0, 10.0, 7.0])
def test_autocorr_with_noise(snr_db):
    fs = 192000.0
    rs = 4800.0
    sig = generate_fsk(
        symbol_rate=rs,
        fs=fs,
        num_symbols=6000,
        order=2,
        snr_db=snr_db,
        seed=42,
    )
    est = estimate_symbol_rate_autocorr(sig.iq, fs)
    assert _rel_err(est, rs) < 0.03


@pytest.mark.parametrize("rs", [800.0, 2400.0, 4800.0, 9600.0])
@pytest.mark.parametrize("method", ["spectral", "autocorr", "auto"])
def test_robustness_sweep_at_low_snr(rs, method):
    fs = 192000.0
    sig = generate_fsk(
        symbol_rate=rs,
        fs=fs,
        num_symbols=6000,
        order=2,
        snr_db=8.0,
        seed=123,
    )
    est = estimate_symbol_rate(sig.iq, fs, method=method)
    assert _rel_err(est, rs) < 0.03


def test_bandlimit_can_be_disabled():
    fs = 100000.0
    rs = 4000.0
    sig = generate_fsk(symbol_rate=rs, fs=fs, num_symbols=4000, seed=4)
    est = estimate_symbol_rate_spectral(sig.iq, fs, bandlimit=False)
    assert _rel_err(est, rs) < 0.02


def test_abs_transform():
    fs = 100000.0
    rs = 5000.0
    sig = generate_fsk(symbol_rate=rs, fs=fs, num_symbols=4000, seed=6)
    est = estimate_symbol_rate_spectral(sig.iq, fs, transform="abs")
    assert _rel_err(est, rs) < 0.02


def test_detect_occupied_band():
    fs = 200000.0
    rs = 4000.0
    sig = generate_fsk(
        symbol_rate=rs, fs=fs, num_symbols=4000, order=2, snr_db=15.0, seed=8
    )
    center, bw = detect_occupied_band(sig.iq, fs)
    # Binary FSK with h=1 occupies roughly 2*Rs; allow a generous band.
    assert abs(center) < rs
    assert rs < bw < 12 * rs


def test_auto_method_agreement_boosts_confidence():
    fs = 100000.0
    rs = 2500.0
    sig = generate_fsk(symbol_rate=rs, fs=fs, num_symbols=5000, seed=5)
    spec = estimate_symbol_rate(sig.iq, fs, method="spectral")
    auto = estimate_symbol_rate(sig.iq, fs, method="auto")
    assert _rel_err(auto, rs) < 0.02
    assert auto.method.startswith("auto")
    # When estimators agree the combined confidence exceeds the spectral one.
    assert auto.confidence >= spec.confidence


def test_result_float_and_repr():
    fs = 50000.0
    rs = 2000.0
    sig = generate_fsk(symbol_rate=rs, fs=fs, num_symbols=3000, seed=9)
    est = estimate_symbol_rate_spectral(sig.iq, fs)
    assert isinstance(est, SymbolRateEstimate)
    assert isinstance(float(est), float)
    assert "SymbolRateEstimate" in repr(est)


def test_search_band_restriction():
    fs = 200000.0
    rs = 10000.0
    sig = generate_fsk(symbol_rate=rs, fs=fs, num_symbols=4000, seed=2)
    est = estimate_symbol_rate_spectral(
        sig.iq, fs, rs_min=5000.0, rs_max=15000.0
    )
    assert _rel_err(est, rs) < 0.02


def test_invalid_method_raises():
    sig = generate_fsk(symbol_rate=1000.0, fs=20000.0, num_symbols=500, seed=1)
    with pytest.raises(ValueError):
        estimate_symbol_rate(sig.iq, 20000.0, method="nonsense")


def test_invalid_transform_raises():
    sig = generate_fsk(symbol_rate=1000.0, fs=20000.0, num_symbols=500, seed=1)
    with pytest.raises(ValueError):
        estimate_symbol_rate_spectral(sig.iq, 20000.0, transform="bogus")


def test_short_signal_raises():
    with pytest.raises(ValueError):
        estimate_symbol_rate_spectral(np.ones(4, dtype=complex), 1000.0)
