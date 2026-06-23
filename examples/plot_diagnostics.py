"""Visualise the internals of the spectral symbol-rate estimator.

Generates a noisy FSK signal, estimates its symbol rate and plots the
transition-signal spectrum with the detected fundamental and its harmonics.

Requires matplotlib (``pip install matplotlib``). Run with::

    python examples/plot_diagnostics.py
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from fsk_baud import estimate_symbol_rate_spectral, generate_fsk  # noqa: E402


def main() -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:  # pragma: no cover
        raise SystemExit("matplotlib is required: pip install matplotlib")

    fs = 192_000.0
    rs = 4800.0
    sig = generate_fsk(
        symbol_rate=rs, fs=fs, num_symbols=6000, order=2, snr_db=10.0, seed=1
    )
    est = estimate_symbol_rate_spectral(sig.iq, fs)

    freqs = est.diagnostics["freqs"]
    spec = est.diagnostics["spectrum"]
    rs_min, rs_max = est.diagnostics["search_band"]

    band = (freqs >= 0) & (freqs <= min(rs_max * 2.2, fs / 2))
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.plot(freqs[band], spec[band], lw=0.8, color="steelblue")
    for k in range(1, 6):
        ax.axvline(est.symbol_rate * k, color="crimson", ls="--", lw=0.8,
                   alpha=0.8, label="estimate & harmonics" if k == 1 else None)
    ax.axvline(rs, color="green", ls=":", lw=1.2, label="true Rs")
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel("|transition spectrum|")
    ax.set_title(
        f"Blind FSK symbol-rate estimate: {est.symbol_rate:.1f} Hz "
        f"(true {rs:.0f} Hz, SNR 10 dB)"
    )
    ax.legend()
    fig.tight_layout()
    out = "examples/diagnostics.png"
    fig.savefig(out, dpi=120)
    print(f"saved {out}; estimate={est.symbol_rate:.2f} Hz "
          f"(true {rs:.0f} Hz, confidence {est.confidence:.1f})")


if __name__ == "__main__":
    main()
