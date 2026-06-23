"""Minimal end-to-end demonstration of blind FSK symbol-rate estimation.

Run with::

    python examples/demo.py
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from fsk_baud import estimate_symbol_rate, generate_fsk  # noqa: E402


def main() -> None:
    fs = 192_000.0
    configs = [
        dict(symbol_rate=2400.0, order=2, snr_db=None),
        dict(symbol_rate=4800.0, order=2, snr_db=10.0),
        dict(symbol_rate=9600.0, order=4, snr_db=12.0),
        dict(symbol_rate=1200.0, order=2, snr_db=8.0),
    ]

    print(f"{'true Rs':>9}  {'order':>5}  {'SNR':>5}  "
          f"{'estimate':>10}  {'error %':>8}  method")
    print("-" * 60)
    for cfg in configs:
        sig = generate_fsk(
            fs=fs, num_symbols=6000, carrier_offset=1500.0, seed=0, **cfg
        )
        est = estimate_symbol_rate(sig.iq, fs, method="auto")
        err = (est.symbol_rate - sig.symbol_rate) / sig.symbol_rate * 100.0
        snr = "clean" if cfg["snr_db"] is None else f"{cfg['snr_db']:.0f}dB"
        print(f"{sig.symbol_rate:9.0f}  {cfg['order']:5d}  {snr:>5}  "
              f"{est.symbol_rate:10.2f}  {err:+8.3f}  {est.method}")


if __name__ == "__main__":
    main()
