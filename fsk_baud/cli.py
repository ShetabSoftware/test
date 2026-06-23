"""Command-line interface for blind FSK symbol-rate estimation.

Examples
--------
Estimate the symbol rate of IQ samples stored in a NumPy ``.npy`` file::

    python -m fsk_baud estimate --input capture.npy --fs 1e6

Estimate from interleaved float32 I/Q (``I0 Q0 I1 Q1 ...``)::

    python -m fsk_baud estimate --input capture.cf32 --format cf32 --fs 1e6

Run a self-contained demo that synthesises a signal and estimates its rate::

    python -m fsk_baud demo --fs 1e6 --symbol-rate 9600 --snr 10
"""

from __future__ import annotations

import argparse
import sys
from typing import Optional

import numpy as np

from .estimation import estimate_symbol_rate
from .signal_gen import generate_fsk


def _load_iq(path: str, fmt: str) -> np.ndarray:
    """Load complex baseband samples from a file in the requested format."""
    fmt = fmt.lower()
    if fmt == "npy":
        data = np.load(path)
        return np.asarray(data).reshape(-1)
    dtype_map = {
        "cf32": np.float32,
        "cf64": np.float64,
        "ci16": np.int16,
        "ci8": np.int8,
    }
    if fmt not in dtype_map:
        raise ValueError(f"unsupported input format: {fmt!r}")
    raw = np.fromfile(path, dtype=dtype_map[fmt]).astype(np.float64)
    if raw.size % 2 != 0:
        raw = raw[:-1]
    return raw[0::2] + 1j * raw[1::2]


def _add_band_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--fs", type=float, required=True, help="Sample rate (Hz).")
    p.add_argument(
        "--rs-min",
        type=float,
        default=None,
        help="Lower bound of the symbol-rate search band (Hz).",
    )
    p.add_argument(
        "--rs-max",
        type=float,
        default=None,
        help="Upper bound of the symbol-rate search band (Hz).",
    )
    p.add_argument(
        "--method",
        choices=["spectral", "autocorr", "auto"],
        default="auto",
        help="Estimation method (default: auto).",
    )
    p.add_argument(
        "--transform",
        choices=["square", "abs"],
        default="square",
        help="Nonlinearity applied to the differenced instantaneous frequency.",
    )


def _print_result(result, fs: float, true_rate: Optional[float] = None) -> None:
    print(f"Estimated symbol rate : {result.symbol_rate:.6g} Hz")
    print(f"Method                : {result.method}")
    print(f"Confidence            : {result.confidence:.4g}")
    sps = fs / result.symbol_rate if result.symbol_rate else float("nan")
    print(f"Samples per symbol    : {sps:.4g}")
    if true_rate is not None:
        err = (result.symbol_rate - true_rate) / true_rate * 100.0
        print(f"True symbol rate      : {true_rate:.6g} Hz")
        print(f"Relative error        : {err:+.3f} %")


def _cmd_estimate(args: argparse.Namespace) -> int:
    iq = _load_iq(args.input, args.format)
    if iq.size < 16:
        print("error: not enough samples in input", file=sys.stderr)
        return 2
    result = estimate_symbol_rate(
        iq,
        fs=args.fs,
        method=args.method,
        rs_min=args.rs_min,
        rs_max=args.rs_max,
        transform=args.transform,
    )
    _print_result(result, args.fs)
    return 0


def _cmd_demo(args: argparse.Namespace) -> int:
    sig = generate_fsk(
        symbol_rate=args.symbol_rate,
        fs=args.fs,
        num_symbols=args.num_symbols,
        order=args.order,
        modulation_index=args.modulation_index,
        carrier_offset=args.carrier_offset,
        snr_db=args.snr,
        seed=args.seed,
    )
    result = estimate_symbol_rate(
        sig.iq,
        fs=args.fs,
        method=args.method,
        rs_min=args.rs_min,
        rs_max=args.rs_max,
        transform=args.transform,
    )
    _print_result(result, args.fs, true_rate=sig.symbol_rate)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="fsk_baud",
        description="Blind symbol-rate (baud-rate) estimation for FSK signals.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    est = sub.add_parser("estimate", help="Estimate the rate of recorded IQ.")
    est.add_argument("--input", "-i", required=True, help="Input IQ file.")
    est.add_argument(
        "--format",
        "-f",
        default="npy",
        choices=["npy", "cf32", "cf64", "ci16", "ci8"],
        help="Input file format (default: npy).",
    )
    _add_band_args(est)
    est.set_defaults(func=_cmd_estimate)

    demo = sub.add_parser("demo", help="Synthesise a signal and estimate it.")
    demo.add_argument("--symbol-rate", type=float, default=9600.0)
    demo.add_argument("--num-symbols", type=int, default=4000)
    demo.add_argument("--order", type=int, default=2)
    demo.add_argument("--modulation-index", type=float, default=1.0)
    demo.add_argument("--carrier-offset", type=float, default=0.0)
    demo.add_argument("--snr", type=float, default=None, help="SNR in dB.")
    demo.add_argument("--seed", type=int, default=None)
    _add_band_args(demo)
    demo.set_defaults(func=_cmd_demo)

    return parser


def main(argv: Optional[list] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
