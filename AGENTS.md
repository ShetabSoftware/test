# AGENTS.md

## Cursor Cloud specific instructions

`fsk-baud` is a pure-Python DSP library + CLI for blind FSK symbol-rate estimation
(see `README.md`). There are **no servers, databases, or environment variables** —
"running" the product means running the CLI/examples and the test suite. A
standalone MATLAB/Octave port lives in `matlab/` (optional; MATLAB/Octave is not
installed here).

- **Python env**: the system Python 3.12 is PEP 668 "externally managed", so
  dependencies are installed into a project virtualenv at `.venv/` (created by the
  startup update script; `python3-venv` is already provisioned on the VM). Run
  everything via `.venv/bin/python ...` (or `source .venv/bin/activate` first).
  `.venv/` is gitignored.
- **Test**: `.venv/bin/python -m pytest -q` (config in `pyproject.toml`,
  `testpaths=["tests"]`).
- **Lint**: no linter/formatter is configured for this repo — there is nothing to run.
- **Run / smoke test**: `.venv/bin/python -m fsk_baud demo --fs 1e6 --symbol-rate 9600 --order 2 --snr 10`
  (or the installed `fsk-baud` console script). `estimate` reads a capture file:
  `.venv/bin/python -m fsk_baud estimate --input capture.npy --fs 1e6`.
- **Build**: `.venv/bin/python -m build` (produces sdist + wheel).
- **Plotting example**: `examples/plot_diagnostics.py` needs `matplotlib`; on this
  headless VM run it with `MPLBACKEND=Agg` so it saves `examples/diagnostics.png`
  instead of opening a window.
