# GNSS Anti-Spoofing Array Processor

An FPGA-oriented redesign of the multi-antenna GPS anti-spoofing method of
Daneshmand, Jafarnia-Jahromi, Broumandan and Lachapelle (ION GNSS 2012),
*"A Low-Complexity GPS Anti-Spoofing Method Using a Multi-Antenna Array"*.

This repository contains a **verified algorithm reference model** intended to
serve as the golden model for an RTL implementation, plus the engineering
analysis behind the design decisions.

---

## Documents

| | |
|---|---|
| [`docs/00-executive-summary.md`](docs/00-executive-summary.md) | The findings and recommendations in one place |
| [`docs/01-array-architecture.md`](docs/01-array-architecture.md) | **Task 1** — 3 vs 4 antennas, geometry, DOF, measured evidence |
| [`docs/02-matlab-review-and-redesign.md`](docs/02-matlab-review-and-redesign.md) | **Task 2** — review of the original code and the redesign |
| [`docs/03-fpga-implementation.md`](docs/03-fpga-implementation.md) | **Task 3** — block diagram, clocking, word lengths, resources |
| [`docs/04-rf-frontend.md`](docs/04-rf-frontend.md) | **Task 4** — AD9361 vs a discrete front end |
| [`docs/05-algorithm-survey.md`](docs/05-algorithm-survey.md) | **Task 5** — MVDR, LCMV, MUSIC, STAP, OSNMA, RAIM, ML |
| [`docs/06-product-architecture.md`](docs/06-product-architecture.md) | **Task 6** — productisation, calibration, test, failure modes |

---

## The three changes that matter most

1. **Add a detector.** The published pipeline nulls unconditionally, so with
   no spoofer present it steers a null into the strongest authentic
   satellite. An anti-spoofing device spends nearly all of its operating
   hours not under attack, so this failure mode dominates its expected
   value. The eigenvalue detector costs nothing — the eigenvalues come out
   of the same decomposition that produces the weights.

2. **Replace the γ/β estimator with the whitened covariance eigenvector.**
   The paper's estimator turns out to be exactly *one power-iteration step*
   applied to the whitened covariance, started from the reference element.
   Running the iteration to convergence costs 3 extra correlators at N = 4
   and buys 9.5 dB of null depth at a 2 ms dwell, removes ~58 BRAM36 of
   epoch-delay memory, and yields the detector, rank estimation and R⁻¹
   for free.

3. **Never round inside the covariance accumulator.** Rounding sample
   products before accumulation injects a *deterministic, dwell-length-
   independent* rank-one error aligned with the boresight steering vector —
   a phantom source at zenith, which is where the satellites are. It caps
   null depth at 20·log₁₀(2^−FL_prod) no matter how long you integrate, and
   it is invisible in floating-point simulation. Exact accumulation in a
   48-bit DSP48 register removes it at zero cost.

---

## Reference model

```matlab
run('matlab/asp_startup.m')

asp_run_all('tests')    % regression suite,           ~1 min
asp_run_all()           % + short studies,            ~5 min
asp_run_all('full')     % + full Monte Carlo studies, ~40 min
```

No toolboxes. Runs unmodified on MATLAB R2018b+ and GNU Octave 7+.

```
matlab/
  config/    single source of truth; powers as C/N0, word lengths derived
  fx/        fixed-point primitives, bit-exact DSP48 covariance model
  model/     ICD-verified C/A codes, geometries, channel and RF impairments
  core/      THE ALGORITHM REFERENCE — maps 1:1 to RTL blocks
  analysis/  metrics, closed-form theory, covariance-domain scene model
  verify/    regression tests
  studies/   the Monte Carlo evidence behind every number in the docs
  export/    bit-exact RTL co-simulation vectors
```

### Design rules the reference obeys

* Deterministic — fixed iteration counts everywhere, no tolerance loops.
* No dividers or square roots on the sample path; all normalisation is
  power-of-two block floating point.
* No matrix inversion anywhere; MVDR and LCMV come from the eigen-decomposition.
* Norm-preserving numerics (unitary Jacobi rotations), so no dynamic-range growth.
* Exact accumulation: full-width products, one rounding at the output.
* Causal — weights from dwell *k* are applied to dwell *k+1*.
* Hardware in the loop: AGC, quantisation, saturation and dwell validity are
  inside the processing loop, not applied afterwards.

---

## Key measured results

4-element Y array (3 on a ring at 120° + 1 centre, 0.45 λ radius), SAPR 5.5 dB,
C/N₀ 45 dB-Hz, full waveform model with real Gold codes, navigation data,
code Doppler and channel mismatch.

| | paper γ/β | this model |
|---|---|---|
| Null depth @ 2 ms | −15.8 dB | **−25.3 dB** |
| Null depth @ 20 ms | −22.2 dB | −26.8 dB |
| First answer available | 2 ms | **1 ms** |
| Epoch delay memory | ~58 BRAM36 | **0** |
| Spoofing detector | none | **λ₁/λ̄ + MDL, free** |
| Min detectable SAPR @ 1 ms | n/a | **+0.3 dB** |
| Dividers on the sample path | 3 | **0** |

**Element count**, all geometries sized to a common 0.5 λ radome:

| | N = 3 | N = 4 (Y) |
|---|---|---|
| Array gain, 1 null | 2.99 dB | **4.68 dB** |
| Array gain, 2 nulls | −0.04 dB | **+2.54 dB** |
| Satellites left worse than one antenna, 1 null | 18.6% | **7.0%** |
| Satellites left worse than one antenna, 2 nulls | 59.9% | **22.3%** |
| P(detect) at 0 dB SAPR, 1 ms | 0.447 | **0.882** |
| AD9361 count | 2 (one RX wasted) | 2 (none wasted) |

**Post-correlation refinement** obeys, to better than 0.1 dB over a 20 dB sweep:

```
    stage-2 null depth  =  -(SAPR + 27.0) dB
```

where 27.0 is the C/A cross-correlation bound (23.9 dB) plus ~3.1 dB.

---

## Provenance

Every number in the documents comes from a script in `matlab/verify/` or
`matlab/studies/`. Where a measurement contradicted an initial hypothesis —
and several did — the document records the correction rather than the
hypothesis. Three worth flagging:

* A specular **ground bounce costs nothing**. A planar array cannot
  distinguish elevation +θ from −θ (measured coherence 1.000000), so the
  reflection shares the direct path's steering vector and is co-nulled for
  free. The rank-2 cases that *do* cost a degree of freedom are arrivals at a
  different **azimuth**.
* **Navigation data bits do not damage** the paper's β statistic, because the
  affected factor is common to all elements and cancels in the scale-invariant
  projector. Bit synchronisation is not a prerequisite for the array processing.
* The paper's **magnitude/phase split is essential, not stylistic**. Using the
  raw covariance column directly gives ρ = 0.74 and a −6.5 dB null, because
  R₁₁ is dominated by noise power while R_i1 contains only the spoofer term.

---

## Reference

Daneshmand, S., A. Jafarnia-Jahromi, A. Broumandan, G. Lachapelle (2012),
*A Low-Complexity GPS Anti-Spoofing Method Using a Multi-Antenna Array*,
ION GNSS 2012, Session B3, Nashville TN, 18–21 September 2012.
