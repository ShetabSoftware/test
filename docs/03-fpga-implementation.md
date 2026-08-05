# Task 3 — FPGA Implementation Strategy

Target: Zynq-7020/7030 or Zynq UltraScale+ ZU3EG, with two AD9361s sharing one LO. Everything below is sized against the reference model in `matlab/`, which is the golden model for RTL verification — not a document that will drift from it.

**Headline sizing result: the entire array processor is about 20 DSP48 slices, zero BRAM, and 15 µs of once-per-millisecond compute.** It is 2% of a Zynq-7020. The FPGA cost of this algorithm is not a design constraint, and any architecture discussion that treats it as one is optimising the wrong thing. The real constraints are the *receiver* correlator bank, the channel equalisers, and the discipline required to avoid the numerical traps in §7.

---

## 1. Block diagram

```
  4x ANTENNA ─┬─ LNA (0.8 dB NF, 30 dB) ─ SAW ─┬── AD9361 #0 RXA ──┐
              │                                └── AD9361 #0 RXB ──┤  LVDS
              │  (external LO fanout, common)                      │  DDR
              ├─ LNA ─ SAW ────────────────────┬── AD9361 #1 RXA ──┤
              └─ LNA ─ SAW ────────────────────┴── AD9361 #1 RXB ──┘
                                                                    │
  ══════════════════════════ FPGA PL ═══════════════════════════════╪═══════
                                                                    ▼
    CLK DOMAIN A: data-clock (DDR LVDS, 61.44 MHz)     ┌──────────────────┐
                                                        │ axi_ad9361 x2    │
                                                        │ deserialise      │
                                                        └────────┬─────────┘
                                                                 │ 4x IQ @ 30.72 MSPS
    CLK DOMAIN B: fabric, 163.68 MHz  ◄── async FIFO ────────────┘
                                                                 │
    ┌────────────────────────────────────────────────────────────▼──────────┐
    │ DIGITAL FRONT END, per channel, 4 instances                           │
    │  DC notch ─ channel equaliser FIR (16 tap, prod. cal) ─ decim /2 ─┐   │
    │  saturation counter ──────────────────────────────────────────┐   │   │
    └──────────────────────────────────────────────────────────────┬┴───┴───┘
                                                4x IQ @ 16.368 MSPS │ s16.15
                 ┌──────────────────────────────────────────────────┤
                 │                                                  │
                 ▼                                                  ▼
    ┌────────────────────────────┐                    ┌───────────────────────┐
    │ COVARIANCE ENGINE          │                    │ BEAMFORMER BANK       │
    │  10 accumulators, s48      │                    │  v_m = f_m^H r        │
    │  exact products, no round  │                    │  12 beams, 12 DSP     │
    │  3 DSP48, 0 BRAM           │                    │  4 cycle latency      │
    └─────────────┬──────────────┘                    └───────────┬───────────┘
                  │ R, once per ms                                │
                  ▼                                               │
    ┌────────────────────────────────────────────┐                │
    │ WEIGHT ENGINE  (1 kHz, 15 us, time-shared) │                │
    │  1. trace normalise, BFP exponent          │                │
    │  2. whiten: 4x CORDIC inverse-sqrt         │                │
    │  3. Jacobi EVD, 6 sweeps x 6 rotations     │                │
    │  4. detector: lambda1/mean(lambda2..N)     │                │
    │  5. weights: f = h - Q(Q^H h)              │                │
    │  ~4 DSP48 + 1 CORDIC                       │────────────────┘
    └──────────────┬─────────────────────────────┘
                   │ eigenvalues, rank                            │
    ═══════════════▼══════════════════════════════════════════════▼══════════
    PS (ARM) or MicroBlaze, 1 kHz:                        CORRELATOR BANK
      MDL rank decision (needs log)                       48-192 channels
      threshold adaptation, hysteresis                    E/P/L, per PRN
      mode/state machine, health, logging                       │
      stage-2 post-correlation weights ◄────────────────────────┘
```

---

## 2. Clock domains

| Domain | Rate | Contents | Crossing |
|---|---|---|---|
| A — AD9361 DATA_CLK | 61.44 MHz DDR | LVDS deserialiser, one per chip | async FIFO into B |
| B — fabric | **163.68 MHz** = 10 × 16.368 | everything on the sample path | — |
| C — AXI/control | 100 MHz | registers, DMA, PS interface | AXI4-Lite CDC |
| D — PS | 667 MHz+ | ARM, mode logic, MDL, logging | AXI |

**Choose the fabric clock as an exact integer multiple of the sample rate.** 163.68 MHz = 10 × 16.368 MHz gives a time-multiplexing ratio of exactly 10, so every time-shared resource has a fixed, statically-scheduled slot assignment with no rate-matching FIFO, no fractional accumulator, and no jitter in the TDM phase. A "nice round" 200 MHz clock would give a ratio of 12.22 and force elastic buffering into every shared block for no benefit. This is the single most useful clocking decision in the design.

**Both AD9361s must run from one reference and one LO.** Distribute a single TCXO/OCXO reference and use the external-LO input on both chips, driven from one synthesiser through a matched-length splitter. Then LO phase noise is *common mode* across all four channels and cannot degrade null depth at all — a common complex scalar multiplies every channel equally and leaves the spatial covariance structure untouched. With independent internal PLLs, the differential phase noise caps the null at roughly 20·log₁₀(σ_φ): 0.5° RMS → −41 dB, 1° → −35 dB. See Task 4.

Use the AD9361 multi-chip sync (MCS) to align the baseband dividers, and **verify sample alignment on every boot**. A one-sample offset between chips at 16.368 MHz is 61 ns of skew — utterly fatal to a narrowband beamformer. Inject a common pilot through a splitter at power-up and cross-correlate; if the lag is non-zero, re-run MCS. This must be an automatic, logged, boot-time self-test, not a bring-up procedure.

---

## 3. Sampling rate

**Recommend fs = 16.368 MHz = 16 × 1.023 MHz.**

| candidate | samples/chip | samples/ms | verdict |
|---|---|---|---|
| 10.000 MHz (paper, your scripts) | 9.775 | 10000 | non-integer chips; forces a code NCO for no benefit |
| 4.092 MHz | 4 | 4092 | integer; adequate for nulling, thin for code tracking |
| **16.368 MHz** | **16** | **16368** | **integer both ways; 8× the C/A main lobe** |
| 30.72 MHz | 30.03 | 30720 | non-integer; an ADI-convenient rate, not a GNSS one |

Integer samples per chip *and* per code period removes the code-NCO phase-accumulator ambiguity from the reference model and makes RTL golden vectors trivially comparable to it. That is worth real money during verification.

**But there is a bandwidth trade in the other direction, and it is not the usual one.** Wider processing bandwidth stresses the narrowband beamformer assumption. Two mechanisms, with very different magnitudes:

* *Aperture dispersion.* The wavefront crosses the array in D/c = 0.82λ/c = 0.52 ns. Over a band B the resulting null floor is ≈ 20·log₁₀(2π(B/√12)·τ_rms). At 16.368 MHz with τ_rms ≈ 0.26 ns this is **−46 dB** — comfortably below the estimator's own −24 dB, so *not* the binding constraint.
* *Inter-channel group-delay skew.* Same formula, but τ is the cable/filter/transceiver mismatch. 100 ps → **−55 dB**; 1 ns → **−35 dB**. This *is* binding if uncalibrated, and it is the reason for the channel equaliser in §4.

I originally expected a third mechanism — decorrelation of a multipath replica across the band — to dominate. Measurement says it does not, because the relevant decorrelation is set by the **signal** bandwidth (1.023 MHz, a 977 ns correlation width), not the sampling rate. Widening the ADC does not make it worse.

---

## 4. Digital front end

Per channel, at 16.368 MHz:

| block | purpose | cost |
|---|---|---|
| DC notch | AD9361 is zero-IF; DC sits on the C/A spectral peak | 1 CIC-based HPF, ~0 DSP |
| **channel equaliser, 16-tap complex FIR** | **flatten inter-channel amplitude and group-delay mismatch across the band** | 16 complex mult TDM'd → 6 DSP for all 4 channels |
| halfband decimator /2 | 30.72 → 16.368 (retune BBPLL to 32.736) | 3 DSP for all 4 channels |
| saturation counter | invalidates the dwell | LUTs |

**The channel equaliser is not optional and it is the highest-value block in the DFE.** Without it, inter-channel group-delay skew caps the null at −35 dB for 1 ns of mismatch, which is easy to accumulate across two separate transceiver chips with independently calibrated baseband filters. With a 16-tap complex equaliser per channel, loaded from a production calibration table, the residual is limited by the calibration accuracy rather than by the hardware.

Calibrate it by injecting a common wideband pilot (or noise) through a matched splitter into all four inputs at manufacturing test, measuring H_i(f) per channel, and solving for taps that equalise all channels to the reference channel. Store per-unit in eFUSE or QSPI. Re-verify at each gain index — see §8.

Use *identical* coefficients for the decimator on every channel so group delay matches exactly by construction. Never use a channel-specific filter anywhere on the sample path.

---

## 5. Covariance engine

The core observation that makes this cheap: only the upper triangle is needed, and the DSP48's native widths fit it exactly.

* Products: 16 × 16 = **32 bits**, held exactly.
* Coherent growth over K = 16368: **+14 bits**.
* Guard: **+2 bits**.
* Total: **48 bits — precisely the DSP48 P register.**

Measured on the bit-exact model: 42 bits used, 6 spare. So one DSP48 per covariance entry, with the accumulator *inside the slice*, no fabric adder, **and no rounding anywhere in the loop**. §7 explains why the last point is the important one.

Resources, N = 4, TDM ratio 10:

| entry type | count | real multiplies | DSP48 |
|---|---|---|---|
| diagonal (real, \|x\|²) | 4 | 8 | 1 |
| off-diagonal (complex) | 6 | 18 (Karatsuba) or 24 | 2 |
| **total** | **10** | **26** | **3** |

Use the 3-multiply complex product,

```
k1 = c(a+b),  k2 = a(d−c),  k3 = b(c+d)
real = k1 − k3,  imag = k1 + k2
```

which maps directly onto the DSP48 pre-adder — the pre-adder does the (a+b), (d−c), (c+d) for free.

Dwell control: accumulate for K samples, latch to a shadow register, zero, continue. Latching (not stalling) keeps the pipeline gapless. A dwell is marked **invalid** if a saturation event or an AGC gain step occurred within it; the weight engine then holds the previous weights rather than consuming a poisoned estimate.

---

## 6. Weight engine — 1 kHz, ~15 µs, fully shared

Runs once per dwell. Everything here is time-shared and the resource cost is negligible.

**Step 1 — trace normalise.** Shift R so trace(R) = 2^k. Carries the exponent separately (block floating point), so the mantissa is not spent on the absolute signal level. Cost: leading-zero count + barrel shift.

**Step 2 — whiten.** D⁻¹RD⁻¹ with D = diag(√R_ii). Four CORDIC inverse-square-roots plus 16 scalings, once per millisecond. This is what preserves the calibration-free property under an eigen-based estimator — see Task 2, M2.

**Step 3 — Jacobi eigen-decomposition.** Fixed 6 sweeps × 6 rotations for N = 4. Each rotation is:

* 1 CORDIC **vectoring** to find the phase that makes the off-diagonal real;
* 1 CORDIC vectoring to find the rotation angle θ = ½·atan2(2|c|, a−b);
* CORDIC **rotations** applied to two columns of A, two rows of A, and two columns of U.

Why Jacobi and nothing else:

* Every operation is a unitary similarity transform, so ‖A‖_F is invariant **exactly** — verified to 1.05×10⁻¹⁵. **Zero dynamic-range growth**, so one word length for the whole block with no rescaling and no overflow analysis beyond the input.
* Unconditionally stable regardless of conditioning — which matters, because R here *is* badly conditioned in the sense that counts: at 5.5 dB SAPR the useful signal is a 25% perturbation of the identity.
* Pure CORDIC. No multipliers at all if you want none.
* **Fixed latency.** No tolerance test, no data-dependent iteration count. Verified: 4 sweeps already reach 4×10⁻¹⁶ at N = 4; 6 is chosen for margin at N = 8.
* Natural parallelism: for N = 4 the three disjoint pairings {(1,2),(3,4)}, {(1,3),(2,4)}, {(1,4),(2,3)} can each be applied simultaneously, so one sweep is **three** parallel double-rotation steps. A sweep-parallel implementation cuts the 15 µs to ~5 µs if you ever need it. You will not.

Fixed-point: 18-bit rotation coefficients track floating point to sin θ < 5×10⁻⁵ (−86 dB), which is **60 dB below** the statistical floor of the estimate. The EVD is nowhere near being the bottleneck.

**Step 4 — detection.** λ₁/mean(λ₂…λ_N) against a threshold. Trivial: one divide (or a compare against threshold × mean, avoiding the divide entirely) once per millisecond.

**Step 5 — weights.** f = h − Q(Qᴴh) with Q from modified Gram-Schmidt on the rank-p signature. Never form P. 2·N·rank multiplies instead of N². For rank 1 at N = 4 that is 8 instead of 16.

**No normalisation.** Power-of-two block-floating-point scaling only (`fx_bfp_scale`): leading-zero count plus barrel shift, ~30 LUTs, one cycle, zero added error. This removes all three inverse-square-roots that the original scripts have on this path.

**Timing budget**

| step | cycles @ 163.68 MHz |
|---|---|
| trace normalise | ~20 |
| whiten (4 CORDIC rsqrt + 16 scale) | ~200 |
| Jacobi, 36 rotations × ~70 cycles | ~2520 |
| detect | ~30 |
| weights | ~120 |
| **total** | **~2890 cycles = 17.7 µs** |

**1.8% duty cycle against a 1 ms dwell.** One instance serves the whole system.

---

## 7. Fixed-point plan, and the trap that must not be missed

| stage | format | derivation |
|---|---|---|
| AD9361 output | s12.11 | native |
| DFE / datapath sample | **s16.15** | 12 ADC bits + ~2 decimation + 2 guard; fits DSP48 A port and Intel 18×19 |
| equaliser taps | s18.16 | DSP48 B port |
| covariance product | s32 exact | 16×16 |
| **covariance accumulator** | **s48** | 32 + ⌈log₂K⌉ + 2 = DSP48 P register exactly |
| EVD datapath | s32.30 + BFP exponent | must resolve λ₁/λ_N up to 10⁶ for 60 dB J/N |
| CORDIC rotation coefficients | s18.16 | −86 dB, 60 dB below the statistical floor |
| **beamformer weights** | **s18.16 (Q1.16)** | null floor 20log₁₀(2⁻¹⁶) − 1.76 = **−98 dB** (measured −97.4) |
| beamformer accumulator | s37 | 34 + 2 + 1 |
| beamformer output | s16.13 | to the correlator |

**AGC operating point: hold the composite (thermal + interference) at 14 dB below full scale.** That gives 5σ of crest headroom, i.e. a per-component clip probability near 6×10⁻⁷.

### 7.1 Never round inside the covariance accumulator

This is the one implementation detail that silently destroys the system, and it does not show up in floating-point simulation.

Rounding sample products before the accumulator injects an error into R that is:

1. **deterministic** — a fixed matrix, not noise;
2. **rank one, aligned with the boresight steering vector** — measured alignment 1.000 with the all-ones vector when the RTL computes both triangles and Hermitises, 0.906 when it computes the upper triangle and mirrors; a random direction gives 0.5. The reason is immediate once stated: a constant bias *c* on every product puts the same *c* in every entry of R, and ones(N) is rank one with the boresight steering vector as its eigenvector. **A truncating implementation synthesises a phantom source at zenith and steers the null into the sky.**
3. **independent of dwell length** — measured ‖E‖_F identical to three significant figures across a 16× change in K.

Property 3 is the trap. Everything else improves as 1/√K, so a designer who validates at 1 ms expects improvement at 100 ms. With product rounding the floor is flat.

| product fractional bits retained | 10 | 14 | 18 | 22 | 26 | full |
|---|---|---|---|---|---|---|
| null-depth floor | −25.4 | −48.9 | −73.0 | −97.1 | −122.0 | −156 dB |

6.02 dB per retained bit. **The fix is free**: full-width products into the 48-bit P register, no rounding. If a narrower accumulator is ever forced, use convergent (round-half-to-even) rounding, which restores 1/√K behaviour — measured 0.28× over a 16× change in K, against 0.25 predicted.

`verify/test_fixedpoint.m` is the regression test. Run it against the RTL.

### 7.2 Converter word length sets your anti-jam capability

After nulling, the interferer is gone but the quantisation noise it forced the AGC to expose is not — quantisation noise is spatially white and survives every beamformer:

$$L = 10\log_{10}\!\left(1+\frac{1/12}{\sigma_\text{th}^2}\right),\qquad \sigma_\text{th}=2^{WL-1}10^{-B/20}10^{-\mathrm{J/N}/20}$$

| word length | 8 | 10 | **12 (AD9361)** | 14 | **16** |
|---|---|---|---|---|---|
| usable J/N at 1 dB loss | 33 | 45 | **57 dB** | 69 | **81 dB** |

Derate the AD9361 to ~51 dB for its ~10.5-bit ENOB. Adequate for spoofing (SAPR 0–30 dB); marginal for jamming. Task 4.

### 7.3 Saturation is not "distortion"

Clipping is a memoryless nonlinearity that redistributes a strong interferer across the whole array manifold. The interferer stops being rank one and **stops being nullable**. Any dwell containing saturation must be *discarded*, not tolerated. Count saturations per channel per dwell, invalidate the dwell, hold weights, and expose the count as a health telemetry point.

---

## 8. What stays in software

| function | where | why |
|---|---|---|
| Streaming covariance, beamforming, correlators | **PL** | rate-driven, trivially parallel, would need >10 GMAC/s in software |
| Jacobi EVD, projector, weight solve | **PL** | 1 kHz but latency-critical and cheap in CORDIC |
| Eigenvalue detection statistic + threshold compare | **PL** | one compare, gates the weights in the same cycle |
| **MDL/AIC rank decision** | **PS** | needs logarithms; 1 kHz; no latency pressure |
| Threshold adaptation, hysteresis, dwell-time policy | **PS** | policy, not arithmetic; must be field-updatable |
| Acquisition search management | **PS** | irregular control flow |
| Tracking loop filters | **PS** or PL | 50–1000 Hz; PS is fine and far easier to tune |
| Stage-2 post-correlation weight computation | **PS** | operates on correlator outputs at 1 kHz, needs flexible per-PRN bookkeeping |
| PVT, RAIM, integrity monitoring | **PS** | 1–10 Hz |
| AD9361 calibration sequencing, gain-index management | **PS** | SPI-rate |
| Health, logging, forensic snapshot capture | **PS** + DDR | |

**Rule of thumb that has served me well:** anything at the sample rate goes in fabric; anything at the dwell rate goes in fabric only if it is on the weight-update critical path; anything that needs a logarithm, a decision policy, or a field update goes in the PS.

One deliberate exception. The **MDL rank estimate** could go in fabric, but keeping it in software means the rank policy — including hysteresis, and refusing to spend a second null unless the second eigenvalue has been elevated for several consecutive dwells — can be tuned in the field against real threats. Rank decisions are where an adaptive nulling system most often misbehaves, and you want that knob in software.

---

## 9. Parallelisation

| opportunity | speedup | worth it? |
|---|---|---|
| Covariance entries fully parallel | 10× | yes, and free — it is 3 DSPs either way |
| Per-PRN beams fully parallel | 12× | yes, one DSP each |
| Jacobi disjoint index pairs | 3× (N=4) | **no** — 17.7 µs is already 1.8% duty |
| Multiple dwells in flight | 2× | **no** — adds weight age, which is a first-order effect |
| Correlator bank across PRNs | 48–192× | yes; this is the real resource consumer |

**Where not to parallelise.** Pipelining two covariance dwells to halve latency would add a full dwell of weight age. Weight age caps the null at 20 log₁₀(2π(D/λ)Δθ/√3) on a rotating platform: −46 dB at 100°/s, −34 dB at 400°/s, −26 dB at 1000°/s with a 1 ms dwell. Doubling it costs 6 dB of that budget to save 17 µs out of 1000. Do not.

---

## 10. Latency and throughput

**Latency, antenna to receiver input**

| stage | latency |
|---|---|
| AD9361 analog + digital filter chain | ~1.5 µs |
| LVDS deserialise + CDC FIFO | 0.3 µs |
| DC notch | 0.1 µs |
| channel equaliser (16 tap) | 0.5 µs |
| halfband decimator | 0.4 µs |
| beamformer | 0.02 µs |
| **total** | **≈ 2.8 µs** |

Two things follow that matter commercially:

* 2.8 µs is **840 m of apparent range**. It is *common to all channels*, so it cancels in the position solution — but it does **not** cancel in a timing product. If you sell a timing receiver, this delay must be characterised, temperature-compensated and specified, and its *variation* (not its value) is what limits your timing accuracy.
* The delay must be **deterministic across power cycles**. An elastic CDC FIFO whose fill level varies at reset makes it non-deterministic. Use a fixed-latency CDC with a deterministic reset sequence, and verify the total delay at boot against the pilot injection of §2.

**Weight update latency:** 1 ms dwell + 17.7 µs = 1.018 ms.

**Throughput:** 16.368 MSPS × 4 channels sustained, no backpressure, no gaps. The covariance engine latches to a shadow register rather than stalling.

**Convergence at attack onset:** the detector fires on the first dwell whose statistic clears threshold. With hysteresis (recommend 2 of 3 consecutive dwells before engaging the null, 10 consecutive below before releasing), attack response is **2–3 ms**. Compare to the tens of seconds a spoofer needs to walk a receiver off its true position — three orders of margin, which is where it should be.

---

## 11. Resource summary

Zynq-7020 (85k logic cells, 220 DSP48E1, 140 BRAM36):

| block | LUT | FF | DSP48 | BRAM36 |
|---|---|---|---|---|
| 2× axi_ad9361 interface | 4 000 | 6 000 | 0 | 4 |
| DFE (4 ch: notch, equaliser, decimator) | 3 500 | 5 000 | 9 | 0 |
| Covariance engine | 900 | 1 400 | 3 | 0 |
| Weight engine (CORDIC, Jacobi, detector) | 4 500 | 5 500 | 4 | 1 |
| Beamformer bank (12 beams) | 2 200 | 3 000 | 12 | 0 |
| **Array processor subtotal** | **15 100** | **20 900** | **28** | **5** |
| | | | **13%** | **4%** |
| Correlator bank, 48 ch (4 streams × 12 PRN) | 12 000 | 18 000 | 96 | 12 |
| DMA, AXI, control | 5 000 | 7 000 | 0 | 8 |
| **Total** | **32 100** | **45 900** | **124** | **25** |
| | 60% | 43% | **56%** | 18% |

Comfortable on a 7020 for a 12-PRN single-frequency product. For a full multi-constellation receiver the correlator bank dominates and you want a 7030/7045 or a ZU3EG — but note that the *anti-spoofing* function is 13% of the DSPs either way. **Size the FPGA for the receiver, not for the array processor.**

---

## 12. Verification strategy

1. **Reference model is golden.** `matlab/` produces bit-exact expected outputs. `export/asp_export_vectors.m` writes stimulus and expected-response files for RTL co-simulation.
2. **Bit-exact comparison at every boundary**: DFE output, covariance accumulator contents, eigenvalues, eigenvectors, weights, beamformer output. Not "close" — identical, because every block in `core/` is specified to be reproducible in integer arithmetic.
3. **The four regression tests that must run in CI**: `test_ca_code` (ICD conformance), `test_evd_herm` (against `eig()` and against the fixed-point datapath), `test_fixedpoint` (the accumulator trap), `test_cov_model` (statistical model agreement).
4. **Hardware-in-the-loop**: replay the recorded IF from `study_*` scenarios through the AD9361 in loopback, or through a Spirent/Orolia simulator combined per the paper's own test setup (Figure 8) — a splitter feeding all channels for the spoofer and an over-the-air array for the authentic signals. That test setup is the correct one and is worth reproducing exactly.
5. **Boot-time self-test**: pilot injection, inter-channel sample alignment, group-delay verification, gain-index phase table validation.
6. **Negative testing** — the part most teams skip and the part that matters here: no spoofer present (verify bypass), spoofer at zenith (verify graceful degeneracy handling), saturation (verify dwell invalidation), AGC step mid-dwell (verify weight hold), one channel dead (verify degradation to N−1), all channels dead (verify safe output).

---

## 13. Things that will bite you

| | Issue | Mitigation |
|---|---|---|
| 1 | **AD9361 per-channel AGC** moves gain indices independently. Each move is an amplitude *and* a several-degree phase step. The array signature changes mid-dwell. | Force MGC or a single common gain index. Model it: set `cfg.fe.commonAGC = false` and watch the null collapse. |
| 2 | **Gain-index phase steps.** Even a common gain change shifts phase. | Per-gain-index phase-correction LUT from production calibration; invalidate the dwell containing the step. |
| 3 | **AD9361 tracking calibrations** (QEC, DC) run per-channel and asynchronously, changing the effective channel response mid-dwell. | Run RF/BB cal once at init, then **freeze tracking cal**. Re-run only on a commanded band/temperature change, and invalidate dwells around it. |
| 4 | **Zero-IF DC offset** sits on the C/A spectral peak. | Tune the LO 1–2 MHz off L1 (low-IF) and digitally downconvert, or accept ZIF with a narrow digital notch. Low-IF is cleaner. |
| 5 | **I/Q image.** A strong interferer's image lands at the mirror frequency with a *different* spatial signature (QEC differs per channel) and therefore cannot be nulled. IRR sets a hard floor on interference suppression. | Budget ≥ 60 dB IRR; verify over temperature with frozen QEC. |
| 6 | **Weight-update phase discontinuities.** Every weight update changes the output carrier phase by a per-satellite amount. At 1 kHz this destroys carrier tracking and makes RTK impossible. | See Task 6 — either constrain fᴴf_prev to be real positive (phase-continuous update), or move combining post-correlation where the phase is tracked per satellite. **This is the classic CRPA problem and it is not addressed anywhere in the paper or your scripts.** |
| 7 | **Inter-chip sample misalignment** after MCS. | Boot-time pilot cross-correlation, automatic re-sync, logged. |
| 8 | **Truncation in the covariance accumulator.** | §7.1. Regression test in CI. |
