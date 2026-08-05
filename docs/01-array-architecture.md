# Task 1 — Antenna Array Architecture: 3 Elements vs 4

**Verdict up front.** Go to four elements. Not primarily for the 1.7 dB of extra array gain — that is the least interesting benefit — but because four elements are the smallest array that can survive the *realistic* threat, which is a spoofer plus one more spatially distinct source (a reflection off a building, a second emitter, or an accompanying jammer). Against that threat a three-element array has **zero** adaptive degrees of freedom left, delivers exactly single-antenna array gain (measured −0.04 dB), and leaves 60% of the constellation worse off than a single passive antenna. Four elements delivers +2.5 dB and 22%.

And the hardware is free. The AD9361 is a 2×2 transceiver, so three antennas requires two chips with one receive path wasted. Four antennas requires the same two chips with nothing wasted. **The RF bill of materials for 3 and 4 elements is identical.** There is no version of this trade study where three wins.

Below is the full analysis, with the measured numbers from `matlab/studies/study_array_size.m` and `study_detector.m`.

---

## 1. Assumptions stated explicitly

| Assumption | Value | Why |
|---|---|---|
| Carrier | GPS L1 C/A, 1575.42 MHz, λ = 190.3 mm | Paper's target |
| Authentic C/N₀ | 45 dB-Hz per PRN, 9 satellites | Paper's simulation |
| SAPR | 5.5 dB per PRN (spoof −153 dBW vs authentic −158.5 dBW) | Paper's simulation |
| Spoofer elevation | 15° | A terrestrial spoofer is a terrestrial object. The paper places it at 45°, which is optimistic — see §7 |
| Radome | 0.5 λ radius (95 mm, a 190 mm puck) | Typical commercial CRPA envelope |
| Grating margin | ≤ −3 dB peak off-mainlobe array-factor response | See §4 |
| Dwell | 1 ms (K = 16368 at fs = 16.368 MHz) | Paper's K = 10000 |

Every geometry below is sized to the **same radome and the same grating margin**. This matters: comparing a 3-element and a 4-element array "both at λ/2 spacing" compares two arrays with different apertures and different ambiguity margins, then attributes the whole difference to the element count.

---

## 2. Degrees of freedom — the argument that actually decides it

An N-element narrowband beamformer has N complex weights, hence **N − 1 independent degrees of freedom** after the overall scale (which is unobservable — see the note on normalisation in Task 2). The ledger for an LCMV formulation with unit gain toward the wanted satellite and hard nulls on the interference is:

```
                                    N = 3      N = 4
  desired-signal constraint            1          1
  spoofer null                         1          1
  ------------------------------------------------------
  remaining adaptive DOF               1          2
  
  ... with a second emitter or reflector (rank-2 spoofing subspace):
  desired-signal constraint            1          1
  spoofer + second-arrival nulls       2          2
  ------------------------------------------------------
  remaining adaptive DOF               0          1
```

Going from 3 to 4 **doubles** the spare adaptive degrees of freedom against a single spoofer, and is the difference between zero and one against a rank-2 spoofing subspace. Zero spare DOF means the beamformer is fully determined by its constraints: it cannot respond to anything it was not explicitly told about.

### 2.1 What makes the spoofing subspace rank 2 — and what does not

I initially assumed the answer was the spoofer's ground bounce. **It is not, and the reason is a genuinely useful property of the planar geometry that is worth knowing.**

For a planar array with every element at z = 0, the response depends only on the direction cosines (cos E·cos A, cos E·sin A). Cosine is even, so **elevation +θ and −θ produce an identical steering vector.** Measured coherence between (az 45°, el +15°) and (az 45°, el −15°): **1.000000** — not approximately, exactly. The array physically cannot distinguish up from down.

A specular ground reflection arrives at exactly the mirror elevation and the same azimuth. It therefore shares the direct path's steering vector at every frequency, and the same rank-one projector that removes the direct path removes the reflection **as a free side effect** — no extra degree of freedom, no wideband penalty. Measured: rank-1 and rank-2 nulls perform identically when the second arrival is a mirror-elevation ground bounce, and spending the second null merely wastes array gain.

This cuts both ways, and the downside belongs in your product documentation: **a planar CRPA provides no spatial multipath rejection for the *authentic* signals either.** It cannot separate a satellite at +θ from its own ground reflection. That job stays with the receiver's code and carrier multipath mitigation, or requires a vertical baseline (a non-planar array).

The cases that genuinely consume a degree of freedom are arrivals at a **different azimuth**:

* a reflection off a building, vehicle, mast or the platform's own structure;
* a second spoofing transmitter;
* an accompanying jammer — and jam-then-spoof is a common attack pattern, because forcing reacquisition is how you get a receiver to accept counterfeit signals quickly.

Measured coherence at 90° of azimuth separation is 0.010, i.e. essentially orthogonal. The rank-2 rows below use that case.

### 2.1 What the DOF ledger costs you, quantitatively

For an array response **a** drawn uniformly on the complex sphere with ‖a‖² = N, and a rank-p orthogonal projector,

$$\frac{\|\mathbf{P}\mathbf{a}\|^2}{N} \sim \mathrm{Beta}(N-p,\;p)$$

Two consequences follow immediately.

**Mean array gain after nulling** (the power-maximisation case, paper eq. 23):

$$\mathbb{E}\left[\|\mathbf{P}\mathbf{a}\|^2\right] = N - p$$

so the array gain over a single antenna is exactly **10 log₁₀(N − p) dB**:

| | N = 3 | N = 4 | Δ |
|---|---|---|---|
| p = 1 (spoofer only), theory | 3.01 dB | 4.77 dB | **+1.76 dB** |
| p = 1, **measured** (y4) | **2.99 dB** | **4.68 dB** | **+1.69 dB** |
| p = 2 (spoofer + azimuthally separate source), theory | 0.00 dB | 3.01 dB | **+3.01 dB** |
| p = 2, **measured** (y4) | **−0.04 dB** | **+2.54 dB** | **+2.58 dB** |

The p = 2 row says it plainly: **with two nulls to place, a 3-element array delivers exactly single-antenna array gain.** −0.04 dB measured against 0.00 dB predicted. Every element beyond the first has been consumed by the constraints, and the array has stopped being an array.

**Outage probability** — the fraction of satellites left with *less* gain than a single passive antenna:

$$P_\text{out} = I_{1/N}(N-p,\;p)$$

| | N = 3 | N = 4 (y4) | N = 7 |
|---|---|---|---|
| p = 1, theory | 11.1% | 1.6% | 0.02% |
| p = 1, **measured** | **18.6%** | **7.0%** | 3.9% |
| p = 2, theory | 55.6% | 15.6% | 0.03% |
| p = 2, **measured** | **59.9%** | **22.3%** | 7.9% |

I regard the p = 2 measured row as the most important number in this trade study. **A 3-element array facing a spoofer plus one more emitter leaves 60% of the constellation worse off than if you had thrown the array away and used one antenna.** With 9 satellites in view that is 5.4 satellites degraded. You will not hold a position fix, let alone a RAIM-protected one.

Note also how much worse the measurements are than the isotropic theory. The theory assumes uncorrelated steering vectors; a physical 0.8 λ aperture has strongly correlated ones. **The isotropic model that most array-processing texts hand you is optimistic here, and it is more optimistic the smaller N is.** That asymmetry works against three elements.

---

## 3. Spatial filtering, null steering and beamforming accuracy

### 3.1 Null depth

Two error mechanisms set the achievable null, and they behave completely differently.

**Noise-limited term.** For a rank-one source in white noise estimated from K snapshots, the classical subspace perturbation result (Kaveh & Barabell, *IEEE T-ASSP* 1986) gives

$$\sin^2\theta \;\approx\; \frac{N-1}{K}\cdot\frac{\lambda_1\lambda_n}{(\lambda_1-\lambda_n)^2},\qquad \lambda_1 = \sigma^2 + N P_s,\;\lambda_n=\sigma^2$$

At low SNR (σ² ≫ N P_s) this reduces to (N−1)σ⁴/(K N² P_s²), so the **ratio between element counts is (N−1)/N²**: going 3 → 4 improves the noise-limited null by 10 log₁₀[(2/9)/(3/16)] = **0.75 dB**. Modest, and in the right direction — the array gain N against the coherent spoofer beats the extra N−1 noise-subspace dimensions.

**Bias-limited term.** The covariance also contains the authentic constellation, Σₘ pₘ **aₘaₘ**ᴴ, which tilts the principal eigenvector by

$$\|\delta\mathbf{u}\| \;\approx\; \frac{\sqrt{(N-1)N_\text{auth}}}{N}\cdot\frac{p_a}{P_s}$$

This term **does not depend on K at all.** Integrating longer cannot remove it. Scaling gives √(N−1)/N = 0.471 (N = 3) vs 0.433 (N = 4), another 0.7 dB in favour of four.

**Measured combined null depth**, 1 ms dwell, radome-normalised geometries:

| geometry | N | max baseline | null depth |
|---|---|---|---|
| tri3 | 3 | 0.42 λ | −20.6 dB |
| sq4 (2×2) | 4 | 0.53 λ | −21.3 dB |
| **y4 (3+centre)** | **4** | **0.82 λ** | **−22.2 dB** |
| circ7 | 7 | 0.96 λ | −24.4 dB |
| circ8 | 8 | 1.00 λ | −24.5 dB |

3 → 4 buys **1.6 dB** of null depth, consistent with the 0.75 + 0.7 dB predicted. This is a real but secondary benefit.

### 3.2 The bias floor is the headline result, and it is independent of N

Measured on the full waveform model (`study_estimators.m`), the eigen-based estimator's null depth versus dwell:

| dwell | 1 ms | 2 ms | 5 ms | 10 ms | 20 ms |
|---|---|---|---|---|---|
| null depth | −23.7 | −25.3 | −26.3 | −26.7 | −26.8 dB |

Over that 20× range a purely noise-limited estimator would have gained 13 dB. It gained 3.1 dB. **The pre-correlation estimator is bias limited beyond about 5 ms.**

This has a direct architectural consequence, developed in Task 2 and Task 5: if you want a deeper null than ~27 dB, more antennas and more integration will not get you there. You have to change the *domain* — estimate the spoofing signature from despread correlator outputs, where code and Doppler isolation suppress the authentic contribution by a further 25–40 dB. Measured stage-2 performance reaches 35–45 dB.

### 3.3 Null width — the aperture argument, and why the usual spacing rule is wrong

Null depth is what a paper reports. **Null *width* is what damages your constellation.** The angular extent of a projection null is set by the array factor's mainlobe, which is set by the aperture:

| geometry | N | max baseline | first null of the array factor (direction cosine) |
|---|---|---|---|
| tri3 | 3 | 0.42 λ | **1.46** |
| sq4 | 4 | 0.53 λ | 1.49 |
| **y4** | **4** | **0.82 λ** | **0.95** |
| circ7 | 7 | 0.96 λ | 0.91 |
| circ8 | 8 | 1.00 λ | 0.77 |

The visible sky spans a unit disk in direction cosines. **A value above 1 means the null's mainlobe covers essentially the entire sky** — every satellite is partially attenuated, and the only question is by how much. Both the 3-element triangle and the 2×2 square are in that regime. The Y-array is not.

---

## 4. Array geometry — a specific recommendation that differs from your design

You are proposing a 2×2 square. **Use three elements on a ring at 120° plus one at the centre instead.** Same four elements, same radome, same electronics. The reason is a lattice property, and it is worth spelling out because the usual rule of thumb hides it.

The standard rule — "keep every baseline under λ/2" — is inherited from linear arrays and is **wrong for planar arrays**. For a planar array the response is

$$a_i(u,v) = \exp\!\left(j2\pi (x_i u + y_i v)/\lambda\right), \qquad u^2+v^2 \le 1$$

and two directions are indistinguishable iff their direction-cosine difference **D** satisfies (**p**ᵢ − **p**ⱼ)·**D**/λ ∈ ℤ for *every* element pair. The set of such **D** is the **dual lattice** of the lattice generated by the baseline vectors, scaled by λ. Ambiguity therefore appears when the shortest non-zero dual vector falls inside the visible disk |**D**| ≤ 2 — a property of the *lattice*, not of the longest baseline.

* A **square** lattice of side d has shortest dual vector λ/d, so it aliases at d = λ/2. A 2×2 array is capped at 95 mm sides and a 0.53 λ maximum baseline.
* A **triangular** lattice of side d has shortest dual vector 2λ/(d√3), so it aliases only at d = 1.155 × (λ/2) — **15.5% more spacing**. And 3-on-a-ring-plus-centre generates exactly a triangular lattice.

This is the familiar hexagonal-sampling advantage from 2-D signal processing, applied to the aperture instead of to an image. Measured consequences at a fixed 0.5 λ radome and −3 dB grating margin:

| | 2×2 square | Y (3 + centre) | gain |
|---|---|---|---|
| max baseline | 0.53 λ | 0.82 λ | **+55%** |
| array-factor null width | 1.49 | 0.95 | **36% narrower** |
| null depth | −21.3 dB | −22.2 dB | +0.9 dB |
| P(satellite worse than 1 antenna), rank 1 | 12.8% | **6.9%** | **1.9× better** |
| P(satellite worse than 1 antenna), rank 2 | 40.4% | **31.6%** | 1.3× better |

**Halving the outage probability for free is a better return than anything else in this document.**

The honest counter-argument: the centre element has three near neighbours while each outer element has one, so their embedded element patterns and mutual coupling differ. This matters for *manifold calibration* (Task 6) but **not** for the projection algorithm, which never uses the manifold — the algorithm sees whatever linear mixing the array imposes and nulls the resulting vector. If you later add MVDR or DOA estimation (Task 5), budget for a per-element calibration table rather than assuming element symmetry.

Two further geometry notes:

* Verify the ambiguity numerically for *your* final layout, including any deliberate asymmetry, using `asp_array_ambiguity.m`. Do not trust the rule of thumb.
* Element *rotation*: use sequential rotation of the patch feeds (0°/120°/240° for the Y) to equalise axial ratio across the array. This is standard GNSS antenna practice and costs nothing in layout.

---

## 5. Detection capability — the benefit nobody quotes

Neither the paper nor any of your three scripts contains a spoofing *detector*; the pipeline nulls unconditionally. This is the most serious defect in the design as it stands (see Task 2), and fixing it makes element count matter in a new way.

The eigenvalue test statistic λ₁/mean(λ₂…λ_N) distinguishes spoofing from authentic signals not by power but by **spatial rank**:

* spoofer: Σₖ pₖ **bb**ᴴ = (N_spoof · p_s)**bb**ᴴ → rank 1, all power in one eigenvalue
* authentic: Σₘ pₘ **aₘaₘ**ᴴ ≈ (N_auth · p_a)**I** → approximately isotropic, raises the floor without creating spread

This is a much stronger discriminant than "the spoofer is louder", which is what the paper leans on — and it is why detection works even at 0 dB SAPR, where the two have equal total power.

Measured, P_FA = 10⁻³, H₀ = a clean sky with 9 satellites at 45 dB-Hz:

**Minimum SAPR for 90% detection**

| geometry | 1 ms | 5 ms | 20 ms |
|---|---|---|---|
| tri3 (N=3) | +1.7 dB | −0.2 dB | −0.3 dB |
| y4 (N=4) | +0.3 dB | −0.7 dB | −1.4 dB |
| circ7 (N=7) | −0.5 dB | −1.6 dB | −2.1 dB |

**Detection probability at 0 dB SAPR, 1 ms** — the paper's own "spoofing is a threat" threshold:

| N = 3 | N = 4 | N = 7 |
|---|---|---|
| 0.447 | **0.882** | 0.998 |

Four elements **doubles the detection probability at exactly the power level where the threat begins**. Expressed as sensitivity, 3 → 4 buys 1.4 dB, which is the array gain N against a coherent source appearing directly in the test statistic.

Two things to notice in that table. First, sensitivity **saturates** near −1 to −2 dB SAPR regardless of dwell: the authentic constellation's own spatial structure sets an irreducible floor on the H₀ statistic. You cannot detect a spoofer weaker than the constellation from the raw spatial covariance, no matter how long you integrate. Second, the H₀ statistic sits measurably **above** the white-noise asymptotic (1 + √(N/K))² at long dwells — 1.021 vs 1.007 for the Y-array at 20 ms. **A product that sets its threshold from the white-noise formula will false-alarm on a clean sky.** Calibrate the threshold against a real constellation.

---

## 6. Cost side of the ledger

### 6.1 Computational complexity — it is linear in N, and that is the point

The streaming (per-sample) cost of the whole method is O(N), not O(N³). Real multipliers per sample, using the 3-multiply complex product that maps onto the DSP48 pre-adder:

| block | N = 3 | N = 4 | ratio |
|---|---|---|---|
| paper's γ + β correlators (2N−1 complex) | 20 | 28 | 1.40 |
| **covariance, N(N+1)/2 entries** | **18** | **32** | 1.78 |
| beamformer, 1 beam | 9 | 12 | 1.33 |
| beamformer, 12 PRN beams | 108 | 144 | 1.33 |

Note that at N = 3 the full covariance is *cheaper* than the paper's γ/β pair, and at N = 4 it costs 14% more than γ/β while eliminating the epoch delay memory entirely.

The O(N³) work — eigen-decomposition, weight solve — runs **once per millisecond**, not per sample. For N = 4, one cyclic Jacobi sweep is 6 rotations and 6 sweeps is 36; at ~75 cycles per pipelined CORDIC rotation that is 2700 cycles ≈ 13.5 µs at 200 MHz, against a 1 ms dwell. **1.4% duty cycle.** Going from N = 3 (3 rotations/sweep) to N = 4 (6 rotations/sweep) doubles a number that is already negligible.

### 6.2 FPGA resources

At fs = 16.368 MHz against a 200 MHz fabric clock the time-multiplexing ratio is 12.2, so one DSP48 serves ~12 multiplies.

| block | N = 3 | N = 4 |
|---|---|---|
| covariance engine | 2 DSP | 3 DSP |
| single output beam | 1 DSP | 1 DSP |
| 12 per-PRN beams | 9 DSP | 12 DSP |
| Jacobi EVD (shared, 1 kHz) | ~4 DSP | ~4 DSP |
| **array processor total** | **~16 DSP** | **~20 DSP** |

Four DSP48s. On a Zynq-7020 (220 DSP48E1) that is **1.8% of the device**. Anyone arguing against the fourth antenna on FPGA-resource grounds has not costed it.

### 6.3 Memory — where the real saving is, and it is not about N

| | N = 3 | N = 4 |
|---|---|---|
| Paper's β delay line (K·N complex, 16-bit I/Q, K = 16368) | 1.57 Mbit ≈ 44 BRAM36 | 2.10 Mbit ≈ 58 BRAM36 |
| **Covariance accumulators (N(N+1)/2 × 48 bit)** | **288 bit (registers)** | **480 bit (registers)** |

The fourth antenna costs 14 BRAM36 *if you keep the paper's β statistic*. Replacing β with the covariance diagonal — which is strictly more accurate, see Task 2 — removes **all** of it. The memory question is not "3 or 4 antennas", it is "β or covariance", and the answer to that is independent of N.

### 6.4 Latency

| stage | latency |
|---|---|
| decimation FIR (64 taps @ 16.368 MHz) | 2.0 µs |
| channel equaliser (16 taps) | 0.5 µs |
| beamformer pipeline | 0.25 µs |
| **sample-path total** | **~2.8 µs** |
| covariance dwell + EVD + weight solve | 1 ms + 15 µs |

The sample path is unaffected by N to within a fraction of a clock cycle. Weight *age* is ~1.015 ms in both cases. On a rotating platform, weight age caps the null at 20 log₁₀(2π(D/λ)Δθ/√3): with D = 0.82 λ that is −46 dB at 100°/s, −34 dB at 400°/s, −26 dB at 1000°/s. All below the estimator's own −24 dB accuracy up to several hundred °/s, so **dwell length is set by estimator statistics, not by platform dynamics**, until about 1000°/s.

### 6.5 Numerical stability and conditioning

Better at N = 4, not worse — the opposite of the usual intuition that bigger matrices are harder.

* The 3-element covariance at 5.5 dB SAPR has eigenvalue spread λ₁/λ₃ ≈ 1.19; the 4-element has ≈ 1.25. **Both are appallingly conditioned in the sense that matters** — the useful signal is a 20% perturbation of the identity — and the extra array gain N makes 4 slightly better, not worse.
* The Jacobi EVD is a sequence of unitary similarity transforms, so ‖A‖_F is invariant *exactly*. There is no dynamic-range growth, and the word length is set by the eigenvalue spread you must resolve, not by N. Verified: 18-bit rotations track floating point to sin θ < 5×10⁻⁵ (−86 dB), which is 60 dB below the statistical floor.
* Nothing in the recommended datapath forms R⁻¹, a Cholesky factor, or a matrix inverse of any kind. Weights come from the eigen-decomposition. Condition number never enters.

### 6.6 Calibration complexity

**Identical for 3 and 4 elements**, because the projection algorithm never uses the array manifold. It nulls whatever vector the spoofer produces at the array output, including the effects of gain mismatch, phase mismatch, cable skew *and mutual coupling* — any fixed invertible linear mixing **M** leaves "the direction the spoofer arrives from" a well-defined vector **MCb̄** in the measured space. `asp_scenario.m` models coupling explicitly so this can be tested rather than asserted.

What *is* affected by element count, and only weakly:

* **Frequency-dependent** mismatch is not absorbed. Inter-channel group-delay skew Δt limits the null to roughly 20 log₁₀(2π(B/√12)Δt): 100 ps of skew → −55 dB, 1 ns → −35 dB. This scales with the *worst* channel pair, so a fourth channel adds one more opportunity to be the worst one. Budget a per-channel FIR equaliser regardless (Task 3).
* If you later adopt MVDR/LCMV or any DOA method (Task 5), you need a calibrated manifold, and calibration effort scales as N. That is an argument for deciding the element count *now* rather than after the chamber campaign.

---

## 7. Challenges to your framing

**Your spoofer is in the wrong place.** The paper puts it at 45° elevation and your scripts inherit that. Terrestrial spoofers sit at low elevation, which puts them angularly close to exactly the low-elevation satellites your projection null will damage, and which is also where nearby structures produce azimuthally-offset reflections. Both effects make three elements look worse than the paper's simulation suggests. My default scenario uses 15°.

**A 2×2 square is not the natural 4-element layout.** See §4 — the Y-array halves your outage probability for free.

**Do not build the product around the "no calibration required" property.** It is a genuine and elegant strength of the algorithm, and it is what makes two independently-locked AD9361s viable (Task 4). But it is a property of *one* algorithm. The moment you want MVDR, LCMV, DOA-based spoof discrimination, or a beam pointed at a satellite you have almanac knowledge of, you need a calibrated manifold. Calibration is a one-off chamber campaign per design plus a cheap per-unit trim, and it unlocks a lot. Plan for it. Treat calibration-free operation as the *graceful degradation mode*, not the design centre.

**Four is a waypoint, not the destination.** The AD9361's 2×2 structure makes the natural SKU ladder 2 / 4 / 8 elements, not 4 / 7. Eight elements (four AD9361s, or two ADRV9026s) gives 6 spare adaptive DOF, 8.3 dB array gain, 2.7% outage, and −24.5 dB null. If a high-end variant is on the roadmap, design the digital architecture parametric in N now — the reference model in `matlab/` already is — so the step is a re-parameterisation rather than a rewrite.

---

## 8. Summary scorecard

| Criterion | N = 3 (triangle) | N = 4 (Y) | Verdict |
|---|---|---|---|
| Spare adaptive DOF, 1 null | 1 | 2 | **4 wins, 2×** |
| Spare adaptive DOF, 2 nulls | 0 | 1 | **4 wins, decisive** |
| Array gain, 1 null | 3.0 dB | 4.7 dB | 4 wins, +1.7 dB |
| Array gain, 2 nulls | −0.04 dB | +2.5 dB | **4 wins, +2.6 dB** |
| Satellites left worse than 1 antenna, 1 null | 18.6% | 7.0% | **4 wins, 2.7×** |
| Satellites left worse than 1 antenna, 2 nulls | 59.9% | 22.3% | **4 wins, 2.7×** |
| Null depth | −20.6 dB | −22.3 dB | 4 wins, +1.7 dB |
| Null width (direction cosine) | 1.46 | 0.95 | **4 wins, 36% narrower** |
| Min detectable SAPR, 1 ms | +1.7 dB | +0.3 dB | 4 wins, 1.4 dB |
| P(detect) at 0 dB SAPR, 1 ms | 0.447 | 0.882 | **4 wins, 2×** |
| Streaming multipliers | 18 | 32 | 3 wins, +14 mult |
| FPGA DSP48 | ~16 | ~20 | 3 wins, 4 DSP (1.8% of a Zynq-7020) |
| BRAM | equal (both zero with the covariance estimator) | equal | tie |
| Sample-path latency | equal | equal | tie |
| Numerical conditioning | slightly worse | slightly better | 4 wins marginally |
| Calibration effort | equal | equal | tie |
| **RF BOM (AD9361 count)** | **2 chips, 1 RX wasted** | **2 chips, 0 wasted** | **tie in cost, 4 wins in value** |
| Antenna/radome cost | 3 elements | 4 elements | 3 wins, one patch + one LNA |

The entire cost of the fourth element is **one patch antenna, one LNA, one SAW filter, four DSP48 slices, and 14 more streaming multiplies**. Call it $12 of BOM and 2% of a small FPGA. Against that, it roughly halves the number of satellites you destroy, doubles your detection probability at the threat threshold, and is the difference between coping and not coping with a second emitter.

**Move to four. Use the Y geometry, not the square. Design the digital chain parametric in N so eight is a configuration change.**
