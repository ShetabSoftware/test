# Executive Summary

## The short version

**Go to four antennas.** The RF hardware cost is zero — the AD9361 is a 2×2 transceiver, so three antennas needs two chips with one receive path wasted and four needs the same two chips with nothing wasted. The digital cost is four DSP48 slices, about 2% of a small Zynq. Against that, four elements roughly halve the number of satellites the null destroys, double the detection probability at the threat threshold, and are the difference between coping and not coping with a second emitter. **Use a Y layout (three on a ring plus one centre), not a 2×2 square** — a triangular baseline lattice tolerates more spacing before aliasing, which buys 55% more aperture inside the same radome and halves the satellite outage, for free.

**The design as it stands has one blocking defect: there is no spoofing detector.** The pipeline nulls unconditionally, so in a clean environment it steers a null into the strongest authentic satellite. Since an anti-spoofing device spends nearly all its operating hours *not* under attack, that failure mode dominates its expected value. Your own File 1 identifies this in a comment; it needs implementing. It is cheap — the eigenvalues come out of the same decomposition that produces the weights.

**The estimator can be improved by 9.5 dB for three extra multipliers.** The paper's γ/β construction is exactly *one power-iteration step* applied to the whitened spatial covariance, started from the reference element. Running the iteration to convergence is the whole improvement, and it also removes ~58 BRAM36 of delay memory, produces the detector for free, enables rank-2 nulling, and gives R⁻¹ for MVDR at no extra cost.

**Keep the AD9361 for the first product, but share the LO.** Independent PLLs cap null depth around −35 to −41 dB through differential phase noise alone; a shared LO makes that term exactly zero. The AD9361's real limitation is its 12-bit converters, which cap usable J/N at ~57 dB — fine for spoofing, marginal for jamming. Migrate to a discrete front end with a quad 16-bit ADC (81 dB J/N) for the volume product, once the digital architecture is frozen.

---

## Findings by severity

### Blocking

| | Finding | Fix |
|---|---|---|
| **B1** | **No spoofing detector.** Nulls unconditionally; damages the receiver whenever there is no threat. | Eigenvalue statistic λ₁/mean(λ₂…λ_N) + MDL rank, gating the weight update. Measured minimum detectable SAPR +0.3 dB at 1 ms, P_FA = 10⁻³. |
| **B2** | **Operating point wrong by 87 dB** in `GPS_SPOOFER_V9` (signals 60 dB above the noise; truth is 27 dB below). Every null-depth number from that script is measured in a regime that does not exist. | Parameterise power as C/N₀ in dB-Hz and derive sample SNR from fs. |
| **B3** | **Non-causal.** Weights estimated from the whole record, then applied to it. | Weights from dwell *k* applied to dwell *k+1*. |
| **B4** | **`error()` on the algorithmic path**, and the degenerate case is reachable. | Graceful degradation to the quiescent beam, then to a single antenna. Never crash. |
| **B5** | **The fixed-point study does not model the accumulator** — the one place fixed-point covariance estimation goes wrong. Its own comment concedes this. | Bit-exact DSP48 model; see the trap below. |
| **B6** | **Word lengths asserted, not derived**, and the study never asks what converter width actually determines (J/N headroom). | All word lengths derived, with the derivations, in `config/asp_fx_plan.m`. |

### The fixed-point trap, because it is invisible in floating point

Rounding sample products *before* the covariance accumulator injects an error into R that is:

* **deterministic** — a fixed matrix, not noise;
* **rank one, aligned with the boresight steering vector** (measured alignment 1.000 with the all-ones vector; a random direction gives 0.5). Boresight is where the satellites are. **A truncating implementation synthesises a phantom source at zenith and steers the null into the sky.**
* **independent of dwell length** — measured identical to three significant figures across a 16× change in K.

The last property is what makes it dangerous: everything else in this estimator improves as 1/√K, so a designer who validates at 1 ms reasonably expects improvement at 100 ms. The floor is flat.

| product fractional bits retained | 10 | 14 | 18 | 22 | 26 | full |
|---|---|---|---|---|---|---|
| null-depth floor | −25.4 | −48.9 | −73.0 | −97.1 | −122.0 | −156 dB |

**The fix is free.** A DSP48 holds the exact 32-bit product of two 16-bit words and accumulates in its 48-bit P register with no rounding. Measured: 42 bits used, 6 spare.

---

## Task 1 — three antennas or four

All geometries sized to a common 0.5 λ radome and a common −3 dB grating margin (1200–1500 scenes each).

| | N = 3 (triangle) | N = 4 (Y) | verdict |
|---|---|---|---|
| Spare adaptive DOF, 1 null | 1 | 2 | 4 wins, 2× |
| Spare adaptive DOF, 2 nulls | **0** | 1 | **4 wins, decisive** |
| Array gain, 1 null | 2.99 dB | 4.68 dB | +1.69 dB |
| Array gain, 2 nulls | **−0.04 dB** | **+2.54 dB** | **+2.58 dB** |
| Satellites worse than one antenna, 1 null | 18.6% | **7.0%** | **2.7×** |
| Satellites worse than one antenna, 2 nulls | 59.9% | **22.3%** | **2.7×** |
| Null depth | −20.6 dB | −22.3 dB | +1.7 dB |
| Null width (direction cosine) | 1.46 | 0.95 | **36% narrower** |
| P(detect) at 0 dB SAPR, 1 ms | 0.447 | **0.882** | **2×** |
| Streaming multipliers | 18 | 32 | 3 wins, +14 |
| FPGA DSP48 | ~16 | ~20 | 3 wins, 4 slices |
| **AD9361 count** | **2, one RX wasted** | **2, none wasted** | **tie in cost** |

The rank-2 row is the one that decides it. **With two nulls to place, a 3-element array delivers exactly single-antenna array gain** (−0.04 dB measured against 0.00 dB predicted) and leaves 60% of the constellation worse off than a single passive antenna. Every element beyond the first has been consumed by the constraints, and the array has stopped being an array.

**Geometry:** use three elements on a ring at 120° plus one at the centre. Ambiguity is a property of the baseline *lattice*, not of the longest baseline — the common "keep every baseline under λ/2" rule is inherited from linear arrays and costs real aperture on a planar one. A triangular lattice aliases at 1.155 × (λ/2) where a square lattice aliases at λ/2, so inside a 0.5 λ radome the Y reaches a 0.82 λ maximum baseline against the square's 0.53 λ. Measured consequence: outage falls from 12.8% to 7.0%.

---

## Task 2 — measured before and after

Full waveform model: real Gold codes, navigation data, code Doppler, channel mismatch, 4-element Y array, SAPR 5.5 dB, C/N₀ 45 dB-Hz.

| Metric | paper γ/β | redesign (whitened EVD) |
|---|---|---|
| Null depth @ 1 ms | not available (needs 2 epochs) | **−23.7 dB** |
| Null depth @ 2 ms | −15.8 dB | **−25.3 dB** (+9.5) |
| Null depth @ 20 ms | −22.2 dB | −26.8 dB (+4.6) |
| Epoch delay memory | ~58 BRAM36 | **0** |
| Streaming multiplies (N=4) | 28 | 32 (+14%) |
| Spoofing detector | none | **free** |
| Rank-2 capable | no | yes |
| Dividers on the sample path | 3 | **0** |
| Causal | no | yes |

**Where it saturates, and what that implies.** The pre-correlation estimator improves only 3.1 dB from a 1 ms to a 20 ms dwell, where a noise-limited estimator would have gained 13 dB. The shortfall is a *bias* from the authentic constellation's own contribution to the covariance, which no amount of integration removes. If you want a deeper null, stop integrating and change domain: post-correlation refinement obeys

$$\text{stage-2 null depth} = -(\mathrm{SAPR} + 27.0)\ \text{dB}$$

to better than 0.1 dB across a 20 dB sweep, where 27.0 is the C/A cross-correlation bound (23.9 dB) plus ~3.1 dB. Honest accounting: at 5.5 dB SAPR that is worth about 6 dB over stage 1's saturation point, not 20 dB. **Its real value is per-PRN attribution** — telling the receiver *which* measurements are counterfeit so it can exclude them and bound the remaining error. That is what turns a gadget into a certifiable component, and it is absent from the paper.

---

## Tasks 3–6 in one paragraph each

**FPGA.** The whole array processor is ~20 DSP48, zero BRAM, and 17.7 µs of once-per-millisecond compute — 2% of a Zynq-7020. Size the device for the *receiver* correlator bank, not for the nulling algorithm. Choose the fabric clock as an exact integer multiple of the sample rate (163.68 MHz = 10 × 16.368) so every time-shared block has a static schedule with no elastic buffering. Use cyclic Jacobi for the eigen-decomposition: it is norm-preserving, so there is no dynamic-range growth, no conditioning concern, a fixed latency, and a direct CORDIC mapping. Keep the MDL rank decision in software so the rank policy stays field-tunable.

**RF front end.** Rank the requirements by *matching*, not by performance — a mediocre receiver with excellent channel matching will out-null an excellent receiver with poor matching, every time. Non-negotiables: one LO, one clock with boot-verified sample alignment, common AGC, a per-gain-index phase table, frozen tracking calibrations, and a digital per-channel equaliser. The AD9361's calibration-free compatibility with this algorithm is a genuine architectural synergy — its inter-chip phase offset changes at every re-sync, and the projection method does not care. Migrate to a discrete front end with a quad 16-bit ADC for volume; evaluate the ADRV9026 (4 RX, one LO, 16-bit, one package) before committing to two AD9361s.

**Algorithms.** Spatial processing is the only technique that works against an arbitrarily sophisticated adversary with a single transmitter, because it exploits physics the attacker cannot control. Everything else is corroboration. Build five layers: AGC monitoring, pre-correlation spatial detection and nulling, post-correlation attribution, measurement-domain consistency plus Galileo OSNMA where available, and RAIM on the *surviving* measurements. Do not market RAIM as anti-spoofing — it assumes independent faults and a spoofer generates a self-consistent constellation, so it sees no residual at all. Skip ICA and end-to-end machine learning; the physics is known and the closed-form estimators are better and certifiable.

**Product.** Decide early whether you are selling an in-line module or an integrated receiver, because the architectural commitments differ and the in-line version cannot have per-satellite power maximisation — which is where *all* of the array gain in this method comes from. Ship the calibration inside a serialised antenna assembly. Handle carrier-phase continuity across weight updates, or you have no RTK and no precise timing product. Test the no-threat case as rigorously as the threat case.

---

## Corrections the measurements forced

Recorded because the negative results are as useful as the positive ones, and because two of them removed constraints I had assumed were real.

* **A specular ground bounce costs nothing.** I expected it to force rank-2 nulling. A planar array cannot distinguish elevation +θ from −θ — measured coherence 1.000000, exactly — so the reflection shares the direct path's steering vector and is co-nulled for free, at every frequency, with no wideband penalty. The rank-2 cases that *do* consume a degree of freedom are arrivals at a different **azimuth**: a building reflection, a second spoofer, a jammer. This also means a planar CRPA gives you **no** spatial multipath rejection for the authentic signals, which belongs in your datasheet.
* **Rank 2 is not automatically better than rank 1.** A second arrival 6 dB down lifts the second eigenvalue only ~0.06 above the noise floor in a 5 ms dwell, giving an eigenvector with ~26° of error — nulling it is measurably *worse* than leaving it alone. The rank decision must be driven by resolvability (MDL on the eigenvalues), never by prior knowledge that a second source exists. This is why the rank decision belongs in software, where its hysteresis policy stays tunable.
* **Navigation data bits do not damage the paper's β statistic.** They invert part of the lagged sum, but only through a factor common to every element, which cancels in the scale-invariant projector. Measured ρ = 0.99668 with bits against 0.99702 without. Bit synchronisation is not a prerequisite for the array processing.
* **The paper's magnitude/phase split is essential, not stylistic.** Taking the raw covariance column directly gives ρ = 0.74 and a −6.5 dB null, because R₁₁ is dominated by noise power while R_i1 contains only the spoofer term. Preserving that structure in your code was right.

---

## Priority order

1. **Eigenvalue detection and gating.** Smallest change, largest impact — without it the product has negative value most of the time.
2. **Whitened covariance eigenvector** replacing γ/β. 9.5 dB at 2 ms, removes 58 BRAM, enables everything below.
3. **Exact covariance accumulation.** Free, and prevents a failure that floating-point simulation cannot show you.
4. **Four elements, Y geometry.** Free in RF, ~$12 of BOM, 2% of the FPGA.
5. **Shared LO and boot-time alignment self-test.** Removes the two silent hardware failure modes.
6. **Post-correlation attribution.** The commercially differentiating feature.
7. **Array calibration and LCMV.** Unlocks the optimal beamformer.
