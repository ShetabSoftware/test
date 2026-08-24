# Task 6 — Product-Level Design Recommendations

Written as the chief architect would: what the paper does not discuss, what industry does as a matter of course, and what will hurt you if you skip it.

---

## 1. Decide what you are selling — it changes almost everything

Two products are hiding inside "a GPS anti-spoofing device", and they have incompatible architectures. Choosing late is expensive.

### SKU-A — in-line module: 4 antennas in, one RF/IF stream out

Sits between the antennas and an *unmodified* customer receiver. Large addressable market (retrofit), no receiver development, easy to sell.

**But it is fundamentally limited, and you should understand the limits before committing.** With one output you can only apply *one* weight vector, so:

* **No per-satellite power maximisation.** The paper's eq. (23) needs one output per PRN. With a single output you are stuck with eq. (14) and a fixed **h**, whose measured mean authentic gain is ≈ 0 dB — no better than a single antenna (measured: +0.79 dB for tri3, −0.61 dB for y4). **All of the array gain in this method comes from the power-maximisation unit,** and SKU-A cannot have it.
* **No attribution.** You cannot tell the downstream receiver which PRNs are counterfeit, so you cannot support an integrity claim.
* **A fixed h nulls a curve on the sky.** With one output and a rank-1 projector at N = 4, the residual is a 3-D subspace and a fixed **h** has a null surface in it. Some satellite direction is always badly attenuated.
* **Carrier phase discontinuity.** See §4. This is the killer for any RTK or timing customer.

### SKU-B — integrated receiver: 4 antennas in, PVT out

Post-correlation combining, per-satellite optimal weights, per-PRN attribution, integrity output. Everything in Task 5's L3 layer becomes available. Higher development cost, much higher value, and the only version that can make a defensible integrity claim.

**Recommendation: build SKU-B, and derive SKU-A from it** as a reduced mode of the same hardware and FPGA image. Do not build SKU-A first and try to grow it — the architectural commitments (single output, no correlators, no receiver state) are exactly the ones you must undo.

If you must ship SKU-A first for revenue, at least design the FPGA with the correlator bank present and the multi-beam datapath in place, even if the first firmware does not use them.

---

## 2. System architecture

```
  ANTENNA ASSEMBLY (sealed, calibrated as a unit, serialised)
    4x patch, Y layout, sequential rotation
    4x LNA + SAW at the element
    temperature sensor
    calibration EEPROM (per-unit equaliser taps, gain-index phase table)
         │  4x phase-matched coax, same reel, same length
  ───────┼──────────────────────────────────────────────────────────
  MAIN UNIT
    RF: 2x AD9361, SHARED external LO, common gain index
    FPGA: DFE -> covariance -> EVD -> detect -> weights -> beamform
                                                        -> correlators
    PS:  MDL rank, mode logic, stage-2 weights, PVT, RAIM, integrity
    DDR: rolling raw-sample buffer for forensic capture
    IF:  Ethernet/CAN/serial, PPS in/out, NMEA + proprietary integrity msg
```

**The antenna assembly must be a calibrated, serialised unit with its own EEPROM.** This is the single most important productisation decision on this list. The calibration that matters — inter-channel amplitude and group-delay response — lives in the antenna, LNAs, filters and cables, not in the main unit. If the assembly is separable and field-swappable without its calibration data, every field swap silently degrades null depth and nobody finds out. Ship the calibration *with* the antenna, in the antenna.

---

## 3. Array, antenna and RF

Covered in Tasks 1 and 4. The decisions that matter:

* **Y geometry (3 on a ring at 120° + 1 centre), ring radius 0.45–0.50 λ.** Triangular baseline lattice tolerates 15.5% more spacing than a square before aliasing, giving a 0.82 λ maximum baseline inside a 0.5 λ radome — 55% more aperture than a 2×2 square, a 36% narrower null, and half the satellite outage. Task 1 §4.
* **Sequential feed rotation** (0°/120°/240°) to equalise axial ratio across the array.
* **Continuous ground plane, ≥ 1.5 λ if mechanics allow.** Sets the low-elevation pattern and front-to-back ratio — which matters because your threat is at low elevation.
* **The centre element has three near neighbours and the outer elements have one**, so their embedded patterns differ. Irrelevant to the projection algorithm (it never uses the manifold); relevant the moment you add MVDR or DOA. Budget a per-element manifold table.
* **Shared LO. Not negotiable.** Independent PLLs cap null depth at −35 to −41 dB through differential phase noise alone; a shared LO makes that term exactly zero.
* **Phase-matched cables**: same reel, same length, same routing, same bend radius, matched to < 100 ps. 1 ns of differential group delay caps the null at −35 dB.
* **A planar array cannot distinguish elevation +θ from −θ.** This is free ground-bounce rejection for the spoofer (Task 1 §2.1) and it means **no** spatial multipath rejection for the authentic signals. Put that in the datasheet honestly; a customer expecting a CRPA to fix their multipath will be disappointed and will blame you.

---

## 4. The problem nobody in the paper mentions: carrier-phase continuity

**This will kill an RTK or timing product and it is not obvious.**

Every weight update changes the output carrier phase of satellite *m* by

$$\Delta\phi_m = \arg(\mathbf{f}_\text{new}^H\mathbf{a}_m) - \arg(\mathbf{f}_\text{old}^H\mathbf{a}_m)$$

which is **different for every satellite** because it depends on **a**ₘ. At a 1 kHz update rate with an estimator whose eigenvector has a random phase perturbation each dwell, this injects per-satellite phase noise directly into the carrier tracking loops. Consequences, in increasing order of severity:

1. Degraded carrier-phase measurement noise → worse velocity.
2. Cycle slips → RTK ambiguity resolution fails.
3. Loss of carrier lock at low C/N₀.

Three fixes, in increasing order of goodness:

**(a) Phase-continuous weight update.** Constrain each new weight vector so that **f**_new **f**_old is real and positive — a one-line rotation that removes the *common* phase jump. It does **not** remove the per-satellite differential, but it removes the largest term and costs nothing.

**(b) Slow, smoothed updates with hysteresis.** Update at 10–50 Hz rather than 1 kHz once converged, and apply a first-order smoothing filter to the weights. Trades tracking of platform dynamics for phase stability. Reasonable for static and low-dynamics installations.

**(c) Post-correlation combining (SKU-B).** Combine *after* despreading, per satellite. The tracking loop then sees a signal whose phase you control explicitly, and the weight change becomes a known, correctable rotation rather than an unknown disturbance. **This is the correct answer** and it is another reason SKU-B is the right product.

For SKU-A, (a) + (b) is the best available, and you should specify the resulting carrier-phase performance honestly rather than let a customer discover it.

---

## 5. Synchronisation and clock distribution

| item | requirement | why |
|---|---|---|
| Reference | one TCXO ±0.5 ppm; OCXO if timing is a feature | common to both transceivers |
| LO | one synthesiser, symmetric on-PCB Wilkinson split | common-mode phase noise |
| Sample clock | derived from the same reference; **verified alignment every boot** | 1 sample = 61 ns = fatal |
| MCS | run at every init; verify, do not assume | AD9361 divider start phase varies per lock |
| Boot self-test | inject a common pilot through a splitter, cross-correlate all pairs, log lag and phase, auto-resync on failure | catches the failure mode that otherwise looks like "poor performance" |
| PPS | discipline from the *authenticated* solution only | a spoofed PPS output is worse than none |

**The boot-time alignment self-test is the highest-value test in the product.** Inter-chip misalignment is silent: the device runs, produces a position, and simply nulls badly. Without an explicit test you will ship units that appear to work and do not. Make it automatic, make it logged, and make the log field-readable.

---

## 6. Calibration

Three layers, three different lifetimes.

**Layer 1 — design calibration (once per design, anechoic chamber).**
Measure the full array manifold including mutual coupling and radome, over azimuth, elevation and frequency. Produces the nominal manifold table. Needed only if you adopt MVDR/LCMV/DOA, but do it anyway: **you cannot write a credible specification without knowing your own embedded element patterns.**

**Layer 2 — production calibration (per unit, minutes, conducted).**
Inject a common wideband pilot through a matched splitter into all four inputs. Measure H_i(f) per channel and solve for the 16-tap complex equaliser that flattens all channels onto the reference. Also build the per-gain-index phase-correction table. Store in the antenna assembly's EEPROM. **This is what converts "matched to 1 ns" into "matched to 50 ps", and it is the difference between a −35 dB null and a −55 dB one.**

**Layer 3 — runtime tracking (continuous, free).**
The algorithm re-estimates the spatial signature every millisecond, so slow thermal drift is tracked automatically. This is the paper's genuine gift and it is why the product can survive with modest production calibration.

**Do not build the product around the calibration-free property.** It is elegant and it is what makes two independently-locked AD9361s viable — but it is a property of *one* algorithm. MVDR, LCMV, DOA attribution and manifold-based spoof discrimination all need calibration, and they are all on the roadmap. Treat calibration-free operation as the **graceful degradation mode**, not the design centre.

---

## 7. PCB, thermal and mechanical

**PCB**
* Separate RF, digital and power sections with continuous ground; stitch vias at ≤ λ/20 at the highest digital harmonic of concern.
* **Route the four receive chains as symmetrically as the layout allows** — same layer, same via count, same reference plane, same length. Symmetry is worth more than absolute length matching because it makes the channels *drift together*.
* Keep the LO distribution symmetric and short; treat the splitter as a controlled-impedance structure.
* Guard the FPGA's switching harmonics away from L1: 1575.42 MHz is a harmonic of many convenient clock rates. **Choose every clock in the system so that no harmonic lands within ±20 MHz of L1** — including the DDR clock, the Ethernet PHY, and switching-regulator frequencies. This is a schematic-review checklist item and it is very expensive to fix in a spin.
* Use LDOs, not switchers, for the LNA and mixer rails, or a switcher with a spread-spectrum mode *disabled* (spread spectrum smears a spur into a noise floor across the GNSS band, which is worse).

**Thermal**
* Drift within a 1 ms dwell is nil, and slow drift is tracked. What matters is **differential** drift: if one channel sits next to the FPGA and another does not, they diverge over minutes and stale the production calibration table.
* Mitigation: symmetric mechanical layout, a temperature sensor per channel (or at least one per RF section), and a calibration table indexed by temperature.
* The AD9361s are the main heat source in the RF section (~1 W each). Place them symmetrically with respect to the four channels, or thermally isolate the RF chains from them.

**Mechanical**
* Antenna assembly sealed, potted where practical, with a defined radome dielectric — the radome is part of the manifold.
* Connector strain relief on the four coax runs; a strained connector is a phase drift.
* Serialise the antenna assembly and bind its calibration data to the serial number.

---

## 8. FPGA and converter selection

| | recommendation | reason |
|---|---|---|
| Prototype | Zynq-7020 (ZC702 / Zedboard / Antsdr) + 2× AD9361 | array processor is 13% of DSPs; ADI stack is mature |
| Product, 12-PRN single frequency | Zynq-7020/7030 | correlator bank dominates, not the array processor |
| Product, multi-constellation | Zynq-7045 or ZU3EG | ZU3EG for the better DSP/W |
| Converter, phase 1 | AD9361, 12-bit | 57 dB usable J/N — adequate for spoofing |
| Converter, phase 2 | discrete + AD9653 quad 16-bit | 81 dB usable J/N; one quad ADC removes the alignment failure mode entirely |
| Alternative single-chip | ADRV9026 (4RX, 16-bit, one LO) | collapses every inter-chip problem into one package; evaluate it before committing to 2× AD9361 |

**Size the FPGA for the receiver, not for the anti-spoofing function.** The array processor is ~20–28 DSP48. The correlator bank is 96–400. Anyone sizing the device from the nulling algorithm will under-specify by an order of magnitude.

---

## 9. Manufacturing and test

**Per-unit production test sequence**

1. Power, rails, current signature.
2. AD9361 SPI/register readback; MCS; **inter-channel sample alignment** (pilot injection).
3. Per-channel gain/phase/group-delay sweep → compute and store equaliser taps.
4. Per-gain-index phase table.
5. Noise-figure spot check per channel.
6. **Conducted spoofing test**: split simulator output into all four inputs at a known SAPR, verify detection and null depth against a limit line.
7. **Conducted no-spoof test**: verify the detector does *not* fire and the nulling is bypassed. This is the test that catches B1 regressions, and it is the one most likely to be omitted.
8. Live-sky sanity: acquire, track, fix.
9. Write calibration data + serial to EEPROM; read back and verify.

Target 5–10 minutes per unit. Steps 3, 6 and 7 are the ones with real yield impact.

**Type test (per design)**
* Anechoic chamber with a rotating positioner: manifold, null depth vs spoofer azimuth/elevation, satellite outage statistics.
* Temperature cycling with continuous null-depth monitoring — the single most informative environmental test, because it directly exercises differential drift.
* Vibration (phase stability under vibration is a real failure mode for cabled arrays).
* EMC, including self-jamming from the unit's own digital sections.

---

## 10. Verification, benchmarking and specification

**Verification pyramid**
1. Unit tests on the reference model (`matlab/verify/`) — run in CI.
2. RTL vs reference model, bit-exact, using `export/asp_export_vectors.m`.
3. Hardware-in-the-loop with recorded IF.
4. Conducted spoofing per the paper's Figure 8 setup (rooftop array + split simulator + conductive combining). That test setup is correct and worth reproducing exactly.
5. Chamber.
6. Field.

**The specification numbers to publish**, with the definitions that make them meaningful:

| parameter | definition | expected |
|---|---|---|
| Null depth | gain toward the spoofer, dB **re one antenna element** | −22 dB @ 1 ms, −27 dB @ 20 ms (stage 1) |
| | | −(SAPR + 27) dB (stage 2) |
| Min detectable SAPR | at P_D = 0.9, P_FA = 10⁻³, against a **live constellation** | +0.3 dB @ 1 ms; −1.4 dB @ 20 ms |
| Detection latency | dwells to declare, with hysteresis | 2–3 ms |
| Authentic array gain | mean over the visible sky, dB re one element | +4.7 dB (power max, rank 1) |
| Satellite outage | P(gain < single antenna) | 7% (rank 1), 22% (rank 2) |
| Usable J/N | 1 dB post-null loss | 51–57 dB (12-bit) |
| Group delay | antenna to receiver, and its **stability** | 2.8 µs; stability is the spec that matters |
| Max platform rotation | for < 3 dB null degradation | ~400 °/s @ 1 ms dwell |

**Define null depth relative to a single antenna element, not relative to the quiescent beam.** The latter is not a property of the beamformer (Task 2, M4) and a competitor using the honest definition will look worse than you while being better.

**Benchmark against a null hypothesis you can defend.** Publish the *no-spoofer* case alongside every performance number: "under attack we deliver X; with no attack present we degrade the receiver by Y". If Y is not ≈ 0, you do not have a product.

---

## 11. Failure modes and how the system must behave

The design rule: **every failure degrades toward a plain single-antenna receiver, and every degradation is reported.** Never crash, never silently continue in a compromised state.

| failure | detection | response |
|---|---|---|
| One channel dead (LNA, cable, ADC) | per-channel power monitor | drop to N−1 elements, report degraded |
| All channels dead | power monitor | fail-safe, report |
| Inter-chip sample misalignment | boot pilot test | auto-resync; if it fails, single-antenna mode + alarm |
| ADC saturation | per-dwell counter | invalidate dwell, hold weights, report J/N |
| AGC gain step mid-dwell | gain-change flag | invalidate dwell, hold weights |
| Degenerate weights (‖f‖ → 0, spoofer aligned with **h**) | norm check | fall back to quiescent beam; **never `error()`** |
| Estimator not converged | eigenvalue sanity | hold previous weights |
| Detector false alarm | hysteresis (2-of-3 to engage, 10 to release) | avoids flapping |
| Calibration EEPROM absent/corrupt | CRC | run uncalibrated, report reduced null-depth spec |
| Temperature out of range | sensor | report, widen calibration uncertainty |
| **Spoofer at a satellite bearing** | rank/attribution | do not null blindly — exclude the measurement (SKU-B) |

**On the last row.** If a spoofer happens to lie along a real satellite's bearing, nulling it destroys a genuine measurement. The correct response is not a deeper null but *attribution*: identify the counterfeit PRNs post-correlation and exclude them. This is another argument for SKU-B and it is the case where a naive nulling-only product does actual harm.

---

## 12. Reliability, security and scalability

**Reliability.** MTBF is dominated by the RF section and the connectors, not the FPGA. Four coax connectors on a moving platform are four failure points; consider an integrated antenna-plus-frontend assembly with a single digital cable to eliminate them. Derate the LNAs for the actual thermal environment (rooftop enclosures reach 70 °C+).

**Security, which the paper does not touch and a defence or infrastructure customer will ask about first.**
* **Signed firmware and secure boot.** An anti-spoofing device with unsigned firmware is an attack surface, not a defence.
* **Tamper-evident calibration storage.** Corrupting the equaliser table silently degrades nulling — an elegant attack.
* **Authenticated telemetry.** If the device reports "no spoofing detected", that report must itself be authenticated, or it is worthless.
* **Forensic capture.** A rolling DDR buffer that snapshots raw samples on detection is worth a great deal commercially: it turns each incident into evidence, and it turns your fleet into a threat-intelligence source.

**Scalability.** Design the digital chain parametric in N *now* — the reference model already is. The natural SKU ladder follows the transceiver's 2×2 structure: **2 / 4 / 8 elements**, not 4 / 7. Eight elements gives 6 spare adaptive DOF, +8.3 dB array gain, 2.7% outage and −24.5 dB null, from four AD9361s or two ADRV9026s. Adding constellations and frequencies costs correlators and RF, not algorithm.

---

## 13. Ten things industry does that the paper does not discuss

1. **Detect before you mitigate.** Unconditional nulling has negative expected value.
2. **Calibrate the array anyway**, even with a calibration-free algorithm. It unlocks MVDR, LCMV and attribution.
3. **Ship the calibration inside the antenna assembly**, bound to a serial number.
4. **One LO, one clock, verified at every boot.**
5. **Common AGC, and a per-gain-index phase table.** Independent AGC quietly destroys the array.
6. **Digital per-channel equalisation.** The difference between a −35 dB and a −55 dB null.
7. **Phase-continuous weight updates**, or post-correlation combining. Otherwise no RTK, no precise timing.
8. **Attribution, not just suppression.** Exclude bad measurements and publish a protection level.
9. **Test the no-threat case as rigorously as the threat case.**
10. **Forensic capture on detection.** Cheap to add, disproportionately valuable commercially.

---

## 14. Suggested sequence

**Foundation.** Adopt the covariance/EVD estimator and the detector (Tasks 2 B1, M2) in the reference model. These are the two changes with the highest ratio of impact to effort, and everything else depends on them. Validate against `verify/` and `studies/`.

**Hardware bring-up.** 2× AD9361 with shared LO on an existing carrier (Zedboard/ZC706 + FMCOMMS5, or an ADI-based SDR with an external LO mod). Bring up MCS, the boot alignment self-test, common AGC and frozen calibrations. Expect the majority of surprises here, not in the algorithm.

**Digital chain.** DFE with channel equaliser, covariance engine, Jacobi EVD, detector, beamformer. Bit-exact against the reference model at every boundary.

**Conducted test campaign.** Reproduce the paper's Figure 8 setup. Sweep SAPR, spoofer bearing, and the no-threat case. This is where the specification numbers come from.

**Receiver integration (SKU-B).** Correlator bank, post-correlation attribution, per-satellite combining, RAIM on surviving measurements, integrity message.

**Productisation.** Antenna assembly with EEPROM, production test sequence, chamber campaign, environmental, EMC, secure boot.

The dependency that matters: **the antenna assembly design and its calibration procedure gate the specification**, so start the chamber work early. It is the long pole, and it is the one most teams start last.
