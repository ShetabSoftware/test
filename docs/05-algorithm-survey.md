# Task 5 — Alternative and Modern Anti-Spoofing Approaches

## 0. Framing: the paper's method is a *mitigator*, and you need three things

The 2012 method is a good spoofing **mitigator** with a genuinely elegant property (no calibration) and a real weakness (it saturates around −27 dB and it cannot tell you *which* measurements are counterfeit). A commercial product needs three functions and they are not the same problem:

| function | question answered | the paper |
|---|---|---|
| **Detection** | am I under attack? | **absent** |
| **Mitigation** | can I still navigate? | yes, this is the contribution |
| **Attribution / integrity** | which measurements do I trust, and can I bound my error? | absent |

Attribution is what regulated markets actually buy. A box that says "spoofing suppressed" is a nice-to-have; a box that says "PRNs 4, 7 and 13 are counterfeit, excluded, and the remaining solution has a 12 m protection level at 10⁻⁷" is a certifiable component. Build toward the second.

A second framing point that shapes everything below. **Spatial processing is the only technique on this list that works against an arbitrarily sophisticated adversary using a single transmitter**, because it exploits physics the attacker cannot control — one transmit antenna means one direction of arrival. Every non-spatial technique listed here is defeatable by a sufficiently capable spoofer. So: spatial processing is the backbone; everything else is corroboration. Do not invert that.

---

## 1. Spatial techniques

### 1.1 Orthogonal projection / blind null steering — the paper's method

**How.** Estimate the dominant spatial signature, project it out. No manifold needed.

* **Advantages.** No calibration. No DOA search. Pre-correlation, so it protects acquisition — which matters because acquisition is when a receiver is most vulnerable. O(N) streaming complexity. Works from a cold start.
* **Disadvantages.** Bias-limited at ~−27 dB by the authentic constellation. Rank-1 by construction as published. No detection, no attribution. Assumes the spoofer dominates the spatial energy — false at SAPR < 0 dB and false for a matched-power drag-off attack.
* **Complexity.** O(N) per sample; O(N³) per dwell.
* **FPGA suitability.** Excellent — see Task 3. ~20 DSP48.
* **Difficulty.** Low. This is the easiest thing on the list.
* **Expected performance.** −22 to −27 dB, detection down to ~0 dB SAPR.

**Verdict: keep as stage 1.** It is the right bootstrap and the right fallback.

### 1.2 Eigen-decomposition of the spatial covariance

Already the recommended replacement for the paper's estimator (Task 2, M2). Everything the paper achieves plus detection, rank estimation, R⁻¹ for free, and no epoch delay memory, for 3 extra correlators at N = 4.

* **Complexity.** N(N+1)/2 correlators + one 4×4 Hermitian EVD at 1 kHz.
* **FPGA suitability.** Excellent. Cyclic Jacobi is norm-preserving, so no dynamic-range growth and no conditioning concerns. 17.7 µs at 163.68 MHz.
* **Verdict: adopt. This is the single highest-value change in the whole review.**

### 1.3 MVDR (Capon)

$$\mathbf{w} = \frac{R^{-1}\mathbf{a}}{\mathbf{a}^HR^{-1}\mathbf{a}}$$

* **Advantages.** Strictly dominates projection when the manifold is known: it suppresses *all* interference including residual spoofer, multipath and unmodelled sources, not merely the one direction you explicitly nulled. Given the EVD you already computed, R⁻¹ = UΛ⁻¹Uᴴ is **free**.
* **Disadvantages.** Requires a calibrated manifold and a known desired direction — i.e. attitude plus almanac plus a chamber campaign. Notoriously sensitive to steering-vector error: a small pointing error causes *signal self-nulling*, which is catastrophic and looks exactly like a spoofing success.
* **Complexity.** Free if you have the EVD.
* **FPGA suitability.** Good, with the important caveat that **diagonal loading is not optional.** Load to roughly the level of the steering-vector uncertainty (Carlson, *IEEE T-AES* 1988). Better still, use a robust variant — worst-case optimisation over an uncertainty ball (Vorobyov, Gershman & Luo, *IEEE T-SP* 2003) — which is a small extra cost and removes the self-nulling failure mode.
* **Verdict: phase 2, once the array is calibrated.** Not for the first product.

### 1.4 LCMV — the right formulation once you have a manifold

$$\mathbf{w} = R^{-1}C(C^HR^{-1}C)^{-1}\mathbf{g},\qquad C = [\mathbf{a}_\text{sat},\;\mathbf{b}_\text{spoof}],\;\mathbf{g}=[1,0]^T$$

This is the natural formulation for this product because it expresses exactly what you want: unit gain toward the satellite, a hard null on the spoofer, and minimum output power with everything left over. It also makes the DOF ledger explicit and is where the 3-vs-4 element argument becomes sharpest (Task 1, §2).

* **Advantages.** Combines everything; hard constraints; degrades gracefully as you add constraints.
* **Disadvantages.** Manifold required. Each constraint costs a DOF. Constraint matrix conditioning must be watched when two constrained directions are close.
* **FPGA suitability.** Good — a 2×2 or 3×3 solve at 1 kHz.
* **Verdict: the target architecture for a calibrated product.**

### 1.5 MUSIC / ESPRIT

* **Advantages.** MUSIC gives super-resolution DOA. ESPRIT is closed-form (no search) and cheap if the array has a shift-invariant structure. DOA of every source is powerful *attribution* evidence: if 9 PRNs all arrive from one bearing, that is a spoofing declaration you can put in a log with a confidence number.
* **Disadvantages.** MUSIC needs a spectral search over the 2-D manifold (expensive) and a **precisely calibrated manifold** — its resolution advantage evaporates with a few degrees of calibration error. Both need more sensors than sources; here there are 10–20 sources and 4 sensors, so **neither works pre-correlation.** Both need a noise subspace, which does not exist when the sources fill the array.
* **Where they *do* work: post-correlation.** After despreading, each PRN is isolated and you have one source per snapshot. Then DOA estimation is well posed and cheap.
* **Complexity.** MUSIC O(N² · grid); ESPRIT O(N³).
* **FPGA suitability.** Poor for MUSIC (2-D search). ESPRIT is feasible but needs a doublet array geometry that constrains your layout.
* **Verdict: post-correlation only, for attribution rather than mitigation.** Do not attempt pre-correlation MUSIC with 4 elements and 18 sources.

### 1.6 STAP — space-time adaptive processing

Replace each element's single complex weight with an M-tap FIR: NM degrees of freedom instead of N.

* **Advantages.** The standard solution for **wideband** nulling. A single spatial weight can only null a source perfectly at one frequency; taps null it across the band. Essential when nulling a wideband jammer whose null depth would otherwise be limited by the aperture and by inter-channel group-delay skew. Also nulls *delayed replicas* with taps instead of with spatial degrees of freedom — which is the right answer to a resolvable multipath reflection, and preserves spatial DOF for the satellites.
* **Disadvantages.** M× the weights, M× the covariance dimension (NM × NM — a 4×5 STAP has a 20×20 covariance, so the EVD goes from 4×4 to 20×20, roughly 100× the arithmetic). Introduces a **frequency-dependent group delay** into the signal path, which biases pseudoranges differently per satellite — a serious issue for a timing or RTK product and a well-known headache in military CRPAs. And each tap costs a degree of freedom that could have been spatial.
* **Complexity.** O((NM)³) per update, O(NM) per sample.
* **FPGA suitability.** Still very tractable at N = 4, M = 3–5: 20×20 Jacobi is ~190 rotations/sweep, ~1 ms at 163 MHz — borderline for a 1 kHz update, so use a 10 ms update or a sweep-parallel implementation.
* **Verdict: not for spoofing; yes if anti-jam is on the roadmap.** For pure spoofing at C/A bandwidth, spatial-only is sufficient and STAP's pseudorange bias is a real cost. Design the datapath so STAP can be added (a tapped delay line ahead of the covariance engine) without an architectural change.

### 1.7 Blind source separation / ICA

* **Advantages.** No manifold, no training. Can separate more sources than a rank-1 projector.
* **Disadvantages.** GNSS signals are 27 dB below the noise pre-despreading, so the non-Gaussianity that ICA exploits is invisible. Convergence is not guaranteed, iteration counts are data-dependent (fatal for a fixed-latency pipeline), and there is a permutation/scaling ambiguity that you must resolve *anyway* using something else.
* **FPGA suitability.** Poor — iterative, data-dependent latency, needs nonlinearities.
* **Verdict: no.** This is a solution looking for a problem here. The covariance approach gets the same answer deterministically.

---

## 2. Non-spatial techniques worth combining

None of these can stand alone against a competent adversary. All of them are cheap, and several are nearly free once you have the receiver. Use them as corroborating evidence in a fused detector (§4).

| technique | detects | cost | defeated by |
|---|---|---|---|
| **C/N₀ + AGC monitoring** | crude power anomalies; AGC gain drop is a strong jamming/spoofing indicator | ~zero | matched-power spoofer |
| **Clock-bias / drift consistency** | receiver clock jumps when it switches to counterfeit signals | ~zero | slow drag-off |
| **Doppler vs position-rate consistency** | spoofed Dopplers inconsistent with the claimed motion | low | good simulator |
| **PRN code-phase / Doppler cross-check against almanac** | signals that should not be visible | low | attacker with the almanac |
| **Signal Quality Monitoring (multi-correlator)** | correlation-peak distortion when authentic and counterfeit overlap | moderate (extra taps) | non-overlapping attack |
| **Carrier-phase double difference across two antennas** | all counterfeit PRNs share a spatial signature | low if you already have the array | nothing simple — this is strong |
| **RAIM / ARAIM** | inconsistent pseudoranges | low | a self-consistent spoofer (which is easy to build) |
| **IMU / odometry coupling** | position solution diverging from inertial | moderate | slow drag-off within IMU drift |
| **Multi-constellation cross-check** | most spoofers do GPS L1 C/A only | low | multi-constellation spoofer |
| **Multi-frequency (L1/L2/L5)** | most spoofers are single-frequency; ionospheric divergence is checkable | high (RF cost) | multi-frequency spoofer |
| **Galileo OSNMA** | **cryptographically authenticated navigation data** | low (software) | replay within the TESLA latency |
| **Chips-Message Robust Authentication (Chimera)** / GPS spreading-code authentication | cryptographic, at the spreading-code level | needs key distribution | not much |

### Two of these deserve emphasis

**Galileo OSNMA is operational and free.** Open Service Navigation Message Authentication reached its operational service declaration in 2023 and authenticates Galileo navigation data via a TESLA hash chain. It is the only *cryptographic* assurance available today on an open civil signal, it costs software only, and it defeats every spoofer that generates its own navigation data. Its limitation is that it does not authenticate the *ranging code*, so a meaconing/replay attack within the TESLA disclosure latency still works — which is exactly the gap that spatial processing fills. **OSNMA and a CRPA are close to complementary, and a product that has both has a genuinely strong story.** If your product will ever see Galileo, budget for OSNMA.

**RAIM's assumptions are false under spoofing, and you should say so in your documentation.** Classical RAIM assumes independent, single-satellite faults. A spoofer generates a *self-consistent* constellation, so all pseudoranges agree perfectly and RAIM sees no residual at all. RAIM is not an anti-spoofing technique and marketing it as one is misleading. Where it *does* help is after your array processing has excluded the counterfeit measurements: then RAIM/ARAIM computes the protection level on what remains, which is what makes an integrity claim quantitative. **Use RAIM downstream of attribution, never as a substitute for it.**

---

## 3. Machine learning — practically justified?

**Mostly no, with two narrow exceptions.**

The general case against: the physics here is fully known, and the optimal estimators are derivable in closed form. A learned detector would have to be trained on spoofing attacks, which are rare, diverse and adversarial — you would be fitting to your own simulator. It would then be a black box in a device whose entire value proposition is trustworthiness, at exactly the point where a regulator asks "why did it decide that?" And its false-alarm rate would be uncharacterisable, whereas the eigenvalue detector's threshold has an analytic H₀ distribution you can calibrate and defend.

Two exceptions where it is genuinely justified:

1. **Multi-feature fusion.** You will have 10–15 heterogeneous detection features (§2) with correlated, non-Gaussian statistics. Hand-tuning that decision boundary is exactly the thing a small, *interpretable* model does better than a human — a logistic regression or a shallow gradient-boosted tree over a fixed feature vector, with the features themselves being physics-derived and individually explainable. This keeps the physics in the features and the ML in the arbitration, which is the defensible split.
2. **Anomaly detection on the RF environment** for *forensics and threat intelligence*, not for real-time decisions. Clustering recorded spectra and spatial signatures across a fleet tells you what your customers are actually seeing. That is product intelligence, and it does not need to be trustworthy in the same way.

Everything else — CNNs on correlator outputs, learned beamformers — is worse than the closed-form answer, harder to certify, and orders of magnitude more expensive in the FPGA.

---

## 4. What I would actually build

A five-layer defence, ordered by latency, with each layer covering the one below's blind spot.

```
 L1  RF / ADC          AGC and saturation monitoring          continuous
                       -> jamming and gross power anomaly
 L2  Pre-correlation   covariance + eigen detection + null    1 ms
     SPATIAL           -> "a dominant spatial source exists"; suppress it
                       -> protects ACQUISITION; works cold-start
 L3  Post-correlation  per-PRN spatial signatures             20-100 ms
     SPATIAL           -> "PRNs 4,7,13 share a bearing" -> ATTRIBUTION
                       -> per-satellite optimal combining
                       -> deeper null: -(SAPR+27) dB
 L4  Measurement       clock/Doppler/C-N0 consistency,        1 s
     domain            multi-constellation, OSNMA
                       -> corroboration, cryptographic where available
 L5  Solution domain   RAIM/ARAIM on the SURVIVING            1 s
                       measurements, IMU coupling
                       -> protection level, integrity claim
```

**L2 and L3 are the backbone** because they exploit physics the attacker cannot fake with one transmitter. L4 and L5 are corroboration and are individually defeatable. L1 is nearly free and catches the crude cases that are, in practice, most of them.

**The most important interface in this diagram is L3 → L5.** Attribution turns "we suppressed something" into "these measurements are excluded and here is the protection level on the rest". That is the difference between a gadget and a certifiable component, and it is entirely absent from the paper.

### Priority order for engineering effort

1. **Eigenvalue detection + gating** (Task 2, B1). Without it the product has negative value most of the time. Small change, largest impact.
2. **Covariance/EVD estimator** replacing γ/β (Task 2, M2). 9.5 dB at 2 ms, removes 58 BRAM, enables everything else.
3. **Post-correlation attribution** (L3). The commercially differentiating feature.
4. **Rank-2 capability** (Task 2, M3). Needed the moment there is a second emitter.
5. **Array calibration + LCMV** (§1.4). Unlocks the optimal beamformer.
6. **OSNMA** if Galileo is in scope. Cheap, cryptographic, complementary.
7. **STAP** only if anti-jam becomes a requirement — and note the pseudorange-bias cost.

### What I would *not* build

* Pre-correlation MUSIC/ESPRIT with 4 elements and 18 sources — ill-posed, not merely hard.
* ICA/BSS — non-deterministic latency, no advantage over the covariance.
* A learned end-to-end detector — uncertifiable, and worse than the closed-form answer.
* RAIM marketed as anti-spoofing — it does not work against a self-consistent spoofer, and claiming otherwise will not survive a serious customer's technical review.

---

## 5. Threat cases your current design does not cover

Worth being explicit, because these define the roadmap and because a customer will ask.

| threat | why the current design struggles | what to add |
|---|---|---|
| **Matched-power drag-off** — starts aligned with the authentic signals, walks slowly | No dominant spatial energy at onset. Detection saturates at ~−1 dB SAPR regardless of dwell | L3 attribution; clock/IMU consistency; carrier-phase double difference |
| **Multiple distributed spoofers** — several transmitters, several bearings | Consumes DOF fast; 4 elements handles 2 sources, not 4 | More elements (8); the SKU ladder in Task 1 |
| **Meaconing / replay** — rebroadcast of genuine signals | Signals are authentic, so cryptography does not help; only the *direction* is wrong | Spatial processing is the primary defence here, which is a strong argument for the whole approach |
| **Spoofer co-located with a satellite bearing** | Null lands on a real satellite | Rank-1 null plus L3 attribution to exclude, rather than relying on nulling |
| **Jam-then-spoof** | Jamming forces reacquisition; needs J/N headroom you may not have at 12 bits | 14/16-bit converters (Task 4); STAP |
| **Spoofer on a moving platform** | Signature moves; 1 ms dwell tracks a few hundred °/s | Adequate as designed; verify against your platform's dynamics |
