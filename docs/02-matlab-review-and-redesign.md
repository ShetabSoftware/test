# Task 2 — MATLAB Review and Redesign

Three files were supplied:

| # | File | Role |
|---|---|---|
| 1 | `GPS_PowerMax_KnownDoppler_Simulation.m` | Known-Doppler acquisition + eigen SSV + null steering |
| 2 | `GPS_SPOOFER_V9_2502.m` | Paper's γ/β estimator + null steering + beam pattern |
| 3 | `FPGA_FixedPoint_NullSteering_QuantTest.m` | Fixed-point quantisation study |

The redesigned reference model is in `matlab/`. It runs unmodified on MATLAB R2018b+ and GNU Octave 7+ with **no toolboxes**, and every claim below is backed by a script in `matlab/verify/` or `matlab/studies/` that you can run.

---

## 0. What is right — stated first, because it matters

I want to be clear about what I am *not* criticising, because several of your choices are correct and a couple are better than the paper.

* **Your C/A code generator (File 1) is correct.** I verified all 32 PRNs against IS-GPS-200 Table 3-Ia (first-ten-chips octal), plus code balance (512/511), the three-valued autocorrelation {1023, 63, −1, −65}, and the −23.9 dB cross-correlation bound. `verify/test_ca_code.m`. The LFSR structure, tap table, ±1 mapping and initialisation are all right. I reused it.
* **The projection algebra is correct.** `f = h - y*((y'*h)/(y'*y))` is exactly eq. (14), and your inline FPGA note that `v = Σ conj(f_i)·r_i` because MATLAB's `f'*r` conjugates the weights is correct and is precisely the kind of thing that silently steers the null to the mirror direction when it is got wrong.
* **The paper's eq. (18) is correct and you implemented it faithfully.** `h = Pa/‖Pa‖` genuinely maximises output *SNR*, not merely output power: since `Ph = h` for any h in the range of P, the noise gain σ²‖Ph‖² is unchanged, so maximising signal power under ‖h‖ = 1 maximises SNR. Some reviewers get this wrong.
* **Applying a single global AGC scalar across all channels** (Files 1 and 3) is the correct hardware model. Per-channel AGC is a genuine failure mode and you avoided it, whether or not that was deliberate.
* **Modelling the spoofer as replaying the same PRN set** with different delays and Dopplers is right, and File 2 asserts it explicitly.
* **The comment in File 1** — *"If no spoof is present, power maximization locks onto authentic GPS. A real system must use spoof detection and bypass nulling in GPS-only"* — is the single most important engineering observation in all three files. It identifies the design's most serious defect. It just needs to be implemented rather than noted.

Now the problems.

---

## 1. Findings by severity

### BLOCKING — would ship a broken product

**B1. There is no spoofing detector. The pipeline nulls unconditionally.**

All three scripts compute weights and apply them regardless of whether a spoofer exists. In a clean environment the estimator's principal direction is simply the strongest authentic satellite, so the device steers a null into it.

This inverts the product's value proposition. An anti-spoofing box spends essentially all of its operating hours *not* under attack. Expected value is

```
P(attack) × benefit  −  P(no attack) × damage
```

and with P(attack) small, the second term dominates. Measured: with no spoofer present and nulling forced on, the mean authentic gain drops by 4–5 dB and 18.6% (N=3) / 6.9% (N=4) of satellites end up worse than a single passive antenna — permanently, for no reason.

*Fix:* `core/asp_detect.m`. Eigenvalue statistic λ₁/mean(λ₂…λ_N) with an MDL rank estimate, gating the weight update in `asp_process.m`. Measured minimum detectable SAPR is +0.3 dB at 1 ms and −1.4 dB at 20 ms on a 4-element array, P_FA = 10⁻³. The detector is essentially free: the eigenvalues come out of the same decomposition that produces the weights.

**B2. The operating point in File 2 is wrong by 87 dB.**

```matlab
P_auth_dBW  = -80.5;   % "nominal GPS L1 C/A power per authentic PRN"
sigma2 = k_B*Tk*B;     % B = 2e6  ->  -140.9 dBW
```

That places every GPS signal **60 dB above** the thermal noise floor. The correct pre-despreading sample SNR is

$$\mathrm{SNR}_\text{sample} = \frac{C/N_0}{f_s} = \frac{10^{4.5}}{16.368\times10^6} = -27.1\ \text{dB}$$

i.e. **27 dB below** the noise. An 87 dB error in the one parameter the estimator's entire behaviour depends on. Every null depth, SSV correlation and beam pattern from that script is measured in a regime that does not exist, where the problem is trivially easy.

The comment `% paper simulation value per spoofing PRN` on `P_spoof_dBW = -50.0` is also wrong: that is SAPR = +30.5 dB, where the paper uses −153 vs −158.5 dBW, i.e. **5.5 dB**. Files 1 and 3 use `noise_rms = 0.50` with `auth_amp = 1.0`, which is +6 dB — better, but still 33 dB adrift.

*Fix:* `config/asp_config.m` parameterises powers as **C/N₀ in dB-Hz** and derives sample SNR from fs. This is the only parameterisation that stays correct when the sample rate changes, and it makes the operating point impossible to get wrong by accident.

**B3. The processing is non-causal.**

All three scripts estimate the SSV from the entire record and then apply the resulting weights to that same record, including the samples used for the estimate. This is not implementable, and it flatters the result: it removes the estimator transient entirely and hides the coupling between weight-update rate and platform dynamics.

*Fix:* `core/asp_process.m` applies weights from dwell *k* to dwell *k+1*. Weight age is then a first-order design parameter: with D = 0.82 λ, the residual phase after a rotation Δθ caps the null at 20 log₁₀(2π(D/λ)Δθ/√3) — about −46 dB at 100°/s, −34 dB at 400°/s, −26 dB at 1000°/s with a 1 ms dwell. Useful conclusion: **dwell length is set by estimator statistics, not platform dynamics, up to roughly 1000°/s.**

**B4. `error()` on the algorithmic path.**

```matlab
if norm(f) < eps
    error('Null weight calculation failed: zero norm.');
end
```

A shipping product must degrade, never crash. And the degenerate case is *reachable*: `f = h − y(yᴴh)/(yᴴy)` vanishes when **y** ∥ **h**, i.e. when the spoofer's signature aligns with the quiescent beam. With `h = ones/√N` that is a spoofer at zenith. Unlikely for a terrestrial spoofer, but "unlikely" is not a specification, and near-degenerate cases produce a small ‖f‖ that then gets normalised, amplifying quantisation noise into the beam.

*Fix:* every failure path in the redesign sets a status and falls back — to the quiescent beam, then to a single antenna. `asp_process.m` also holds the previous weights whenever a dwell is invalidated by an AGC step or by ADC saturation, on the principle that a slightly stale estimate always beats a poisoned one.

**B5. The fixed-point study does not test the thing that breaks.**

File 3 quantises the input samples, `y_est`, `f` and the beamformer output — and runs the estimator itself in double precision. Its own comment concedes it: *"Correlation accumulators are not saturated in this MATLAB model."* The accumulator is exactly where fixed-point covariance estimation goes wrong.

Measured on a bit-exact model (`verify/test_fixedpoint.m`), rounding the sample products *before* the covariance accumulator injects an error into R that is:

1. **deterministic** — a fixed matrix, not noise;
2. **rank one and aligned with the boresight steering vector** (measured alignment 1.000 with all-ones when the RTL Hermitises, 0.906 when it mirrors the upper triangle; a random direction gives 0.5). Boresight is where the satellites are. **A truncating implementation synthesises a phantom source at zenith and steers the null into the sky;**
3. **independent of dwell length** — measured ‖E‖_F identical to three digits across a 16× change in K.

Property 3 is what makes it dangerous. Everything else in this estimator improves as 1/√K, so a designer who validates at 1 ms reasonably expects improvement at 100 ms. With product rounding the floor is flat and nothing in a floating-point simulation predicts it.

| product fractional bits retained | 10 | 14 | 18 | 22 | 26 | full (30) |
|---|---|---|---|---|---|---|
| resulting null-depth floor | −25.4 | −48.9 | −73.0 | −97.1 | −122.0 | −156 dB |

6.02 dB per retained bit, exactly as theory predicts.

**The fix is free.** A DSP48 holds the exact 32-bit product of two 16-bit words and accumulates in its 48-bit P register with *no rounding at all*. At K = 16368 the accumulator needs 42 bits — measured, with 6 bits spare in the 48-bit budget. Rounding inside the accumulation loop is a design choice, and it is the wrong one. If a narrower accumulator is ever unavoidable, convergent (round-half-to-even) rounding converts the flat floor back into a 1/√K error: measured 0.28× over a 16× change in K, against the 0.25 that 1/√K predicts.

**B6. The word-length plan is asserted, not derived.**

```matlab
fmt.sample.WL = 18;  fmt.sample.FL = 17;   % from a 12-bit AD9361
fmt.weight.WL = 24;  fmt.weight.FL = 22;
```

Where do 18 and 24 come from? The AD9361 delivers 12 bits; decimation adds ~1–2 bits of processing gain. And weight quantisation has a closed-form effect on null depth:

$$\text{null floor} = 20\log_{10}(2^{-\mathrm{FL}_w}) - 1.76\ \text{dB}$$

FL = 12 gives −75 dB, FL = 16 gives −98 dB. Measured −97.4 dB against −98.1 dB predicted. Since the estimator itself only reaches −24 dB, **anything at or above 12 fractional bits is irrelevant**; Q1.22 in a 24-bit word buys nothing and costs a DSP cascade on Intel parts where the native multiplier is 18×19.

More importantly, the study never asks the question that converter word length actually answers: **how much interference can this receiver null before quantisation noise dominates?** After nulling, the interferer is gone but the quantisation noise it forced the AGC to expose is not — quantisation noise is spatially white and survives every beamformer. With the AGC holding the composite at backoff *B* below full scale:

$$L = 10\log_{10}\!\left(1 + \frac{1/12}{\sigma_\text{th}^2}\right),\qquad \sigma_\text{th} = 2^{WL-1}\,10^{-B/20}\,10^{-\mathrm{J/N}/20}$$

| word length | 8 | 10 | **12 (AD9361)** | 14 | **16** |
|---|---|---|---|---|---|
| usable J/N at 1 dB loss | 33 dB | 45 dB | **57 dB** | 69 dB | **81 dB** |

Derate the AD9361 to ~51 dB for its ~10.5-bit ENOB. This is the number that decides Task 4, and it is absent from the study. `analysis/asp_adc_jn_limit.m`, `config/asp_fx_plan.m`.

### MAJOR — correct, but significantly suboptimal

**M1. β is redundant *and* inferior. Delete it.**

The paper's eq. (7) takes the *phase* from γᵢ = Σ rᵢ r₁* and the *magnitude* from √|βᵢ|, βᵢ = Σ rᵢ(n) rᵢ*(n−T). The magnitude/phase split is essential — I confirm this below — but β is the wrong source for the magnitude.

From the paper's own eq. (11), βᵢ → |Cᵢ|²·d with

$$d = K\sum_k p_k^s e^{j2\pi f_k^s T} + K\sum_m p_m^a e^{j2\pi f_m^a T}$$

With |f| up to 5 kHz and T = 1 ms the phasors wrap many times, so **d is a random walk over ~18 emitters**, of typical magnitude √18·p rather than 18·p, and occasionally near zero. Measured across 40 scenes on the full waveform model: |β| varies by **13.7×**. The covariance diagonal, which estimates the same |Cᵢ|², varies by a few percent.

Why the diagonal works: the channel mismatch is *post-LNA* (cables, filters, mixer, transceiver gain), so it scales signal and noise identically, and R_ii = |Cᵢ|²(P_total + σ₀²) is proportional to |Cᵢ|² with the bracket common to all elements — and the projector is scale invariant, so the common factor cancels.

Three costs of keeping β:

* **Memory:** one full code period of delay per channel. At K = 16368 and 16-bit I/Q that is 2.10 Mbit ≈ **58 BRAM36** on a 7-series part (44 for N = 3), to compute a statistic the covariance already contains.
* **Accuracy:** measured null depth, full waveform model, 4-element Y array:

| dwell | paper γ/β | `column` (raw R column) | **`gamma` (diagonal + phase)** | **`evd`** |
|---|---|---|---|---|
| 1 ms | — (needs 2 epochs) | −6.5 | **−21.8** | **−23.7** |
| 2 ms | −15.8 | −6.5 | −22.7 | **−25.3** |
| 5 ms | −19.5 | −6.5 | −23.2 | **−26.3** |
| 10 ms | −21.5 | −6.5 | −23.4 | **−26.7** |
| 20 ms | −22.2 | −6.6 | −23.2 | **−26.8** |

* **Latency:** β needs two epochs, so the paper's estimator cannot produce an answer in the first millisecond. The covariance can.

*Aside, and a compliment:* the `column` row shows the paper's magnitude/phase split is **not stylistic**. R₁₁ = Σ|r₁|² is dominated by *noise* power while R_i1 (i ≠ 1) contains only the spoofer cross-term, so the raw column's first entry is larger than the rest by σ²/P_s ≈ 12 dB. Using it directly gives ρ = 0.74 and a −6.5 dB null. Splitting magnitude from phase is what makes the estimator work at all, and preserving that structure in your code was right.

**M2. Use the full covariance and its eigenvector. This is the central redesign.**

Here is the observation that unifies everything. Let D = diag(√R₁₁,…,√R_NN) and R_w = D⁻¹RD⁻¹ be the whitened (unit-diagonal) covariance. Then

> **the paper's estimator is exactly D × (the reference column of R_w), which is one step of power iteration started from e_ref.**

Replacing "one step" with "the converged answer" *is* the improvement:

$$\hat{\mathbf{y}}_\text{paper} \sim D\cdot R_w \mathbf{e}_\text{ref} \qquad\longrightarrow\qquad \hat{\mathbf{y}}_\text{evd} = D\cdot u_1(R_w)$$

Whitening is not decoration. With post-LNA mismatch, R = C R₀ Cᴴ, and the principal eigenvector of R is *not* C**b̄** unless C is a scalar times a unitary. Whitening by the measured diagonal turns C into a pure-phase diagonal, which *is* unitary, so u₁(R_w) is the phase-only steering vector and D·u₁ recovers C**b̄** in the measured domain. **This is what preserves the paper's calibration-free property under an eigen-based estimator** — plain `eig(R)` does not have it.

Cost, N = 4: 10 correlators instead of 7 (32 real multiplies vs 28, +14%). At N = 3 the covariance is actually *cheaper* than γ/β (18 vs 20). What it buys:

| | paper γ/β | covariance + EVD |
|---|---|---|
| Detection statistic | none | all N eigenvalues, free |
| Rank > 1 nulling | no | yes |
| R⁻¹ for MVDR/LCMV | no | free from the EVD |
| Diagonal loading | n/a | free (shift Λ) |
| Epoch delay memory | 44–58 BRAM36 | none |
| First answer available | 2 ms | 1 ms |
| Null depth @ 2 ms | −15.8 dB | −25.3 dB |

**M3. Rank-1 only — and the ground bounce is the single largest unmodelled effect in the design.**

A terrestrial spoofer at 100 m with a 2 m antenna height produces a specular reflection ~20 ns later at −3 to −10 dB, from a *different* elevation. None of the three scripts models it. It turns out to matter more than anything else.

Twenty nanoseconds is 0.02 C/A chips, so the bounce is **not resolvable** at C/A bandwidth — it adds coherently inside a single correlation cell. That does not make it harmless; it makes the spoofing source **partially coherent**. The covariance of a two-ray source is

$$R_s = P_s\left[\mathbf{bb}^H + \alpha^2\mathbf{b_m b_m}^H + \alpha\rho\,(\mathbf{b b_m}^H + \mathbf{b_m b}^H)\right],\qquad \rho = \mathrm{sinc}(B\,\Delta t)$$

which is **rank two**, with a second eigenvalue that grows as the processing bandwidth widens and ρ falls away from 1. A rank-one projector cannot remove it, so the achievable null is capped at λ₂/λ₁. For α = −6 dB and Δt = 20 ns:

| processing bandwidth | ρ | rank-1 null cap |
|---|---|---|
| 4.1 MHz | 0.989 | −24.4 dB |
| 16.4 MHz | 0.804 | **−11.9 dB** |

The 16 MHz figure is 12 dB *worse* than what the estimator itself achieves, so it — not the estimator — would set system performance. **This is the one limit in the whole design that gets worse as you improve the front end**, which is the direction every other consideration pushes you. It is a genuine architectural fork, developed in Task 3 and Task 5:

1. narrow the processing bandwidth toward the C/A main lobe, at the cost of code-tracking resolution;
2. spend a spatial degree of freedom on a rank-2 null — which a 3-element array does not have to spare (Task 1, §2);
3. go to **space-time adaptive processing**, where per-element tapped delay lines null the delayed replica with *taps* instead of with spatial degrees of freedom. This is what production wideband anti-jam CRPAs do, and this measurement is the reason why.

`asp_scenario.m` models it; `asp_detect.m` estimates the rank via MDL; `asp_weights('project', Y, h)` takes a Y of any rank; `studies/study_multipath.m` measures the effect against the closed-form cap.

*A negative result worth recording:* I expected the **post-correlation** covariance to reveal the bounce through its λ₁/λ₂ ratio. It does not — measured 4103 without a bounce and 9189 with one, i.e. no drop at all. The reason is the same coherence: at 0.02 chips the despread snapshot sees one composite vector **b** + α**b**ₘ, not two. Detecting the bounce requires *delay* resolution — extra correlator taps or wider bandwidth — and its effect on nulling is a wideband effect that only appears across the band.

**M4. The null-depth metric is not a property of the beamformer.**

```matlab
NullDepth_dB = 10*log10(abs(f'*b_vector)^2 / abs(h'*b_vector)^2)
```

This measures the null relative to whatever gain the *arbitrary* quiescent vector **h** happened to have toward the spoofer. If **h** is unlucky and already has 6 dB of loss that way, the reported null depth is 6 dB better than the truth.

*Fix:* `analysis/asp_pattern_metrics.m` uses one consistent definition throughout,

$$G(\mathbf{a}) = \frac{|\mathbf{f}^H\mathbf{a}|^2}{\|\mathbf{f}\|^2}$$

which equals 1 for a single antenna element (since |aᵢ| = 1) and has unit mean over directions with E[**aa**ᴴ] = I. So dB values read directly as "gain over one antenna", and a beamformer cannot show positive average gain — it can only redistribute. Both metrics are reported so the numbers stay comparable to your scripts.

**M5. Three separate `norm()` normalisations, all removable.**

`y_est/norm(y_est)`, `f/norm(f)` and `q_m/norm(q_m)` are each an inverse square root — a CORDIC or Newton-Raphson block in RTL with its own latency, pipeline and corner cases. None is necessary:

* the projector I − **yy**ᴴ/(**y**ᴴ**y**) is invariant to the scale of **y**;
* the beamformer output feeds a correlator, C/N₀ estimator and tracking loop, all invariant to a constant complex gain;
* likewise **q**ₘ.

What *does* matter is that the weights sit at the top of their word. `fx/fx_bfp_scale.m` does that with a leading-zero count and a barrel shift: ~30 LUTs, one cycle, no divider, no convergence question, and **exactly zero** additional error. File 2 also computes `den_y = real(y_est'*y_est)` *after* normalising `y_est` to unit norm, so that division is by 1.0.

**M6. Do not form the projector.** `P = I − yyᴴ/(yᴴy)` then `f = Ph` costs N² storage and N² multiplies per beam. The factored form `f = h − Q(Qᴴh)` costs 2·N·rank: 8 multiplies instead of 16 at N = 4, rank 1, and the gap widens with N. `asp_weights('project', …)` orthonormalises Y by modified Gram-Schmidt first, which is also far better conditioned than forming YᴴY and inverting it.

**M7. File 1's snapshot handling loses information.**

```matlab
z = z*exp(-1j*angle(z(ref_idx)));      % de-rotate
snapshots(:,m) = z/(norm(z) + eps);    % normalise
R = snapshots*snapshots';
```

Two problems. The de-rotation is **worse than useless**: `z*z'` is invariant to a common phase, so the rotation cannot change R — but `angle(z(ref_idx))` is a *noisy* estimate, and multiplying by it injects reference-channel noise into every other element for no benefit. And normalising each snapshot to unit norm **discards the SNR weighting**, so a weak, noise-dominated PRN peak contributes to R exactly as much as a strong clean one. Accumulate `z*z'` unnormalised, or weight by peak power.

**M8. File 1's estimator needs known Doppler, which discards the paper's main advantage.**

The paper's whole selling point is that the SSV is obtained *before* despreading, at low complexity and with no 2-D search. File 1's estimator acquires each PRN with known Doppler first. That is a legitimate and in fact *better* estimator — it is essentially the stage-2 refinement I recommend — but it cannot run at cold start, cannot run during reacquisition (exactly when an attack is most effective), and File 1 provides no fallback. The right answer is not to choose: run both, as two stages with different jobs (§3).

Also, File 1 selects the global max over delay per PRN with **no consistency check**. If the spoofer does not transmit some PRN, that PRN's global max is the *authentic* peak and it contaminates the SSV estimate. The check is cheap and is itself an excellent detector: if all selected peaks share a spatial direction, λ₁/Σλ ≈ 1 and they are spoofed; if not, they are not.

**M9. Two flaws in the paper's power-maximisation unit that will bite when you implement it.**

*Doppler aliasing.* Eq. (20) accumulates e^(−j2πf_m T w) over snapshots w spaced by T = 1 ms, i.e. **sampling the Doppler phasor at 1 kHz**. GPS Dopplers span ±5 kHz. Two satellites whose Dopplers differ by a multiple of 1 kHz are indistinguishable, and with ~10 satellites the chance of a collision is not small. *Fix:* wipe off each PRN's carrier before accumulation (which a receiver does anyway), or use sub-epoch snapshots, or work from the receiver's own prompt correlators where the residual Doppler is near zero by construction.

*Spoofer leakage in the despreading reference.* Eq. (19) correlates the *projected* vector **x** against the conjugate of one period of the **raw reference antenna** r₁ — which still contains the spoofer at full strength. The spoof×spoof cross term is suppressed only by the single null depth (~25 dB), against a spoofer whose total power exceeds the authentic by 15 dB. *Fix:* use the projected reference x₁, so the residual is squared (~−50 dB). One-line change, ~25 dB improvement.

**M10. `f_m = P·ĥ_m` is a no-op.** With ĥₘ = **q**ₘ/‖**q**ₘ‖ and **q**ₘ = P(·) already in the range of P, P**q**ₘ = **q**ₘ, so eq. (23) reduces to fₘ = ĥₘ. Harmless, and re-projecting is arguably a useful cheap error-cleanup if the projector was recomputed in between — but it should be a deliberate choice with a comment, not an unnoticed redundancy.

**M11. Nothing physical is modelled.** No front-end filter, IQ imbalance, DC offset, LO phase noise, inter-channel group delay, mutual coupling, AGC dynamics, or ADC saturation. For an "algorithm reference for an FPGA implementation" these are not optional: **inter-channel group-delay skew is the impairment that actually limits null depth in hardware.** 100 ps of skew caps the null at −55 dB; 1 ns caps it at −35 dB, which is *above* what the algorithm delivers. A model that cannot express this cannot be used to justify a null-depth specification.

### MODERATE — code quality that matters at product scale

| | Finding |
|---|---|
| C1 | `rng()` re-seeded with fixed constants *inside* functions. Monte Carlo is impossible, and in File 1 both the GPS-only and GPS+spoof cases share an identical noise realisation. Good for a controlled A/B; fatal for statistics. |
| C2 | `code_fft = fft(codes(m,:))` recomputed inside the epoch loop. |
| C3 | `corr_store = zeros(N_ant,K,num_epochs)` — 3.2 MB per PRN, allocated in a loop, to find one maximum. A running max needs O(1). |
| C4 | Beam-pattern loop is 361×91 with a matrix-vector product inside. Vectorise. |
| C5 | Large commented-out blocks (LPF, file I/O, ADC export) that encode real intent but cannot be tested. Promote to tested functions — see `export/asp_export_vectors.m`. |
| C6 | `spoof_off = 0` set and never used. |
| C7 | Missing semicolon on `NullDepth_dB`; `caxis` is deprecated in favour of `clim`. |
| C8 | Delays applied by `circshift` of an integer number of samples, so every correlation peak lands exactly on a sample. Real delays are continuous; sample-aligned delays are optimistic. |
| C9 | No code Doppler. Negligible over 1 ms (1.6×10⁻³ chips), 1.6 chips over 1 s — i.e. total decorrelation. Silently misleads anyone who later lengthens the dwell. |
| C10 | No navigation data bits. (Worth noting: I *expected* this to break the paper's β and tested it — it does not, because the affected factor d is common-mode and cancels in the scale-invariant projector. Measured ρ = 0.99668 with bits vs 0.99702 without. A useful negative result: bit synchronisation is **not** a prerequisite for the array processing.) |
| C11 | `auth_el = rand(1,N)*90` is uniform in elevation *angle*, which over-populates zenith by 1/cos(el) and understates how often a satellite sits near a low-elevation spoofer — the case the projection null damages. Use uniform in solid angle. |
| C12 | Validation correlates only the first epoch and uses genie delay/Doppler, so it tests neither acquisition nor the estimator's effect on it. |
| C13 | Three files with independently hard-coded, mutually inconsistent parameters (fs, N_ant, powers, seeds, word lengths). This is how the 87 dB error in B2 survived. |

---

## 2. The redesigned reference model

```
matlab/
  asp_startup.m            path setup
  asp_run_all.m            full regression + studies
  config/
    asp_config.m           SINGLE source of truth; powers as C/N0
    asp_fx_plan.m          word lengths DERIVED, with the derivations
  fx/                      fixed-point primitives (no toolbox)
    fx_fmt, fx_quant, fx_round_convergent
    fx_cov_accum           bit-exact DSP48 model
    fx_bfp_scale           power-of-two normalisation
  model/                   simulation only, NOT for synthesis
    asp_ca_code            ICD-verified Gold codes
    asp_array_geometry     + baseline / lattice metadata
    asp_steering           carrier phase AND the matching delay
    asp_scenario           emitters, channel, impairments
    asp_rx_generate        block-streaming waveform generator
    asp_agc_adc            AGC + quantisation + saturation
  core/                    THE ALGORITHM REFERENCE (maps to RTL)
    asp_ssv_paper          baseline, for comparison
    asp_ssv_from_cov       column / gamma / evd estimators
    asp_evd_herm           fixed-sweep cyclic Jacobi
    asp_detect             eigenvalue detection + MDL rank
    asp_weights            projection / MVDR / LCMV / MRC
    asp_beamform           streaming, bit-exact fixed point
    asp_postcorr_ssv       stage-2 refinement
    asp_process            causal streaming pipeline
  analysis/                metrics and theory
  verify/                  regression tests
  studies/                 the Monte Carlo evidence
  export/                  RTL golden vectors
```

### Design rules the reference obeys

1. **Deterministic.** Fixed Jacobi sweep count, no tolerance-based loops, no data-dependent iteration counts anywhere in `core/`. A hardware block with data-dependent latency is a timing-closure hazard and a real-time scheduling hazard.
2. **Streaming and block-based.** Nothing allocates O(total samples). The generator produces one dwell at a time, which is how the FPGA consumes data and what makes causal weight application natural.
3. **No dividers, no square roots on the sample path.** Every normalisation is a power-of-two shift. The only reciprocal-square-roots are N per millisecond (the whitening), which is a 1 kHz event.
4. **No matrix inverse anywhere.** Projection is factored. MVDR and LCMV come from the eigen-decomposition. Condition number never enters the datapath.
5. **Norm-preserving numerics.** Jacobi rotations are unitary, so ‖A‖_F is invariant *exactly* (verified to 1.05×10⁻¹⁵). No dynamic-range growth, no intermediate rescaling.
6. **Exact accumulation.** Full-width products into a wide accumulator, one rounding at the output. Bit-exact, verified.
7. **Hardware in the loop.** AGC, quantisation, saturation and dwell validity are *inside* the processing loop, not applied afterwards to an already-computed floating-point answer.
8. **Every claim tested.** `verify/` is a regression suite, not documentation.

### What maps to what

| Paper | Original scripts | Redesign | Change |
|---|---|---|---|
| eq. (8) γᵢ | `phase_acc` | R off-diagonals | now part of the covariance |
| eq. (8) βᵢ | `beta_acc` + K-sample delay line | R diagonal | 58 BRAM36 → 0, +8 dB statistic |
| eq. (7) y | `sqrt(abs(beta)).*exp(1j*angle(phase))` | `asp_ssv_from_cov` | one power-iteration step → converged eigenvector |
| eq. (13) P | explicit `I − yyᴴ/(yᴴy)` | factored, never formed | N² → 2N·rank multiplies |
| eq. (14) f | `h − y(yᴴh)/(yᴴy)`, normalised | same, BFP-scaled | divider removed |
| eq. (19)-(22) | known-Doppler acquisition | `asp_postcorr_ssv` | fixes aliasing + reference leakage |
| eq. (23) fₘ | `P·ĥₘ` | `ĥₘ` (identity noted) | redundant matvec removed |
| — | *(absent)* | `asp_detect` | **the missing block** |

---

## 3. The recommended architecture: two stages that fail differently

The most valuable thing the measurements say is *where the pre-correlation estimator stops improving*, because that tells you what a longer dwell cannot fix.

Measured EVD null depth: −23.7 dB at 1 ms → −26.8 dB at 20 ms. A noise-limited estimator would have gained 13 dB over that range; it gained 3.1. The shortfall is the authentic constellation's contribution Σₘ pₘ**aₘaₘ**ᴴ to the covariance, which tilts the principal eigenvector by an amount set by p_auth/P_spoof and **not** by K.

So: **stop integrating and change domain.**

```
   4x RF  ──►  DFE  ──┬──────────────────────────────────────►  beamformer ──► receiver
                      │                                              ▲
                      ▼                                              │
              covariance + Jacobi EVD                                │
                      │                                              │
                      ├──► detector (gate + rank) ──► weights ───────┤
                      │                                              │
                      ▼                                              │
              STAGE 1: open loop, ~1 ms, no receiver state           │
              22-25 dB, works at cold start                          │
                                                                     │
              STAGE 2: closed loop, ~100 ms                          │
              despread snapshots ──► covariance ──► EVD ─────────────┘
              35-45 dB, needs tracking state
```

**Stage 1** (pre-correlation) needs no receiver state, so it works from a cold start and during reacquisition — exactly when an attack is most effective. It reaches −23.7 dB at 1 ms and saturates near −26.8 dB, which is enough to drop the spoofer below the authentic signals so the receiver can acquire the *real* peaks. Bias limited by the authentic constellation.

**Stage 2** (post-correlation) refines the signature from despread snapshots, raising the wanted term by the full 36–42 dB despreading gain. It is **also bias limited** — I expected it to become noise limited and the measurement says otherwise — but by a much lower floor set by the C/A code itself. Measured across a 20 dB sweep (10 scenes each, 20 ms coherent):

| SAPR | 0 dB | 5.5 dB | 12 dB | 20 dB |
|---|---|---|---|---|
| stage-2 null depth | −26.96 | −32.47 | −38.96 | −46.92 dB |

$$\boxed{\text{stage-2 null depth} = -(\mathrm{SAPR} + 27.0)\ \text{dB}}$$

to better than 0.1 dB across the whole sweep. The 27.0 is the C/A cross-correlation bound of 23.9 dB (65/1023) plus ~3.1 dB, because 65/1023 is a worst-case peak and the RMS over random code phases sits below it. The 1 dB-per-dB slope is the authentic contribution falling relative to the spoofer.

**Honest accounting of what stage 2 is worth.** At 5.5 dB SAPR it buys about **6 dB** over stage 1's saturation point, not the 20 dB one might hope for. The deeper null is a bonus. Its real value is what stage 1 cannot do at *any* dwell length:

* **Per-PRN attribution.** Comparing each PRN's despread spatial snapshot against the estimated spoofing signature tells you *which* PRNs are counterfeit. Stage 1 sees only an aggregate. This is what lets the receiver exclude specific measurements rather than discarding the whole solution — and it is the basis of any credible integrity claim (Task 5, RAIM).
* **Per-satellite optimal combining with no Doppler ambiguity**, because the residual Doppler after the tracking loop is near zero by construction. This is the clean fix for M9.
* **Rank and structure** of the spoofing source.

And (*) tells you where to look if you ever need to go deeper than −(SAPR+27): the limiting resource is the **code**, not the array and not the dwell. That points at the modernised pilots (L1C and L5, 10230 chips, roughly 10 dB better cross-correlation) or at successive interference cancellation.

The two stages fail for different reasons, so combining them is coverage, not redundancy. Your File 1 was already reaching for stage 2; it needs stage 1 underneath it as the bootstrap and fallback.

---

## 4. Measured before/after

Full waveform model, 4-element Y array, SAPR 5.5 dB, C/N₀ 45 dB-Hz, real C/A codes, navigation data, code Doppler, channel mismatch.

| Metric | Original (paper γ/β) | Redesign (whitened EVD) | Redesign + stage 2 |
|---|---|---|---|
| Null depth @ 1 ms | not available (needs 2 epochs) | **−23.7 dB** | — |
| Null depth @ 2 ms | −15.8 dB | **−25.3 dB** (+9.5) | — |
| Null depth @ 20 ms | −22.2 dB | **−26.8 dB** (+4.6) | see `study_postcorr.m` |
| SSV correlation @ 20 ms | 0.99759 | **0.99912** | — |
| Delay memory | 58 BRAM36 | **0** | 0 |
| Streaming multiplies (N=4) | 28 | 32 (+14%) | +despreader |
| Spoofing detector | none | **λ₁/λ̄ + MDL, free** | rank from λ₁/λ₂ |
| Min detectable SAPR @ 1 ms | n/a | **+0.3 dB** | — |
| Rank-2 capable | no | **yes** | yes |
| Dividers on the sample path | 3 | **0** | 0 |
| Matrix inversions | 0 | **0** | 0 |
| Causal | no | **yes** | yes |
| Crashes on degenerate input | yes (`error()`) | **no** (graceful degradation) | no |

---

## 5. How to run it

```matlab
run('matlab/asp_startup.m')
asp_run_all                 % regression suite + short studies
asp_run_all('full')         % + the long Monte Carlo studies

% or individually
test_ca_code; test_evd_herm; test_fixedpoint; test_cov_model
study_array_size(1500)
study_estimators(40)
study_detector(600)
study_postcorr(20)
```

`export/asp_export_vectors.m` writes stimulus and expected-response files for RTL co-simulation, so the reference model is the golden model rather than a separate document that drifts from it.
