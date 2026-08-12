# Follow-up Engineering Decisions

Answers to four specific questions. Rationale only where it changes a number.

---

## 1. 12-bit DAC output format

Implemented in `stage10_txdac` of `matlab/golden/asp_golden_model.m`.

| | |
|---|---|
| Format | **signed, s12.11 (Q0.11)** — 1 sign bit, 0 integer bits, 11 fractional |
| Integer range | **[−2048, +2047]**, two's complement (AD9361 TX data port) |
| Real range | [−1.0, +0.99951], LSB = 2⁻¹¹ = 4.883×10⁻⁴ |
| Target level | **256 LSB RMS per component** = −18.1 dBFS |
| Crest headroom | 8σ, clip probability < 10⁻¹⁵ per sample |
| Gain control | **power-of-two barrel shift** with hysteresis at 181 / 362 LSB |
| Rounding | **convergent** (round-half-to-even) |
| Overflow | **hard saturation**, never wrap |

**Why 18 dB of backoff on TX when RX uses 14 dB.** The two sides solve opposite problems. On RX, backoff buys headroom for an interferer you have not removed yet. On TX the interferer is already nulled, so there is nothing to leave room for, and backing off further costs only quantisation noise — 10·log₁₀(1 + (1/12)/256²) = 5×10⁻⁶ dB, i.e. nothing. Clipping, by contrast, is a memoryless nonlinearity that would intermodulate the residual spoofer back into the band and partially undo the nulling. On TX, clipping margin is free; buy it.

**Why a shift and not a multiplier.** Amplitude steps on TX are harmless — unlike RX, where a per-channel gain step corrupts the array mid-dwell. So a coarse power-of-two gain is sufficient, and a barrel shifter plus a leading-zero count adds exactly zero error where a fractional multiplier would add a rounding.

**Why convergent rounding here.** Truncation would place a −0.5 LSB DC offset on the DAC, which becomes a spur at the TX LO and can disturb the downstream receiver's DC and AGC logic. It costs one OR gate. (This is a much weaker requirement than in the covariance accumulator, where truncation bias is a correctness bug rather than a cosmetic one.)

```matlab
% tx_scale.vhd
sh = F_BEAM - F_DAC + shift;                 % net right shift
y  = shiftRoundSat(real(v), sh, 12) + 1i*shiftRoundSat(imag(v), sh, 12);
% shiftRoundSat = arithmetic shift, round-half-to-even, clamp to [-2048, 2047]

% AGC update, once per dwell, hysteretic
rms = sqrt(sum(real(v).^2 + imag(v).^2)/(2*K)) / 2^sh;
if     rms > 362, shift = shift + 1;
elseif rms < 181, shift = shift - 1;  end
```

First dwell uses fast acquisition (jump straight to the right shift) rather than stepping one bit per dwell; otherwise the first several milliseconds clip, and clipped dwells have to be discarded.

Measured over 6 ms: 229.5 LSB RMS, **0 clipped samples of 49104**.

---

## 2. LO offset versus DC notch

**Use A: offset the LO and do a digital DDC.** Not a close call, and the deciding argument is one that neither option's usual framing mentions.

**The decider.** LO leakage appears as a *per-channel* DC offset. Four channels with four different DC offsets form a rank-one term **dd**ᴴ in the spatial covariance — which is exactly the signature of a spatially coherent source. It is therefore not merely an artefact to be cleaned up; it is a **phantom emitter that the eigen-detector will find and null**. With no real spoofer present, it is the *dominant* eigenvalue and the system nulls a direction determined by nothing but analogue mismatch.

That reframes the choice. Option B has to notch DC *before the covariance*, deeply, at the exact frequency where the C/A spectrum peaks, per channel, while the AD9361's own DC-tracking calibration is asynchronously moving the thing being notched. Option A moves it to a frequency where an ordinary decimation filter removes it by 70 dB with no signal loss and no interaction with anything.

Secondary advantages of A, each real:

* I/Q image becomes **separable by filtering** instead of folding onto the signal.
* 1/f noise corner and even-order distortion products land outside the signal band.
* No dependence on AD9361 DC-tracking calibration, which must be frozen anyway (a tracking cal that runs mid-dwell changes the channel response and corrupts the covariance).

Cost of A: one shared NCO and four complex multipliers. That is 5–7 DSP48 with time multiplexing.

### Recommended DDC

| | |
|---|---|
| RX LO | **1577.466 MHz** = L1 + 2.046 MHz |
| Offset | **2.046 MHz = f_ADC/16 exactly** |
| f_ADC | **32.736 MHz** (32 × 1.023) |
| NCO | 16-entry ROM indexed by a 4-bit counter — **no phase accumulator**, hence no phase-truncation spurs |
| NCO format | s16.14 (Q1.14) so +1.0 is representable; ±16384 |
| Mixer | one complex multiply per channel, 3 DSP48 each with the pre-adder |
| Stage 1 | halfband ÷2 → **16.368 MHz** (16 samples/chip, 16368 per ms), 11 taps, 3 multipliers |
| Stage 2 | 63-tap symmetric FIR, passband ±1.2 MHz, stopband from 2.046 MHz, ≥70 dB |
| Working rate | 16.368 MHz |
| TX | LO = **1573.374 MHz** (L1 − 2.046 MHz), digital upconversion by the same NCO, so TX LO leakage lands 2 MHz off L1 |

**Choose the offset as an exact binary fraction of the sample rate.** f_ADC/16 makes the NCO a 16-entry table with no accumulator and no truncation spur. f_ADC/4 would make the mixer multiplier-free ({1, j, −1, −j}) but needs an 8.184 MHz offset and a wider analogue passband; f_ADC/16 is the better balance.

**One NCO drives all four mixers.** A common complex rotation applied to every element is invisible to the array algorithm — the covariance, the projector and the beamformer are all invariant to a common complex scalar — so a shared NCO contributes **exactly zero** inter-channel mismatch. Four independent NCOs would not.

No CIC anywhere. Decimation is only ÷2, and CIC passband droop plus its non-linear phase would have to be equalised per channel. Halfband FIRs are linear phase with identical group delay by construction, which is what keeps the four channels aligned.

Fixed point: ADC s12.11 → mixer s16.15 (product s28, shift 10, round, saturate) → halfband s16.15 → FIR s16.15. Coefficients s18.17, DC gain forced to exactly 2¹⁷ so the filters are unity-gain and no scaling drift accumulates.

---

## 3. Building the receiver

**Yes — and for adaptive, multi-directional attacks it is not an enhancement, it is the enabling capability.** Pre-correlation spatial processing gives you aggregate statistics over a mixture. It cannot tell you *which* PRNs are counterfeit, and against several transmitters it runs out of degrees of freedom before it runs out of threats.

### The discriminant that makes multi-directional attacks tractable

An attacker with T transmitters spoofing M PRNs is constrained by physics: **the M despread spatial signatures collapse into T clusters.** Authentic signals give M distinct signatures because they come from M satellites.

So the test is *number of distinct spatial clusters versus number of tracked PRNs* — not power, not direction, not any property the attacker controls cheaply. Adding transmitters costs the attacker linearly and only raises T; it never makes T equal M. Ten PRNs from four transmitters is still four clusters against ten.

This test needs per-PRN spatial signatures, which exist only after despreading. That is the argument for the receiver, and it is sufficient on its own.

### What despreading adds, quantitatively

| Quantity | Pre-correlation | Post-correlation |
|---|---|---|
| Per-PRN spatial signature | no | yes |
| Effective SNR of the spatial estimate | −26 dB/sample | +36–42 dB processing gain |
| Which PRNs are spoofed | no | yes |
| Simultaneous authentic + counterfeit for one PRN | invisible | two peaks, different (τ, f_d), different signature |
| Achievable null | −24 dB, saturates | −(SAPR + 27) dB |

### Algorithms, ranked by value for this threat model

| Method | What it catches | FPGA suitability | Verdict |
|---|---|---|---|
| **Per-PRN spatial clustering** (signature + MDL/k-means on cluster count) | multi-directional spoofing, any power | correlators in PL, clustering in PS at 1–10 Hz | **build first** |
| **Cross-antenna carrier-phase single difference** | all counterfeit PRNs share one differential phase; authentic ones do not | nearly free once you correlate per antenna | **build second** — strongest cost/benefit on the list |
| **Dual-peak / CAF ridge search** | overlapping authentic + counterfeit, drag-off onset | needs the acquisition engine anyway | build third |
| **LCMV with per-PRN constraints** | beamform toward authentic while nulling identified spoof directions | R⁻¹ free from the existing EVD; needs a calibrated manifold | after calibration |
| MUSIC / ESPRIT post-correlation | actual DOA per PRN, for logging and attribution | well posed only post-correlation (one source per snapshot); MUSIC's 2-D search belongs in PS | optional |
| Successive interference cancellation | strong spoofer masking a weak authentic signal | heavy; PL regeneration + PS control | later |
| RAIM / ARAIM on surviving measurements | protection level after exclusion | PS, 1 Hz | required for an integrity claim |

Pre-correlation MUSIC/ESPRIT remains ill-posed: 4 sensors, 18 sources, no noise subspace. Do not attempt it.

### The degree-of-freedom reality, and what follows from it

With N = 4 you have 3 usable nulls. A four-transmitter attack exhausts them and leaves nothing to beamform with. **Against multi-directional attacks the answer is not deeper nulling — it is attribution and exclusion**: identify the counterfeit PRNs, drop those measurements, and beamform toward the ones that survive. Nulling then only has to handle whichever one or two emitters are strong enough to raise the noise floor.

If multi-transmitter attacks are genuinely in your threat model, this is the strongest argument for **8 elements** (7 DOF, 4 AD9361s or 2 ADRV9026s). The digital chain is already parametric in N.

Split of work: correlator bank, per-PRN snapshot accumulation and per-PRN beamforming in the PL; clustering, rank decisions, DOA, exclusion logic, RAIM and all hysteresis policy in the PS at 1–10 Hz, where they stay field-tunable.

---

## 4. Golden reference model

`matlab/golden/asp_golden_model.m` — one file, no toolboxes, MATLAB R2018b+ or Octave 7+.

```matlab
G = asp_golden_model('durationMs',6, 'outDir','./gold_vectors');
```

Every signal is a MATLAB double holding an **exact integer** — the raw two's-complement register contents of the corresponding VHDL signal. Nothing in the datapath is a fraction; scaling lives in the comments and constants, never in the data. So a MATLAB value and a `std_logic_vector` are the same number and can be diffed with **no tolerance**.

**Block map — one local function per VHDL entity**

| MATLAB | VHDL | Format out |
|---|---|---|
| `stage1_ddc` | `ddc_mixer.vhd` | s16.15 @ 32.736 MHz ×4 |
| `stage2_hbdec` | `hb_decim2.vhd` | s16.15 @ 16.368 MHz ×4 |
| `stage3_fir` | `fir_shape.vhd` | s16.15 @ 16.368 MHz ×4 |
| `stage4_cov` | `cov_accum.vhd` | s48 × 10 entries @ 1 kHz |
| `stage5_whiten` | `whiten.vhd` | s32.26 @ 1 kHz |
| `stage6_evd` | `jacobi_evd.vhd` | s20.16 (U), s32.26 (λ) |
| `stage7_detect` | `detect.vhd` | integer compare, **no divider** |
| `stage8_weights` | `weight_calc.vhd` | s18.16 × 4 complex |
| `stage9_beamform` | `beamformer.vhd` | s16.15 @ 16.368 MHz |
| `stage10_txdac` | `tx_scale.vhd` | **s12.11** @ 16.368 MHz |

**Rules obeyed, each because breaking it makes the VHDL un-diffable**

* every multiply is followed by an explicit shift, round and saturate;
* no divisions, square roots, sin, cos or atan on the datapath — all transcendentals are CORDIC or Newton iterations with a **fixed** iteration count, modelled bit for bit;
* no `filter()`, `fft()`, `eig()`, `inv()` or `norm()` anywhere in the chain;
* fixed iteration counts everywhere, no tolerance-based loops;
* vectorised expressions appear only where the hardware also accumulates at full precision and rounds once (FIR taps, covariance accumulation).

**Transcendentals, specified exactly**

* CORDIC, 16 iterations, angle scaled by 2¹⁶, gain compensated by the seed 1/K = 39797 (Q16). Two vectoring operations give the rotation angles; one rotation operation gives cos/sin in Q1.16, which are then applied with multipliers — deriving the twiddle once per rotation avoids imposing the CORDIC gain on the data at every pass.
* Reciprocal square root: even-shift range reduction to [0.25, 1), chord seed, **4 Newton iterations**, fixed. Range reduction is what makes the iteration count data-independent and hence the latency deterministic.
* Integer square root: 24-step non-restoring array.

**Verification output.** Each stage writes its input and output to `sNN_<name>_{in,out}.txt` as signed decimal integers, one per line, complex interleaved I,Q, multi-channel column-major. `MANIFEST.txt` records the format, the rates, and every constant the VHDL package needs (NCO table, filter taps, CORDIC parameters, thresholds).

**Measured, 6 ms, 4-element Y array, SAPR 5.5 dB**

| | |
|---|---|
| Determinism | **bit-identical across runs** at every stage |
| Null depth, fixed point | −24.53 dB |
| Null depth, double precision, same data | −24.44 dB |
| **Fixed-point loss** | **0.09 dB** |
| Worst eigenvalue relative error vs `eig()` | 2.4×10⁻³ |
| Worst weight angle vs double-precision solution | 0.066° |
| Covariance accumulator | 42 of 48 bits used |
| DAC | 229.5 LSB RMS, 0 clipped samples of 49104 |

The self-check is deliberately a comparison against **double precision on the same data**, not against the model's own earlier output. A bit-exact model that computes the wrong thing is worse than no model, because it makes the VHDL wrong too.
