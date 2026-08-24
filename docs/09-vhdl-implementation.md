# VHDL implementation for XC7Z020

The finalized MATLAB model (`matlab/golden/asp_golden_model.m`) is the sole
functional reference. Every RTL block in `rtl/` is co-simulated against that
model's own dumped stage vectors with **no tolerance** — both sides are
integers, so any difference is a bug rather than noise.

---

## 1. Verification status

Run `./rtl/sim/run_sim.sh`. Every number below is a measured result from a clean
build over 4 ms of the reference scenario, not an intention. **13 of 13
testbenches pass; about 2.2 million values are bit-identical.**

| Block | Stage | Values compared | Result |
|---|---|---|---|
| `asp_ddc_mixer` | 1 | 1,047,552 | bit-identical |
| `asp_hb_decim2` | 2 | 400,000 | bit-identical |
| `asp_fir_shape` | 3 | 523,528 | bit-identical |
| `asp_cov_accum` | 4 | 128 (4 dwells, 48-bit) | bit-identical + Hermitian verified |
| `asp_rsqrt` / `asp_isqrt` | — | 371 + 377 cases | bit-identical |
| `asp_cordic` | — | 544 vectoring + 408 rotation | bit-identical |
| `asp_whiten` | 5 | 128 | bit-identical, diagonal real |
| `asp_jacobi_evd` | 6 | 144 (eigenvectors **and** eigenvalues) | bit-identical, trace preserved |
| `asp_detect` | 7 | 16 | bit-identical |
| `asp_weight_calc` | 8 | 32 + rank-0 path | bit-identical |
| `asp_beamformer` | 9 | 130,944 | bit-identical |
| `asp_tx_scale` | 10 | 98,208 (dwell 2 onward) | bit-identical |
| `asp_datapath` | 1–10 | 24 weights, from raw ADC samples | bit-identical, detector fired on every dwell |

The last row is the one that matters most. `tb_asp_datapath` runs several million
clocks from raw AD9361 samples and checks the **weight sequence** against the
model. The weights are a function of every stage from the mixer through the EVD —
mixer, halfband, FIR, covariance, whitening, reciprocal-sqrt, integer sqrt,
CORDIC, 36 Jacobi rotations, the sort, the detector and the Gram-Schmidt
projection. If any one of them drifted by a single LSB the eigenvector would
rotate and the weights would not match.

Two deviations from the model exist, both forced, both documented in the RTL
where they bite:

- **`asp_fir_shape`** computes the causal convolution and discards 31 warm-up
  outputs per channel; the model uses the centred form. After the discard the
  two are sample-aligned.
- **`asp_tx_scale`** cannot reproduce dwell 1. The model's fast acquisition
  measures a dwell's RMS and applies the result to *that same dwell*, which is
  not causal. Reproducing it would mean buffering 16,368 complex samples
  (≈15 BRAM36) to improve one millisecond after reset. The RTL scales dwell 1
  with a register and performs the identical jump at the end of it, so dwell 2
  onward is exact.

---

## 2. PL / PS partitioning

**In the PL** — anything at sample rate, or that must be deterministic:

| Stage | Block | Rate | Why it cannot be software |
|---|---|---|---|
| 1 | DDC mixer | 130.944 MS/s | sample rate |
| 2 | Halfband ÷2 | 130.944 MS/s | sample rate |
| 3 | Shaping FIR | 65.472 MS/s | 4190 M MAC/s |
| 4 | Covariance | 65.472 MS/s in | 32 MAC per sample vector |
| 5–8 | Whiten, EVD, detect, weights | 1 kHz | must complete inside one dwell with bounded latency |
| 9 | Beamformer | 16.368 MS/s out | sample rate |
| 10 | TX scale / DAC | 16.368 MS/s | sample rate |

Stages 5–8 run at 1 kHz and would fit in software on raw throughput. They are
in the PL anyway because the weight update has to land on a dwell boundary with
**bounded** latency, and a Linux userspace process cannot promise that. They
occupy about 4% of a dwell in the PL, which is why they are also heavily
resource-shared rather than parallel.

**In the PS** — policy, not arithmetic:

- **MDL rank estimation and detection hysteresis.** The PL exports both sides of
  the detector's cross-multiplied comparison (`DET_LHS`, `DET_RHS`) so the PS can
  apply any policy at 1 kHz without needing a divider. The M-of-N vote before
  committing to a null lives here — and it is load-bearing: the detector
  threshold trades false-alarm margin against sensitivity, and the hysteresis is
  what absorbs an isolated tail excursion.
- **AD9361 configuration over SPI**, calibration, and **freezing the DC tracking
  calibration**. A tracking cal that runs mid-dwell changes the channel response
  and corrupts the covariance.
- **Threshold re-calibration** whenever `fs`, the FIR, the dwell length or the
  satellite count changes.
- Logging, health monitoring, operator interface.

Not yet in either: the post-correlation stage (per-PRN spatial clustering). It
needs correlators in the PL and clustering in the PS. The hooks are the
AXI4-Stream capture port and the dwell interrupt.

---

## 3. Architecture decisions and their consequences

### 3.1 One clock domain at 130.944 MHz

`clk_dsp = 130.944 MHz = 4 × FS_ADC = 8 × FS_WORK`, derived from the AD9361
sample clock — **not** from the PS PLL.

Every rate in the datapath is an exact integer ratio of the ADC clock, so every
rate change is a clock enable rather than a FIFO, and **there is no true
clock-domain crossing anywhere between the ADC pins and the DAC pins**. That
removes an entire class of metastability, reconvergence and gray-code bugs
rather than managing them.

The price: the PS AXI clock is asynchronous to the datapath. That is dealt with
once, by an AXI Clock Converter — a reviewed IP instance instead of a
synchroniser per register.

### 3.2 Channel-serial (TDM) datapath

The four antennas traverse the **same** logic rather than four copies of it.

This is not primarily a resource decision, though it is a large one — it folds
4 complex multipliers into 1 in the mixer and 256 into 32 in the FIR, which is
what makes the design fit an XC7Z020 at all. The more important property is that
the four antennas **cannot acquire different group delays**, because they pass
through the same registers. Inter-channel delay mismatch is the one error the
array processing downstream can neither tolerate nor detect.

### 3.3 Weight latency is two dwells, not one

The model applies dwell *k*'s estimate to dwell *k+1*. Hardware cannot: the
covariance for dwell *k* only closes at the end of dwell *k*, and the estimator
chain then needs

```
whiten ~270 + EVD ~4700 + detect 2 + weights ~350 ≈ 5400 clocks
```

which is about 4% of a dwell. The new weights are ready *partway into* dwell
*k+1*.

Switching weights mid-dwell would put a discontinuity inside a coherent
integration period, so the weight register is loaded only at a dwell boundary
and dwell *k*'s estimate is applied from the start of dwell *k+2*.

The weights are therefore 1 ms staler than the model's. The spatial signature of
a spoofer moves on the timescale of platform and satellite motion — seconds — so
this changes the null depth by far less than the dwell-to-dwell estimation
noise. It is not a free choice; it is what the estimator's own latency costs.

### 3.4 No backpressure

Deliberately. Every rate is an exact integer ratio, so the pipeline is
rate-matched by construction and a `tready` would never deassert. The one place
elasticity *is* needed — the burst out of the decimator — has an 8-deep FIFO
inside `asp_fir_shape` with a **sticky overflow flag** wired to a status
register, so a rate-plan error shows up as a bit rather than as silently
corrupted data.

---

## 4. Fixed-point specification

Every value in the datapath is a raw two's-complement integer. `W` is total bits
including sign; `F` records where the binary point sits, and is used only to
derive shifts.

| Signal | Format | Range | Set by |
|---|---|---|---|
| ADC input | s12.11 | ±2048 | AD9361 RX port |
| NCO twiddle | s16.14 | ±16384 | Q1.14 so +1.0 is representable |
| Mixer output | s16.15 | — | `>> (F_ADC+F_NCO−F_MIX) = >>10` |
| FIR coefficients | s18.17 | DC gain exactly 2¹⁷ | unity gain, no scaling drift |
| Working sample | s16.15 | — | `>>17` after each filter |
| Covariance accumulator | s48 | peak 2⁴⁵ | DSP48E1 P register |
| Whitened matrix | s32.26 | diagonal ≈ 2²⁶ | see below |
| CORDIC rotations | s18.16 | — | Q1.16 |
| Eigenvector accumulator | s20 | — | |
| Weights | s18.16 | — | DSP48E1 B port |
| Beamformer output | s16.15 | — | `>> (F_DAT+F_WGT−F_BEAM) = >>16` |
| DAC output | s12.11 | ±2048 | AD9361 TX port |

### Three rounding rules, and why they are distinct

`asp_pkg` carries three separate primitives because substituting one for another
passes most vectors and fails on exact ties:

- `conv_round_shr` — round half to **even**. The datapath default. Zero mean
  error, unlike truncation (−0.5 LSB).
- `round_away_shr` — round half **away from zero**. This is MATLAB `round()`, and
  it is what the model uses inside `bfpScale` and the reciprocal-sqrt range
  reduction. The block-float shift is a power of two, so exact ties occur
  constantly there.
- `trunc_shr` — truncate **toward zero**. This is MATLAB `fix()`, used inside the
  Newton iteration and the CORDIC. It differs from an arithmetic shift for
  negative operands, every iteration, and across 16 CORDIC iterations that
  compounds into a visibly wrong angle.

### Why nothing is rounded inside the covariance accumulator

This is correctness, not precision. A constant rounding bias *c* on every product
puts the same *c* into every entry of **R**. `ones(N)` is rank one with the
boresight steering vector as its eigenvector, so the estimator would invent a
source at zenith — exactly where the satellites are — and null it. Because the
error is a fixed matrix rather than a noise term it does **not** shrink with
dwell length, so longer integration never reveals it, and a floating-point
simulation never shows it at all.

Products are 32 bits and held exactly. Accumulating K = 16,368 adds
⌈log₂K⌉ = 14 bits; two products per real entry adds one more; 47 bits with a
guard is 48 — exactly the DSP48E1 P register.

### Whitened diagonal

The whitened diagonal lands at 2²⁶ to about **3 parts in 10⁵**, not exactly: Y
carries 17 bits and each of the two factors reproduces unity only to ~10⁻⁶.
Measured spread over the reference dwells is −1250…+1912 LSB. That is what keeps
`trace(Rw) ≈ N·2^F_EVD`, which is the bound the EVD's word-length proof rests on.

---

## 5. Block reference

Latency figures are deterministic unless stated; nothing in the datapath has a
data-dependent path except the reciprocal-sqrt normalisation loop (at most two
extra iterations) and the EVD's final sort.

| Block | Latency | Throughput | DSP48 | Notes |
|---|---|---|---|---|
| `asp_ddc_mixer` | 4 clk | 1 ch/clk | 4 | 16-entry NCO ROM, no phase accumulator |
| `asp_hb_decim2` | 5 clk | 1 ch/2 clk | **0** | all coefficients are power-of-two sums |
| `asp_fir_shape` | 12 clk | 1 ch/2 clk | 32 | symmetric, DSP48 pre-adder, 8-deep input FIFO |
| `asp_cov_accum` | 4 + 16 clk | 1 dwell/ms | 16 | 2:1 folded, P-register accumulate |
| `asp_rsqrt` | ~40 clk | 6/dwell | 2 | range reduction + 4 Newton steps |
| `asp_isqrt` | 26 clk | 4/dwell | 0 | restoring, exact by construction |
| `asp_cordic` | 19–21 clk | 4/rotation | 2 | shared vectoring/rotation |
| `asp_whiten` | ~270 clk | 1/dwell | 3 | |
| `asp_jacobi_evd` | ~4700 clk | 1/dwell | 8 | 36 rotations, fixed |
| `asp_detect` | 2 clk | 1/dwell | **0** | 5/4 threshold ⇒ shifts only |
| `asp_weight_calc` | ~350 clk | 1/dwell | 6 | |
| `asp_beamformer` | 4 clk | 1 sample/8 clk | 4 | |
| `asp_tx_scale` | 3 clk | 1 sample/8 clk | 2 | no sqrt, no log |

### Design points worth restating

**`asp_hb_decim2` uses zero DSP48.** Every non-zero coefficient is a sum of at
most two powers of two (−2¹², 2¹⁵+2¹², 2¹⁶). Elaboration assertions tie that
claim to the generated coefficient table, so a model re-tune breaks the build
instead of the filter. This is the strong reason it is custom RTL and not FIR
Compiler, which would spend multipliers on it.

**`asp_detect` uses zero DSP48 and no divider.** The threshold is 1280/1024 = 5/4
exactly, so the cross-multiplied test `lam1·3·1024 > 1280·tail` degenerates to
`lam1·3·4 > 5·tail`, and both 4× and 5× are shifts plus one add.

**`asp_jacobi_evd` sorts by insertion, not selection.** MATLAB's `sort()` is
stable and selection sort is not, so equal eigenvalues would return their
eigenvector columns in a different order than the model. Ties have probability
zero in real data, which is exactly why that would survive every test and then
differ on some captured dwell.

**`asp_weight_calc` rank 0 returns the quiescent beam untouched** — it does not
block-float it, which would double it (h = 2¹⁵ sits one bit below the top of an
18-bit word). The chain is invariant to a common real scale, so a doubled
quiescent beam is invisible to every measurement the array can make. Rank 0 is
the state the system is in for essentially all of its operating hours.

---

## 6. Required Vivado IP cores

Only infrastructure is IP. The algorithm is custom RTL because bit-exactness
against the golden model is the acceptance criterion, and a black-box
accumulation order cannot be diffed.

### 6.1 Clocking Wizard — **required**

```
IP NAME: clk_wiz
PRIMITIVE                   = MMCM
PRIM_IN_FREQ                = 245.760      (AD9361 DATA_CLK)
CLKOUT1_REQUESTED_OUT_FREQ  = 130.944      -> clk_dsp
CLKOUT2_REQUESTED_OUT_FREQ  = 32.736       (optional, ADC-rate enable)
USE_RESET = true, RESET_TYPE = ACTIVE_LOW
USE_LOCKED = true
Jitter filter: minimise output jitter
CONNECTIONS:
  clk_in1  <- AD9361 DATA_CLK via IBUFDS + BUFG
  clk_out1 -> asp_top.clk_dsp
  locked   -> ANDed into aresetn
```

**Check the M/D the wizard picks reports zero frequency error.** 245.760 /
130.944 is not an integer ratio. If the AD9361 is clocked from a 40 MHz
reference, derive 130.944 = 40 × 3.2736 instead. The one thing that must not
happen is `clk_dsp` being asynchronous to the ADC data — the entire enable-based
rate plan assumes it is not.

### 6.2 AXI Clock Converter — **required**

```
IP NAME: axi_clock_converter
PROTOCOL   = AXI4LITE
ADDR_WIDTH = 12
DATA_WIDTH = 32
ASYNC_CLK  = 1
CONNECTIONS:
  S_AXI  <- PS M_AXI_GP0 (through AXI Interconnect), s_axi_aclk = FCLK_CLK0
  M_AXI  -> asp_top.s_axi_*,                          m_axi_aclk = clk_dsp
```

This is the only asynchronous boundary in the design. Doing it here, once, is
why `asp_axi_lite_regs` contains no synchronisers.

### 6.3 AXI AD9361 — **required**, two instances

```
IP NAME: axi_ad9361 (Analog Devices)
MODE_1R1T = 0  (2R2T: each instance supplies TWO antennas)
DAC_DDS_DISABLE = 1
ADC_DATAPATH_DISABLE = 0
DELAY_REFCLK_FREQUENCY = 200
CONNECTIONS:
  device 0 -> rx_re[1:0], rx_im[1:0]
  device 1 -> rx_re[3:2], rx_im[3:2]
  both     -> one common rx_valid strobe at FS_ADC
```

**Both devices must share one LO and be MCS synchronised.** An independent LO per
pair puts an uncalibrated phase offset between antenna pairs that the array
cannot detect — it would appear as a fixed, plausible steering vector.

RX LO = **1577.466 MHz** (L1 + FS_ADC/16). TX LO = **1573.374 MHz**
(L1 − FS_ADC/16), so TX LO leakage lands 2 MHz off L1.

### 6.4 Zynq7 Processing System — **required**

```
IP NAME: processing_system7
FCLK_CLK0             = 100 MHz        (AXI-Lite only; NOT the datapath clock)
M_AXI_GP0             = enabled
S_AXI_HP0             = enabled        (only if the DMA capture path is built)
IRQ_F2P               = enabled, 1 bit  <- asp_top.irq
FCLK_RESET0_N         = enabled
```

### 6.5 AXI Interconnect — **required**

```
IP NAME: axi_interconnect
Number of Slave Interfaces  = 1   (PS M_AXI_GP0)
Number of Master Interfaces = 2   (asp control, axi_dma control)
```

### 6.6 AXI Direct Memory Access — **optional** (capture path only)

```
IP NAME: axi_dma
Enable Scatter Gather        = 0
Write (S2MM) only, Read Channel = 0
Width of Buffer Length Register = 26
Memory Map Data Width = 64, Stream Data Width = 32
Max Burst Size = 256
CONNECTIONS:
  S_AXIS_S2MM <- asp_top.m_axis_*
  M_AXI_S2MM  -> PS S_AXI_HP0
  s2mm_introut -> PS IRQF2P
```

Not required for the anti-spoofing function — it exists so captured beamformed
data can be logged, and as the hook for the future post-correlation stage.

### 6.7 FIR Compiler — **not used**, documented as an alternative

```
IP NAME: fir_compiler 7.2
Filter Type = Single_Rate, Number of Channels = 8
Input Sample Frequency = 16.368, Clock Frequency = 130.944
Coefficient Width = 18, Data Width = 16
Quantization = Integer_Coefficients
Output Rounding Mode = Convergent_Rounding_to_Even, Output Width = 16
Coefficient File = fir_shape.coe   (export FIR_COEF from asp_coef_pkg)
```

Rejected because bit-exactness against the golden model is the acceptance
criterion and the IP's internal accumulation order is not user-visible. It is a
legitimate substitute if `tb_asp_fir_shape` still passes with it — that test is
the arbiter, not this note.

---

## 7. Block design

```
                 ┌──────────────────────────┐
                 │  Zynq7 Processing System │
                 │  FCLK_CLK0 100 MHz       │
                 └──┬────────┬─────────┬────┘
             M_AXI_GP0   S_AXI_HP0   IRQ_F2P
                    │        │         │
          ┌─────────▼──┐     │         │
          │ AXI        │     │         │
          │ Interconn. │     │         │
          └──┬──────┬──┘     │         │
             │      │        │         │
   ┌─────────▼──┐ ┌─▼──────┐ │         │
   │ AXI Clock  │ │axi_dma │─┘         │
   │ Converter  │ │ (opt)  │◄──┐       │
   └─────┬──────┘ └────────┘   │       │
         │ clk_dsp             │ AXIS  │
   ┌─────▼───────────────────────────┐ │
   │            asp_top              ├─┘ irq
   │  ┌───────────────────────────┐  │
   │  │ asp_axi_lite_regs         │  │
   │  └───────────┬───────────────┘  │
   │  ┌───────────▼───────────────┐  │
   │  │ asp_datapath (stages 1-10)│  │
   │  └───┬───────────────────┬───┘  │
   └──────┼───────────────────┼──────┘
      rx_* │                  │ tx_*
   ┌───────▼──────┐    ┌──────▼───────┐
   │ axi_ad9361 0 │    │ axi_ad9361   │
   │ axi_ad9361 1 │    │ (TX path)    │
   └──────┬───────┘    └──────┬───────┘
          │  DATA_CLK         │
   ┌──────▼────────────┐      │
   │ Clocking Wizard   ├──────┘
   │ -> clk_dsp        │
   └───────────────────┘
```

Address map: assign `asp_top` a 4 kB window on `M_AXI_GP0` (e.g.
`0x43C0_0000`). The register map is in the header of `asp_axi_lite_regs.vhd`.

---

## 8. Timing and implementation

### Resource estimate (XC7Z020)

| Resource | Used | Available | % |
|---|---|---|---|
| DSP48E1 | ~85 | 220 | 39% |
| BRAM36 | **0** | 140 | 0% |
| FF | ~19,000 | 106,400 | 18% |
| LUT | ~13,000 | 53,200 | 24% |

There is no BRAM in the datapath at all, so placement never competes with a DMA
buffer or a PS cache path for block RAM columns.

### Reset strategy

One **synchronous, active-high** reset for the datapath, asserted asynchronously
and released synchronously through a two-flop synchroniser in `asp_top`.

Synchronous because Xilinx SRL primitives have no reset input: an asynchronous
reset would force the FIR and halfband delay lines into ~8000 flip-flops instead
of SRL16s. Datapath pipeline registers that self-flush are **not** reset — only
control state and valid pipelines are. That keeps the reset net small, which
matters more than it sounds: a global high-fanout reset is a routing and timing
problem, not a safety feature.

### Paths that were split for timing closure

Three places previously stacked a long combinational cone into one FSM
state. They are now spread across states (still **not** multicycle
exceptions — the FSMs advance every clock):

1. `asp_weight_calc / S_BF_PK` + `S_BF_PK_CH` — peak magnitude one channel
   (two abs + two compares) per clock instead of eight 48-bit magnitudes
   in one cycle.
2. `asp_whiten / S_NORM` + `S_NORM_SH` — max of four diagonals, then
   `ceil_log2` on the registered peak.
3. `asp_tx_scale / p_agc` — fast acquisition is a 5-step binary search
   (one constant-threshold compare per clock) plus hysteresis, finishing
   inside the 8-clock FS_WORK sample gap after the dwell tick.

### CDC inventory

| Crossing | Mechanism |
|---|---|
| PS AXI ↔ clk_dsp | AXI Clock Converter IP (the only one) |
| `aresetn` → clk_dsp | two-flop reset synchroniser, async assert / sync release |
| `irq` → PS | level held until the ISR clears `IRQ_STATUS`; PS IRQF2P inputs are level-sensitive and synchronised internally |

Everything else is one synchronous domain by construction.

---

## 9. Software requirements (PS)

Minimum to bring the system up (skeleton in `sw/`):

1. **Configure both AD9361 devices** over SPI: RX LO 1577.466 MHz, TX LO
   1573.374 MHz, RX bandwidth ~10 MHz, sample rate to give FS_ADC = 32.736 MHz.
   Run the full calibration sequence, then **freeze DC tracking**.
2. **MCS-synchronise the two devices** so all four antennas share a time base.
3. Write `SHIFT_INIT` (0x58) from a power measurement made during calibration, so
   dwell 0 does not clip.
4. Write `CONTROL` (0x04) bit 0 to enable.
5. Enable the dwell interrupt: `IRQ_ENABLE` (0x60) bit 0.

Per-dwell ISR (1 kHz) — see `sw/src/asp_driver.c` / `asp_dwell_policy`:

1. Read `STATUS`, `LAMBDA0..3`, `DET_LHS/RHS`.
2. Run the MDL test and the M-of-N hysteresis vote.
3. Set or clear `CONTROL` bit 2 (`rank2_en`) accordingly.
4. Clear `IRQ_STATUS` bit 0 (write 1).

A 1 kHz ISR reading ~15 registers is a few microseconds of Cortex-A9 time; it
does not need a real-time kernel, but it does need the interrupt not to be
starved by a userspace logger.

Vivado project bootstrap: `vivado/create_project.tcl` and
`vivado/bd_checklist.tcl`.

---

## 10. Assumptions, limitations, open items

**Assumptions**

- XC7Z020-1 or faster. The design targets 130.944 MHz; the three formerly
  long combinational cones (§8) are now multi-state.
- Two AD9361 devices sharing one LO, MCS-synchronised. Four antennas.
- `clk_dsp` is derived from the AD9361 sample clock. If it is not, the
  enable-based rate plan is invalid and the design needs FIFOs at every rate
  change.
- The AD9361 interface supplies four antennas in parallel with one common valid
  strobe at FS_ADC. Whether that comes from the ADI IP or a custom SelectIO front
  end does not matter to `asp_top`.

**Limitations**

- **Weight latency is two dwells, not one** (§3.3). Unavoidable given the
  estimator's latency; quantified above.
- **Dwell 0 of `asp_tx_scale` differs from the model** (§1). Unavoidable without
  a full-dwell buffer.
- `RANK2_ENABLE` defaults off. Rank-2 nulling is implemented and verified in the
  RTL but is gated by the PS, because a second arrival only 6 dB down gives an
  eigenvector with tens of degrees of error and nulling it measures *worse* than
  leaving it alone.
- No post-correlation stage. Pre-correlation spatial processing runs out of
  headroom below SAPR ≈ 0 dB; per-PRN clustering is the answer to
  multi-directional attacks and it is not built here.
- The design has not been through Vivado synthesis or place-and-route in this
  environment — no Xilinx toolchain is available here. Resource and timing
  figures are engineering estimates from the structures actually written, not
  tool reports. **Treat the first implementation run as the confirmation step.**
  Bootstrap with `vivado/create_project.tcl`.

**Open items**

- Confirm the MMCM can produce 130.944 MHz from the chosen AD9361 reference with
  zero frequency error (§6.1).
- Decide whether the capture DMA path is built for production or only for bring-up.
- The `set_clock_groups` constraint in `asp_timing.xdc` names `clk_fpga_0`; check
  that matches the PS clock name in the actual block design.

---

## 11. Per-module verification reference

Behaviour, ports, formats and latency are in each file's header. What follows is
the part that is easy to leave implicit: the **edge cases each module has to
survive** and **how each one is actually tested**.

Every module shares the same clock/reset contract: single `clk`, synchronous
active-high `rst`, asserted asynchronously and released synchronously by the
synchroniser in `asp_top`. Datapath pipeline registers are not reset; control
state and valid pipelines are.

| Module | Edge cases that matter | Verification |
|---|---|---|
| `asp_ddc_mixer` | NCO phase wrap at 16; all four TDM slots; a missing `s_valid` must repeat a sample rather than permute the antennas | `tb_asp_ddc_mixer`, 523,776 values; TDM channel order asserted every slot |
| `asp_hb_decim2` | decimation **phase** — keeping odd instead of even samples still looks correct in a spectrum; zero history at reset must match the model's zero-padded convolution | `tb_asp_hb_decim2`, 261,888 values |
| `asp_fir_shape` | the 4-on/4-off input burst from the decimator; FIFO overflow; the 31-sample warm-up discard | `tb_asp_fir_shape` drives the bursty pattern deliberately, asserts the overflow flag stayed low, and asserts a minimum comparison count so an early exit cannot look like a pass |
| `asp_cov_accum` | dwell boundary (last sample must not be lost); accumulator width at peak input; Hermitian structure | `tb_asp_cov_accum` checks 48-bit values exactly **and** asserts symmetry, antisymmetry, real diagonal and non-negative diagonal on every dwell |
| `asp_rsqrt` | powers of two and their neighbours, where the range reduction changes *k*; `a = 1`, where it shifts **left**; `a = 0`; the top of the range (2⁴⁰) that the weight stage produces | `tb_asp_rsqrt`, 371 cases, boundary values included explicitly rather than relying on random draws |
| `asp_isqrt` | `a = 0`; perfect squares and the values either side | `tb_asp_rsqrt`, 377 cases |
| `asp_cordic` | all four quadrants and both axes (the left-half-plane pre-rotation is where hand-written CORDICs fail); the ±π/2 fold in rotation mode; magnitudes to 2²⁹ | `tb_asp_cordic`, 544 vectoring + 408 rotation cases |
| `asp_whiten` | all-zero dwell / dead antenna (the clamp exists so it saturates rather than wraps); negative off-diagonals, which round half **away from zero** | `tb_asp_whiten`, plus assertions that the diagonal is exactly real and within 2¹² LSB of 2²⁶ |
| `asp_jacobi_evd` | equal eigenvalues (sort stability); trace preservation, which fails if a rotation is applied on one side only | `tb_asp_jacobi_evd` checks eigenvectors **and** eigenvalues, asserts descending order and trace |
| `asp_detect` | statistic exactly at the threshold; rank 2 gated off | `tb_asp_detect` compares both operands of the cross-multiplied test, not just the flag — the flag alone would pass with a threshold wrong by 2× |
| `asp_weight_calc` | **rank 0**, the normal operating state; a degenerate (numerically zero) eigenvector column; block-float saturation at the `−W` clamp | `tb_asp_weight_calc` checks the weights exactly and separately asserts that rank 0 returns the quiescent beam untouched |
| `asp_beamformer` | conjugate on the **weights** (getting it backwards passes a broadside test and fails everywhere else); accumulator clear at channel 0 | `tb_asp_beamformer`, 130,944 values, driven with the causal weight sequence |
| `asp_tx_scale` | clipping; AGC hysteresis band edges; dwell boundary power accumulation | `tb_asp_tx_scale`, 98,208 values from dwell 2 onward, with dwell 1 excluded by design and the reason stated in the testbench |
| `asp_datapath` | the dwell handshake between blocks; the two-dwell weight schedule; telemetry only valid from the second tick | `tb_asp_datapath` runs several million clocks from raw ADC samples and checks the weight sequence bit-exactly — the strongest end-to-end check available, since the weights depend on every stage from the mixer through the EVD |
| `asp_axi_lite_regs` | independent AW/W arrival order; write-1-to-clear racing a dwell tick | not co-simulated; recommend the Xilinx AXI VIP in AXI4-Lite protocol-checker mode, or a directed read/write testbench against the register map in the file header |
| `asp_top` | reset release; `enable` gating | recommend bring-up on hardware with the ILA described below |

### Recommended hardware bring-up

The RTL is verified against the model in simulation; what simulation cannot check
is the AD9361 interface and the clocking. In order:

1. **Clock first.** Confirm `clk_dsp` is exactly 130.944 MHz and MMCM `locked` is
   stable. Everything downstream assumes the integer rate plan.
2. **ILA on the mixer input**, triggered on `rx_valid`. Confirm four antennas
   arrive in the expected slot order and that the data is not all zeros or all
   ones — a mis-wired AD9361 interface usually shows as one antenna stuck.
3. **Read `DWELL_CNT`** (0x0C) and confirm it increments at 1 kHz. If it does
   not, the sample rate is wrong, not the algorithm.
4. **Point the array at a clean sky** and confirm `STATUS.detected` stays low.
   This is the single most valuable hardware test, and it is the one that fails
   silently: a false alarm produces no symptom except an authentic satellite
   being nulled.
5. **Inject a spoofer** (a second signal generator through a splitter with a
   phase offset) and confirm `STATUS.detected` goes high and the weights in
   0x30–0x4C move away from the quiescent value.
6. **Check `CLIP_CNT`** (0x54) stays at zero once `SHIFT_INIT` is set correctly.
