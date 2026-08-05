# Task 4 — RF Front-End Architecture

## 0. A clarification first, because it changes the answer

Your Option B is described as:

> Low-frequency DAC → analog reconstruction filtering → mixer/upconverter → LO → RF chain up to GPS L1 → custom RF front end

That is a **transmit** chain. DAC, reconstruction filter and upconverter to L1 describe a signal *generator* — which, given that one of your files is named `GPS_SPOOFER_V9`, may well be deliberate: you need a spoofer to test against.

An anti-spoofing *receiver* needs the reverse: antenna → LNA → filter → downconverter → IF filter → ADC. I answer both, because you need both and they are different products with different requirements:

* **§1–§7** compare AD9361 versus a discrete **receive** front end. This is the anti-spoofing product.
* **§8** covers the transmit chain as a **test asset**, where the requirements are genuinely different and mostly stricter.

---

## 1. What this application actually demands of a front end

Ranked by how much they constrain the design. Note how little of this list is about conventional RF figures of merit.

| # | Requirement | Why it dominates |
|---|---|---|
| 1 | **Inter-channel phase and amplitude stability over temperature and time** | The null depth *is* the channel matching. Everything else is secondary. |
| 2 | **Common LO** | With independent LOs, differential phase noise caps the null near −35 to −41 dB. With a shared LO the phase noise is common mode and contributes **exactly zero**. |
| 3 | **Deterministic inter-channel sample alignment** | One sample of skew at 16 MHz is 61 ns — fatal. |
| 4 | **Common AGC** | Independent per-channel gain steps change the array signature mid-dwell. |
| 5 | **Converter word length** | Sets usable J/N: 57 dB at 12 bits, 81 dB at 16 bits. |
| 6 | **Flat inter-channel group delay across the band** | 1 ns of skew caps the null at −35 dB. |
| 7 | **Image rejection** | An interferer's I/Q image has a *different* spatial signature and cannot be nulled. |
| 8 | Noise figure | Set by the external LNA; the transceiver's own NF is irrelevant. |
| 9 | Absolute phase noise | Common mode with a shared LO. Matters for carrier tracking, not for nulling. |

**The striking thing about this list is that items 1–4 and 6 are all about *matching*, not about performance.** A mediocre four-channel receiver with excellent matching will out-null an excellent four-channel receiver with poor matching, every time. This is the single most important framing for the architecture decision, and it is why "which chip has the better datasheet" is the wrong question.

---

## 2. Option A — AD9361

### What it does well here

* **Two RX per chip, so four antennas is exactly two chips with nothing wasted.** (Three antennas is also two chips, with one path idle — see Task 1.)
* External LO input on both chips, driven from one synthesiser. ADI's own FMCOMMS5 reference design (2× AD9361, shared LO, MCS) exists precisely to do this, so the hard part is de-risked and the errata are known.
* MCS aligns the baseband dividers and digital clocks across chips.
* Enormous flexibility: 70 MHz–6 GHz, up to 56 MHz bandwidth. If you later add L2/L5/E5 or want to survey the interference environment out of band, the same hardware does it.
* Mature ADI HDL/Linux stack (`axi_ad9361`, IIO), which is worth months of schedule.
* Excellent prototype-to-product continuity: the same silicon in the lab and in the field.

### What it does badly here, in order of severity

**12-bit ADCs cap usable J/N at ~57 dB, ~51 dB at realistic ENOB.** Fine for spoofing (SAPR 0–30 dB). Marginal the moment a jammer is present — and real attacks are frequently jam-then-spoof, because forcing reacquisition is how you get a receiver to accept counterfeit signals quickly.

**Per-channel AGC, QEC and DC-offset tracking calibrations run independently and asynchronously.** Every one of them changes a channel's amplitude/phase response at a time the array processor does not control. All must be forced into manual/frozen operation, which means giving up the automation that is a large part of the AD9361's appeal, and then re-implementing gain control yourself with a per-gain-index phase-correction table.

**Zero-IF.** DC offset lands on the C/A spectral peak, and the I/Q image of a strong interferer appears at the mirror frequency with a *different* spatial signature (because QEC differs per channel) — so it **cannot be nulled**. Image rejection therefore sets a hard floor on interference suppression. ~65 dB after calibration is fine; degraded to 40 dB by temperature drift with frozen QEC, it becomes the limit.

**Inter-chip baseband filter mismatch.** The two chips' analog and digital filters are calibrated independently. Amplitude and group-delay mismatch across the band of a few tenths of a dB and a few hundred picoseconds is normal, which caps the null in the −35 to −45 dB region without digital equalisation.

**Power.** ~700 mW–1.2 W per AD9361 in 2×2 receive. Two chips is 1.5–2.5 W before the FPGA. For a vehicle installation this is nothing; for a handheld or a UAV it is a real constraint.

**Cost.** Roughly $80–130 each in volume, so $160–260 for four channels, plus a substantial support burden: multiple LDO rails, a large BOM of passives, and a non-trivial PCB.

### Verdict on Option A

**AD9361 is the right choice for the prototype and for the first shipping product, and this is not a compromise.** The reason is specific: the algorithm's *calibration-free* property is precisely what makes two independently-locked transceivers viable. Their inter-chip phase offset is constant after a given sync event but changes on every re-sync — and the projection method does not care, because it re-estimates the spatial signature every millisecond from the data. A manifold-based algorithm (MVDR toward a known direction, MUSIC, AOA discrimination) would need per-boot calibration to use two AD9361s at all.

That is a genuine architectural synergy between the paper's algorithm and this silicon, and it is worth stating plainly because it is not obvious.

---

## 3. Option B — discrete GNSS receive front end

The credible version is not "a mixer and an LO". It is:

```
  4× [ patch → LNA (0.8 dB NF, 30 dB) → SAW (L1, ±20 MHz) → LNA ]
        │
        ├─► 4× downconverter IC in ANALOG IF OUTPUT mode
        │      (MAX2771 or equivalent, external LO input)
        │        ▲
        │        └── one LMX2594 / ADF4351 synthesiser
        │            → matched-length 1:4 splitter (Wilkinson, on PCB)
        │
        └─► 4× IF SAW / LC filter
              │
              └─► ONE quad ADC: AD9653 (16-bit, 125 MSPS, 4 channels)
                     │  single chip ⇒ one clock, one sampling instant,
                     │  inherently sample-aligned by construction
                     └─► LVDS → FPGA
```

### Advantages, and they are substantial

| | |
|---|---|
| **16-bit converters** | 81 dB usable J/N versus 57 dB. **This is the single biggest technical difference between the options.** |
| **One quad ADC** | All four channels share one clock and one sampling instant. Inter-channel sample alignment is guaranteed by construction — no MCS, no boot-time verification, no failure mode. |
| **Genuinely shared LO** | Not an external-LO mode with internal dividers whose start phase varies; one synthesiser through a symmetric on-PCB splitter. |
| **Low-IF or real IF** | No zero-IF DC offset on the spectral peak. Image rejection becomes a filter problem, not a matched-quadrature problem, which is far more stable over temperature. |
| **No hidden calibrations** | Nothing runs an autonomous cal loop behind your back. Everything is deterministic. |
| **GNSS-optimised filtering** | SAWs sized for L1, not general-purpose ±28 MHz. Better out-of-band rejection means better resilience to LTE/harmonic interference. |
| **Power** | ~100–150 mW per channel versus ~350–600 mW. Roughly 3× better. |
| **Cost at volume** | ~$15–25 per channel versus ~$45–65. |

**Trap to avoid:** many GNSS front-end ICs (NT1065, MAX2769 in default mode) output 2- or 3-bit quantised samples. **Those parts are unusable for nulling.** Two-bit quantisation caps interference suppression around 15–20 dB regardless of what the array does. You must use the analog-IF output mode and your own 14/16-bit ADC. This is the most common way a "GNSS-optimised front end" turns out to be a dead end.

### Disadvantages

* Substantially more RF NRE: LO distribution, IF filtering, layout, shielding, EMC.
* A synthesiser, splitter and clock tree you now own and must qualify over temperature.
* No vendor reference design; ADI's HDL/Linux stack does not apply.
* Fixed to the bands you designed for.
* Longer schedule to first light, and the failure modes are analog ones that are harder to debug than register writes.
* Four separate downconverter ICs means four separate device-to-device variations — you have traded *inter-chip* mismatch for *inter-device* mismatch. It is smaller and more stable, but it does not vanish, and you still need the digital channel equaliser.

---

## 4. Head-to-head

| Criterion | A: 2× AD9361 | B: discrete + quad ADC | Winner |
|---|---|---|---|
| Usable J/N (converter limited) | 51–57 dB | **81 dB** | **B, by 24 dB** |
| Noise figure (system) | 1.2 dB | 1.2 dB | tie (external LNA sets it) |
| Phase noise impact on nulling | zero with shared LO | zero | tie |
| Phase noise impact on carrier tracking | good | **better** (dedicated synth) | B, marginally |
| Inter-channel phase stability | good after equalisation | **better** | B |
| Inter-channel sample alignment | MCS + boot verification | **guaranteed by one quad ADC** | **B, removes a failure mode** |
| Image rejection | 55–70 dB, QEC dependent, drifts | **filter-set, stable** | B |
| Spurious | internal PLL/divider spurs to manage | fewer, but LO splitter leakage | tie |
| Multi-channel coherence | good | **excellent** | B |
| Autonomous calibration interference | **significant** (AGC/QEC/DC) | none | **B** |
| Band flexibility | **70 MHz–6 GHz** | fixed | **A, decisively** |
| BOM, 4 ch, volume | $160–260 | **$60–100** | B |
| PCB complexity | moderate (2 chips, many rails) | **high** (RF layout, LO tree) | A |
| Manufacturability / test | **easier** (digital cal) | more RF test steps | A |
| Software/HDL maturity | **mature ADI stack** | you own it | **A, worth months** |
| Schedule to first light | **weeks** | months | **A** |
| Engineering risk | **low** | moderate-high | **A** |
| Scalability to 8 channels | 4 chips | 2 quad ADCs + 8 RF chains | tie |
| Power, 4 ch | 1.5–2.5 W | **0.4–0.6 W** | B |

---

## 5. Recommendation

**Phase 1 — prototype and first product: Option A.** Two AD9361s with a **shared external LO** (mandatory, not optional), MCS, forced manual gain with a common gain index, frozen tracking calibrations, and a digital per-channel channel equaliser. This gets you to a demonstrable, sellable anti-spoofing product with low RF risk and a mature software stack. It is fully adequate against spoofing, which is the product's stated purpose.

**Phase 2 — volume product: Option B**, once the algorithm, the FPGA architecture and the calibration procedure are frozen and you know what the field actually throws at you. The 24 dB of extra J/N, the 3× power reduction, the 2–3× BOM reduction and the elimination of the sample-alignment failure mode are all real — but none of them is worth taking analog risk *before* the digital architecture is settled.

**Consider a third option seriously.** If power and cost allow, a single **ADRV9026** (4 RX in one package, one shared internal LO, 16-bit converters) collapses almost every problem above into one chip: no inter-chip anything, no MCS, no LO distribution, 16 bits. It is expensive (~$300–400) and power-hungry (~4 W), so it is not a handheld part — but for a vehicle, base-station or infrastructure-timing product it is arguably the *correct* answer and it deserves an evaluation before you commit to two AD9361s. It also scales to 8 elements as two chips.

### Non-negotiables, whichever option you choose

1. **One LO for all channels.** Independent PLLs cap the null at −35 to −41 dB through differential phase noise alone. With a shared LO this term is exactly zero.
2. **One clock, and verified sample alignment at every boot.** Inject a common pilot through a splitter, cross-correlate, log the result, re-sync automatically if non-zero.
3. **Common AGC, or manual gain.** Never let per-channel AGC run.
4. **Per-gain-index phase-correction table** from production calibration, and invalidate any dwell containing a gain change.
5. **Digital per-channel equaliser** (16-tap complex FIR) with per-unit production calibration. This is what converts "matched to 1 ns" into "matched to 50 ps".
6. **≥ 14-bit converters** if anti-jam is on the roadmap at all.
7. **Freeze autonomous tracking calibrations** during array processing dwells.

---

## 6. RF chain details worth specifying now

| Item | Recommendation | Reason |
|---|---|---|
| Antenna | Ceramic patch, RHCP, ~25 mm, sequential feed rotation across the array | Axial ratio consistency matters more than peak gain in a CRPA |
| Ground plane | Continuous, ≥ 1.5λ diameter if the mechanics allow | Sets the low-elevation pattern and the front-to-back ratio |
| LNA | 0.6–0.8 dB NF, 30 dB gain, integrated at the element | NF is set here; nothing downstream can recover it |
| Pre-filter | SAW immediately after the first LNA stage | Protects against out-of-band overload before any gain stage can be driven into compression |
| Cabling | **Phase-matched, same reel, same length, same routing, same bend radius** | 1 ns of differential delay caps the null at −35 dB. Match to < 100 ps |
| LO distribution | On-PCB Wilkinson, symmetric layout, equal trace lengths | Any LO phase imbalance is a fixed offset, absorbed by the algorithm — but keep it stable over temperature |
| Reference | TCXO ±0.5 ppm minimum; OCXO if timing is a product feature | |
| Bias tees | Identical parts, same lot, per channel | Common-mode gain differences are absorbed; *drift* differences are not |

**Note on cable matching.** The algorithm absorbs an arbitrary *fixed* linear channel mismatch, including mutual coupling — so a constant phase offset between channels costs nothing. What it does not absorb is **frequency-dependent** mismatch (group-delay skew) and **time-varying** mismatch (thermal drift within a dwell). Specify and test those, not absolute phase match.

---

## 7. Thermal and drift

The array signature is re-estimated every millisecond, so slow thermal drift is tracked automatically and is a non-issue. What matters instead:

* **Drift *within* a dwell.** 1 ms. Nothing thermal moves that fast. Non-issue.
* **Differential thermal gradients across the array.** If one channel's chain sits next to the FPGA and another does not, they drift differently over minutes. Still tracked by the re-estimation, but it makes the *production calibration table* stale. Mitigation: symmetric mechanical layout, a temperature sensor per channel, and a calibration table indexed by temperature.
* **Gain-index phase steps** are the fast, discontinuous drift that actually hurts. See §5.4.

---

## 8. The transmit chain, as a test asset

If Option B was meant as a signal generator, the requirements are different and mostly stricter.

You need a controlled spoofer to test against, and buying one is not straightforward: transmitting GNSS signals over the air is illegal essentially everywhere. The paper's own test setup (its Figure 8) is the correct pattern and worth reproducing exactly:

> a rooftop antenna array receives the authentic signals; a hardware simulator generates the spoofing constellation; the simulator output is split N ways and **conductively combined** with each antenna's signal.

This is exactly right, because it reproduces the property the algorithm exploits — all spoofing PRNs arriving through **one** propagation channel, hence one spatial signature — while the authentic signals arrive from genuinely different directions. And it does it inside a cable, so nothing radiates.

Requirements for the generator, in priority order:

1. **A single, self-consistent constellation from one clock.** All PRNs must share one carrier phase reference. Independent per-PRN phases would make the spatial energy add incoherently and understate the threat — which is precisely the property the whole method depends on.
2. **Calibrated, repeatable output power** with a programmable attenuator, so SAPR can be swept in 1 dB steps. Every performance curve in `matlab/studies/` is a function of SAPR; without a calibrated sweep you cannot reproduce them.
3. **Phase-matched N-way split.** The splitter defines the "spoofer direction" in the test. Its phase imbalance *is* the spoofing spatial signature — which is fine and controllable, and lets you place the synthetic spoofer at an arbitrary apparent direction by inserting known phase offsets. This is a genuinely useful capability: you can sweep the spoofer's apparent direction across the sky without moving anything.
4. Realistic scenarios: matched-power, over-powered, and the hard case — a **drag-off** attack that starts aligned with the authentic signals and walks away slowly.

For the generator itself: a commercial simulator (Spirent, Orolia/Safran, Syntony) is the right answer if the budget exists, because a calibrated reference removes an entire category of "is it the test or the DUT" ambiguity. If you build it, a second AD9361 in transmit mode is far simpler than a discrete DAC-plus-upconverter chain and is more than good enough for a spoofer — spoofer signal quality requirements are much looser than receiver requirements, because you are deliberately generating an impairment.

**One caution.** A test setup built entirely on conductive combining will never exercise the effects that dominate in the field: mutual coupling under a real radome, the antenna's actual embedded element patterns, ground-bounce geometry, and platform scattering. Budget for an anechoic-chamber campaign with a rotating positioner before you claim any specification. Conducted testing verifies the algorithm; only chamber testing verifies the *product*.
