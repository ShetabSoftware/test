# Using the Code

## 1. Getting it

The code is on a **branch**, not on `main`. That is why it does not appear in your default checkout.

| | |
|---|---|
| Repository | `github.com/ShetabSoftware/test` |
| Branch | `cursor/gnss-antispoof-array-redesign-4a0a` |
| Pull request | #3 |

```bash
git fetch origin
git checkout cursor/gnss-antispoof-array-redesign-4a0a
ls matlab/
```

If you would rather not use git, the PR page has a **"Files changed"** tab, and
`Code → Download ZIP` on the branch page gives you the whole tree.

To merge it into `main` once you are happy with it:

```bash
git checkout main
git merge cursor/gnss-antispoof-array-redesign-4a0a
```

---

## 2. Requirements

**None beyond MATLAB or Octave.** No toolboxes — not Fixed-Point Designer, not
Signal Processing, not Phased Array System Toolbox. This is deliberate: an
algorithm reference that an FPGA team has to license three toolboxes to run is
an algorithm reference nobody runs.

* MATLAB R2018b or newer, **or**
* GNU Octave 7 or newer (`sudo apt install octave` on Debian/Ubuntu)

Everything in this repository was developed and verified under Octave 8.4, and
the code avoids `RandStream`, `contains`, `validateattributes` and `arguments`
blocks so that it behaves identically on both.

---

## 3. First run

```matlab
cd <repo>/matlab
run('asp_startup.m')     % puts config/ fx/ model/ core/ analysis/ verify/ studies/ export/ on the path

asp_run_all('tests')     % regression suite only          ~30 s
asp_run_all()            % + short studies                ~5 min
asp_run_all('full')      % + full Monte Carlo studies     ~40 min
```

From a shell, without opening the GUI:

```bash
cd matlab
octave-cli --no-gui -q --eval "asp_startup; asp_run_all('tests')"
```

`asp_run_all('tests')` should end with:

```
 ALL REGRESSION TESTS PASSED   (26.8 s)
```

If it does not, stop and read the failure — every assertion is a claim made
somewhere in the documents, so a failure means a document is now wrong.

---

## 4. What to read, in what order

| Read this | To understand |
|---|---|
| `config/asp_config.m` | every system parameter, in one place, with the reasoning |
| `core/asp_process.m` | the whole pipeline end to end — start here for the algorithm |
| `core/asp_ssv_from_cov.m` | the estimator, and why it supersedes the paper's |
| `core/asp_detect.m` | the missing block: detection and rank |
| `verify/test_fixedpoint.m` | the numerical trap that costs you 25 dB if you miss it |
| `studies/study_array_size.m` | the evidence behind the 3-vs-4 decision |

The comment blocks are the design rationale. They are longer than usual on
purpose: each one states what the code does, what the alternative was, and what
was measured to choose between them. If you disagree with a decision, the
comment tells you which measurement to re-run.

---

## 5. Directory map

```
matlab/
  asp_startup.m     path setup — run this first, every session
  asp_run_all.m     regression suite + studies

  config/           SINGLE SOURCE OF TRUTH
    asp_config.m      all system parameters; powers as C/N0, not amplitudes
    asp_fx_plan.m     word lengths, DERIVED with the derivations shown

  fx/               fixed-point primitives (replaces Fixed-Point Designer)
    fx_fmt            format descriptor (word length, fraction, rounding)
    fx_quant          quantise with saturation reporting
    fx_cov_accum      BIT-EXACT DSP48 covariance accumulator model
    fx_bfp_scale      power-of-two normalisation (replaces every norm())
    fx_round_convergent

  model/            SIMULATION ONLY — none of this goes in the FPGA
    asp_ca_code       ICD-verified GPS C/A Gold codes
    asp_array_geometry  tri3 sq4 circ4 y4 circ7 circ8 + lattice metadata
    asp_steering      array response AND the matching propagation delay
    asp_scenario      emitters, channel mismatch, coupling, multipath
    asp_rx_generate   block-streaming waveform generator
    asp_agc_adc       AGC + ADC quantisation + saturation counting

  core/             THE ALGORITHM REFERENCE — maps 1:1 onto RTL blocks
    asp_process       causal streaming pipeline (the top level)
    asp_ssv_from_cov  spoofing signature: 'column' | 'gamma' | 'evd'
    asp_ssv_paper     the 2012 estimator, kept as the comparison baseline
    asp_evd_herm      fixed-sweep cyclic Jacobi Hermitian eigen-decomposition
    asp_detect        eigenvalue detection + MDL rank estimation
    asp_weights       'project' | 'mvdr' | 'lcmv' | 'mrc'
    asp_beamform      streaming beamformer, float or bit-exact fixed point
    asp_postcorr_ssv  stage-2 refinement from despread snapshots

  analysis/         metrics and closed-form theory
  verify/           regression tests — run these in CI
  studies/          the Monte Carlo evidence behind every number in docs/
  export/           bit-exact RTL co-simulation vectors
```

---

## 6. Recipes

Every one of these has been run exactly as written.

### Run one scenario and look at the result

```matlab
cfg = asp_config();
scn = asp_scenario(cfg, 'durationMs', 10);
out = asp_process(scn, 'durationMs', 10, 'verbose', true);
```

You will see one line per millisecond:

```
blk   1  det=1 rank=1 stat=1.9398  null= -25.31 dB  authMean= -0.42 dB  rho=0.9986
```

`det` is the detector decision, `rank` the number of nulls placed, `stat` the
eigenvalue statistic, `null` the gain toward the spoofer in dB relative to one
antenna element, `rho` the SSV estimate accuracy.

### Change the array

```matlab
cfg = asp_config('geometry','tri3', 'nAnt',3, 'elementSpacingLambda',0.45);   % the paper's array
cfg = asp_config('geometry','y4',   'nAnt',4, 'elementSpacingLambda',0.45);   % recommended
cfg = asp_config('geometry','circ8','nAnt',8, 'elementSpacingLambda',0.50);   % 8-element SKU
```

Available: `tri3`, `sq4`, `circ4`, `y4`, `circ7`, `circ8`. Check any layout for
grating ambiguity before committing to it:

```matlab
[pos, meta] = asp_array_geometry('y4', cfg.lambda, 0.45);
meta.maxBaselineLambda      % 0.779
meta.gratingLevelDB         % -3.59  (below -3 dB is safe)
meta.nullWidthU             % 0.995  (above ~1 means the null covers the whole sky)
```

### Sweep spoofing power

```matlab
for sapr = [0 10 20]
    cfg = asp_config('fs', 4*1.023e6, 'saprDB', sapr);
    scn = asp_scenario(cfg, 'seed', 1, 'durationMs', 3);
    out = asp_process(scn, 'durationMs', 2);
    fprintf('SAPR %4.1f dB -> null %7.2f dB\n', sapr, out.block(end).metrics.nullGainDB);
end
```

```
SAPR  0.0 dB -> null  -14.65 dB
SAPR 10.0 dB -> null  -31.08 dB
SAPR 20.0 dB -> null  -44.34 dB
```

The null tracks the threat — deeper attacks get deeper nulls — which is the
behaviour you want and a genuine strength of the method.

### Compare estimators, including the paper's

```matlab
cfg = asp_config('fs', 4*1.023e6);
scn = asp_scenario(cfg, 'seed', 7, 'durationMs', 4);
x   = asp_rx_generate(scn, 0, cfg.K*3, []);
R   = (x*x')/size(x,2);
h   = ones(cfg.nAnt,1)/sqrt(cfg.nAnt);

for m = {'column','gamma','evd'}
    y = asp_ssv_from_cov(R, m{1}, struct('rank',1,'jacobiSweeps',6));
    f = asp_weights('project', y, h);
    fprintf('%-8s rho=%.4f  null=%7.2f dB\n', m{1}, ...
        asp_ssv_correlation(y, scn.bTrue), ...
        10*log10(abs(f'*scn.bTrue)^2/real(f'*f)));
end

y = asp_ssv_paper(x, cfg.K, 1);          % the 2012 method
```

**Caution:** a single seed is one realisation and the spread between dwells is
several dB. Do not draw conclusions from one run — use `study_estimators(40)`,
which averages over scenes, before believing any ordering.

### Float versus fixed point

```matlab
o1 = asp_process(scn, 'durationMs', 3);
o2 = asp_process(scn, 'durationMs', 3, 'fixedPoint', true);
```

With no interference present these agree closely, because the AGC keeps thermal
noise at ~400 LSB and quantisation noise sits 60 dB below it. The word length
only starts to matter when a jammer forces the AGC to back off — see
`analysis/asp_adc_jn_limit.m`.

### Deliberately break things (this is the useful part)

```matlab
% Per-channel AGC instead of a common gain index — watch the null collapse
cfg = asp_config('fe.commonAGC', false);

% Uncalibrated inter-channel group delay
cfg = asp_config('fe.delayMismatchSec', 1e-9);

% Independent LOs instead of one shared synthesiser
cfg = asp_config('fe.sharedLO', false, 'fe.loPhaseNoiseRmsDeg', 1.0);

% Mutual coupling — demonstrates the algorithm is manifold-agnostic
cfg = asp_config('fe.mutualCoupling', 0.25);

% No spoofer at all — the case that must NOT trigger nulling
scn = asp_scenario(cfg, 'spoofEnabled', false);
```

These are how you justify a hardware specification. `fe.delayMismatchSec` in
particular converts a cable-matching tolerance into a null-depth number.

### Generate RTL co-simulation vectors

```matlab
asp_export_vectors('/path/to/vectors', asp_config(), 8);
```

Writes ADC stimulus, expected covariance accumulator contents (**raw
integers**), eigenvalues, weights and beamformer output, plus a `manifest.txt`
recording the formats and the config that produced them.

Compare **bit-exact**, not with a tolerance. The failure this design is most
exposed to — a rounding bias in the covariance accumulator — is a few LSBs that
any sensible tolerance passes, while capping null depth at −25 dB.

### Run the studies individually

```matlab
study_array_size(1500)   % 3 vs 4 vs 7 vs 8 elements       ~2 min
study_estimators(40)     % paper vs gamma vs evd           ~5 min
study_detector(600)      % detection ROC vs SAPR           ~5 min
study_postcorr(20)       % two-stage architecture          ~4 min
study_multipath(14)      % ground bounce vs reflectors    ~20 min
```

The first argument is the trial count. Reduce it for a quick look; the numbers
quoted in the documents used the values shown.

---

## 7. Using your own recorded data

Replace the generator, keep everything else. `asp_process` consumes an
`nAnt × nSamples` complex matrix, so:

```matlab
cfg = asp_config('fs', <your sample rate>, 'nAnt', 4);
K   = cfg.K;

% x is nAnt x nSamples, complex baseband, one row per antenna
R = (x(:,1:K) * x(:,1:K)') / K;

[y, dbg] = asp_ssv_from_cov(R, 'evd', ...
    struct('rank', cfg.est.maxNullRank, 'jacobiSweeps', 6));
det = asp_detect(dbg.lam, K, cfg.est.detectThreshold, cfg.est.maxNullRank);

if det.rank > 0
    f = asp_weights('project', y(:,1:det.rank), ones(cfg.nAnt,1)/sqrt(cfg.nAnt));
else
    f = ones(cfg.nAnt,1)/sqrt(cfg.nAnt);      % no threat: pass the quiescent beam
end

v = asp_beamform(x, f);                        % your spoofing-suppressed stream
```

Two things to get right:

**Scaling.** The model normalises thermal noise to unit variance per antenna.
If your data is raw ADC integers, either scale it or override
`cfg.fx.agcBackoffDB`. The estimator itself is scale invariant, but the
detector threshold is not, so recalibrate:

```matlab
thr = asp_detect_threshold(cfg, 1e-3);
```

**Threshold validity.** `cfg.est.detectThreshold` is calibrated for the `cfg`
that produced it. It depends on N, on dwell length **and on fs**, because the
authentic per-sample SNR is (C/N₀)/fs. Change any of those and recalibrate — a
threshold carried across configurations is a latent bug that only appears later.

---

## 8. Gotchas

| | |
|---|---|
| **Run `asp_startup` first**, every session. Otherwise you get "undefined function". |
| **`asp_config()` takes ~1 s** because it calibrates the detector threshold by Monte Carlo. Build it once and reuse it in loops. |
| **The default `fs` is 16.368 MHz**, which is realistic but slow to simulate. Studies use `asp_config('fs', 4*1.023e6)` for speed; the operating point is preserved because power is specified as C/N₀. |
| **Single-seed results are noisy** — several dB of spread between dwells. Use the studies for any comparison. |
| **`scn.bTrue` and `scn.aTrue` are ground truth for scoring only.** No function in `core/` ever reads them. If you add a feature, keep it that way. |
| **Override syntax is dotted**: `asp_config('fe.sharedLO', false)`, not a nested struct. |
| **Changing `geometry` also changes `nAnt`** — the geometry name is authoritative. Pass both to be explicit. |

---

## 9. Where to start, given what you want to do

**"I want to check the 3-vs-4 decision myself."**
`study_array_size(1500)`, then read `docs/01-array-architecture.md` alongside it.

**"I want to start the RTL."**
Read `docs/03-fpga-implementation.md` for the block diagram and word lengths,
`core/asp_process.m` for the dataflow, then `asp_export_vectors` for the golden
vectors. Put `verify/test_fixedpoint.m` in CI on day one.

**"I want to justify a hardware specification."**
Sweep the `cfg.fe.*` impairments (§6) and read the null depth out. That turns
cable-matching and LO-sharing decisions into dB.

**"I want to see what was wrong with my original scripts."**
`docs/02-matlab-review-and-redesign.md`, findings ordered by severity, with the
measurement backing each one.

**"I want the one-page version for a colleague."**
`docs/00-executive-summary.md`.
