function fx = asp_fx_plan(cfg)
%ASP_FX_PLAN Derive the fixed-point word-length plan from the system config.
%
%   fx = ASP_FX_PLAN(cfg)
%
%   Every word length below is DERIVED, not guessed.  The derivations are
%   reproduced in the comments and are re-checked at run time by
%   verify/test_fx_plan.m.
%
%   Target device assumption: Xilinx 7-series / Zynq-7000 (DSP48E1, 25x18
%   signed multiplier, 48-bit accumulator, 25-bit pre-adder) or UltraScale
%   (DSP48E2, 27x18).  Intel equivalents (18x19 in Cyclone V, 27x27 in
%   Arria/Stratix variable-precision mode) also accept these widths.
%
%   ------------------------------------------------------------------
%   STAGE 1: ADC / RF interface
%   ------------------------------------------------------------------
%   AD9361 delivers 12-bit signed I and Q.  Nothing downstream can create
%   information that is not in those 12 bits, so the ADC word length caps
%   the achievable interference-to-noise headroom.  See
%   analysis/asp_adc_jn_limit.m for the derivation of
%
%       SNR loss after nulling = 10*log10(1 + (q^2/12)/sigma_thermal^2)
%
%   which gives ~55-60 dB of usable J/N for a 12-bit converter operating at
%   14 dB of AGC backoff, versus ~80 dB for 16-bit.
%
%   ------------------------------------------------------------------
%   STAGE 2: digital front end (equaliser + decimator)
%   ------------------------------------------------------------------
%   16 bits.  This is 12 ADC bits plus the ~2 bits of processing gain from
%   decimation plus 2 bits of guard.  16 also happens to sit inside the
%   DSP48E1 25-bit A port with room for the sign extension needed by the
%   pre-adder, and inside the Intel 18x19 multiplier.
%
%   ------------------------------------------------------------------
%   STAGE 3: covariance accumulator
%   ------------------------------------------------------------------
%   product width  = 2*16                                    = 32 bits
%   coherent growth= ceil(log2(K))  with K = samples per dwell
%   guard          = 2 bits
%   For K = 16368 that is 32 + 14 + 2 = 48 bits, which is EXACTLY the DSP48
%   P register.  This is not a coincidence worth ignoring: it means one
%   DSP48 per covariance entry with the accumulator in the slice, no fabric
%   adder, no rounding anywhere inside the loop, and therefore no
%   accumulated rounding bias (see FX_FMT for why that matters).
%
%   ------------------------------------------------------------------
%   STAGE 4: eigen-decomposition
%   ------------------------------------------------------------------
%   The Jacobi method applies only unitary similarity transforms, so the
%   Frobenius norm of the working matrix is invariant and there is no
%   dynamic range growth at all.  This is the reason to prefer it over
%   anything based on explicit matrix inversion.  The word length is set by
%   the eigenvalue spread that must be resolved:
%
%       lambda_1/lambda_N ~ 1 + N*J/N_thermal
%
%   For 60 dB of J/N that is 10^6, i.e. 20 bits of range, plus ~6 bits to
%   keep the smallest eigenvalue meaningful.  26 bits minimum; 32 chosen.
%   The matrix is trace-normalised into block floating point first, so the
%   exponent is carried separately and does not consume mantissa.
%
%   ------------------------------------------------------------------
%   STAGE 5: beamformer weights
%   ------------------------------------------------------------------
%   Weight quantisation sets a hard floor on null depth.  With weights
%   perturbed by e (uniform, variance q^2/12 per real component, q = 2^-FLw)
%   the residual gain toward the nulled direction is
%
%       E|e^H b|^2 / ||f||^2 = N*(q^2/6) / ||f||^2
%
%   With block-floating-point scaling ||f||^2 ~ N/4, so the floor is
%   (2/3)*q^2, i.e.
%
%       null floor [dB] = 20*log10(q) - 1.76
%
%   FL = 12 -> -75 dB, FL = 16 -> -100 dB.  Anything at or above 12
%   fractional bits is far below the ~25 dB the estimator itself can
%   deliver, so 16 is chosen purely to land on the DSP48 18-bit B port with
%   a sign bit and one guard bit.  The original scripts used Q1.22 in a
%   24-bit word, which buys nothing and wastes a whole DSP cascade on
%   Intel parts.

nAntBits = ceil(log2(max(cfg.nAnt, 2)));

if isfield(cfg, 'K') && ~isempty(cfg.K)
    K = cfg.K;
else
    K = round(cfg.fs * cfg.Tcode);
end
if isfield(cfg, 'est') && isfield(cfg.est, 'coherentMs')
    K = K * cfg.est.coherentMs;
end

% --- Stage 1
fx.adc = fx_fmt(12, 11);                       % Q0.11, AD9361 native

% --- Stage 2
fx.sample = fx_fmt(16, 15);                    % Q0.15 datapath sample
fx.eqCoef = fx_fmt(18, 16);                    % channel equaliser taps

% --- Stage 3
fx.covProductBits = 2 * fx.sample.wl;                       % 32
fx.covGrowthBits  = ceil(log2(K));                          % 14 at K=16368
fx.covGuardBits   = 2;
fx.covAccWl       = fx.covProductBits + fx.covGrowthBits + fx.covGuardBits;

% --- Stage 4
fx.evd     = fx_fmt(32, 30);                   % trace-normalised, BFP exponent carried apart
fx.evdRot  = fx_fmt(18, 16);                   % CORDIC rotation sin/cos
fx.eigVal  = fx_fmt(32, 30);

% --- Stage 5
fx.weight  = fx_fmt(18, 16);                   % Q1.16, DSP48 B port

% --- Stage 6: beamformer
fx.beamProductBits = fx.sample.wl + fx.weight.wl;           % 34
fx.beamAccWl       = fx.beamProductBits + nAntBits + 1;     % 34 + 2 + 1 = 37
fx.beamOut         = fx_fmt(16, 13);           % back to 16 bits for the correlator

% --- Derived predictions, checked by verify/test_fx_plan.m
fx.predictedWeightNullFloorDB = 20*log10(2^(-fx.weight.fl)) - 1.76;
fx.predictedAdcJnLimitDB      = asp_adc_jn_limit(fx.adc.wl, 14, 1.0);
fx.predictedSampleJnLimitDB   = asp_adc_jn_limit(fx.sample.wl, 14, 1.0);

% --- AGC operating point
% Composite (thermal + interference) RMS is held at this backoff below full
% scale.  14 dB gives 5 sigma of crest headroom for a Gaussian composite,
% i.e. a per-component clip probability near 6e-7.  Lower backoff clips on
% interference peaks; higher backoff throws away J/N headroom one-for-one.
fx.agcBackoffDB = 14;

end
