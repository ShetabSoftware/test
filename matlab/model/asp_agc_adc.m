function [xq, agc] = asp_agc_adc(x, cfg, agc)
%ASP_AGC_ADC Automatic gain control followed by ADC quantisation.
%
%   [XQ, AGC] = ASP_AGC_ADC(X, CFG, AGC)
%
%   Holds the COMPOSITE (thermal + interference) per-component RMS at
%   CFG.fx.agcBackoffDB below full scale, then quantises to CFG.fx.adc.
%   AGC state is carried between calls so the loop has realistic memory.
%
%   TWO HARDWARE RULES ENCODED HERE
%   -------------------------------
%   1. ONE gain for all channels.  The AD9361 (and every other multichannel
%      transceiver) will happily run an independent AGC per channel.  If you
%      let it, each channel's gain index moves at its own time, and every
%      move is both an amplitude step and - because the gain is realised by
%      switching attenuator and amplifier stages - a PHASE step of several
%      degrees.  The array's spatial signature then changes mid-dwell, the
%      covariance estimate is a mixture of two different arrays, and the
%      null collapses.  Set cfg.fe.commonAGC = false to watch this happen.
%
%   2. Gain updates are quantised to the hardware's gain step and are rate
%      limited.  A gain change invalidates the covariance dwell that
%      straddles it; the pipeline flags that so the estimate can be held
%      rather than corrupted.
%
%   AGC.gain          current linear gain
%   AGC.changed       true if the gain moved during this block
%   AGC.satFraction   fraction of quantiser outputs that saturated

if nargin < 3 || isempty(agc)
    agc.gain      = 1;
    agc.stepDB    = 0.5;      % AD9361 gain table granularity, approx.
    agc.maxRateDB = 3;        % max change per block
    agc.changed   = false;
end

fmt      = cfg.fx.adc;
targetRms = fmt.maxval * 10^(-cfg.fx.agcBackoffDB/20);

% Per-component RMS: I and Q each carry half the complex power.
if cfg.fe.commonAGC
    curRms = sqrt(mean(abs(x(:)).^2)/2);
    desired = targetRms / max(curRms, eps);
    agc.gain = applyGainLimits(agc.gain, desired, agc);
    gains = repmat(agc.gain, size(x,1), 1);
    agc.changed = abs(20*log10(desired/agc.gain)) > agc.stepDB;
else
    curRms = sqrt(mean(abs(x).^2, 2)/2);
    desired = targetRms ./ max(curRms, eps);
    if isscalar(agc.gain)
        agc.gain = repmat(agc.gain, size(x,1), 1);
    end
    for i = 1:size(x,1)
        agc.gain(i) = applyGainLimits(agc.gain(i), desired(i), agc);
    end
    gains = agc.gain;
    agc.changed = true;
end

xg = bsxfun(@times, gains(:), x);

[xq, info] = fx_quant(xg, fmt);
agc.satFraction = info.satFraction;
agc.appliedGain = gains;

end

% -------------------------------------------------------------------------
function g = applyGainLimits(g, desired, agc)
deltaDB = 20*log10(desired/g);
deltaDB = max(min(deltaDB, agc.maxRateDB), -agc.maxRateDB);
% Quantise to the hardware gain step.
deltaDB = round(deltaDB/agc.stepDB)*agc.stepDB;
g = g * 10^(deltaDB/20);
end
