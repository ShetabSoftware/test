function cfg = asp_config(varargin)
%ASP_CONFIG Single source of truth for every system parameter.
%
%   cfg = ASP_CONFIG()
%   cfg = ASP_CONFIG('fs', 4.092e6, 'nAnt', 3, ...)
%
%   Every module in this reference model takes its numbers from here.  The
%   original scripts hard-coded fs, K, N_ant, powers, seeds and word lengths
%   independently in three files, which is how they ended up disagreeing
%   with each other and with the paper.
%
%   POWER PARAMETERISATION
%   ----------------------
%   Powers are specified as C/N0 in dB-Hz, NOT as absolute dBW and NOT as
%   arbitrary amplitudes.  This is the only parameterisation that stays
%   correct when the sample rate changes, because the pre-despreading sample
%   SNR is
%
%       SNR_sample = (C/N0) / fs                                        (1)
%
%   for complex baseband sampling at rate fs with a front end that passes
%   the full fs of noise bandwidth.  At C/N0 = 45 dB-Hz and fs = 16.368 MHz
%   this is -27.1 dB, i.e. every GPS signal sits 27 dB BELOW the thermal
%   noise in each sample.  Getting this number right is not a detail: the
%   entire behaviour of the spatial covariance estimator is governed by it.
%   (The original GPS_SPOOFER_V9 script placed the signals 60 dB ABOVE the
%   noise, an operating-point error of 87 dB, which makes every reported
%   null depth in that script meaningless.)

% ---------------------------------------------------------------- RF plan
cfg.fc      = 1575.42e6;      % GPS L1 centre frequency [Hz]
cfg.c       = 299792458;      % speed of light [m/s]
cfg.lambda  = cfg.c / cfg.fc; % 0.190293 m

% Sample rate: 16 x 1.023 MHz gives an integer number of samples per C/A
% chip (16) and per code period (16368), which removes the code NCO phase
% accumulator ambiguity from the reference model and makes the RTL golden
% vectors trivially comparable.  4.092e6 (4 samples/chip) is the low-rate
% option used by the Monte Carlo studies for speed.
cfg.fs      = 16 * 1.023e6;

cfg.Tcode   = 1e-3;           % C/A code period [s]
cfg.chipRate = 1.023e6;

% ---------------------------------------------------------------- Array
cfg.nAnt        = 4;
% 'y4' = three elements on a ring at 120 deg plus one at the centre.
%
% This beats the 2x2 square for the same four elements and the same radome,
% and the reason is a lattice property rather than anything about the
% element count.  The baseline set of the Y generates a TRIANGULAR lattice,
% which aliases only when the spacing reaches 2*lambda/sqrt(3) = 1.155*(lambda/2),
% whereas the square lattice of a 2x2 array aliases at lambda/2 exactly.
% Inside a 0.5 lambda radome and holding the grating response at -3 dB, that
% buys a maximum baseline of 0.82 lambda instead of 0.53 lambda, a 37%
% narrower spatial null, and it cuts the fraction of satellites left worse
% off than a single antenna from 12.8% to 6.9% (studies/study_array_size.m).
% The cost is that the centre element has three near neighbours while the
% outer elements have one, so their embedded patterns and mutual coupling
% differ - which matters for manifold calibration but not for the
% projection algorithm, since it never uses the manifold.
cfg.geometry    = 'y4';
cfg.elementSpacingLambda = 0.45;   % ring radius in wavelengths
cfg.refElement  = 1;

% ---------------------------------------------------------------- Signals
cfg.authCN0dBHz   = 45;       % per authentic PRN, paper value
cfg.nAuth         = 9;
cfg.authPrn       = [2 5 10 12 21 25 29 30 31];
cfg.saprDB        = 5.5;      % spoof-to-authentic power ratio, per PRN
cfg.nSpoof        = 9;
cfg.spoofPrn      = [2 5 10 12 21 25 29 30 31];
cfg.spoofAzEl     = [45 15];  % [azimuth elevation] deg.  A terrestrial
                              % spoofer is LOW, not at 45 deg elevation as
                              % in the paper's simulation.
cfg.elevationMaskDeg = 5;

% Second coherent arrival from the spoofer (a reflection off a building,
% vehicle or structure).  This is what makes the spoofing subspace rank 2.
%
% IMPORTANT AND NON-OBVIOUS: a SPECULAR GROUND BOUNCE IS NOT SUCH A CASE.
% For a planar array with every element at z = 0, the response depends only
% on the direction cosines (cos(el)cos(az), cos(el)sin(az)), and cos is
% even, so elevation +theta and -theta produce an IDENTICAL steering vector.
% Measured coherence between (az 45, el +15) and (az 45, el -15) is
% 1.000000.  A specular ground reflection arrives at exactly the mirror
% elevation and the same azimuth, so it is AUTOMATICALLY CO-NULLED by the
% same rank-one projector that removes the direct path.  No extra degree of
% freedom is needed and no wideband penalty appears, because the two
% arrivals share one steering vector at every frequency.
%
% This is a real and favourable property of the planar geometry, and it
% cuts the other way too: the array cannot spatially separate a satellite
% at +theta from its own ground reflection either, so a planar CRPA gives
% you no multipath rejection against ground bounce of the AUTHENTIC
% signals.  Rejecting that needs a vertical baseline (a non-planar array)
% or the receiver's own code/carrier multipath mitigation.
%
% The rank-2 cases that DO consume a degree of freedom are arrivals at a
% different AZIMUTH: a reflector off to one side, a second spoofer, or an
% accompanying jammer.  The default below is a building reflection at
% azimuth 135 deg, which is essentially orthogonal to the direct path
% (measured coherence 0.010).
cfg.spoofMultipath.enable   = false;
cfg.spoofMultipath.relDB    = -6;
cfg.spoofMultipath.delaySec = 60e-9;      % ~18 m of excess path
cfg.spoofMultipath.azEl     = [135 20];

% ---------------------------------------------------------------- Front end
cfg.fe.gainMismatchDB     = 0.5;   % per-channel amplitude mismatch, +/- dB
cfg.fe.phaseMismatchDeg   = 5;     % per-channel phase mismatch, +/- deg
cfg.fe.delayMismatchSec   = 100e-12;  % per-channel group delay mismatch, +/- s
cfg.fe.mutualCoupling     = 0;     % nearest-neighbour coupling coefficient (linear)
cfg.fe.sharedLO           = true;  % true: common LO, phase noise is common mode
cfg.fe.loPhaseNoiseRmsDeg = 0.4;   % per-channel RMS residual phase noise
cfg.fe.iqImbalanceDB      = 0;     % amplitude imbalance
cfg.fe.iqImbalanceDeg     = 0;     % quadrature error
cfg.fe.dcOffsetDBFS       = -Inf;  % per-channel DC offset
cfg.fe.commonAGC          = true;  % MANDATORY in hardware; false models the
                                   % failure mode of per-channel AGC

% ---------------------------------------------------------------- Estimator
cfg.est.coherentMs      = 1;     % covariance accumulation length [ms]
cfg.est.jacobiSweeps    = 6;     % FIXED sweep count => deterministic latency
cfg.est.maxNullRank     = 2;     % max number of spatial nulls to place
cfg.est.diagonalLoadDB  = -20;   % MVDR diagonal loading, relative to trace/N
cfg.est.mode            = 'evd'; % 'evd' | 'paper' | 'gamma' | 'column'
cfg.est.pfa             = 1e-3;  % target false-alarm probability
% cfg.est.detectThreshold is DERIVED below, not hard-coded.  It depends on
% N, on the dwell length AND on fs (because the authentic per-sample SNR is
% (C/N0)/fs), so any constant written here would be silently wrong the
% moment one of those changed.  See analysis/asp_detect_threshold.m.
cfg.est.detectThreshold = [];

% ---------------------------------------------------------------- Numerics
cfg.fx = asp_fx_plan(cfg);

% ---------------------------------------------------------------- Repeatability
cfg.seed = 20120926;   % ION GNSS 2012 presentation date, for luck

% ---------------------------------------------------------------- Overrides
for k = 1:2:numel(varargin)
    key = varargin{k};
    val = varargin{k+1};
    cfg = asp_setfield_dotted(cfg, key, val);
end

% ---------------------------------------------------------------- Derived
cfg.K = round(cfg.fs * cfg.Tcode);
if abs(cfg.K - cfg.fs*cfg.Tcode) > 1e-9
    error('asp_config:fs', ...
        'fs*Tcode = %.6f is not an integer; choose fs as a multiple of 1 kHz.', ...
        cfg.fs*cfg.Tcode);
end
cfg.samplesPerChip = cfg.fs / cfg.chipRate;

% Per-sample SNR of one authentic PRN, equation (1) above.
cfg.authSnrSample = 10^(cfg.authCN0dBHz/10) / cfg.fs;
cfg.spoofSnrSample = cfg.authSnrSample * 10^(cfg.saprDB/10);

% Re-derive the fixed-point plan if the user overrode anything it depends on.
cfg.fx = asp_fx_plan(cfg);

% Calibrate the detector threshold against a realistic H0 (a clean sky with
% the authentic constellation present), unless the caller supplied one.
if isempty(cfg.est.detectThreshold)
    [cfg.est.detectThreshold, cfg.est.detectThresholdInfo] = ...
        asp_detect_threshold(cfg, cfg.est.pfa, 200);
end

end

% -------------------------------------------------------------------------
function s = asp_setfield_dotted(s, key, val)
% Support 'fe.gainMismatchDB' style keys without needing setfield chains.
parts = strsplit(key, '.');
switch numel(parts)
    case 1
        s.(parts{1}) = val;
    case 2
        s.(parts{1}).(parts{2}) = val;
    case 3
        s.(parts{1}).(parts{2}).(parts{3}) = val;
    otherwise
        error('asp_config:key', 'Key "%s" is nested too deeply.', key);
end
end
