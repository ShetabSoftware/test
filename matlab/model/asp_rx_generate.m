function [x, state] = asp_rx_generate(scn, blockIdx, blockLen, state)
%ASP_RX_GENERATE Generate one block of received antenna-array samples.
%
%   [X, STATE] = ASP_RX_GENERATE(SCN, BLOCKIDX, BLOCKLEN, STATE)
%
%   X is nAnt x BLOCKLEN complex baseband samples at scn.cfg.fs, with the
%   thermal noise variance normalised to 1 per antenna.  BLOCKIDX is
%   0-based; blocks are contiguous.
%
%   Block-at-a-time generation, rather than one giant matrix, is deliberate:
%   it bounds memory for long dwells, it matches how the FPGA consumes data,
%   and it lets the pipeline apply CAUSAL weights (weights derived from
%   block n-1 applied to block n), which the original scripts do not do.
%
%   SIGNAL MODEL
%   ------------
%   For emitter s at antenna i,
%
%     r_i(t) = A_s * a_i(s) * d_s(t - T) * c_s(t - T) * exp(j(2 pi f_s t + phi_s))
%     T      = tau0_s + tauGeom_i(s) + tauChan_i
%
%   with a_i(s) = exp(j 2 pi (u_s . p_i)/lambda) and the CODE evaluated at
%   the per-antenna delayed time.  Carrying the geometric delay into the
%   code argument is what makes this a wideband-correct model; the paper and
%   the original scripts keep only the carrier phase term.  The difference
%   is small (0.5 ns across a 27 cm aperture) but it is precisely the
%   physical effect that sets the achievable null floor, so a reference
%   model that omits it cannot be used to justify a null-depth
%   specification.
%
%   Code Doppler is included: the chipping rate is scaled by (1 + fd/fc).
%   Over 1 ms this moves the code by only 1.6e-3 chips, but over a 1 s
%   dwell it is 1.6 chips, i.e. total decorrelation.  Reference models that
%   omit it silently mislead anyone who later increases the dwell.

cfg  = scn.cfg;
nAnt = scn.nAnt;

if nargin < 4 || isempty(state)
    rng(scn.seed + 991);          % noise stream, disjoint from the scenario stream
    state.noiseInit = true;
end

n0 = blockIdx * blockLen;
t  = (n0 + (0:blockLen-1)) / cfg.fs;          % 1 x blockLen

sig = zeros(nAnt, blockLen);

for s = 1:numel(scn.src)
    src = scn.src(s);

    [aVec, tauGeom] = asp_steering(src.az, src.el, scn.pos, cfg.lambda, cfg.c);

    code = asp_ca_code(src.prn);

    % Total per-antenna delay: source code delay + geometric + channel skew.
    tauTot = src.tau0 + tauGeom(:) + scn.tauChan(:);        % nAnt x 1

    % Code phase in chips, per antenna, including code Doppler.
    tRel  = bsxfun(@minus, t, tauTot);                      % nAnt x blockLen
    chips = tRel * cfg.chipRate * (1 + src.fd/cfg.fc);
    idx   = mod(floor(chips), 1023) + 1;
    cSamp = code(idx);
    if blockLen == 1
        cSamp = reshape(cSamp, nAnt, 1);
    end

    % Navigation data at 50 bps.  Uses the reference-element timing; the
    % 0.5 ns of array skew is 2.5e-8 of a bit and cannot straddle an edge in
    % any way that matters.
    bitIdx = floor((t - src.tau0 - src.bitEpoch)/0.02) + 2;
    bitIdx = min(max(bitIdx, 1), numel(src.bits));
    dSamp  = src.bits(bitIdx);                              % 1 x blockLen

    carrier = exp(1i*(2*pi*src.fd*t + src.phi0));           % 1 x blockLen

    common = src.amp * (dSamp .* carrier);                  % 1 x blockLen
    sig = sig + bsxfun(@times, aVec(:), bsxfun(@times, cSamp, common));
end

% ----------------------------------------------------- mutual coupling
% Coupling acts on the antenna terminals, before the LNA.
sig = scn.M * sig;

% ----------------------------------------------------- thermal noise
% Injected at the LNA input, spatially white, unit variance per antenna.
eta = (randn(nAnt, blockLen) + 1i*randn(nAnt, blockLen)) / sqrt(2);
x   = sig + eta;

% --------------------------------- post-LNA amplitude / phase mismatch
% Cables, filters and transceiver channel gain scale signal AND noise
% identically.  This is why the covariance diagonal is a valid per-channel
% amplitude reference (asp_estimate_ssv), and it is a physical statement
% about where the mismatch lives, not a modelling convenience.
x = bsxfun(@times, scn.Cdiag(:), x);

% ----------------------------------------------------- LO phase noise
if cfg.fe.loPhaseNoiseRmsDeg > 0
    sigmaRad = deg2rad(cfg.fe.loPhaseNoiseRmsDeg);
    if cfg.fe.sharedLO
        % One synthesiser feeding every mixer: the phase error is common
        % mode.  A common complex scalar multiplies every channel equally,
        % leaves the spatial covariance structure untouched, and therefore
        % cannot degrade null depth at all.  This is the single strongest
        % argument for distributing one LO rather than running independent
        % PLLs per transceiver.
        ph = sigmaRad * randn(1, blockLen);
        x  = bsxfun(@times, exp(1i*ph), x);
    else
        ph = sigmaRad * randn(nAnt, blockLen);
        x  = x .* exp(1i*ph);
    end
end

% ----------------------------------------------------- IQ imbalance
if cfg.fe.iqImbalanceDB ~= 0 || cfg.fe.iqImbalanceDeg ~= 0
    g   = 10^(cfg.fe.iqImbalanceDB/20);
    phi = deg2rad(cfg.fe.iqImbalanceDeg);
    alpha = (1 + g*exp(-1i*phi))/2;
    beta  = (1 - g*exp( 1i*phi))/2;
    x = alpha*x + beta*conj(x);
end

% ----------------------------------------------------- DC offset
if isfinite(cfg.fe.dcOffsetDBFS)
    dc = 10^(cfg.fe.dcOffsetDBFS/20) * exp(1i*2*pi*(0:nAnt-1).'/nAnt);
    x  = bsxfun(@plus, x, dc);
end

end
