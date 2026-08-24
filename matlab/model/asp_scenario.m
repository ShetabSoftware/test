function scn = asp_scenario(cfg, varargin)
%ASP_SCENARIO Build a repeatable spoofing scenario.
%
%   scn = ASP_SCENARIO(cfg)
%   scn = ASP_SCENARIO(cfg, 'spoofEnabled', false, 'seed', 7, ...)
%
%   Returns a struct describing every emitter in the scene plus the array
%   and channel realisation.  Nothing downstream generates random numbers,
%   so a scenario fully determines a run.
%
%   scn.src is a struct array with one entry per emitter:
%     .prn   PRN identifier
%     .kind  'auth' | 'spoof' | 'spoofmp'
%     .amp   complex baseband amplitude, referred to unit noise variance
%     .az    azimuth [deg]
%     .el    elevation [deg]
%     .tau0  code delay [s]
%     .fd    Doppler [Hz]
%     .phi0  initial carrier phase [rad]
%     .bitEpoch  offset of the 20 ms navigation bit grid [s]
%     .bits  navigation data bits (+/-1) covering the run
%
%   MODELLING CHOICES THAT DIFFER FROM THE ORIGINAL SCRIPTS
%   -------------------------------------------------------
%   1. Satellite directions are drawn uniformly in SOLID ANGLE above the
%      elevation mask, not uniformly in the elevation ANGLE.  Uniform-in-
%      elevation (as in the original scripts, rand*90) over-populates zenith
%      by a factor of 1/cos(el) and therefore understates how often an
%      authentic satellite sits near a low-elevation spoofer, which is
%      exactly the case the projection null damages.
%
%   2. The spoofer is placed at LOW elevation by default (15 deg).  A
%      terrestrial spoofer is a terrestrial object.  Putting it at 45 deg
%      elevation, as the paper's simulation does, is optimistic in two ways:
%      it separates the spoofer from the low-elevation satellites and it
%      removes the ground-bounce geometry entirely.
%
%   3. Navigation data bits are present, at 50 bps with an independent bit
%      grid offset per emitter.  Without them, any estimator that correlates
%      across a 20 ms boundary appears to work when in hardware it does not.
%
%   4. The spoofer optionally has a specular ground reflection, making the
%      spoofing subspace rank 2.

opt.spoofEnabled = true;
opt.seed         = cfg.seed;
opt.authAzEl     = [];      % nAuth x 2, overrides the random draw
opt.durationMs   = 20;

for k = 1:2:numel(varargin)
    if ~isfield(opt, varargin{k})
        error('asp_scenario:opt', 'Unknown option "%s".', varargin{k});
    end
    opt.(varargin{k}) = varargin{k+1};
end

rng(opt.seed);

% ------------------------------------------------------------------ array
[pos, geo] = asp_array_geometry(cfg.geometry, cfg.lambda, cfg.elementSpacingLambda);
nAnt = size(pos,2);
if nAnt ~= cfg.nAnt
    % The geometry name is authoritative; keep the config consistent.
    cfg.nAnt = nAnt;
end

% ---------------------------------------------------------- channel state
% Per-channel amplitude/phase mismatch (post-LNA: cables, filters, mixer,
% transceiver channel gain).  Because it is post-LNA it scales the thermal
% noise identically, which is the property that makes the covariance
% diagonal a valid per-channel amplitude reference (see asp_estimate_ssv).
gainLin  = 10.^(cfg.fe.gainMismatchDB/20 * (2*rand(1,nAnt)-1));
phaseRad = deg2rad(cfg.fe.phaseMismatchDeg * (2*rand(1,nAnt)-1));
Cdiag    = gainLin .* exp(1i*phaseRad);

% Per-channel group-delay mismatch (cable length, filter group delay).  This
% is the impairment that actually limits null depth in hardware and that the
% original scripts do not model at all.
tauChan = cfg.fe.delayMismatchSec * (2*rand(1,nAnt)-1);

% Mutual coupling.  Modelled as a symmetric Toeplitz-like matrix whose
% off-diagonal magnitude decays with element separation.  Crucially, the
% projection algorithm is INDIFFERENT to this: it never uses the array
% manifold, so an arbitrary invertible linear mixing of the channels leaves
% "the direction the spoofer arrives from" a well-defined vector in the
% measured space.  Coupling is modelled so that this claim can be TESTED
% rather than asserted.
M = eye(nAnt);
if cfg.fe.mutualCoupling ~= 0
    dmat = zeros(nAnt);
    for i = 1:nAnt
        for j = 1:nAnt
            dmat(i,j) = norm(pos(:,i)-pos(:,j))/cfg.lambda;
        end
    end
    M = cfg.fe.mutualCoupling .^ (dmat/max(dmat(~eye(nAnt))));
    M(logical(eye(nAnt))) = 1;
    M = M .* exp(1i*pi*dmat);      % coupling is complex
    M(logical(eye(nAnt))) = 1;
end

% ------------------------------------------------------------- directions
if isempty(opt.authAzEl)
    az = rand(1,cfg.nAuth)*360;
    % Uniform in solid angle above the mask: sin(el) uniform in
    % [sin(mask), 1].
    sMin = sind(cfg.elevationMaskDeg);
    el   = asind(sMin + (1-sMin)*rand(1,cfg.nAuth));
else
    az = opt.authAzEl(:,1).';
    el = opt.authAzEl(:,2).';
end

% ----------------------------------------------------------------- powers
% amplitude^2 / noise variance = (C/N0)/fs, with noise variance normalised
% to 1.
authAmp  = sqrt(cfg.authSnrSample);
spoofAmp = sqrt(cfg.spoofSnrSample);

nBits = ceil(opt.durationMs/20) + 2;

src = struct('prn',{},'kind',{},'amp',{},'az',{},'el',{}, ...
             'tau0',{},'fd',{},'phi0',{},'bitEpoch',{},'bits',{});

for m = 1:cfg.nAuth
    src(end+1) = mkSrc(cfg.authPrn(min(m,numel(cfg.authPrn))), 'auth', ...
        authAmp, az(m), el(m), rand*cfg.Tcode, (rand*2-1)*5000, ...
        rand*2*pi, rand*0.02, nBits); %#ok<AGROW>
end

if opt.spoofEnabled
    % A real spoofer synthesises a self-consistent constellation from a
    % single clock, so all its PRNs share one carrier phase reference and
    % one bit grid.  Modelling them with independent phases would make the
    % spatial energy add incoherently and understate the threat.
    spoofPhi  = rand*2*pi;
    spoofBitEpoch = rand*0.02;
    spoofBits = 2*randi([0 1],1,nBits)-1;

    for k = 1:cfg.nSpoof
        s = mkSrc(cfg.spoofPrn(min(k,numel(cfg.spoofPrn))), 'spoof', ...
            spoofAmp, cfg.spoofAzEl(1), cfg.spoofAzEl(2), ...
            rand*cfg.Tcode, (rand*2-1)*5000, spoofPhi, spoofBitEpoch, nBits);
        s.bits = spoofBits;
        src(end+1) = s; %#ok<AGROW>
    end

    if cfg.spoofMultipath.enable
        mpAmp = spoofAmp * 10^(cfg.spoofMultipath.relDB/20);
        for k = 1:cfg.nSpoof
            base = src(cfg.nAuth + k);
            s = base;
            s.kind = 'spoofmp';
            s.amp  = mpAmp;
            s.az   = cfg.spoofMultipath.azEl(1);
            s.el   = cfg.spoofMultipath.azEl(2);
            s.tau0 = base.tau0 + cfg.spoofMultipath.delaySec;
            % A specular ground bounce inverts the phase for the horizontal
            % polarisation component and adds the extra path length.
            s.phi0 = base.phi0 + pi + 2*pi*cfg.fc*cfg.spoofMultipath.delaySec;
            src(end+1) = s; %#ok<AGROW>
        end
    end
end

% ------------------------------------------------------------- assemble
scn.cfg      = cfg;
scn.pos      = pos;
scn.geo      = geo;
scn.nAnt     = nAnt;
scn.Cdiag    = Cdiag;
scn.tauChan  = tauChan;
scn.M        = M;
scn.src      = src;
scn.seed     = opt.seed;
scn.durationMs = opt.durationMs;
scn.spoofEnabled = opt.spoofEnabled;

% True spatial signature vectors, for scoring only.  The algorithm never
% sees these.
scn.aTrue = zeros(nAnt, numel(src));
for s = 1:numel(src)
    scn.aTrue(:,s) = M * (Cdiag(:) .* asp_steering(src(s).az, src(s).el, pos, cfg.lambda, cfg.c));
end
idxSpoof = find(strcmp({src.kind},'spoof'), 1);
if ~isempty(idxSpoof)
    scn.bTrue = scn.aTrue(:, idxSpoof);
else
    scn.bTrue = [];
end
scn.kinds = {src.kind};

end

% -------------------------------------------------------------------------
function s = mkSrc(prn, kind, amp, az, el, tau0, fd, phi0, bitEpoch, nBits)
s.prn      = prn;
s.kind     = kind;
s.amp      = amp;
s.az       = az;
s.el       = el;
s.tau0     = tau0;
s.fd       = fd;
s.phi0     = phi0;
s.bitEpoch = bitEpoch;
s.bits     = 2*randi([0 1],1,nBits)-1;
end
