function m = asp_cov_model(p)
%ASP_COV_MODEL Fast covariance-domain model of the array scene.
%
%   M = ASP_COV_MODEL(P) builds the true spatial covariance of a scene and
%   draws a sample covariance from it, without generating any waveform.
%
%   P fields (all optional, defaults shown):
%     .geometry  'circ4'      array layout, see asp_array_geometry
%     .spacing   0.5          spacing / ring radius in wavelengths
%     .lambda    0.190293     wavelength [m]
%     .nAuth     9            number of authentic satellites
%     .authSnr   1.93e-3      per-satellite power / noise variance
%     .nSpoof    9            number of spoofing PRNs
%     .saprDB    5.5          spoof-to-authentic power ratio per PRN
%     .spoofAzEl [45 15]      spoofer direction [az el] deg
%     .mpRelDB   -Inf         ground-bounce relative power (-Inf = none)
%     .mpAzEl    [45 -12]     ground-bounce direction
%     .jamSnr    0            broadband jammer power / noise variance
%     .jamAzEl   [200 8]      jammer direction
%     .K         16368        snapshots in the covariance dwell
%     .elMask    5            elevation mask [deg]
%     .Cdiag     []           channel mismatch (default: drawn)
%     .gainMismatchDB    0.5
%     .phaseMismatchDeg  5
%     .draw      true         draw a sample covariance (false: return truth)
%
%   M.Rtrue, M.Rhat, M.b, M.A (authentic responses), M.bmp, M.jam, M.pos
%
%   WHY A COVARIANCE-DOMAIN MODEL IS LEGITIMATE HERE
%   ------------------------------------------------
%   Over a dwell of K >~ 10^3 samples, the sample covariance of a sum of
%   independent wideband sources converges to a complex Wishart matrix with
%   mean R and K degrees of freedom, regardless of whether the sources are
%   BPSK or Gaussian: the codes are pseudo-random, the Dopplers and delays
%   are distinct, and the fourth-order corrections scale as 1/K.  Drawing
%   R directly is therefore statistically equivalent to generating and
%   correlating 16368 samples per dwell, and it is roughly 500x faster,
%   which is what makes 10^4-trial Monte Carlo studies practical.
%
%   The equivalence is not assumed: verify/test_cov_model.m checks the
%   waveform model and this model against each other on the same scenario.
%
%   Note that this model correctly reproduces the DETERMINISTIC bias that
%   limits the estimator, because that bias lives in Rtrue (the authentic
%   satellites' contribution to the covariance) and not in the sampling
%   fluctuation.

d = @(f, v) getdef(p, f, v);

geometry = d('geometry','circ4');
spacing  = d('spacing', 0.5);
lambda   = d('lambda', 299792458/1575.42e6);
nAuth    = d('nAuth', 9);
authSnr  = d('authSnr', 10^4.5/(16*1.023e6));
nSpoof   = d('nSpoof', 9);
saprDB   = d('saprDB', 5.5);
spoofAzEl= d('spoofAzEl', [45 15]);
mpRelDB  = d('mpRelDB', -Inf);
mpAzEl   = d('mpAzEl', [45 -12]);
jamSnr   = d('jamSnr', 0);
jamAzEl  = d('jamAzEl', [200 8]);
K        = d('K', 16368);
elMask   = d('elMask', 5);
drawSamp = d('draw', true);
gainMisDB = d('gainMismatchDB', 0.5);
phaseMisDeg = d('phaseMismatchDeg', 5);

pos = asp_array_geometry(geometry, lambda, spacing);
n = size(pos,2);

Cdiag = getdef(p,'Cdiag',[]);
if isempty(Cdiag)
    Cdiag = (10.^(gainMisDB/20*(2*rand(1,n)-1))) .* ...
            exp(1i*deg2rad(phaseMisDeg*(2*rand(1,n)-1)));
end
Cdiag = Cdiag(:);

% --- authentic directions, uniform in solid angle above the mask
authAzEl = getdef(p,'authAzEl',[]);
if isempty(authAzEl)
    az = rand(1,nAuth)*360;
    sMin = sind(elMask);
    el = asind(sMin + (1-sMin)*rand(1,nAuth));
    authAzEl = [az(:) el(:)];
end

A = bsxfun(@times, Cdiag, asp_steering(authAzEl(:,1), authAzEl(:,2), pos, lambda));
b = Cdiag .* asp_steering(spoofAzEl(1), spoofAzEl(2), pos, lambda);

pAuth  = authSnr;
pSpoof = authSnr * 10^(saprDB/10);

Rtrue = eye(n);
Rtrue = Rtrue + pAuth * (A*A');
Rtrue = Rtrue + (nSpoof*pSpoof) * (b*b');

bmp = [];
if isfinite(mpRelDB)
    bmp = Cdiag .* asp_steering(mpAzEl(1), mpAzEl(2), pos, lambda);
    Rtrue = Rtrue + (nSpoof*pSpoof*10^(mpRelDB/10)) * (bmp*bmp');
end

jam = [];
if jamSnr > 0
    jam = Cdiag .* asp_steering(jamAzEl(1), jamAzEl(2), pos, lambda);
    Rtrue = Rtrue + jamSnr * (jam*jam');
end

Rtrue = (Rtrue + Rtrue')/2;

if drawSamp
    L = chol(Rtrue + 1e-12*eye(n), 'lower');
    G = (randn(n,K) + 1i*randn(n,K))/sqrt(2);
    X = L*G;
    Rhat = (X*X')/K;
    Rhat = (Rhat + Rhat')/2;
else
    Rhat = Rtrue;
end

m.Rtrue = Rtrue;
m.Rhat  = Rhat;
m.b     = b;
m.A     = A;
m.bmp   = bmp;
m.jam   = jam;
m.pos   = pos;
m.n     = n;
m.Cdiag = Cdiag;
m.authAzEl = authAzEl;
m.pAuth = pAuth;
m.pSpoof = pSpoof;
m.K = K;

end

% -------------------------------------------------------------------------
function v = getdef(s, f, dflt)
if isfield(s, f) && ~isempty(s.(f))
    v = s.(f);
else
    v = dflt;
end
end
