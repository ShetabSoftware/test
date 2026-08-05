function res = study_array_size(nTrial, verbose)
%STUDY_ARRAY_SIZE Quantitative comparison of array element counts.
%
%   RES = STUDY_ARRAY_SIZE(NTRIAL, VERBOSE)
%
%   Produces the evidence for the 3-vs-4 element decision.  For each
%   geometry it measures, over NTRIAL random scenes:
%
%     rho          SSV estimation accuracy (paper Figure 3's quantity)
%     nullDepth    gain toward the spoofer, dB re one antenna element
%     Gfix         authentic gain with the paper's fixed-h weights (eq. 14)
%     Gmax         authentic gain with per-satellite power maximisation (23)
%     P(G<0dB)     probability a satellite ends up worse than ONE antenna
%     P(dG<-6dB)   probability a satellite loses >6 dB vs the quiescent beam
%     detStat      eigenvalue detection statistic
%     spareDOF     degrees of freedom left after the desired signal and nulls
%
%   THEORY BEING TESTED
%   -------------------
%   With N elements and a rank-p projector, an authentic response a with
%   ||a||^2 = N retains ||P a||^2, and for isotropically distributed a
%
%       ||P a||^2 / N  ~  Beta(N - p, p)                                (1)
%       E[||P a||^2]   =  N - p                                        (2)
%
%   (2) says the power-maximisation array gain is exactly 10*log10(N-p) dB
%   over a single antenna: 3.01 dB for N=3, 4.77 dB for N=4 - an
%   improvement of 10*log10(3/2) = 1.76 dB.
%
%   (1) is the more important one.  The probability that a satellite is left
%   worse off than a single antenna is
%
%       P(||P a||^2 < 1) = (1/N)^(N-p)                                 (3)
%
%   i.e. 1/9 = 11.1% at N=3 and 1/64 = 1.6% at N=4 - a SEVENFOLD reduction,
%   and a far more consequential difference than the 1.76 dB of mean gain.
%   With ten satellites in view, three elements statistically sacrifice one
%   of them; four elements sacrifice one in six scenes.
%
%   With a rank-2 null (spoofer plus its ground bounce, or spoofer plus a
%   jammer) the same formulas give E[gain] = N-2: 0 dB at N=3, i.e. THREE
%   ELEMENTS COLLAPSE TO SINGLE-ANTENNA PERFORMANCE, versus 3.01 dB at N=4.

if nargin < 1 || isempty(nTrial), nTrial = 3000; end
if nargin < 2 || isempty(verbose), verbose = true; end

cfg = asp_config();
geoms = {'tri3','sq4','circ4','y4','circ7','circ8'};

radomeR = 0.5;      % wavelengths; ~9.5 cm at L1, a 19 cm CRPA puck

res = struct([]);

for gi = 1:numel(geoms)
    g = geoms{gi};
    % Size every geometry to the SAME radome and the same ambiguity margin,
    % so the comparison isolates the element count and the lattice.
    [spacing, si] = asp_geometry_spacing_limit(g, cfg.lambda, radomeR, -3.0);
    for rankNull = [1 2]
        r = runGeometry(cfg, g, spacing, nTrial, rankNull);
        r.geometry  = g;
        r.rankNull  = rankNull;
        r.spacing   = spacing;
        r.limitedBy = si.limitedBy;
        r.gratingDB = si.meta.gratingLevelDB;
        r.nullWidthU = si.meta.nullWidthU;
        if isempty(res), res = r; else, res(end+1) = r; end %#ok<AGROW>
    end
end

if verbose
    printTable(res, nTrial, radomeR);
end

end

% -------------------------------------------------------------------------
function r = runGeometry(cfg, geom, spacing, nTrial, rankNull)

[pos, meta] = asp_array_geometry(geom, cfg.lambda, spacing);
n = size(pos,2);

p = struct('geometry', geom, 'spacing', spacing, ...
           'lambda', cfg.lambda, 'nAuth', cfg.nAuth, ...
           'authSnr', cfg.authSnrSample, 'nSpoof', cfg.nSpoof, ...
           'saprDB', cfg.saprDB, 'spoofAzEl', cfg.spoofAzEl, ...
           'K', cfg.K, 'draw', false);

if rankNull == 2
    p.mpRelDB = cfg.spoofMultipath.relDB;
    p.mpAzEl  = cfg.spoofMultipath.azEl;
end

rho = zeros(nTrial,1);
nd  = zeros(nTrial,1);
detS = zeros(nTrial,1);
gFix = []; gMax = []; dGfix = [];

h = ones(n,1)/sqrt(n);

for t = 1:nTrial
    rng(90000 + t);
    m = asp_cov_model(p);
    Rhat = asp_wishart_draw(m.Rtrue, cfg.K);

    ssvOpt = struct('rank', rankNull, 'jacobiSweeps', cfg.est.jacobiSweeps);
    [y, dbg] = asp_ssv_from_cov(Rhat, 'evd', ssvOpt);

    rho(t) = asp_ssv_correlation(y(:,1), m.b);
    d = asp_detect(dbg.lam, cfg.K, cfg.est.detectThreshold, max(n-2,1));
    detS(t) = d.stat;

    f = asp_weights('project', y, h);
    nd(t) = 10*log10(abs(f'*m.b)^2 / real(f'*f));

    % Authentic gains, referenced to a single antenna element.
    for k = 1:size(m.A,2)
        a = m.A(:,k);
        gFix(end+1,1)  = abs(f'*a)^2/real(f'*f); %#ok<AGROW>
        % Power maximisation: MRC inside the projected subspace.
        fm = asp_weights('project', y, a);
        gMax(end+1,1) = abs(fm'*a)^2/real(fm'*fm); %#ok<AGROW>
        dGfix(end+1,1) = gFix(end) / (abs(h'*a)^2/real(h'*h)); %#ok<AGROW>
    end
end

r.n            = n;
r.maxBaseline  = meta.maxBaselineLambda;
r.aperture     = meta.apertureLambda;
r.isotropy     = meta.isotropy;
r.ambiguous    = meta.ambiguous;
r.rhoMean      = mean(rho);
r.rhoP10       = prctile_local(rho, 10);
r.nullMeanDB   = 10*log10(mean(10.^(nd/10)));
r.nullP90DB    = prctile_local(nd, 90);
r.detStatMean  = mean(detS);
r.GfixMeanDB   = 10*log10(mean(gFix));
r.GmaxMeanDB   = 10*log10(mean(gMax));
r.GmaxTheoryDB = 10*log10(max(n - rankNull, eps));
r.pBelowUnity  = mean(gMax < 1);
r.pBelowUnityTheory = (1/n)^(n-rankNull);
r.pFixDrop6    = mean(10*log10(dGfix) < -6);
r.spareDOF     = n - 1 - rankNull;

end

% -------------------------------------------------------------------------
function printTable(res, nTrial, radomeR)

fprintf('\n');
fprintf('====================================================================================\n');
fprintf(' ARRAY SIZE STUDY  -  %d scenes/geometry, SAPR %.1f dB, 1 ms dwell\n', nTrial, 5.5);
fprintf(' All geometries sized to a common radome radius of %.2f lambda (%.1f cm at L1)\n', ...
    radomeR, radomeR*19.03);
fprintf(' and a common grating-response margin of -3 dB.\n');
fprintf('====================================================================================\n');

fprintf('\n--- Geometry (independent of the scene) ---\n');
fprintf('%-8s %3s %8s %7s %8s %9s %9s\n', ...
    'geom','N','spacing','Dmax','grating','nullWidth','limitedBy');
fprintf('%s\n', repmat('-',1,60));
for k = 1:numel(res)
    r = res(k);
    if r.rankNull ~= 1, continue; end
    fprintf('%-8s %3d %8.3f %7.3f %8.2f %9.3f %9s\n', ...
        r.geometry, r.n, r.spacing, r.maxBaseline, r.gratingDB, r.nullWidthU, r.limitedBy);
end

for rankNull = [1 2]
    if rankNull == 1
        fprintf('\n--- Rank-1 null: single point spoofer, no ground bounce ---\n');
    else
        fprintf('\n--- Rank-2 null: spoofer + specular ground bounce at -6 dB ---\n');
    end
    fprintf('%-8s %3s  %6s %8s  %7s %7s %7s  %8s %8s %6s\n', ...
        'geom','N','rho','null dB','Gfix','Gmax','theory', ...
        'P(G<0dB)','isotrop','spare');
    fprintf('%s\n', repmat('-',1,84));
    for k = 1:numel(res)
        r = res(k);
        if r.rankNull ~= rankNull, continue; end
        fprintf('%-8s %3d  %6.4f %8.2f  %7.2f %7.2f %7.2f  %8.4f %8.4f %6d\n', ...
            r.geometry, r.n, r.rhoMean, r.nullMeanDB, ...
            r.GfixMeanDB, r.GmaxMeanDB, r.GmaxTheoryDB, ...
            r.pBelowUnity, r.pBelowUnityTheory, r.spareDOF);
    end
end

fprintf('\nColumns:\n');
fprintf('  spacing   side length (tri3/sq4) or ring radius (circ*/y4), wavelengths\n');
fprintf('  Dmax      max baseline in wavelengths\n');
fprintf('  grating   peak off-mainlobe array-factor response over the visible region, dB\n');
fprintf('  nullWidth first null of the array factor in direction-cosine units.  The visible\n');
fprintf('            sky spans a unit disk, so a value above ~1 means the projection null\n');
fprintf('            covers essentially the WHOLE sky and every satellite is affected.\n');
fprintf('  rho       |<yhat,b>|, the quantity plotted in Figure 3 of the paper\n');
fprintf('  null dB   gain toward the spoofer, dB re one antenna element\n');
fprintf('  Gfix      mean authentic gain with the paper eq. (14) fixed-h weights\n');
fprintf('  Gmax      mean authentic gain with the paper eq. (23) power maximisation\n');
fprintf('  theory    10*log10(N - rank), the isotropic-manifold prediction for Gmax\n');
fprintf('  P(G<0dB)  fraction of satellites left worse off than a SINGLE antenna\n');
fprintf('  isotrop   (1/N)^(N-rank), the isotropic-manifold prediction for P(G<0dB).\n');
fprintf('            Measured values exceed it because a physical small-aperture array\n');
fprintf('            has strongly correlated steering vectors - the isotropic model is\n');
fprintf('            optimistic, and increasingly so as N falls.\n');
fprintf('  spare     adaptive degrees of freedom left after 1 desired + rank nulls\n');
fprintf('\n');

end

function s = ternary(c,a,b)
if c, s = a; else, s = b; end
end

function v = prctile_local(x, p)
x = sort(x(:));
idx = max(1, min(numel(x), round(p/100*numel(x))));
v = x(idx);
end
