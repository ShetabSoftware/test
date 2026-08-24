function res = study_detector(nTrial, verbose)
%STUDY_DETECTOR Detection performance of the eigenvalue spoofing detector.
%
%   RES = STUDY_DETECTOR(NTRIAL, VERBOSE)
%
%   Answers the question the paper never asks: at what spoofing power, and
%   with how much integration, can the system KNOW it is under attack?
%
%   This matters more than null depth for a shipping product.  An
%   anti-spoofing device spends almost all of its operating hours not under
%   attack.  If it nulls unconditionally - as the paper's pipeline and all
%   three provided scripts do - then during those hours it steers a null
%   into the strongest authentic satellite and degrades the receiver it is
%   supposed to protect.  The expected value of the device is then
%
%       P(attack) * benefit  -  P(no attack) * damage
%
%   and with P(attack) small, the second term dominates.  Detection is not
%   a feature; it is the precondition for the product being worth fitting.
%
%   METHOD
%   ------
%   H0: authentic constellation only, no spoofer.  Note this is NOT white
%   noise: 9 satellites at 45 dB-Hz are present and DO perturb the
%   eigenvalues.  Calibrating the threshold against pure white noise would
%   understate the false alarm rate.
%
%   H1: same, plus a spoofer at the swept SAPR.
%
%   Statistic: lambda_1 / mean(lambda_2..lambda_N) on the WHITENED sample
%   covariance.  Threshold set empirically from the H0 distribution at the
%   requested false-alarm probability, which is what a product would do at
%   calibration time.  The asymptotic (1+sqrt(N/K))^2 expression is printed
%   alongside as a cross-check.
%
%   WHY THE STATISTIC WORKS AT ALL AT 0 dB SAPR
%   -------------------------------------------
%   At 0 dB SAPR the spoofer's TOTAL power equals the authentic
%   constellation's total power, so a power-based test cannot separate them.
%   The eigenvalue test still can, because the two differ in SPATIAL RANK:
%
%       spoofer:   sum_k p_k * b b'  =  (N_spoof * p_s) * b b'   -> rank 1
%       authentic: sum_m p_m a_m a_m' ~ (N_auth * p_a) * I       -> rank N
%
%   The authentic sum is approximately isotropic because E[a a'] = I for
%   directions spread over the sky, so it raises the noise FLOOR without
%   creating an eigenvalue SPREAD.  The spoofer puts all of its power into
%   one eigenvalue.  The detector is therefore measuring spatial coherence,
%   not power - which is the correct physical discriminant and a stronger
%   one than the "the spoofer is louder" argument the paper relies on.

if nargin < 1 || isempty(nTrial),  nTrial = 4000; end
if nargin < 2 || isempty(verbose), verbose = true; end

cfg = asp_config();
pfaTarget = 1e-3;

geoms   = {'tri3','y4','circ7'};
saprSweep = -10:2:14;
dwellsMs  = [1 5 20];

res.geoms = geoms;
res.sapr  = saprSweep;
res.dwellsMs = dwellsMs;
res.pd    = nan(numel(geoms), numel(dwellsMs), numel(saprSweep));
res.thr   = nan(numel(geoms), numel(dwellsMs));
res.minDetectableSapr = nan(numel(geoms), numel(dwellsMs));

for gi = 1:numel(geoms)
    spacing = asp_geometry_spacing_limit(geoms{gi}, cfg.lambda, 0.5, -3.0);

    for di = 1:numel(dwellsMs)
        K = cfg.K * dwellsMs(di);

        base = struct('geometry', geoms{gi}, 'spacing', spacing, ...
                      'lambda', cfg.lambda, 'nAuth', cfg.nAuth, ...
                      'authSnr', cfg.authSnrSample, 'nSpoof', cfg.nSpoof, ...
                      'K', K, 'draw', false);

        % --- H0 distribution
        stat0 = zeros(nTrial,1);
        for t = 1:nTrial
            rng(600000 + t);
            p0 = base; p0.nSpoof = 0; p0.saprDB = -200;
            m = asp_cov_model(p0);
            Rh = asp_wishart_draw(m.Rtrue, K);
            stat0(t) = statOf(Rh, cfg);
        end
        thr = quantile_local(stat0, 1 - pfaTarget);
        res.thr(gi,di) = thr;
        res.h0mean(gi,di) = mean(stat0);
        res.h0asym(gi,di) = (1 + sqrt(m.n/K))^2;

        % --- H1 sweep
        for si = 1:numel(saprSweep)
            hit = 0;
            for t = 1:nTrial
                rng(700000 + 1000*si + t);
                p1 = base; p1.saprDB = saprSweep(si);
                m = asp_cov_model(p1);
                Rh = asp_wishart_draw(m.Rtrue, K);
                hit = hit + (statOf(Rh, cfg) > thr);
            end
            res.pd(gi,di,si) = hit/nTrial;
        end

        pdv = squeeze(res.pd(gi,di,:));
        idx = find(pdv >= 0.9, 1, 'first');
        if ~isempty(idx)
            if idx == 1
                res.minDetectableSapr(gi,di) = saprSweep(1);
            else
                % linear interpolation on Pd
                x0 = saprSweep(idx-1); x1 = saprSweep(idx);
                y0 = pdv(idx-1);       y1 = pdv(idx);
                res.minDetectableSapr(gi,di) = x0 + (0.9-y0)*(x1-x0)/(y1-y0);
            end
        end
    end
end

if verbose
    printResults(res, nTrial, pfaTarget);
end

end

% -------------------------------------------------------------------------
function s = statOf(R, cfg)
[~, dbg] = asp_ssv_from_cov(R, 'evd', struct('rank',1,'jacobiSweeps',cfg.est.jacobiSweeps));
lam = dbg.lam;
s = lam(1) / mean(lam(2:end));
end

function printResults(res, nTrial, pfa)
fprintf('\n');
fprintf('=========================================================================\n');
fprintf(' SPOOFING DETECTOR  -  %d trials per point, target Pfa = %.0e\n', nTrial, pfa);
fprintf(' H0 = authentic constellation only (9 sats @ 45 dB-Hz), no spoofer\n');
fprintf('=========================================================================\n');

fprintf('\n--- Decision threshold on lambda_1 / mean(lambda_2..N) ---\n');
fprintf('%-8s %8s %10s %10s %12s\n','geom','dwell','H0 mean','threshold','white-noise');
fprintf('%s\n', repmat('-',1,52));
for gi = 1:numel(res.geoms)
    for di = 1:numel(res.dwellsMs)
        fprintf('%-8s %6d ms %10.4f %10.4f %12.4f\n', res.geoms{gi}, ...
            res.dwellsMs(di), res.h0mean(gi,di), res.thr(gi,di), res.h0asym(gi,di));
    end
end
fprintf('\nThe "white-noise" column is the asymptotic (1+sqrt(N/K))^2 for a pure\n');
fprintf('noise input.  The measured H0 mean sits ABOVE it because the authentic\n');
fprintf('constellation is itself a set of spatially structured sources.  A product\n');
fprintf('that set its threshold from the white-noise formula would false-alarm\n');
fprintf('on a clean sky.\n');

fprintf('\n--- Minimum SAPR for 90%% detection probability, dB ---\n');
fprintf('%-8s', 'geom');
for di = 1:numel(res.dwellsMs), fprintf(' %8d ms', res.dwellsMs(di)); end
fprintf('\n%s\n', repmat('-',1,8+11*numel(res.dwellsMs)));
for gi = 1:numel(res.geoms)
    fprintf('%-8s', res.geoms{gi});
    for di = 1:numel(res.dwellsMs)
        v = res.minDetectableSapr(gi,di);
        if isnan(v), fprintf(' %11s','>max');
        else, fprintf(' %11.1f', v); end
    end
    fprintf('\n');
end

fprintf('\n--- Detection probability vs SAPR (1 ms dwell) ---\n');
fprintf('%-8s', 'SAPR dB');
for gi = 1:numel(res.geoms), fprintf(' %8s', res.geoms{gi}); end
fprintf('\n%s\n', repmat('-',1,8+9*numel(res.geoms)));
for si = 1:numel(res.sapr)
    fprintf('%8.0f', res.sapr(si));
    for gi = 1:numel(res.geoms)
        fprintf(' %8.3f', res.pd(gi,1,si));
    end
    fprintf('\n');
end
fprintf('\n');

end

function v = quantile_local(x, q)
x = sort(x(:));
n = numel(x);
idx = max(1, min(n, ceil(q*n)));
v = x(idx);
end
