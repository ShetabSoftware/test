function [thr, info] = asp_detect_threshold(cfg, pfa, nTrial)
%ASP_DETECT_THRESHOLD Calibrate the spoofing detector threshold.
%
%   [THR, INFO] = ASP_DETECT_THRESHOLD(CFG, PFA, NTRIAL)
%
%   Returns the decision threshold on lambda_1/mean(lambda_2..lambda_N) that
%   achieves a false-alarm probability of PFA (default 1e-3) under H0.
%
%   H0 IS NOT WHITE NOISE, AND THAT IS THE WHOLE POINT
%   --------------------------------------------------
%   The natural instinct is to set the threshold from the asymptotic
%   largest-eigenvalue result for a complex Wishart matrix,
%
%       E[lambda_1] -> (1 + sqrt(N/K))^2                                (1)
%
%   That is wrong here, and it is wrong in the dangerous direction.  The
%   no-attack condition is not an empty sky: it contains 8-12 authentic
%   satellites, which are spatially structured sources.  They raise the
%   test statistic above (1), and unlike the sampling fluctuation their
%   contribution does NOT fall as 1/sqrt(K).  Measured for a 4-element Y
%   array at fs = 16.368 MHz:
%
%       dwell    H0 mean stat   eq. (1) predicts
%        1 ms      1.0313          1.0315
%        5 ms      1.0230          1.0140
%       20 ms      1.0211          1.0070
%
%   At 1 ms the two agree and (1) looks fine.  At 20 ms the authentic
%   contribution is twice the sampling term, and a threshold set from (1)
%   would false-alarm continuously on a clean sky.  Anyone who validates
%   their threshold at a short dwell and then lengthens it will walk
%   straight into this.
%
%   The threshold also depends on fs, because the per-sample SNR of the
%   authentic signals is (C/N0)/fs.  A threshold calibrated at one sample
%   rate is not valid at another.  Hard-coding it - as a constant in a
%   config file - is therefore a latent bug that only appears when someone
%   changes an apparently unrelated parameter.
%
%   So: calibrate by Monte Carlo against a realistic H0, which is what a
%   product does at design time anyway.  The O(N^2) Bartlett draw makes
%   this cheap enough to do inside asp_config.
%
%   INFO.h0Mean, INFO.h0Std, INFO.asymptotic, INFO.margin

if nargin < 2 || isempty(pfa),    pfa = 1e-3;   end
if nargin < 3 || isempty(nTrial), nTrial = 400; end

K = cfg.K * cfg.est.coherentMs;

p = struct('geometry', cfg.geometry, 'spacing', cfg.elementSpacingLambda, ...
           'lambda', cfg.lambda, 'nAuth', cfg.nAuth, ...
           'authSnr', cfg.authSnrSample, 'nSpoof', 0, 'saprDB', -200, ...
           'K', K, 'draw', false, ...
           'gainMismatchDB', cfg.fe.gainMismatchDB, ...
           'phaseMismatchDeg', cfg.fe.phaseMismatchDeg);

stat = zeros(nTrial,1);
n = 0;
for t = 1:nTrial
    rng(770000 + t);
    m = asp_cov_model(p);
    n = m.n;
    Rh = asp_wishart_draw(m.Rtrue, K);
    [~, dbg] = asp_ssv_from_cov(Rh, 'evd', ...
        struct('rank',1,'jacobiSweeps',cfg.est.jacobiSweeps));
    lam = dbg.lam;
    stat(t) = lam(1)/mean(lam(2:end));
end

s = sort(stat);
idx = max(1, min(nTrial, ceil((1-pfa)*nTrial)));
empirical = s(idx);

% With a few hundred trials the empirical 1-pfa quantile is noisy, so use a
% Gaussian tail extrapolation from the bulk and take the larger of the two.
% Erring high costs detection sensitivity; erring low costs false alarms on
% a clean sky, which is the failure the product cannot afford.
mu = mean(stat);
sd = std(stat);
zq = sqrt(2)*erfinv(1 - 2*pfa);
extrapolated = mu + zq*sd;

thr = max(empirical, extrapolated);

info.h0Mean     = mu;
info.h0Std      = sd;
info.empirical  = empirical;
info.extrapolated = extrapolated;
info.asymptotic = (1 + sqrt(n/K))^2;
info.margin     = thr - mu;
info.nTrial     = nTrial;
info.pfa        = pfa;

end
