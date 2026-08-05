function ok = test_cov_model()
%TEST_COV_MODEL Cross-validate the covariance-domain model against waveforms.
%
%   The Monte Carlo studies use asp_cov_model, which draws a Wishart sample
%   covariance instead of generating and correlating C/A waveforms.  That is
%   a 500x speedup, and it is only legitimate if the two agree.  This test
%   checks agreement on the three quantities the studies report:
%
%     1. the eigenvalue spectrum of the sample covariance
%     2. the SSV estimation accuracy (rho against the true b)
%     3. the resulting null depth
%
%   Agreement is checked in distribution, not sample by sample, because the
%   two models cannot share a noise realisation.

ok = true;
fprintf('test_cov_model:\n');

fsLow = 4*1.023e6;
cfg = asp_config('fs', fsLow, 'fe.delayMismatchSec', 0, ...
                 'fe.loPhaseNoiseRmsDeg', 0);
K = cfg.K;
nTrial = 24;

lamWave = zeros(nTrial, cfg.nAnt);
rhoWave = zeros(nTrial,1);
ndWave  = zeros(nTrial,1);

for t = 1:nTrial
    scn = asp_scenario(cfg, 'seed', 4000+t, 'durationMs', 1);
    [x, ~] = asp_rx_generate(scn, 0, K, []);
    R = (x*x')/K;
    [y, dbg] = asp_ssv_from_cov(R, 'evd', struct('rank',1,'jacobiSweeps',6));
    lamWave(t,:) = dbg.lam(:).';
    rhoWave(t) = asp_ssv_correlation(y, scn.bTrue);
    f = asp_weights('project', y, ones(cfg.nAnt,1)/sqrt(cfg.nAnt));
    ndWave(t) = 10*log10(abs(f'*scn.bTrue)^2/real(f'*f));
end

lamCov = zeros(nTrial, cfg.nAnt);
rhoCov = zeros(nTrial,1);
ndCov  = zeros(nTrial,1);

p = struct('geometry', cfg.geometry, 'spacing', cfg.elementSpacingLambda, ...
           'lambda', cfg.lambda, 'nAuth', cfg.nAuth, ...
           'authSnr', cfg.authSnrSample, 'nSpoof', cfg.nSpoof, ...
           'saprDB', cfg.saprDB, 'spoofAzEl', cfg.spoofAzEl, 'K', K, ...
           'gainMismatchDB', cfg.fe.gainMismatchDB, ...
           'phaseMismatchDeg', cfg.fe.phaseMismatchDeg);

for t = 1:nTrial
    rng(5000+t);
    m = asp_cov_model(p);
    [y, dbg] = asp_ssv_from_cov(m.Rhat, 'evd', struct('rank',1,'jacobiSweeps',6));
    lamCov(t,:) = dbg.lam(:).';
    rhoCov(t) = asp_ssv_correlation(y, m.b);
    f = asp_weights('project', y, ones(m.n,1)/sqrt(m.n));
    ndCov(t) = 10*log10(abs(f'*m.b)^2/real(f'*f));
end

% Compare in dB / linear terms.  The waveform eigenvalues carry the channel
% gain mismatch, so compare the WHITENED spectrum shape (ratio to the mean).
shapeWave = mean(bsxfun(@rdivide, lamWave, mean(lamWave,2)), 1);
shapeCov  = mean(bsxfun(@rdivide, lamCov,  mean(lamCov ,2)), 1);

fprintf('        eigen shape waveform: %s\n', sprintf('%.4f ', shapeWave));
fprintf('        eigen shape cov model: %s\n', sprintf('%.4f ', shapeCov));
ok = report(ok, max(abs(shapeWave - shapeCov)) < 0.05, ...
    sprintf('eigenvalue spectra agree to %.4f', max(abs(shapeWave-shapeCov))));

fprintf('        rho  waveform %.4f +/- %.4f  |  cov model %.4f +/- %.4f\n', ...
    mean(rhoWave), std(rhoWave), mean(rhoCov), std(rhoCov));
ok = report(ok, abs(mean(rhoWave)-mean(rhoCov)) < 0.02, ...
    sprintf('mean SSV correlation agrees to %.4f', abs(mean(rhoWave)-mean(rhoCov))));

mnW = 10*log10(mean(10.^(ndWave/10)));
mnC = 10*log10(mean(10.^(ndCov/10)));
fprintf('        null waveform %.2f dB  |  cov model %.2f dB\n', mnW, mnC);
ok = report(ok, abs(mnW-mnC) < 3, ...
    sprintf('mean null depth agrees to %.2f dB', abs(mnW-mnC)));

end

function ok = report(ok, cond, msg)
if cond, fprintf('  PASS  %s\n', msg);
else,    fprintf('  FAIL  %s\n', msg); ok = false; end
end
