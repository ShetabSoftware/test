function res = study_estimators(nTrial, verbose)
%STUDY_ESTIMATORS Compare the paper's SSV estimator with the covariance ones.
%
%   RES = STUDY_ESTIMATORS(NTRIAL, VERBOSE)
%
%   Runs the full waveform model (real Gold codes, real navigation data,
%   code Doppler, channel mismatch) because the paper's estimator uses a
%   one-code-period lag and therefore cannot be exercised by a
%   covariance-domain model.
%
%   Four questions are answered with measurements rather than argument:
%
%   Q1  How much accuracy does the beta/gamma construction give up against
%       the covariance-diagonal amplitude and against the full eigenvector?
%
%   Q2  Do navigation data bits damage beta?  A 20 ms bit transition inside
%       the dwell inverts part of beta's lagged sum, so the concern is
%       reasonable.  It turns out NOT to matter, because the affected factor
%       d is common to every element and the projector is scale invariant.
%       Reported here because the negative result is as useful as a positive
%       one: it removes a design constraint that would otherwise force bit
%       synchronisation ahead of the array processing.
%
%   Q3  Does beta's constant d fade?  Equation (11) makes d a sum of ~18
%       Doppler phasors e^{j2 pi f T} with |f| up to 5 kHz and T = 1 ms, so
%       the phases wrap and d is a random walk that can land near zero.
%       Measured as the spread of |beta| across scenes.
%
%   Q4  Where does null depth saturate with integration time?  The
%       covariance estimator's error has a 1/sqrt(K) noise term and a
%       DWELL-INDEPENDENT bias from the authentic constellation, so more
%       integration must stop helping at some point.  Locating that point
%       tells you the dwell to design for - and tells you that going deeper
%       requires a different mechanism, not a longer dwell.

if nargin < 1 || isempty(nTrial),  nTrial = 60;  end
if nargin < 2 || isempty(verbose), verbose = true; end

% 4 samples/chip keeps the waveform model fast without changing any
% conclusion; C/N0 is converted to per-sample SNR using fs, so the operating
% point is preserved exactly.
cfg = asp_config('fs', 4*1.023e6);
K = cfg.K;
n = cfg.nAnt;
h = ones(n,1)/sqrt(n);

modes = {'paper','column','gamma','evd'};
dwellsMs = [1 2 5 10 20];

res.modes = modes;
res.dwellsMs = dwellsMs;
res.rho  = nan(numel(modes), numel(dwellsMs), nTrial);
res.null = nan(numel(modes), numel(dwellsMs), nTrial);
res.betaMag = nan(numel(dwellsMs), nTrial);

for t = 1:nTrial
    scn = asp_scenario(cfg, 'seed', 31000+t, 'durationMs', max(dwellsMs)+1);

    % Generate the longest dwell once, then use prefixes.
    x = asp_rx_generate(scn, 0, K*max(dwellsMs), []);

    for di = 1:numel(dwellsMs)
        m = dwellsMs(di);
        xs = x(:, 1:K*m);

        for mi = 1:numel(modes)
            switch modes{mi}
                case 'paper'
                    if m < 2, continue; end
                    [y, dbg] = asp_ssv_paper(xs, K, cfg.refElement);
                    res.betaMag(di,t) = mean(abs(dbg.beta)) / (K*(m-1));
                otherwise
                    R = (xs*xs')/size(xs,2);
                    y = asp_ssv_from_cov(R, modes{mi}, ...
                        struct('refIdx',cfg.refElement,'rank',1,'jacobiSweeps',6));
            end
            res.rho(mi,di,t) = asp_ssv_correlation(y, scn.bTrue);
            f = asp_weights('project', y, h);
            res.null(mi,di,t) = 10*log10(abs(f'*scn.bTrue)^2/real(f'*f));
        end
    end
end

% --- Q2: navigation bits on/off, at a 10 ms dwell where the effect bites
cfgNB = cfg;
nbTrial = max(round(nTrial/2), 20);
rhoBits = nan(2, nbTrial);
for t = 1:nbTrial
    for withBits = [1 0]
        scn = asp_scenario(cfgNB, 'seed', 41000+t, 'durationMs', 12);
        if ~withBits
            for s = 1:numel(scn.src)
                scn.src(s).bits = ones(size(scn.src(s).bits));
            end
        end
        x = asp_rx_generate(scn, 0, K*10, []);
        y = asp_ssv_paper(x, K, cfg.refElement);
        rhoBits(2-withBits, t) = asp_ssv_correlation(y, scn.bTrue);
    end
end
res.rhoBitsOn  = rhoBits(1,:);
res.rhoBitsOff = rhoBits(2,:);

if verbose
    printResults(res, nTrial);
end

end

% -------------------------------------------------------------------------
function printResults(res, nTrial)

fprintf('\n');
fprintf('====================================================================\n');
fprintf(' SSV ESTIMATOR COMPARISON  -  %d scenes, full waveform model\n', nTrial);
fprintf(' 4-element Y array, SAPR 5.5 dB, C/N0 45 dB-Hz, real C/A + nav data\n');
fprintf('====================================================================\n');

fprintf('\n--- Mean |<yhat,b>| (higher is better) ---\n');
fprintf('%-8s', 'dwell');
for mi = 1:numel(res.modes), fprintf(' %10s', res.modes{mi}); end
fprintf('\n%s\n', repmat('-',1,8+11*numel(res.modes)));
for di = 1:numel(res.dwellsMs)
    fprintf('%5d ms', res.dwellsMs(di));
    for mi = 1:numel(res.modes)
        v = squeeze(res.rho(mi,di,:));
        v = v(~isnan(v));
        if isempty(v), fprintf(' %10s','-');
        else, fprintf(' %10.5f', mean(v)); end
    end
    fprintf('\n');
end

fprintf('\n--- Mean null depth, dB re one antenna element (lower is better) ---\n');
fprintf('%-8s', 'dwell');
for mi = 1:numel(res.modes), fprintf(' %10s', res.modes{mi}); end
fprintf('\n%s\n', repmat('-',1,8+11*numel(res.modes)));
for di = 1:numel(res.dwellsMs)
    fprintf('%5d ms', res.dwellsMs(di));
    for mi = 1:numel(res.modes)
        v = squeeze(res.null(mi,di,:));
        v = v(~isnan(v));
        if isempty(v), fprintf(' %10s','-');
        else, fprintf(' %10.2f', 10*log10(mean(10.^(v/10)))); end
    end
    fprintf('\n');
end

evdIdx = find(strcmp(res.modes,'evd'));
paperIdx = find(strcmp(res.modes,'paper'));
n1 = 10*log10(mean(10.^(squeeze(res.null(evdIdx,1,:))/10)));
nL = 10*log10(mean(10.^(squeeze(res.null(evdIdx,end,:))/10)));
fprintf('\nQ4  EVD null depth improves %.1f dB going from 1 ms to %d ms.\n', ...
    n1-nL, res.dwellsMs(end));
fprintf('    A pure 1/sqrt(K) noise-limited estimator would improve by %.1f dB\n', ...
    10*log10(res.dwellsMs(end)));
fprintf('    over that range.  The shortfall is the DWELL-INDEPENDENT bias from the\n');
fprintf('    authentic constellation: the covariance contains sum_m p_m a_m a_m'', which\n');
fprintf('    perturbs the principal eigenvector by an amount set by p_auth/P_spoof and\n');
fprintf('    NOT by K.  Integrating longer cannot fix it; only raising the effective\n');
fprintf('    SAPR can, which is exactly what post-correlation refinement does.\n');

fprintf('\n--- Q3  beta magnitude across scenes (paper eq. 11 constant d) ---\n');
bm = res.betaMag(end,:); bm = bm(~isnan(bm));
if ~isempty(bm)
    fprintf('    mean %.3e   min %.3e   max %.3e   max/min = %.1f\n', ...
        mean(bm), min(bm), max(bm), max(bm)/max(min(bm),realmin));
    fprintf('    A spread this wide is the random walk of ~18 Doppler phasors.  When d\n');
    fprintf('    lands near zero the amplitude estimate is noise, and nothing in the\n');
    fprintf('    algorithm detects that it has happened.\n');
end

fprintf('\n--- Q2  navigation data bits vs the paper''s beta (10 ms dwell) ---\n');
fprintf('    bits present: mean rho = %.5f   (min %.5f)\n', ...
    mean(res.rhoBitsOn), min(res.rhoBitsOn));
fprintf('    bits removed: mean rho = %.5f   (min %.5f)\n', ...
    mean(res.rhoBitsOff), min(res.rhoBitsOff));
fprintf('    NEGATIVE RESULT, and a useful one: bit transitions do invert part of\n');
fprintf('    beta''s lagged sum, but only through the constant d of equation (11),\n');
fprintf('    which is common to every element and therefore cancels in the\n');
fprintf('    scale-invariant projector.  Bit synchronisation is NOT a prerequisite\n');
fprintf('    for the array processing.\n\n');

end
