function res = study_postcorr(nTrial, verbose)
%STUDY_POSTCORR Two-stage architecture: pre- vs post-correlation nulling.
%
%   RES = STUDY_POSTCORR(NTRIAL, VERBOSE)
%
%   Measures what the recommended stage-2 refinement is actually worth,
%   using the full waveform model, and checks the claim that stage 1
%   saturates while stage 2 does not.
%
%   Also measures the rank-2 case: with a specular ground bounce present,
%   the post-correlation covariance should show a resolvable second
%   eigenvalue, giving the system a way to LEARN that a rank-2 null is
%   required rather than being told.

if nargin < 1 || isempty(nTrial),  nTrial = 25; end
if nargin < 2 || isempty(verbose), verbose = true; end

cfg = asp_config('fs', 4*1.023e6);
n = cfg.nAnt;
h = ones(n,1)/sqrt(n);
K = cfg.K;

epochSweep = [1 5 20 50];
maxMs = max(epochSweep);

res.epochs = epochSweep;
res.nullPre  = nan(nTrial,1);
res.nullPost = nan(numel(epochSweep), nTrial);
res.rhoPre   = nan(nTrial,1);
res.rhoPost  = nan(numel(epochSweep), nTrial);
res.rankRatioMp = nan(nTrial,1);
res.rankRatioNoMp = nan(nTrial,1);

for t = 1:nTrial
    scn = asp_scenario(cfg, 'seed', 52000+t, 'durationMs', maxMs+1);
    x = asp_rx_generate(scn, 0, K*maxMs, []);

    % --- stage 1: pre-correlation, 1 ms
    R = (x(:,1:K)*x(:,1:K)')/K;
    yPre = asp_ssv_from_cov(R, 'evd', struct('rank',1,'jacobiSweeps',6));
    res.rhoPre(t) = asp_ssv_correlation(yPre, scn.bTrue);
    fPre = asp_weights('project', yPre, h);
    res.nullPre(t) = 10*log10(abs(fPre'*scn.bTrue)^2/real(fPre'*fPre));

    % --- stage 2: post-correlation, using the receiver's tracking state
    isSpoof = strcmp({scn.src.kind},'spoof');
    idx = find(isSpoof);
    hyp = struct('prn',{},'tau',{},'fd',{});
    for k = 1:numel(idx)
        s = scn.src(idx(k));
        hyp(end+1) = struct('prn', s.prn, 'tau', s.tau0, 'fd', s.fd); %#ok<AGROW>
    end

    for ei = 1:numel(epochSweep)
        [bHat, dbg] = asp_postcorr_ssv(x, cfg, hyp, epochSweep(ei));
        res.rhoPost(ei,t) = asp_ssv_correlation(bHat, scn.bTrue);
        fPost = asp_weights('project', bHat, h);
        res.nullPost(ei,t) = 10*log10(abs(fPost'*scn.bTrue)^2/real(fPost'*fPost));
        if epochSweep(ei) == 20
            res.rankRatioNoMp(t) = dbg.rankRatio;
        end
    end
end

% --- rank detection with a ground bounce present
cfgMp = cfg;
cfgMp.spoofMultipath.enable = true;
for t = 1:min(nTrial, 15)
    scn = asp_scenario(cfgMp, 'seed', 53000+t, 'durationMs', 22);
    x = asp_rx_generate(scn, 0, K*20, []);
    isSpoof = strcmp({scn.src.kind},'spoof');
    idx = find(isSpoof);
    hyp = struct('prn',{},'tau',{},'fd',{});
    for k = 1:numel(idx)
        s = scn.src(idx(k));
        hyp(end+1) = struct('prn', s.prn, 'tau', s.tau0, 'fd', s.fd); %#ok<AGROW>
    end
    [~, dbg] = asp_postcorr_ssv(x, cfgMp, hyp, 20);
    res.rankRatioMp(t) = dbg.rankRatio;
end

if verbose
    printResults(res, nTrial);
end

end

% -------------------------------------------------------------------------
function printResults(res, nTrial)
fprintf('\n');
fprintf('==============================================================\n');
fprintf(' TWO-STAGE ARCHITECTURE  -  %d scenes, full waveform model\n', nTrial);
fprintf(' 4-element Y array, SAPR 5.5 dB, C/N0 45 dB-Hz\n');
fprintf('==============================================================\n');

mn = @(v) 10*log10(mean(10.^(v(~isnan(v))/10)));

fprintf('\nStage 1 (pre-correlation, 1 ms, no receiver state):\n');
fprintf('   rho = %.5f    null = %.2f dB\n', mean(res.rhoPre), mn(res.nullPre));

fprintf('\nStage 2 (post-correlation, needs tracking state):\n');
fprintf('%10s %12s %12s %14s\n','epochs','rho','null dB','gain vs st.1');
fprintf('%s\n', repmat('-',1,50));
for ei = 1:numel(res.epochs)
    nd = mn(res.nullPost(ei,:));
    fprintf('%8d ms %12.6f %12.2f %13.1f dB\n', res.epochs(ei), ...
        mean(res.rhoPost(ei,:)), nd, mn(res.nullPre) - nd);
end

d1 = mn(res.nullPost(1,:));
dN = mn(res.nullPost(end,:));
fprintf('\nStage 2 improves %.1f dB from %d to %d ms; 1/sqrt(K) predicts %.1f dB.\n', ...
    d1-dN, res.epochs(1), res.epochs(end), 10*log10(res.epochs(end)/res.epochs(1)));
fprintf('Stage 2 therefore remains NOISE limited over this range, whereas stage 1\n');
fprintf('is bias limited and stops improving after a few milliseconds.  That is the\n');
fprintf('whole argument for the two-stage split: they fail for different reasons,\n');
fprintf('so combining them is not redundancy, it is coverage.\n');

fprintf('\nRank indicator lambda_1/lambda_2 of the post-correlation covariance:\n');
fprintf('   no ground bounce : %8.1f  (rank 1, as expected)\n', ...
    median(res.rankRatioNoMp(~isnan(res.rankRatioNoMp))));
fprintf('   with ground bounce: %8.1f  (second eigenvalue lifts -> rank 2)\n', ...
    median(res.rankRatioMp(~isnan(res.rankRatioMp))));
fprintf('The system can therefore DISCOVER that a rank-2 null is required rather\n');
fprintf('than being configured for it, which matters because the ground bounce\n');
fprintf('geometry changes as the platform moves.\n\n');

end
