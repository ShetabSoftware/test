function out = asp_process(scn, varargin)
%ASP_PROCESS Streaming, causal anti-spoofing pipeline (the FPGA reference).
%
%   OUT = ASP_PROCESS(SCN, 'name', value, ...)
%
%   Options:
%     'estimator'   'evd' (default) | 'gamma' | 'column' | 'paper'
%     'fixedPoint'  false (default) | true
%     'detect'      true (default): gate the nulling on the detector
%     'hMode'       'ones' (default) | 'principal'
%     'dwellMs'     covariance dwell length in ms (default cfg.est.coherentMs)
%     'verbose'     false
%
%   OUT.v            1 x M beamformer output
%   OUT.vQuiescent   1 x M output of the quiescent (non-adaptive) beam
%   OUT.block(k)     per-dwell diagnostics: weights, eigenvalues, detection,
%                    pattern metrics, saturation, timing
%
%   THREE STRUCTURAL DIFFERENCES FROM THE ORIGINAL SCRIPTS
%   ------------------------------------------------------
%   1. CAUSALITY.  Weights computed from dwell k are applied to dwell k+1.
%      The original scripts estimate the SSV from the whole record and then
%      apply the resulting weights to that same record, including the
%      samples used for the estimate.  That is not implementable and it
%      flatters the result: it removes the estimator's transient entirely
%      and hides the interaction between weight-update rate and platform
%      dynamics.  On a rotating platform the weight age is a first-order
%      effect - at 100 deg/s and a 1 ms dwell the signature rotates 0.1 deg
%      between estimate and use, capping the null at about -42 dB; at
%      400 deg/s the cap is about -30 dB.
%
%   2. GATING.  The detector decides whether to null at all, and with what
%      rank.  Without it, the system attacks its own satellites whenever the
%      threat is absent, which is almost all of the time.
%
%   3. HARDWARE IN THE LOOP.  AGC, ADC quantisation, saturation counting and
%      the fixed-point covariance accumulator are inside the loop, not
%      applied afterwards to an already-computed floating-point answer.
%      The distinction matters: an AGC gain step or a saturating sample
%      corrupts the dwell in which it occurs, and the pipeline has to notice
%      and hold the previous weights rather than use a poisoned estimate.

cfg = scn.cfg;

opt.estimator  = cfg.est.mode;
opt.fixedPoint = false;
opt.detect     = true;
opt.hMode      = 'ones';
opt.dwellMs    = cfg.est.coherentMs;
opt.verbose    = false;
opt.durationMs = scn.durationMs;

for k = 1:2:numel(varargin)
    if ~isfield(opt, varargin{k})
        error('asp_process:opt', 'Unknown option "%s".', varargin{k});
    end
    opt.(varargin{k}) = varargin{k+1};
end

K        = cfg.K;
blockLen = K * opt.dwellMs;
nBlocks  = floor(opt.durationMs / opt.dwellMs);
nAnt     = scn.nAnt;

h = ones(nAnt,1)/sqrt(nAnt);

fCur   = h;           % weights in force during the current block
rankCur = 0;
genState = [];
agc      = [];

v          = zeros(1, nBlocks*blockLen);
vQuiescent = zeros(1, nBlocks*blockLen);
blocks = struct([]);

% Paper-mode estimator needs raw epochs, so keep the previous epoch when
% required.  This buffer IS the 58-BRAM cost that the covariance path avoids.
prevEpoch = [];

for kb = 1:nBlocks
    [x, genState] = asp_rx_generate(scn, kb-1, blockLen, genState);

    if opt.fixedPoint
        [xq, agc] = asp_agc_adc(x, cfg, agc);
        agcChanged = agc.changed;
        satFrac    = agc.satFraction;
    else
        xq = x;
        agcChanged = false;
        satFrac = 0;
    end

    % ------------------------------------------------- beamform (causal)
    idx = (kb-1)*blockLen + (1:blockLen);
    if opt.fixedPoint
        [vb, bfInfo] = asp_beamform(xq, fCur, cfg.fx);
        vq           = asp_beamform(xq, h,    cfg.fx);
    else
        [vb, bfInfo] = asp_beamform(xq, fCur);
        vq           = asp_beamform(xq, h);
    end
    v(idx)          = vb;
    vQuiescent(idx) = vq;

    % ------------------------------------------------- estimate for NEXT block
    tEst = tic;
    switch lower(opt.estimator)
        case 'paper'
            if isempty(prevEpoch)
                yEst = [];
                lam  = [];
                det  = struct('detected',false,'rank',0,'stat',NaN,'snrEstDB',NaN,'lam',[]);
            else
                [yEst, ~] = asp_ssv_paper([prevEpoch, xq], K, cfg.refElement);
                lam = [];
                % The paper's estimator produces no eigenvalues, so it
                % cannot gate itself.  Model that honestly: always detected.
                det = struct('detected',true,'rank',1,'stat',NaN,'snrEstDB',NaN,'lam',[]);
            end
            R = [];

        otherwise
            if opt.fixedPoint
                accBlock = fx_cov_accum(xq, cfg.fx.adc, cfg.fx.covAccWl);
                R = accBlock.R;
                covBits = accBlock.count;
            else
                R = (xq * xq') / blockLen;
                covBits = blockLen;
            end

            ssvOpt = struct('refIdx', cfg.refElement, 'rank', cfg.est.maxNullRank, ...
                            'jacobiSweeps', cfg.est.jacobiSweeps);
            [yAll, dbg] = asp_ssv_from_cov(R, opt.estimator, ssvOpt);

            if isempty(dbg.lam)
                % Column/gamma modes carry no spectrum; fall back to a
                % dedicated whitened EVD purely for detection.
                [~, lamOnly] = asp_evd_herm(dbg.Rw, cfg.est.jacobiSweeps);
                lam = lamOnly;
            else
                lam = dbg.lam;
            end

            if opt.detect
                det = asp_detect(lam, blockLen, cfg.est.detectThreshold, cfg.est.maxNullRank);
            else
                det = asp_detect(lam, blockLen, -Inf, cfg.est.maxNullRank);
                det.rank = max(det.rank, 1);
            end

            if det.rank > 0
                yEst = yAll(:, 1:min(det.rank, size(yAll,2)));
            else
                yEst = [];
            end
    end
    estTime = toc(tEst);

    % ------------------------------------------------- quiescent vector
    switch lower(opt.hMode)
        case 'ones'
            hUse = h;
        case 'principal'
            % Maximise total residual power after projection: the leading
            % eigenvector of the projected covariance.  Calibration free.
            if ~isempty(R) && ~isempty(yEst)
                fTmp = asp_weights('project', yEst, eye(nAnt));
                Rp = fTmp' * R * fTmp;
                [Up,~] = asp_evd_herm((Rp+Rp')/2, cfg.est.jacobiSweeps);
                hUse = fTmp * Up(:,1);
            else
                hUse = h;
            end
        otherwise
            error('asp_process:hMode','Unknown hMode "%s".', opt.hMode);
    end

    % ------------------------------------------------- weights for NEXT block
    dwellValid = ~agcChanged && satFrac < 1e-3;
    if dwellValid
        if isempty(yEst)
            fNext = h;             % no threat: pass the quiescent beam
            rankCur = 0;
        else
            fNext = asp_weights('project', yEst, hUse);
            rankCur = size(yEst,2);
        end
        fCur = fNext;
    end
    % If the dwell was invalid, fCur is simply held: a corrupted estimate is
    % strictly worse than a slightly stale one.

    % ------------------------------------------------- diagnostics
    bl.k          = kb;
    bl.f          = fCur;
    bl.rank       = rankCur;
    bl.lam        = lam;
    bl.detected   = det.detected;
    bl.detectStat = det.stat;
    bl.snrEstDB   = det.snrEstDB;
    bl.dwellValid = dwellValid;
    bl.satFraction = satFrac;
    bl.estTimeSec = estTime;
    bl.metrics    = asp_pattern_metrics(fCur, scn, h);
    if ~isempty(yEst)
        bl.ssvCorr = asp_ssv_correlation(yEst(:,1), scn.bTrue);
    else
        bl.ssvCorr = NaN;
    end
    if opt.fixedPoint
        bl.accBitsUsed = bfInfo.accBitsUsed;
    end

    if isempty(blocks)
        blocks = bl;
    else
        blocks(end+1) = bl; %#ok<AGROW>
    end

    if opt.verbose
        fprintf(['  blk %3d  det=%d rank=%d stat=%.4f  null=%7.2f dB  ' ...
                 'authMean=%6.2f dB  rho=%.4f\n'], ...
            kb, det.detected, rankCur, det.stat, bl.metrics.nullGainDB, ...
            bl.metrics.authGainMeanDB, bl.ssvCorr);
    end

    prevEpoch = xq(:, end-K+1:end);
end

out.v          = v;
out.vQuiescent = vQuiescent;
out.block      = blocks;
out.opt        = opt;
out.nBlocks    = nBlocks;
out.blockLen   = blockLen;

end
