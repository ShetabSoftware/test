function outDir = asp_export_vectors(outDir, cfg, nDwell)
%ASP_EXPORT_VECTORS Write RTL co-simulation stimulus and expected responses.
%
%   outDir = ASP_EXPORT_VECTORS(OUTDIR, CFG, NDWELL)
%
%   Writes, for NDWELL consecutive dwells of a spoofing scenario:
%
%     adc_ch<i>_i.txt / adc_ch<i>_q.txt   quantised ADC samples, signed
%                                         integers, one per line
%     cov_expected.txt                    raw 48-bit covariance accumulator
%                                         contents per dwell, as integers
%     eig_expected.txt                    eigenvalues per dwell
%     weights_expected.txt                beamformer weights, as the signed
%                                         integers the RTL must produce
%     beam_expected_i.txt / _q.txt        beamformer output samples
%     manifest.txt                        formats, scaling, and the exact
%                                         config that produced the vectors
%
%   WHY THIS MATTERS MORE THAN IT LOOKS
%   -----------------------------------
%   The reference model in matlab/core is written so that every quantity on
%   this list is reproducible in EXACT integer arithmetic: full-width
%   products, a 48-bit accumulator with no intermediate rounding, unitary
%   Jacobi rotations, and power-of-two weight scaling.  That is not an
%   aesthetic choice - it means the RTL can be verified by BIT-EXACT
%   comparison rather than by "close enough" tolerance checking.
%
%   Tolerance-based verification of a nulling receiver is close to useless,
%   because the failure this design is most exposed to - a rounding bias in
%   the covariance accumulator - produces an error of a few LSBs that any
%   sensible tolerance would pass, while capping null depth at -25 dB
%   (verify/test_fixedpoint.m).  A bit-exact comparison catches it on the
%   first vector.
%
%   The values written here are the RAW INTEGERS, not scaled reals, for the
%   same reason.

if nargin < 1 || isempty(outDir)
    outDir = fullfile(tempdir, 'asp_vectors');
end
if nargin < 2 || isempty(cfg)
    cfg = asp_config('fs', 4*1.023e6);
end
if nargin < 3 || isempty(nDwell)
    nDwell = 4;
end

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

scn = asp_scenario(cfg, 'seed', 20120926, 'durationMs', nDwell + 1);
K = cfg.K;
nAnt = scn.nAnt;

agc = [];
genState = [];
h = ones(nAnt,1)/sqrt(nAnt);
fCur = h;

adcI = zeros(nAnt, K*nDwell);
adcQ = zeros(nAnt, K*nDwell);
beamI = zeros(1, K*nDwell);
beamQ = zeros(1, K*nDwell);

fidCov = fopen(fullfile(outDir,'cov_expected.txt'),'w');
fidEig = fopen(fullfile(outDir,'eig_expected.txt'),'w');
fidW   = fopen(fullfile(outDir,'weights_expected.txt'),'w');

for d = 1:nDwell
    [x, genState] = asp_rx_generate(scn, d-1, K, genState);
    [xq, agc] = asp_agc_adc(x, cfg, agc);

    idx = (d-1)*K + (1:K);
    adcI(:,idx) = round(real(xq) * 2^cfg.fx.adc.fl);
    adcQ(:,idx) = round(imag(xq) * 2^cfg.fx.adc.fl);

    % Beamform with the CURRENT weights (causal), then update.
    [vb, bi] = asp_beamform(xq, fCur, cfg.fx);
    beamI(idx) = round(real(vb) * 2^cfg.fx.beamOut.fl);
    beamQ(idx) = round(imag(vb) * 2^cfg.fx.beamOut.fl);

    acc = fx_cov_accum(xq, cfg.fx.adc, cfg.fx.covAccWl);

    % Raw accumulator contents, upper triangle, exactly as the DSP48 P
    % registers hold them.
    fprintf(fidCov, '# dwell %d\n', d);
    for i = 1:nAnt
        for j = i:nAnt
            fprintf(fidCov, '%d %d %.0f %.0f\n', i-1, j-1, acc.re(i,j), acc.im(i,j));
        end
    end

    [y, dbg] = asp_ssv_from_cov(acc.R, 'evd', struct('refIdx', cfg.refElement, ...
        'rank', cfg.est.maxNullRank, 'jacobiSweeps', cfg.est.jacobiSweeps));
    fprintf(fidEig, '# dwell %d\n', d);
    fprintf(fidEig, '%.12e\n', dbg.lam);

    det = asp_detect(dbg.lam, K, cfg.est.detectThreshold, cfg.est.maxNullRank);
    if det.rank > 0
        fNext = asp_weights('project', y(:,1:min(det.rank,size(y,2))), h);
    else
        fNext = h;
    end
    fCur = fNext;

    fq = fx_quant(fCur, cfg.fx.weight);
    fprintf(fidW, '# dwell %d  detected=%d rank=%d stat=%.6f\n', ...
        d, det.detected, det.rank, det.stat);
    for i = 1:nAnt
        fprintf(fidW, '%.0f %.0f\n', ...
            round(real(fq(i))*2^cfg.fx.weight.fl), ...
            round(imag(fq(i))*2^cfg.fx.weight.fl));
    end
end

fclose(fidCov); fclose(fidEig); fclose(fidW);

for i = 1:nAnt
    writeCol(fullfile(outDir, sprintf('adc_ch%d_i.txt', i-1)), adcI(i,:));
    writeCol(fullfile(outDir, sprintf('adc_ch%d_q.txt', i-1)), adcQ(i,:));
end
writeCol(fullfile(outDir,'beam_expected_i.txt'), beamI);
writeCol(fullfile(outDir,'beam_expected_q.txt'), beamQ);

fid = fopen(fullfile(outDir,'manifest.txt'),'w');
fprintf(fid, 'GNSS anti-spoofing array processor - RTL co-simulation vectors\n');
fprintf(fid, 'generated %s\n\n', datestr(now, 31));
fprintf(fid, 'fs                 %.6f MHz\n', cfg.fs/1e6);
fprintf(fid, 'samples per dwell  %d\n', K);
fprintf(fid, 'dwells             %d\n', nDwell);
fprintf(fid, 'antennas           %d  (%s, spacing %.3f lambda)\n', ...
    nAnt, cfg.geometry, cfg.elementSpacingLambda);
fprintf(fid, '\nFIXED-POINT FORMATS (signed two''s complement)\n');
fprintf(fid, '  adc      s%d.%d   scale 2^%d\n', cfg.fx.adc.wl, cfg.fx.adc.fl, cfg.fx.adc.fl);
fprintf(fid, '  weight   s%d.%d   scale 2^%d\n', cfg.fx.weight.wl, cfg.fx.weight.fl, cfg.fx.weight.fl);
fprintf(fid, '  beamOut  s%d.%d   scale 2^%d\n', cfg.fx.beamOut.wl, cfg.fx.beamOut.fl, cfg.fx.beamOut.fl);
fprintf(fid, '  cov accumulator  s%d, EXACT products, NO rounding in the loop\n', cfg.fx.covAccWl);
fprintf(fid, '\nSCENARIO\n');
fprintf(fid, '  authentic C/N0     %.1f dB-Hz x %d PRN\n', cfg.authCN0dBHz, cfg.nAuth);
fprintf(fid, '  SAPR               %.1f dB x %d PRN\n', cfg.saprDB, cfg.nSpoof);
fprintf(fid, '  spoofer az/el      %.1f / %.1f deg\n', cfg.spoofAzEl(1), cfg.spoofAzEl(2));
fprintf(fid, '  detector threshold %.6f (derived, Pfa = %.0e)\n', ...
    cfg.est.detectThreshold, cfg.est.pfa);
fprintf(fid, '\nVERIFICATION NOTE\n');
fprintf(fid, '  cov_expected.txt holds RAW accumulator integers.  Compare BIT-EXACT.\n');
fprintf(fid, '  A rounding bias of a few LSB here caps null depth at -25 dB while\n');
fprintf(fid, '  passing any tolerance-based check.  See verify/test_fixedpoint.m.\n');
fclose(fid);

fprintf('Wrote RTL vectors to %s\n', outDir);

end

% -------------------------------------------------------------------------
function writeCol(fname, v)
fid = fopen(fname, 'w');
fprintf(fid, '%.0f\n', v);
fclose(fid);
end
