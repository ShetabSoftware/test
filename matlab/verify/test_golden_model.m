function ok = test_golden_model()
%TEST_GOLDEN_MODEL Interface and calibration regression for the VHDL golden model.
%
%   Two classes of defect are covered, and they fail in opposite ways.
%
%   INPUT INTERFACE.  The golden model is the reference the RTL is diffed
%   against, so a run driven from captured data must produce bit-identical
%   results to the same samples supplied in memory.  If the file format and
%   the in-memory format ever disagree, every subsequent RTL mismatch gets
%   blamed on the RTL.  A mistyped option name is checked here too: silently
%   falling back to a default is the worst way for a verification run to be
%   wrong, because the run still completes and prints a plausible number.
%
%   DETECTOR CALIBRATION.  The clean-sky case is the one that fails
%   silently.  A missed detection shows up as a bad null depth; a FALSE
%   ALARM shows up as nothing at all, while the array steers a null into an
%   authentic satellite.  The threshold must therefore be calibrated against
%   the REALISTIC null hypothesis - thermal noise, plus the sample
%   correlation the shaping FIR imposes, plus the authentic constellation,
%   which is itself a structured non-white term in the spatial covariance.
%   Calibrating on white noise puts the threshold below the H0 median and
%   yields a false alarm rate near 60%, which is what this test caught.

ok = true;
fprintf('test_golden_model:\n');

% ------------------------------------------------------------------ 1
% Option-name validation.
bad = { {'duratonMs',2},        'typo in an option name'
        {'durationMs'},         'odd number of arguments'
        {'adc', zeros(3,80000)},'wrong channel count'
        {'adc', zeros(4,1000)}, 'input shorter than one dwell'
        {'adcFile','/nonexistent/nope.txt'}, 'missing input file' };
nCaught = 0;
for k = 1:size(bad,1)
    try
        asp_golden_model(bad{k,1}{:}, 'dump',false, 'verbose',false);
    catch
        nCaught = nCaught + 1;
    end
end
ok = report(ok, nCaught == size(bad,1), ...
    sprintf('bad arguments rejected, not silently defaulted (%d of %d)', ...
    nCaught, size(bad,1)));

% ------------------------------------------------------------------ 2
% The three input paths must agree BIT FOR BIT.
tmp = tempname; mkdir(tmp);
G1 = asp_golden_model('durationMs',2, 'outDir',tmp,  'verbose',false);
G2 = asp_golden_model('adc',G1.adc,   'dump',false,  'verbose',false);
G3 = asp_golden_model('adcFile',fullfile(tmp,'s00_adc_out.txt'), ...
                                      'dump',false,  'verbose',false);

ok = report(ok, isequal(G3.adc, G1.adc), ...
    'vector file round-trips to the same ADC samples');
ok = report(ok, isequal(G2.dac, G1.dac) && isequal(G2.beam, G1.beam), ...
    'in-memory input reproduces the generated run bit-exactly');
ok = report(ok, isequal(G3.dac, G1.dac) && isequal(G3.beam, G1.beam), ...
    'file-driven input reproduces the generated run bit-exactly');
ok = report(ok, ~isfield(G2,'check') && ~isfield(G3,'check'), ...
    'self-check is withheld when there is no ground truth');

G1b = asp_golden_model('durationMs',2, 'dump',false, 'verbose',false);
ok = report(ok, isequal(G1b.dac, G1.dac), ...
    'same seed reproduces the same run (no hidden state)');

% ------------------------------------------------------------------ 3
% Detector calibration against the realistic H0.
%
% Deliberately checked through the FULL chain rather than on synthetic
% covariances: the two effects that move the threshold - FIR sample
% correlation and the authentic constellation - both live in the chain and
% neither appears in a white-noise model.
nSeed = 8;
h0 = []; h1 = [];
for s = 1:nSeed
    A = asp_golden_model('durationMs',3, 'spoof',false, 'seed',9000+s, ...
        'dump',false, 'verbose',false);
    B = asp_golden_model('durationMs',3, 'spoof',true,  'seed',9000+s, ...
        'dump',false, 'verbose',false);
    h0 = [h0, [A.dwell.det]];   %#ok<AGROW>
    h1 = [h1, [B.dwell.det]];   %#ok<AGROW>
end
statOf = @(d) arrayfun(@(z) z.lam(1)*3/sum(z.lam(2:4)), d);
s0 = statOf(h0); s1 = statOf(h1);
fa = mean([h0.detected]); pd = mean([h1.detected]);

fprintf('        H0 statistic  mean %.4f  max %.4f   (%d dwells)\n', ...
    mean(s0), max(s0), numel(s0));
fprintf('        H1 statistic  min  %.4f            (%d dwells)\n', ...
    min(s1), numel(s1));
fprintf('        threshold     %.4f\n', 1280/1024);

ok = report(ok, fa == 0, ...
    sprintf('CLEAN SKY: no false alarm in %d dwells (%.0f%%)', numel(s0), 100*fa));
ok = report(ok, pd == 1, ...
    sprintf('spoofer at SAPR +5.5 dB: detected on %.0f%% of dwells', 100*pd));

% The margin, not just the outcome.  A threshold that merely happens to sit
% between two samples is not calibrated; it has to clear the H0 spread.
% The design figure is 5.2 sigma, measured over 160 dwells; this test runs a
% few dozen, where the sigma estimate is itself noisy, so it asserts the
% weaker bound that still fails loudly if the threshold slips back under the
% H0 distribution.
marginSigma = (1280/1024 - mean(s0)) / std(s0);
fprintf('        threshold sits %.1f sigma above the H0 mean\n', marginSigma);
ok = report(ok, marginSigma > 3, ...
    sprintf('threshold clears the realistic H0 by >3 sigma (%.1f)', marginSigma));
ok = report(ok, min(s1) > 1280/1024, ...
    sprintf('H1 minimum (%.3f) stays above the threshold', min(s1)));

% ------------------------------------------------------------------ 4
% Sensitivity floor: the useful range of the PRE-correlation detector.
sapr0 = asp_golden_model('durationMs',3, 'spoof',true, 'saprDB',0, ...
    'seed',424242, 'dump',false, 'verbose',false);
ok = report(ok, all(arrayfun(@(z) z.det.detected, sapr0.dwell)), ...
    'still detects a spoofer at SAPR 0 dB');

rmdir(tmp,'s');
end


function ok = report(ok, cond, msg)
if cond, fprintf('  PASS  %s\n', msg);
else,    fprintf('  FAIL  %s\n', msg); ok = false; end
end
