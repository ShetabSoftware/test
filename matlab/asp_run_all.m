function ok = asp_run_all(mode)
%ASP_RUN_ALL Regression suite and studies for the anti-spoofing reference model.
%
%   ok = ASP_RUN_ALL()          % regression tests + short studies (~3 min)
%   ok = ASP_RUN_ALL('tests')   % regression tests only (~30 s)
%   ok = ASP_RUN_ALL('full')    % + the long Monte Carlo studies (~30 min)
%
%   Returns true if every regression test passed.  The studies are
%   measurements, not pass/fail, so they do not affect the return value.

if nargin < 1 || isempty(mode)
    mode = 'short';
end

t0 = tic;
ok = true;

fprintf('\n');
fprintf('#############################################################\n');
fprintf('#  GNSS ANTI-SPOOFING ARRAY PROCESSOR - REFERENCE MODEL      #\n');
fprintf('#  regression suite, mode = %-32s#\n', mode);
fprintf('#############################################################\n\n');

tests = {@test_ca_code, @test_evd_herm, @test_fixedpoint, ...
         @test_wishart_draw, @test_cov_model, @test_pipeline};

for k = 1:numel(tests)
    try
        ok = tests{k}() && ok;
    catch err
        fprintf('  ERROR in %s: %s\n', func2str(tests{k}), err.message);
        ok = false;
    end
    fprintf('\n');
end

fprintf('=============================================================\n');
if ok
    fprintf(' ALL REGRESSION TESTS PASSED   (%.1f s)\n', toc(t0));
else
    fprintf(' *** REGRESSION FAILURES ***   (%.1f s)\n', toc(t0));
end
fprintf('=============================================================\n');

if strcmpi(mode, 'tests')
    return;
end

fprintf('\n\n');
fprintf('#############################################################\n');
fprintf('#  STUDIES                                                   #\n');
fprintf('#############################################################\n');

if strcmpi(mode, 'full')
    study_array_size(1500);
    study_estimators(40);
    study_detector(600);
    study_postcorr(20);
    study_multipath(20);
else
    study_array_size(200);
    study_estimators(8);
    study_detector(120);
    study_postcorr(5);
    study_multipath(4);
end

fprintf('\nTotal elapsed: %.1f s\n', toc(t0));

end
