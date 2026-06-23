function runTests()
%RUNTESTS  Self-checking test suite for the MATLAB FSK baud-rate estimators.
%   Runs a set of assertions across modulation orders, symbol rates, carrier
%   offsets and SNRs and reports pass/fail counts. Compatible with both MATLAB
%   and GNU Octave.

    npass = 0;
    nfail = 0;

    relerr = @(est, truth) abs(est - truth) / truth;

    % --- DSP helper checks --------------------------------------------------
    fs = 1000; n = (0:4095); f0 = 123;
    x = exp(1i * 2 * pi * f0 * n / fs);
    [npass, nfail] = check(abs(mean(instantaneousFrequency(x, fs)) - f0) < 1e-6, ...
        'instantaneousFrequency pure tone', npass, nfail);

    yv = [0 1 4 1 0];
    [idx, val] = quadraticPeakInterp(yv, 3);
    [npass, nfail] = check(abs(idx - 3) < 1e-9 && abs(val - 4) < 1e-9, ...
        'quadraticPeakInterp centred', npass, nfail);

    xp = (0:5); yp = -((xp - 2.3) .^ 2) + 10;
    [~, k] = max(yp);
    [idx, val] = quadraticPeakInterp(yp, k);
    [npass, nfail] = check(abs((idx - 1) - 2.3) < 1e-6 && abs(val - 10) < 1e-6, ...
        'quadraticPeakInterp offset', npass, nfail);

    % --- Generator sanity ---------------------------------------------------
    sig = generateFSK(1000, 20000, 200, 'seed', 1);
    [npass, nfail] = check(abs(sig.samplesPerSymbol - 20) < 1e-9, ...
        'generateFSK samplesPerSymbol', npass, nfail);
    [npass, nfail] = check(max(abs(abs(sig.iq) - 1)) < 1e-9, ...
        'generateFSK constant envelope', npass, nfail);

    % --- Noise-free accuracy across rates and methods -----------------------
    fs = 192000;
    for rs = [800 2400 4800 9600]
        sig = generateFSK(rs, fs, 4000, 'order', 2, 'seed', 7);
        for m = {'spectral', 'autocorr', 'auto'}
            est = estimateSymbolRate(sig.iq, fs, 'method', m{1});
            [npass, nfail] = check(relerr(est.symbolRate, rs) < 0.01, ...
                sprintf('clean rs=%d method=%s', rs, m{1}), npass, nfail);
        end
    end

    % --- 4-FSK --------------------------------------------------------------
    sig = generateFSK(5000, 100000, 4000, 'order', 4, 'seed', 3);
    est = estimateSymbolRate(sig.iq, 100000, 'method', 'spectral');
    [npass, nfail] = check(relerr(est.symbolRate, 5000) < 0.02, '4-FSK spectral', ...
        npass, nfail);

    % --- Carrier offset robustness -----------------------------------------
    sig = generateFSK(4000, 100000, 4000, 'order', 2, ...
        'carrierOffset', 7777, 'seed', 11);
    est = estimateSymbolRate(sig.iq, 100000, 'method', 'auto');
    [npass, nfail] = check(relerr(est.symbolRate, 4000) < 0.02, ...
        'carrier offset', npass, nfail);

    % --- Noise robustness ---------------------------------------------------
    fs = 192000;
    for rs = [800 2400 4800 9600]
        for snr = [20 10 8]
            sig = generateFSK(rs, fs, 6000, 'order', 2, 'snrDb', snr, 'seed', 42);
            est = estimateSymbolRate(sig.iq, fs, 'method', 'auto');
            [npass, nfail] = check(relerr(est.symbolRate, rs) < 0.03, ...
                sprintf('rs=%d snr=%ddB', rs, snr), npass, nfail);
        end
    end

    % --- Occupied-band detection -------------------------------------------
    sig = generateFSK(4000, 200000, 4000, 'order', 2, 'snrDb', 15, 'seed', 8);
    [center, bw] = detectOccupiedBand(sig.iq, 200000);
    [npass, nfail] = check(abs(center) < 4000 && bw > 4000 && bw < 12 * 4000, ...
        'detectOccupiedBand', npass, nfail);

    fprintf('\n%d passed, %d failed\n', npass, nfail);
    if nfail > 0
        error('runTests:failures', '%d test(s) failed', nfail);
    end
end

function [npass, nfail] = check(cond, name, npass, nfail)
    if cond
        npass = npass + 1;
        fprintf('  ok   %s\n', name);
    else
        nfail = nfail + 1;
        fprintf('  FAIL %s\n', name);
    end
end
