function est = estimateSymbolRateAutocorr(iq, fs, varargin)
%ESTIMATESYMBOLRATEAUTOCORR  Autocorrelation blind FSK symbol-rate estimate.
%   EST = ESTIMATESYMBOLRATEAUTOCORR(IQ, FS, ...) estimates the symbol rate from
%   the autocorrelation of the transition signal. The symbol period is the
%   smallest lag whose integer multiples are all autocorrelation peaks (which
%   resolves the period/multiple ambiguity and skips any band-pass correlation
%   main lobe); it is then refined against a far harmonic peak for fine
%   resolution.
%
%   Name-value options:
%     'rsMin', 'rsMax'  - symbol-rate search band (Hz); empty => auto defaults.
%     'transform'       - 'square' (default) or 'abs'.
%     'bandlimit'       - true (default) to band-pass to the occupied band.
%     'maxHarmonics'    - number of harmonics used in the comb (default 8).
%     'familyTol'       - tolerance for the tied-score family (default 0.9).
%
%   EST is a struct with fields: symbolRate, method, confidence, diagnostics.

    opts = parseOptions(struct( ...
        'rsMin', [], 'rsMax', [], 'transform', 'square', ...
        'bandlimit', true, 'maxHarmonics', 8, 'familyTol', 0.9), varargin);

    if fs <= 0
        error('estimateSymbolRateAutocorr:fs', 'fs must be positive');
    end

    [y, bandCenter, bandBw, effBw] = preprocessSignal(iq, fs, opts.bandlimit);

    trans = transitionSignal(y, fs, opts.transform);
    n = numel(trans);
    if n < 8
        error('estimateSymbolRateAutocorr:short', ...
            'signal is too short for autocorrelation estimation');
    end

    [rsMin, rsMax] = resolveSearchBand(fs, n, opts.rsMin, opts.rsMax);

    % Skip the correlation main lobe introduced by the band-pass filter; its
    % width is roughly the reciprocal of the pass bandwidth.
    if effBw < fs
        lagHump = ceil(1.3 * fs / effBw);
    else
        lagHump = 1;
    end
    lagMin = max([floor(fs / rsMax), lagHump, 2]);
    lagMax = min(ceil(fs / rsMin), floor(n / 3));
    if lagMax <= lagMin + 2
        error('estimateSymbolRateAutocorr:lag', ...
            'search band yields an empty autocorrelation lag range');
    end

    maxLag = min(n - 1, lagMax * opts.maxHarmonics);
    acf = autocorrelation(trans, maxLag);   % acf(tau+1) = autocorr at lag tau
    zeroLag = acf(1);
    if zeroLag == 0
        zeroLag = 1e-30;
    end
    a = max(acf / zeroLag, 0);

    scores = zeros(1, lagMax + 1);          % scores(tau+1)
    for tau = lagMin:lagMax
        m = min(opts.maxHarmonics, floor(maxLag / tau));
        if m < 2
            continue;
        end
        ks = 1:m;
        scores(tau + 1) = mean(a(ks * tau + 1));
    end
    best = max(scores);
    if best <= 0
        error('estimateSymbolRateAutocorr:periodic', ...
            'no periodic structure found in autocorrelation');
    end
    familyTaus = find(scores >= opts.familyTol * best) - 1;  % lags
    familyTaus = familyTaus(familyTaus >= lagMin);
    tau0 = familyTaus(1);

    % Fine refinement using the highest reliable harmonic peak.
    kmax = min(opts.maxHarmonics, floor(maxLag / tau0));
    refinedPeriod = tau0;
    for k = kmax:-1:1
        center = k * tau0;            % lag
        w = max(2, floor(tau0 / 3));
        lo = max(center - w, 1);
        hi = min(center + w, numel(acf) - 2);
        if hi <= lo
            continue;
        end
        seg = acf((lo:hi) + 1);
        [~, mi] = max(seg);
        localLag = lo + (mi - 1);
        interpIdx = quadraticPeakInterp(acf, localLag + 1);  % 1-based index
        interpLag = interpIdx - 1;
        refinedPeriod = interpLag / k;
        break;
    end

    symbolRate = fs / refinedPeriod;
    confidence = a(tau0 + 1);

    diagnostics = struct( ...
        'acf', acf, 'scores', scores, 'lagRange', [lagMin, lagMax], ...
        'coarsePeriod', tau0, 'refinedPeriod', refinedPeriod, ...
        'occupiedBand', [bandCenter, bandBw]);
    est = struct('symbolRate', symbolRate, 'method', 'autocorr', ...
                 'confidence', confidence, 'diagnostics', diagnostics);
end
