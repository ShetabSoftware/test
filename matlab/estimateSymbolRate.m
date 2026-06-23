function est = estimateSymbolRate(iq, fs, varargin)
%ESTIMATESYMBOLRATE  Blind symbol (baud) rate estimation for FSK signals.
%   EST = ESTIMATESYMBOLRATE(IQ, FS) estimates the symbol rate (Hz) of the
%   complex baseband FSK signal IQ sampled at FS hertz, with no prior knowledge
%   of the modulation parameters (tone count, spacing, modulation index or
%   carrier offset).
%
%   EST = ESTIMATESYMBOLRATE(IQ, FS, 'method', M, ...) selects the estimator:
%     'spectral' - frequency-domain cyclostationary line (default; accurate).
%     'autocorr' - time-domain periodicity of the transition signal.
%     'auto'     - run both and combine / cross-check.
%
%   Additional name-value options are forwarded to the underlying estimators,
%   e.g. 'rsMin', 'rsMax', 'transform' ('square'|'abs'), 'window', 'bandlimit'.
%
%   EST is a struct with fields:
%     symbolRate   - estimated symbol rate in Hz
%     method       - estimator that produced the result
%     confidence   - peak-to-floor ratio of the detected feature
%     diagnostics  - struct of intermediate quantities
%
%   Example:
%     sig = generateFSK(9600, 192000, 6000, 'order', 2, 'snrDb', 10, 'seed', 0);
%     est = estimateSymbolRate(sig.iq, 192000, 'method', 'auto');
%     fprintf('Rs = %.1f Hz\n', est.symbolRate);

    opts = parseOptions(struct( ...
        'method', 'spectral', 'rsMin', [], 'rsMax', [], ...
        'transform', 'square', 'window', 'hann', 'bandlimit', true, ...
        'maxHarmonics', 8, 'familyTol', 0.9), varargin);

    specArgs = {'rsMin', opts.rsMin, 'rsMax', opts.rsMax, ...
                'transform', opts.transform, 'window', opts.window, ...
                'bandlimit', opts.bandlimit};
    acArgs = {'rsMin', opts.rsMin, 'rsMax', opts.rsMax, ...
              'transform', opts.transform, 'bandlimit', opts.bandlimit, ...
              'maxHarmonics', opts.maxHarmonics, 'familyTol', opts.familyTol};

    switch lower(opts.method)
        case 'spectral'
            est = estimateSymbolRateSpectral(iq, fs, specArgs{:});
        case 'autocorr'
            est = estimateSymbolRateAutocorr(iq, fs, acArgs{:});
        case 'auto'
            spec = estimateSymbolRateSpectral(iq, fs, specArgs{:});
            ac = estimateSymbolRateAutocorr(iq, fs, acArgs{:});
            relDiff = abs(spec.symbolRate - ac.symbolRate) / ...
                      max(spec.symbolRate, 1e-30);
            if relDiff < 0.02
                est = struct( ...
                    'symbolRate', 0.5 * (spec.symbolRate + ac.symbolRate), ...
                    'method', 'auto(agree)', ...
                    'confidence', spec.confidence + ac.confidence, ...
                    'diagnostics', struct('spectral', spec, 'autocorr', ac));
            else
                if spec.confidence >= ac.confidence
                    chosen = spec;
                else
                    chosen = ac;
                end
                est = struct( ...
                    'symbolRate', chosen.symbolRate, ...
                    'method', ['auto(' chosen.method ')'], ...
                    'confidence', chosen.confidence, ...
                    'diagnostics', struct('spectral', spec, 'autocorr', ac));
            end
        otherwise
            error('estimateSymbolRate:method', ...
                'unknown method: %s', opts.method);
    end
end
