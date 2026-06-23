function est = estimateSymbolRateSpectral(iq, fs, varargin)
%ESTIMATESYMBOLRATESPECTRAL  Spectral-line blind FSK symbol-rate estimate.
%   EST = ESTIMATESYMBOLRATESPECTRAL(IQ, FS, ...) estimates the symbol rate of
%   the complex baseband FSK signal IQ (sample rate FS, Hz) by locating the
%   fundamental cyclostationary line of the transition signal in a long
%   periodogram.
%
%   Name-value options:
%     'rsMin', 'rsMax'  - symbol-rate search band (Hz); empty => auto defaults.
%     'nfft'            - FFT length; empty => next power of two >= length.
%     'transform'       - 'square' (default) or 'abs'.
%     'window'          - 'hann' (default), 'hamming', 'blackman' or 'none'.
%     'bandlimit'       - true (default) to band-pass to the occupied band.
%
%   EST is a struct with fields: symbolRate, method, confidence, diagnostics.

    opts = parseOptions(struct( ...
        'rsMin', [], 'rsMax', [], 'nfft', [], ...
        'transform', 'square', 'window', 'hann', 'bandlimit', true), varargin);

    if fs <= 0
        error('estimateSymbolRateSpectral:fs', 'fs must be positive');
    end

    [y, bandCenter, bandBw] = preprocessSignal(iq, fs, opts.bandlimit);

    trans = transitionSignal(y, fs, opts.transform);
    n = numel(trans);
    if n < 8
        error('estimateSymbolRateSpectral:short', ...
            'signal is too short for spectral estimation');
    end

    [rsMin, rsMax] = resolveSearchBand(fs, n, opts.rsMin, opts.rsMax);

    w = analysisWindow(opts.window, n);
    if isempty(opts.nfft)
        nfft = 2 ^ nextpow2(n);
    else
        nfft = opts.nfft;
    end

    specFull = abs(fft(trans .* w, nfft));
    spectrum = specFull(1:floor(nfft / 2) + 1);

    jsel = selectFundamental(spectrum, fs, nfft, rsMin, rsMax);
    interpIdx = quadraticPeakInterp(spectrum, jsel);
    symbolRate = (interpIdx - 1) * fs / nfft;

    df = fs / nfft;
    freqs = (0:numel(spectrum) - 1) * df;
    inBand = (freqs >= rsMin) & (freqs <= rsMax);
    med = median(spectrum(inBand));
    if med == 0
        med = 1e-30;
    end
    confidence = spectrum(jsel) / med;

    diagnostics = struct( ...
        'freqs', freqs, 'spectrum', spectrum, ...
        'searchBand', [rsMin, rsMax], 'peakIndex', jsel, 'nfft', nfft, ...
        'occupiedBand', [bandCenter, bandBw]);
    est = struct('symbolRate', symbolRate, 'method', 'spectral', ...
                 'confidence', confidence, 'diagnostics', diagnostics);
end
