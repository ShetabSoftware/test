function sig = generateFSK(symbolRate, fs, numSymbols, varargin)
%GENERATEFSK  Generate a continuous-phase M-FSK signal at complex baseband.
%   SIG = GENERATEFSK(SYMBOLRATE, FS, NUMSYMBOLS, ...) returns a struct with the
%   generated signal and its ground-truth parameters. Mainly intended for
%   testing and demonstrating the blind symbol-rate estimators, but also a
%   usable standalone CPFSK modulator.
%
%   Name-value options:
%     'order'           - modulation order M (number of tones), default 2.
%     'modulationIndex' - CPFSK index h; tone spacing = h*SYMBOLRATE unless
%                         'freqSeparation' is given. Default 1.0.
%     'freqSeparation'  - explicit tone spacing (Hz); overrides the index.
%     'carrierOffset'   - residual carrier offset (Hz), default 0.
%     'snrDb'           - if set, add complex AWGN at this SNR (dB).
%     'symbols'         - explicit symbol sequence (values in 0..M-1).
%     'toneFrequencies' - explicit tone frequencies (Hz); overrides computed.
%     'amplitude'       - linear amplitude of the noise-free signal, default 1.
%     'seed'            - RNG seed for symbols and noise.
%
%   SIG fields: iq, fs, symbolRate, symbols, toneFrequencies, carrierOffset,
%   snrDb, samplesPerSymbol.

    opts = parseOptions(struct( ...
        'order', 2, 'modulationIndex', 1.0, 'freqSeparation', [], ...
        'carrierOffset', 0.0, 'snrDb', [], 'symbols', [], ...
        'toneFrequencies', [], 'amplitude', 1.0, 'seed', []), varargin);

    if symbolRate <= 0, error('generateFSK:rs', 'symbolRate must be positive'); end
    if fs <= 0, error('generateFSK:fs', 'fs must be positive'); end
    if numSymbols < 1, error('generateFSK:ns', 'numSymbols must be >= 1'); end

    if ~isempty(opts.seed)
        try
            rng(opts.seed);
        catch
            rand('state', opts.seed);   %#ok<RAND>
            randn('state', opts.seed);  %#ok<RAND>
        end
    end

    order = opts.order;
    if ~isempty(opts.toneFrequencies)
        tones = opts.toneFrequencies(:).';
        order = numel(tones);
    else
        if order < 2
            error('generateFSK:order', 'FSK order must be >= 2');
        end
        freqSep = opts.freqSeparation;
        if isempty(freqSep)
            freqSep = opts.modulationIndex * symbolRate;
        end
        if freqSep <= 0
            error('generateFSK:sep', 'tone separation must be positive');
        end
        idx = (0:order - 1) - (order - 1) / 2.0;
        tones = idx * freqSep;
    end

    if ~isempty(opts.symbols)
        sym = opts.symbols(:).';
        numSymbols = numel(sym);
        if min(sym) < 0 || max(sym) >= order
            error('generateFSK:sym', 'symbol values must lie in 0..order-1');
        end
    else
        sym = randi(order, 1, numSymbols) - 1;   % values 0..order-1
    end

    sps = fs / symbolRate;
    totalSamples = round(numSymbols * sps);
    if totalSamples < 2
        error('generateFSK:len', ...
            'signal is too short; increase numSymbols or fs');
    end

    sampleIdx = 0:totalSamples - 1;
    symbolOfSample = min(floor(sampleIdx / sps), numSymbols - 1);  % 0-based
    instFreq = tones(sym(symbolOfSample + 1) + 1) + opts.carrierOffset;

    phase = 2 * pi * cumsum(instFreq) / fs;
    iq = opts.amplitude * exp(1i * phase);

    if ~isempty(opts.snrDb)
        iq = addAWGN(iq, opts.snrDb);
    end

    sig = struct( ...
        'iq', iq, 'fs', fs, 'symbolRate', symbolRate, ...
        'symbols', sym, 'toneFrequencies', tones, ...
        'carrierOffset', opts.carrierOffset, 'snrDb', opts.snrDb, ...
        'samplesPerSymbol', sps);
end
