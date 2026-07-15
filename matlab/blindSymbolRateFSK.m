function [Rs_hat, D] = blindSymbolRateFSK(rx, Fs, Rs_min, Rs_max, varargin)
%BLINDSYMBOLRATEFSK  Blind symbol-rate estimation for 2-FSK (HF-channel ready).
%
%   RS_HAT = BLINDSYMBOLRATEFSK(RX, FS, RS_MIN, RS_MAX) estimates the symbol
%   (baud) rate of a complex-baseband 2-FSK signal RX sampled at FS hertz,
%   searching the band [RS_MIN, RS_MAX] Hz, with NO knowledge of the tone
%   frequencies, tone spacing (Delta f) or symbol rate.
%
%   [RS_HAT, D] = BLINDSYMBOLRATEFSK(...) also returns a diagnostics struct D
%   (instantaneous frequency, estimated tone midpoint and spacing, the notched
%   transition spectrum over the search band, and the chosen bin).
%
%   Algorithm (IF hard-limit -> transition clock FFT)
%   -------------------------------------------------
%     1. Instantaneous frequency (IF) via a phase discriminator, then a short
%        moving-average to suppress noise jitter.
%     2. Blind two-level decision: the IF histogram is bimodal (one cluster per
%        tone). A robust midpoint hard-limits the IF to a binary switching
%        waveform. This is amplitude-invariant, so it is immune to Watterson /
%        Rayleigh fading -- only the symbol-edge *timing* matters.
%     3. transitions = |diff(binary)| : unit pulses at the symbol edges, whose
%        fundamental spectral line is the symbol rate.
%     4. FFT of the transition signal over [Rs_min, Rs_max]. On a multipath
%        channel the IF also *beats* at the tone spacing Delta f (when the
%        direct and delayed paths carry different tones), producing a strong
%        spurious line at Delta f and its harmonics -- these are notched using a
%        blind estimate of Delta f from the two IF clusters.
%     5. Peak pick, a sub-harmonic preference test (so a harmonic is not taken
%        for the clock), and parabolic interpolation for sub-bin accuracy.
%
%   Name-value options (defaults match a tuned HF 2-FSK setup)
%     'Nfft'             - FFT length (default 2^18 -> very fine bin spacing).
%     'smoothDiv'        - IF moving-average length = max(5, round(Fs/
%                          (smoothDiv*Rs_nom))), Rs_nom = sqrt(Rs_min*Rs_max).
%                          Default 2.
%     'notchToneSpacing' - true (default): notch Delta f and its harmonics.
%                          IMPORTANT: if the symbol rate can EQUAL the tone
%                          spacing (here Rs = Delta f = 300 Hz), set this false
%                          for that case -- otherwise the true clock line is
%                          notched. Telling a 300-baud clock from a 300-Hz tone
%                          beat is an inherent blind ambiguity.
%     'notchHz'          - half-width of each notch, Hz (default 25).
%     'notchHarmonics'   - number of Delta f harmonics to notch (default 3).
%     'subharmThresh'    - prefer sub-harmonic Rs/h if its line exceeds this
%                          fraction of the peak (default 0.35).
%     'subharmMax'       - highest sub-harmonic divisor tested (default 4).
%
%   No toolboxes required; runs under MATLAB and GNU Octave.
%
%   Example (drop-in for a hand-rolled estimation block):
%     Rs_estimated = blindSymbolRateFSK(rxChanNoisy, Fs, Rs_min, Rs_max);

    opt = struct('Nfft', 2^18, 'smoothDiv', 2, 'notchToneSpacing', true, ...
                 'notchHz', 25, 'notchHarmonics', 3, ...
                 'subharmThresh', 0.35, 'subharmMax', 4);
    opt = localParse(opt, varargin);

    rx = rx(:).';
    if numel(rx) < 64
        error('blindSymbolRateFSK:short', 'need at least 64 samples');
    end

    % --- 1) Instantaneous frequency + smoothing -----------------------------
    inst_phase = unwrap(angle(rx));
    inst_freq = diff(inst_phase) * Fs / (2 * pi);

    Rs_nom = sqrt(Rs_min * Rs_max);                 % blind nominal rate
    win_len = max(5, round(Fs / (opt.smoothDiv * Rs_nom)));
    inst_freq = localMovMean(inst_freq, win_len);

    % --- 2) Blind tone midpoint and spacing (Delta f) -----------------------
    f_lo = localPercentile(inst_freq, 20);
    f_hi = localPercentile(inst_freq, 80);
    mid_f = 0.5 * (f_lo + f_hi);
    delta_f = abs(f_hi - f_lo);

    % --- 3) Hard-limit and detect transitions -------------------------------
    bit_seq = double(inst_freq > mid_f);
    transitions = abs(diff(bit_seq));

    % --- 4) Transition-clock spectrum over the search band ------------------
    Nfft = opt.Nfft;
    Y = abs(fft(transitions - mean(transitions), Nfft));
    f_axis = (0:Nfft - 1) * Fs / Nfft;
    sel = (f_axis >= Rs_min) & (f_axis <= Rs_max);
    f_search = f_axis(sel);
    Y_search = Y(sel);

    % Notch the tone-spacing beat (Delta f) and its low harmonics.
    if opt.notchToneSpacing && delta_f > 0
        for h = 1:opt.notchHarmonics
            Y_search(abs(f_search - h * delta_f) < opt.notchHz) = 0;
        end
    end

    % --- 5) Peak, sub-harmonic preference, parabolic refine -----------------
    [peak_mag, iPk] = max(Y_search);
    Rs_hat = f_search(iPk);

    for h = 2:opt.subharmMax
        f_fund = f_search(iPk) / h;
        if f_fund >= Rs_min && f_fund <= Rs_max
            [~, jH] = min(abs(f_search - f_fund));
            if Y_search(jH) > opt.subharmThresh * peak_mag
                Rs_hat = f_fund;
            end
        end
    end

    [~, iR] = min(abs(f_search - Rs_hat));
    if iR > 1 && iR < numel(Y_search)
        a = Y_search(iR - 1); b = Y_search(iR); c = Y_search(iR + 1);
        den = a - 2 * b + c;
        if den ~= 0
            p = 0.5 * (a - c) / den;
            dfBin = f_search(2) - f_search(1);
            Rs_hat = f_search(iR) + p * dfBin;
        end
    end

    D = struct('inst_freq', inst_freq, 'mid_f', mid_f, 'delta_f', delta_f, ...
               'f_search', f_search, 'Y_search', Y_search, 'peakBin', iR, ...
               'peak_mag', peak_mag, 'win_len', win_len);
end

% =========================================================================
% Local helpers (no toolboxes -- work in MATLAB and Octave)
% =========================================================================
function y = localMovMean(x, w)
    if w <= 1, y = x; return; end
    n = numel(x);
    c = cumsum([0, x]);
    half = floor(w / 2);
    lo = max((1:n) - half, 1);
    hi = min((1:n) + half, n);
    y = (c(hi + 1) - c(lo)) ./ (hi - lo + 1);
end

function v = localPercentile(x, p)
    xs = sort(x(:));
    m = numel(xs);
    if m == 0, v = 0; return; end
    if m == 1, v = xs(1); return; end
    rank = (p / 100) * (m - 1) + 1;
    loi = max(min(floor(rank), m), 1);
    hii = max(min(ceil(rank), m), 1);
    frac = rank - loi;
    v = xs(loi) * (1 - frac) + xs(hii) * frac;
end

function opt = localParse(opt, args)
    if isempty(args), return; end
    if mod(numel(args), 2) ~= 0
        error('blindSymbolRateFSK:pairs', 'options must be name-value pairs');
    end
    for i = 1:2:numel(args)
        name = char(args{i});
        if ~isfield(opt, name)
            error('blindSymbolRateFSK:unknown', 'unknown option: %s', name);
        end
        opt.(name) = args{i + 1};
    end
end
