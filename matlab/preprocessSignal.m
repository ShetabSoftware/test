function [y, center, bw, effBw] = preprocessSignal(iq, fs, doBandlimit)
%PREPROCESSSIGNAL  Optionally band-limit a signal to its occupied band.
%   [Y, CENTER, BW, EFFBW] = PREPROCESSSIGNAL(IQ, FS, DOBANDLIMIT) returns the
%   (optionally filtered) signal Y, the detected band CENTER and width BW, and
%   the effective pass bandwidth EFFBW. If DOBANDLIMIT is false, or the detected
%   band already spans almost the whole spectrum, no filtering is applied.

    if ~doBandlimit
        y = iq(:).';
        center = 0.0;
        bw = fs;
        effBw = fs;
        return;
    end
    [center, bw] = detectOccupiedBand(iq, fs);
    if bw <= 0 || bw >= fs * 0.98
        y = iq(:).';
        effBw = fs;
        return;
    end
    [y, effBw] = bandlimitSignal(iq, fs, center, bw);
end
