function [y, effBw] = bandlimitSignal(iq, fs, center, bw, margin)
%BANDLIMITSIGNAL  Band-pass a signal to its occupied band (brick-wall FFT).
%   [Y, EFFBW] = BANDLIMITSIGNAL(IQ, FS, CENTER, BW) zeroes all spectral
%   content outside [CENTER-H, CENTER+H], where H = BW/2*MARGIN (MARGIN
%   default 1.5). EFFBW is the two-sided pass bandwidth actually applied, which
%   the autocorrelation estimator uses to skip the filter's correlation lobe.

    if nargin < 5 || isempty(margin)
        margin = 1.5;
    end
    x = iq(:).';
    n = numel(x);
    nfft = 2 ^ nextpow2(n);
    X = fftshift(fft(x, nfft));

    half = floor(nfft / 2);
    freqs = ((-half):(nfft - 1 - half)) / nfft * fs;

    h = max(bw / 2 * margin, fs / nfft * 8);
    h = min(h, fs / 2);
    outOfBand = (freqs < center - h) | (freqs > center + h);
    X(outOfBand) = 0;
    y = ifft(ifftshift(X));
    y = y(1:n);
    effBw = 2 * h;
end
