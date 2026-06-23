function [center, bw] = detectOccupiedBand(iq, fs, thresholdMult)
%DETECTOCCUPIEDBAND  Estimate the centre frequency and width of the FSK band.
%   [CENTER, BW] = DETECTOCCUPIEDBAND(IQ, FS) returns the centre frequency and
%   bandwidth (Hz) of the occupied band. A Welch power spectrum is thresholded
%   at THRESHOLDMULT (default 4) times a robust noise-floor estimate (the 25th
%   percentile of the PSD). Even at a few dB of wideband SNR the in-band PSD
%   sits well above the floor for oversampled FSK, making this reliable.

    if nargin < 3 || isempty(thresholdMult)
        thresholdMult = 4.0;
    end
    x = iq(:).';
    n = numel(x);
    nperseg = min(n, max(256, floor(n / 32)));
    nperseg = max(nperseg, 16);
    noverlap = floor(nperseg / 2);
    step = max(nperseg - noverlap, 1);
    w = analysisWindow('hann', nperseg);
    winNorm = sum(w .^ 2);
    if winNorm == 0
        winNorm = 1;
    end

    accum = zeros(1, nperseg);
    count = 0;
    for s = 1:step:(n - nperseg + 1)
        seg = x(s:s + nperseg - 1);
        accum = accum + abs(fft(seg .* w)) .^ 2 / (fs * winNorm);
        count = count + 1;
    end
    if count == 0
        accum = abs(fft(x, nperseg)) .^ 2 / (fs * winNorm);
        count = 1;
    end
    psd = fftshift(accum / count);

    half = floor(nperseg / 2);
    freqs = ((-half):(nperseg - 1 - half)) / nperseg * fs;

    floorVal = percentileValue(psd, 25);
    if floorVal == 0
        floorVal = 1e-30;
    end
    mask = psd > thresholdMult * floorVal;
    if ~any(mask)
        center = 0.0;
        bw = fs;
        return;
    end
    idx = find(mask);
    fLo = freqs(idx(1));
    fHi = freqs(idx(end));
    center = 0.5 * (fLo + fHi);
    bw = fHi - fLo;
end
