function jsel = selectFundamental(spectrum, fs, nfft, rsMin, rsMax, ...
                                  peakFrac, maxHarmonics, familyTol)
%SELECTFUNDAMENTAL  FFT bin (1-based) of the fundamental symbol-rate line.
%   JSEL = SELECTFUNDAMENTAL(SPECTRUM, FS, NFFT, RSMIN, RSMAX) returns the
%   1-based index into SPECTRUM (a one-sided magnitude spectrum, element J at
%   frequency (J-1)*FS/NFFT) of the fundamental cyclostationary line.
%
%   The transition pulse train produces lines at every integer multiple of the
%   symbol rate, frequently with comparable amplitude, so the global maximum is
%   often a harmonic. Each candidate peak is scored by the mean energy of its
%   harmonic comb (noise floor removed); the fundamental and every clean
%   harmonic attain a similar score while spurious peaks and sub-harmonics score
%   poorly. The fundamental is the lowest-frequency member of the tied family.

    if nargin < 6 || isempty(peakFrac),     peakFrac = 0.2;     end
    if nargin < 7 || isempty(maxHarmonics), maxHarmonics = 8;   end
    if nargin < 8 || isempty(familyTol),    familyTol = 0.9;    end

    spectrum = spectrum(:).';
    N = numel(spectrum);
    df = fs / nfft;

    % Work in 0-based bin space (bin b is SPECTRUM(b+1)), mirroring the
    % reference implementation, to keep the harmonic arithmetic exact.
    kMin = max(ceil(rsMin / df), 1);
    kMax = min(floor(rsMax / df), N - 1);
    if kMax <= kMin
        error('selectFundamental:band', ...
            'search band contains no FFT bins; widen the band');
    end

    bandBins = kMin:kMax;                    % 0-based
    band = spectrum(bandBins + 1);
    noiseFloor = median(band);
    positive = max(spectrum - noiseFloor, 0);   % 1-based array

    bandPos = positive(bandBins + 1);
    bandMax = max(bandPos);
    if bandMax <= 0
        [~, mi] = max(band);
        jsel = bandBins(mi) + 1;
        return;
    end

    relMax = localMaxima(bandPos);              % indices into bandPos
    candBins = [];                              % 0-based bins
    for ii = relMax
        if bandPos(ii) >= peakFrac * bandMax
            candBins(end + 1) = kMin + (ii - 1); %#ok<AGROW>
        end
    end
    if isempty(candBins)
        [~, mi] = max(bandPos);
        candBins = kMin + (mi - 1);
    end

    scores = zeros(1, numel(candBins));
    for ci = 1:numel(candBins)
        b0 = candBins(ci);
        f0frac1 = quadraticPeakInterp(positive, b0 + 1);  % 1-based fractional
        f0bin = f0frac1 - 1;                              % 0-based fractional
        nHarm = min(maxHarmonics, floor(kMax / f0bin));
        if nHarm < 2
            scores(ci) = 0;
            continue;
        end
        total = 0.0;
        for k = 1:nHarm
            center0 = round(k * f0bin);
            w = max(2, round(0.004 * center0));
            lo0 = max(center0 - w, 0);
            hi0 = min(center0 + w, N - 1);
            total = total + max(positive((lo0:hi0) + 1));
        end
        scores(ci) = total / nHarm;
    end

    best = max(scores);
    if best <= 0
        [~, mi] = max(positive(candBins + 1));
        jsel = candBins(mi) + 1;
        return;
    end
    family = candBins(scores >= familyTol * best);
    jsel = min(family) + 1;     % lowest-frequency member -> fundamental
end
