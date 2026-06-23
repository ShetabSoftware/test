function [rsMin, rsMax] = resolveSearchBand(fs, n, rsMin, rsMax, minSymbols)
%RESOLVESEARCHBAND  Determine the symbol-rate search band, applying defaults.
%   [RSMIN, RSMAX] = RESOLVESEARCHBAND(FS, N, RSMIN, RSMAX) fills in sensible
%   defaults for empty bounds: RSMAX defaults to FS/4 and RSMIN to the rate at
%   which MINSYMBOLS symbols fit in the N-sample record.

    if nargin < 5 || isempty(minSymbols)
        minSymbols = 16;
    end
    if isempty(rsMax)
        rsMax = fs / 4.0;
    end
    if isempty(rsMin)
        rsMin = max(minSymbols * fs / max(n, 1), fs / 1e6);
    end
    if rsMin <= 0
        error('resolveSearchBand:rsMin', 'rsMin must be positive');
    end
    if rsMax <= rsMin
        error('resolveSearchBand:order', 'rsMax must be greater than rsMin');
    end
    if rsMax >= fs / 2.0
        rsMax = fs / 2.0 * 0.999;
    end
end
