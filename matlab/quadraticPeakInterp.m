function [idx, val] = quadraticPeakInterp(y, k)
%QUADRATICPEAKINTERP  Parabolic refinement of a discrete peak location.
%   [IDX, VAL] = QUADRATICPEAKINTERP(Y, K) fits a parabola through the three
%   samples centred on the (1-based) index K and returns the interpolated
%   fractional index IDX and the interpolated peak value VAL.

    y = y(:).';
    n = numel(y);
    if k <= 1 || k >= n
        idx = k;
        val = y(k);
        return;
    end
    ym1 = y(k - 1);
    y0  = y(k);
    yp1 = y(k + 1);
    denom = ym1 - 2 * y0 + yp1;
    if denom == 0
        idx = k;
        val = y0;
        return;
    end
    delta = 0.5 * (ym1 - yp1) / denom;
    delta = max(min(delta, 1), -1);
    idx = k + delta;
    val = y0 - 0.25 * (ym1 - yp1) * delta;
end
