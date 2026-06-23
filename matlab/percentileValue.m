function v = percentileValue(x, p)
%PERCENTILEVALUE  Linear-interpolated percentile (no Statistics toolbox needed).
%   V = PERCENTILEVALUE(X, P) returns the P-th percentile (0..100) of X using
%   the same linear interpolation convention as NumPy's percentile.

    xs = sort(x(:));
    n = numel(xs);
    if n == 0
        v = 0;
        return;
    end
    if n == 1
        v = xs(1);
        return;
    end
    rank = (p / 100) * (n - 1) + 1;   % 1-based fractional rank
    lo = floor(rank);
    hi = ceil(rank);
    lo = max(min(lo, n), 1);
    hi = max(min(hi, n), 1);
    frac = rank - lo;
    v = xs(lo) * (1 - frac) + xs(hi) * frac;
end
