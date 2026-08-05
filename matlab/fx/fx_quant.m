function [y, info] = fx_quant(x, fmt)
%FX_QUANT Quantise real or complex data to a fixed-point format.
%
%   [Y, INFO] = FX_QUANT(X, FMT)
%
%   Y is returned as a double whose value is exactly representable in FMT,
%   i.e. Y*2^FMT.fl is an integer.  Complex X is quantised component-wise,
%   which is what a hardware I/Q datapath does (each of I and Q is its own
%   two's-complement word).
%
%   INFO.satFraction  fraction of scalar components (I and Q counted
%                     separately) that hit the saturation limit
%   INFO.satCount     number of saturating components
%   INFO.n            number of scalar components examined
%
%   Saturation is reported rather than silently applied because in a nulling
%   receiver a clipped sample is not merely distorted: clipping is a
%   memoryless nonlinearity that redistributes a strong interferer across
%   the whole array manifold.  The interferer stops being rank-one and stops
%   being nullable.  Any epoch containing saturation must be discarded from
%   the covariance estimate, not merely tolerated.

if ~isreal(x)
    [re, iRe] = fx_quant(real(x), fmt);
    [im, iIm] = fx_quant(imag(x), fmt);
    y = re + 1i*im;
    info.satCount = iRe.satCount + iIm.satCount;
    info.n        = iRe.n + iIm.n;
    info.satFraction = info.satCount / max(info.n, 1);
    return;
end

scaled = x * 2^fmt.fl;

switch fmt.round
    case 'convergent'
        q = fx_round_convergent(scaled);
    case 'nearest'
        q = round(scaled);
    case 'floor'
        q = floor(scaled);
    case 'trunc'
        q = fix(scaled);
    otherwise
        error('fx_quant:round', 'Unknown rounding mode "%s".', fmt.round);
end

switch fmt.overflow
    case 'saturate'
        sat = (q > fmt.maxint) | (q < fmt.minint);
        q   = min(max(q, fmt.minint), fmt.maxint);
    case 'wrap'
        sat  = (q > fmt.maxint) | (q < fmt.minint);
        span = 2^fmt.wl;
        q    = mod(q - fmt.minint, span) + fmt.minint;
    otherwise
        error('fx_quant:overflow', 'Unknown overflow mode "%s".', fmt.overflow);
end

y = q * fmt.lsb;

info.satCount    = sum(sat(:));
info.n           = numel(q);
info.satFraction = info.satCount / max(info.n, 1);

end
