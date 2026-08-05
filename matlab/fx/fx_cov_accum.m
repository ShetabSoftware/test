function [acc, info] = fx_cov_accum(x, sampleFmt, accWl, acc)
%FX_COV_ACCUM Bit-exact streaming spatial covariance accumulation.
%
%   [ACC, INFO] = FX_COV_ACCUM(X, SAMPLEFMT, ACCWL)
%   [ACC, INFO] = FX_COV_ACCUM(X, SAMPLEFMT, ACCWL, ACC)  % continue an accumulation
%
%   X          N x M matrix of already-quantised complex samples (values must
%              be exactly representable in SAMPLEFMT).
%   SAMPLEFMT  fx_fmt descriptor of the sample words.
%   ACCWL      accumulator word length in bits (48 maps to one Xilinx DSP48
%              P register, 64 maps to a two-DSP cascade or fabric adder).
%   ACC        optional accumulator state from a previous call.
%
%   ACC.re / ACC.im are N x N integer-valued double matrices holding the raw
%   accumulator contents.  ACC.count is the number of samples accumulated.
%   ACC.R is the scaled Hermitian covariance in real units,
%       R = (acc.re + 1i*acc.im) * 2^(-2*SAMPLEFMT.fl) / count.
%
%   MODEL FIDELITY
%   --------------
%   This is bit-exact, not approximate.  Each product of two WL-bit words is
%   held exactly in 2*WL bits and is added into the accumulator with NO
%   intermediate rounding, exactly as a DSP48 does when its P port is fed
%   back into the adder.  Because doubles carry a 53-bit mantissa and the
%   accumulator never exceeds 2*WL + ceil(log2(M)) bits, the arithmetic here
%   is exact for the word lengths this design uses (16-bit samples,
%   M <= 2^20 => 52 bits).  An assertion enforces that.
%
%   The absence of intermediate rounding is the whole point.  Rounding each
%   product to the accumulator format before adding would inject a per-
%   product bias of up to 0.5 LSB into every entry of R identically, which
%   is a rank-one perturbation along the all-ones vector, i.e. a phantom
%   boresight source.  See FX_FMT for the magnitude of that effect.

[n, m] = size(x);

if nargin < 4 || isempty(acc)
    acc.re    = zeros(n, n);
    acc.im    = zeros(n, n);
    acc.count = 0;
    acc.n     = n;
elseif acc.n ~= n
    error('fx_cov_accum:dim', 'Channel count changed mid-accumulation.');
end

% Convert to raw integer sample words.
xi = x * 2^sampleFmt.fl;
if any(abs(xi(:) - round(xi(:))) > 1e-6)
    error('fx_cov_accum:notQuantised', ...
        'Input is not exactly representable in the declared sample format.');
end
xi = round(real(xi)) + 1i*round(imag(xi));

% Exact outer-product accumulation.  x*x' is the Hermitian rank-M update.
p = xi * xi';

acc.re    = acc.re + real(p);
acc.im    = acc.im + imag(p);
acc.count = acc.count + m;

maxMag  = max(max(abs(acc.re(:))), max(abs(acc.im(:))));
bitsUsed = 1 + ceil(log2(max(maxMag, 1) + 1));

if bitsUsed > 52
    error('fx_cov_accum:exactness', ...
        'Accumulator needs %d bits; the double-precision model is no longer exact.', bitsUsed);
end

info.bitsUsed  = bitsUsed;
info.accWl     = accWl;
info.overflow  = bitsUsed > accWl;
info.headroom  = accWl - bitsUsed;

if info.overflow
    warning('fx_cov_accum:overflow', ...
        'Covariance accumulator needs %d bits but only %d were budgeted.', bitsUsed, accWl);
end

acc.R = (acc.re + 1i*acc.im) * 2^(-2*sampleFmt.fl) / max(acc.count, 1);
acc.R = (acc.R + acc.R')/2;   % enforce exact Hermitian symmetry (free in HW: only the upper triangle is stored)

end
