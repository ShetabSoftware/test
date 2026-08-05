function [v, info] = asp_beamform(x, f, fx)
%ASP_BEAMFORM Apply beamformer weights, optionally in bit-exact fixed point.
%
%   [V, INFO] = ASP_BEAMFORM(X, F)          floating point
%   [V, INFO] = ASP_BEAMFORM(X, F, FX)      fixed point, FX from asp_fx_plan
%
%   X   nAnt x M complex samples
%   F   nAnt x nBeam complex weights
%   V   nBeam x M beamformer outputs, v = F' * X
%
%   The conjugate is on the WEIGHTS, matching the paper's f^H r convention.
%   In RTL this means the stored coefficient is conj(f_i) and the datapath
%   computes sum_i conj(f_i) * r_i[n]; getting this backwards conjugates the
%   whole spatial response and steers the null to the mirror direction,
%   which is a bug that passes a broadside test and fails everywhere else.
%
%   FIXED-POINT MODEL
%   -----------------
%   Weights are quantised to FX.weight (Q1.16, matching the DSP48E1 18-bit B
%   port).  Products are held exactly.  Accumulation is exact across the
%   nAnt terms.  Only the final output is rounded to FX.beamOut.  This is
%   what the hardware does and it means the only quantisation error that
%   reaches the correlator is one rounding at the output, not nAnt of them.
%
%   Resource note: one complex multiply is 4 real multiplies, or 3 using the
%   Gauss/Karatsuba identity
%       (a+jb)(c+jd): k1 = c(a+b), k2 = a(d-c), k3 = b(c+d)
%       real = k1 - k3,  imag = k1 + k2
%   The 3-multiply form maps directly onto the DSP48 pre-adder, so an
%   N = 4 single-beam beamformer costs 12 DSP48s at the sample rate, or 2
%   DSP48s if time-multiplexed 8:1 against a 160 MHz fabric clock.

if nargin < 3 || isempty(fx)
    v = f' * x;
    info.satFraction = 0;
    info.fixedPoint = false;
    return;
end

fq = fx_quant(f, fx.weight);

acc = fq' * x;                   % exact products and exact accumulation

[v, qi] = fx_quant(acc, fx.beamOut);

info.satFraction = qi.satFraction;
info.fixedPoint  = true;
info.weightsQ    = fq;

% Bit-growth check against the budgeted accumulator.
peak = max(abs([real(acc(:)); imag(acc(:))]));
if peak > 0
    bits = 1 + ceil(log2(peak * 2^(fx.sample.fl + fx.weight.fl) + 1));
else
    bits = 1;
end
info.accBitsUsed = bits;
info.accBitsBudget = fx.beamAccWl;

end
