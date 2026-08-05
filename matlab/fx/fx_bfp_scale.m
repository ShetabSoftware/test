function [y, shift] = fx_bfp_scale(x, targetIwlHeadroom)
%FX_BFP_SCALE Block-floating-point normalisation by a power of two.
%
%   [Y, SHIFT] = FX_BFP_SCALE(X)
%   [Y, SHIFT] = FX_BFP_SCALE(X, HEADROOM)
%
%   Scales X by 2^SHIFT so that max(abs(Y(:))) lies in [0.5, 1) (or in
%   [0.5, 1)*2^-HEADROOM when HEADROOM > 0).  Returns the shift so the caller
%   can track the exponent.
%
%   WHY THIS REPLACES EVERY norm()-BASED NORMALISATION
%   --------------------------------------------------
%   Three separate places in the original scripts divide by a Euclidean norm:
%   y_est/norm(y_est), f/norm(f) and q_m/norm(q_m).  Each is an inverse
%   square root, which in RTL is a CORDIC or a Newton-Raphson iteration with
%   its own latency, its own pipeline, and its own overflow corner cases.
%
%   None of them is necessary:
%
%     * The orthogonal projector P = I - y*y'/(y'*y) is invariant to the
%       scale of y.  Normalising y is pure waste.
%     * The beamformer output v = f'*r feeds a correlator followed by a
%       C/N0 estimator and a tracking loop, all of which are invariant to a
%       constant complex gain.  Only the RELATIVE weights matter.
%     * q_m/norm(q_m) likewise only sets an overall output scale.
%
%   What DOES matter is that the weights occupy the top of their word so
%   quantisation noise stays small relative to them.  A power-of-two
%   normalisation achieves that with a barrel shifter and a leading-zero
%   count: about 30 LUTs, one cycle, no divider, no convergence question,
%   and it is exact (introduces zero additional error).

if nargin < 2 || isempty(targetIwlHeadroom)
    targetIwlHeadroom = 0;
end

peak = max(abs(x(:)));
if peak == 0
    y = x;
    shift = 0;
    return;
end

% Choose the integer shift that puts the peak just below 2^-headroom.
shift = -ceil(log2(peak)) - targetIwlHeadroom;
y = x * 2^shift;

% Guard against the exact-power-of-two boundary landing at 1.0.
if max(abs(y(:))) >= 2^(-targetIwlHeadroom)
    shift = shift - 1;
    y = x * 2^shift;
end

end
