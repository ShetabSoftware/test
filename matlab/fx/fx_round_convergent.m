function q = fx_round_convergent(x)
%FX_ROUND_CONVERGENT Round half to even (banker's rounding, IEEE default).
%
%   q = FX_ROUND_CONVERGENT(x) rounds each element of x to the nearest
%   integer, resolving exact .5 cases toward the even integer.
%
%   Unlike round() (half away from zero) and floor() (truncation), this has
%   zero mean error on both symmetric and one-sided data, so it does not
%   inject a rank-one bias into a long covariance accumulation.  In RTL this
%   is the standard "round half to even" logic: add 2^(s-1) - 1 + lsb_out to
%   the value before the right shift by s, where lsb_out is bit s of the
%   input.  It costs one extra OR gate over round-half-up.

f = floor(x);
r = x - f;

q = f;
up   = r > 0.5;
tie  = (r == 0.5);

q(up) = f(up) + 1;
% Ties go to the even neighbour.
qTie = f(tie);
oddTie = mod(qTie, 2) ~= 0;
qTie(oddTie) = qTie(oddTie) + 1;
q(tie) = qTie;

end
