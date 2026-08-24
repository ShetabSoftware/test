function pOut = asp_projection_outage(N, p, thresholdLinear)
%ASP_PROJECTION_OUTAGE Probability that a projected signal falls below a gain.
%
%   POUT = ASP_PROJECTION_OUTAGE(N, P, THRESHOLDLINEAR)
%
%   For an array response vector a drawn uniformly on the complex sphere
%   with ||a||^2 = N, and an orthogonal projector of rank N-P (i.e. P nulls
%   placed), the retained power obeys
%
%       ||Pa||^2 / N  ~  Beta(N-P, P)                                   (1)
%
%   so the probability that the satellite is left with less gain than
%   THRESHOLDLINEAR (relative to a single antenna element, whose gain is 1)
%   is the regularised incomplete beta function
%
%       POUT = I_{thr/N}(N-P, P)                                        (2)
%
%   Default THRESHOLDLINEAR = 1, i.e. "worse than one antenna".
%
%   For integer parameters (2) has the closed binomial form
%
%       I_x(a,b) = sum_{j=a}^{a+b-1} C(a+b-1, j) x^j (1-x)^(a+b-1-j)     (3)
%
%   which is what is evaluated here - no toolbox needed.
%
%   Worked values at threshold 1:
%       P = 1 null:   N=3 -> 11.1%    N=4 -> 1.6%     N=7 -> 0.02%
%       P = 2 nulls:  N=3 -> 55.6%    N=4 -> 15.6%    N=7 -> 0.03%
%
%   The rank-2 row is the one that decides the element count.  A terrestrial
%   spoofer with a specular ground bounce is a rank-2 source, and against it
%   a 3-element array leaves more than half of the constellation worse off
%   than a single passive antenna.
%
%   These are OPTIMISTIC.  They assume isotropically distributed responses;
%   a physical array of 0.8 lambda aperture has strongly correlated steering
%   vectors across the visible sky, and the measured outage in
%   studies/study_array_size.m is roughly 1.5-2x higher.  The ORDERING and
%   the ratios between element counts survive, which is what the decision
%   rests on.

if nargin < 3 || isempty(thresholdLinear)
    thresholdLinear = 1;
end

if p <= 0
    pOut = 0;
    return;
end
if p >= N
    pOut = 1;
    return;
end

a = N - p;
b = p;
x = min(max(thresholdLinear / N, 0), 1);

m = a + b - 1;
pOut = 0;
for j = a:m
    pOut = pOut + nchoosek(m, j) * x^j * (1-x)^(m-j);
end

end
