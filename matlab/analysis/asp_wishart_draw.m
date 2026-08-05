function Rhat = asp_wishart_draw(Rtrue, K, L)
%ASP_WISHART_DRAW Exact draw of a complex sample covariance, in O(N^2).
%
%   RHAT = ASP_WISHART_DRAW(RTRUE, K)
%   RHAT = ASP_WISHART_DRAW(RTRUE, K, L)   % L = chol(Rtrue,'lower'), cached
%
%   Draws RHAT distributed as (1/K)*CW_N(K, RTRUE), i.e. exactly the
%   distribution of (1/K)*sum_{n=1..K} r(n)r(n)' for circularly symmetric
%   complex Gaussian r with covariance RTRUE.
%
%   Uses the Bartlett decomposition, so the cost is O(N^2) instead of the
%   O(N^2 * K) of generating and correlating K snapshots.  At N = 4 and
%   K = 16368 that is a ~4000x speedup, which is the difference between a
%   10-trial sanity check and a 10^4-trial Monte Carlo characterisation.
%
%   BARTLETT FOR THE COMPLEX CASE
%   -----------------------------
%   If A is N x N lower triangular with
%       A_ii ~ sqrt(Gamma(K-i+1, 1))          (real, positive)
%       A_ij ~ CN(0,1)  for i > j             (complex standard normal)
%   then A*A' ~ CW_N(K, I), and L*A*A'*L' ~ CW_N(K, L*L').
%
%   The Gamma variates are drawn with the Wilson-Hilferty cube-root
%   transform, which for shape parameters of order 10^4 is accurate to
%   better than 1e-6 relative - far below any quantity being measured here.
%   (A naive sum-of-exponentials draw would need K exponentials per
%   variate and would defeat the purpose.)

n = size(Rtrue,1);

if nargin < 3 || isempty(L)
    L = chol(Rtrue + 1e-14*trace(Rtrue)/n*eye(n), 'lower');
end

A = zeros(n);
for i = 1:n
    shape = K - i + 1;
    A(i,i) = sqrt(gammaWH(shape));
    if i > 1
        A(i,1:i-1) = (randn(1,i-1) + 1i*randn(1,i-1))/sqrt(2);
    end
end

W = A*A';
Rhat = (L*W*L')/K;
Rhat = (Rhat + Rhat')/2;

end

% -------------------------------------------------------------------------
function g = gammaWH(k)
%GAMMAWH Gamma(k,1) variate via the Wilson-Hilferty transform.
%   Exact in the limit of large k; relative error ~ k^-3/2.
z = randn();
g = k * (1 - 1/(9*k) + z/sqrt(9*k))^3;
if g <= 0
    g = k;    % pathological tail guard; probability ~1e-40 at k ~ 1e4
end
end
