function det = asp_detect(lam, K, threshold, maxRank)
%ASP_DETECT Eigenvalue-based detection of a dominant spatial source.
%
%   DET = ASP_DETECT(LAM, K, THRESHOLD, MAXRANK)
%
%   LAM        eigenvalues of the WHITENED sample covariance, descending
%   K          number of samples in the covariance dwell
%   THRESHOLD  decision threshold on the test statistic (see below)
%   MAXRANK    largest number of nulls the beamformer may place
%
%   DET.stat        lambda_1 / mean(lambda_2..lambda_N)
%   DET.detected    stat > threshold
%   DET.rank        estimated number of dominant sources (0..MAXRANK)
%   DET.mdl         MDL criterion value per candidate rank
%   DET.snrEstDB    estimated interference-to-noise ratio of the dominant
%                   source, referred to the per-element noise
%   DET.h0Mean      asymptotic mean of lambda_1 under H0, for context
%
%   THE MISSING BLOCK
%   -----------------
%   Neither the paper nor any of the three provided scripts contains a
%   spoofing detector, yet both scripts contain comments noting that the
%   algorithm nulls the strongest authentic satellite when no spoofer is
%   present.  An anti-spoofing device that degrades the receiver whenever it
%   is NOT under attack has negative expected value: attacks are rare, so
%   the failure mode dominates the operating hours.  The detector must gate
%   the nulling, and it must be part of the same computation, not a bolted-on
%   afterthought.
%
%   TEST STATISTIC
%   --------------
%   Under H0 (spatially white input) the whitened sample covariance is
%   complex Wishart; the largest eigenvalue concentrates near
%
%       E[lambda_1] -> (1 + sqrt(N/K))^2                                (1)
%
%   with Tracy-Widom fluctuations of order K^(-2/3).  For N = 4 and
%   K = 16368, (1) gives 1.031, so the H0 statistic sits just above unity
%   with a very tight spread.  Under H1 a spoofer with total received power
%   P_s per element contributes
%
%       lambda_1 -> 1 + N*P_s/(sigma^2 + N_auth*p_a)                    (2)
%
%   The N in (2) is the array gain against the coherent spoofer, and it is
%   the term that improves with element count: the SAME spoofer is 1.25 dB
%   more detectable on 4 elements than on 3.
%
%   Note what the denominator of (2) says.  The authentic constellation
%   raises the effective noise floor by N_auth*p_a but, being spread across
%   many directions, contributes almost NOTHING to lambda_1 - the sum of
%   many isotropically distributed rank-one terms is approximately a
%   multiple of the identity.  That asymmetry, coherent-in-space versus
%   spread-in-space, is the real physical basis of the whole method, and it
%   is stronger than the "spoofer is louder" argument the paper leans on.
%
%   RANK ESTIMATION
%   ---------------
%   MDL (Wax and Kailath, IEEE T-ASSP 1985) on the whitened eigenvalues.
%   MDL is consistent (AIC is not) and under-estimates rather than
%   over-estimates, which is the safe direction here: an over-estimated rank
%   spends degrees of freedom that the array does not have to spare.

lam = real(lam(:));
lam = sort(lam, 'descend');
n = numel(lam);

if nargin < 4 || isempty(maxRank)
    maxRank = max(n-2, 1);
end

noiseFloor = mean(lam(2:end));
det.stat   = lam(1) / max(noiseFloor, eps);
det.detected = det.stat > threshold;
det.h0Mean = (1 + sqrt(n/max(K,1)))^2;

% Interference-to-noise of the dominant source, referred to one element.
det.snrEstDB = 10*log10(max(lam(1) - noiseFloor, eps) / max(noiseFloor, eps));

% --- MDL rank estimation on the whitened eigenvalues
mdl = zeros(1, n);
for k = 0:n-1
    tail = lam(k+1:end);
    nt = numel(tail);
    gm = exp(mean(log(max(tail, eps))));
    am = mean(tail);
    ll = -K*nt*log(max(gm,eps)/max(am,eps));
    pen = 0.5*k*(2*n - k)*log(max(K,2));
    mdl(k+1) = ll + pen;
end
[~, kBest] = min(mdl);
det.mdl  = mdl;
det.rank = min(kBest - 1, maxRank);

if ~det.detected
    det.rank = 0;
end

% Guard: never spend more than N-2 degrees of freedom on nulls.  With N
% elements, p nulls leave N-p dimensions; the expected post-null array gain
% for a random satellite direction is exactly (N-p), so p = N-1 collapses
% the array to a single element and p = N annihilates everything.
det.rank = min(det.rank, max(n-2, 0));
det.lam  = lam;

end
