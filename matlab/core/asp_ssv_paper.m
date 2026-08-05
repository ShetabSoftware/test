function [y, dbg] = asp_ssv_paper(x, K, refIdx)
%ASP_SSV_PAPER Baseline: the spoofing SSV estimator of Daneshmand et al. (2012).
%
%   [Y, DBG] = ASP_SSV_PAPER(X, K, REFIDX)
%
%   X       nAnt x M samples, M must be an integer multiple of K
%   K       samples per PRN code period (one epoch)
%   REFIDX  reference antenna index (default 1)
%
%   Implements equations (7)-(12) of the paper verbatim, for use as the
%   comparison baseline only:
%
%       gamma_i = sum_n r_i(n) * conj(r_1(n))            i = 2..N     (8a)
%       beta_i  = sum_n r_i(n) * conj(r_i(n - K))        i = 1..N     (8b)
%       y_i     = sqrt(|beta_i|) * exp(j*angle(gamma_i))              (7)
%
%   THE ALGEBRA IS CORRECT.  THE STATISTICS ARE NOT THE BEST AVAILABLE.
%   ------------------------------------------------------------------
%   The construction is sound: gamma uses SPATIAL decorrelation of the noise
%   (E[eta_i conj(eta_1)] = 0 for i != 1) to get a noise-free phase, and
%   beta uses TEMPORAL decorrelation (E[eta(n) conj(eta(n-K))] = 0) to get
%   an amplitude, since gamma_1 = sum |r_1|^2 is swamped by noise power and
%   cannot supply the reference amplitude.  Three problems follow.
%
%   (a) beta is the WEAKER of the two amplitude statistics, by a large
%       margin.  From the paper's own equation (11),
%           beta_i -> |C_i|^2 * d,   d = K*sum_k p_k^s e^{j2 pi f_k^s T}
%                                      + K*sum_m p_m^a e^{j2 pi f_m^a T}
%       The Doppler phasors e^{j2 pi f T} with |f| up to 5 kHz and T = 1 ms
%       wrap many times, so d is a RANDOM WALK over the ~18 emitters, of
%       typical magnitude sqrt(18)*p rather than 18*p - and occasionally
%       near zero, at which point the amplitude estimate is pure noise.
%       Meanwhile gamma_i -> K*P_s*C_1^*|C_i| is a COHERENT sum over all
%       spoofing PRNs, larger by roughly N_spoof/sqrt(N_total) ~ 8 dB, and
%       it has no such fading mode.
%
%   (b) beta is DATA-BIT SENSITIVE.  It correlates samples one code period
%       apart, so any 20 ms navigation bit transition inverts the sign of
%       the contribution from every emitter whose bit flipped.  Accumulating
%       beta beyond ~10 ms therefore partially cancels.  gamma correlates
%       samples at the SAME instant, so the data bit appears on both factors
%       and cancels identically - gamma can be accumulated indefinitely.
%
%   (c) beta needs a full code period of delay memory per channel: at
%       fs = 16.368 MHz and 16-bit I/Q that is 4 x 16368 x 32 bits = 2.1 Mbit,
%       about 58 BRAM36 on a 7-series part, purely to compute a statistic
%       that the covariance diagonal already provides for free.
%
%   The replacement is in ASP_SSV_FROM_COV.  This function exists so that
%   the replacement can be measured against the original rather than merely
%   asserted to be better.

if nargin < 3 || isempty(refIdx)
    refIdx = 1;
end

[n, m] = size(x);
nEpoch = floor(m/K);
if nEpoch < 2
    error('asp_ssv_paper:len', ...
        'Need at least 2 code periods; beta is a one-epoch-lagged product.');
end

gammaAcc = zeros(n,1);
betaAcc  = zeros(n,1);

for e = 2:nEpoch
    cur  = x(:, (e-1)*K + (1:K));
    prev = x(:, (e-2)*K + (1:K));

    gammaAcc = gammaAcc + sum(cur .* conj(repmat(cur(refIdx,:), n, 1)), 2);
    betaAcc  = betaAcc  + sum(cur .* conj(prev), 2);
end

y = sqrt(abs(betaAcc)) .* exp(1i*angle(gammaAcc));

dbg.gamma = gammaAcc;
dbg.beta  = betaAcc;
dbg.nEpochUsed = nEpoch - 1;

% gamma(refIdx) is sum|r_ref|^2, a positive real, so its phase is
% identically zero.  The paper's construction is consistent because it uses
% only the phase - but note that the estimator is spending a full-rate
% complex correlator on a quantity known a priori to be real.
end
