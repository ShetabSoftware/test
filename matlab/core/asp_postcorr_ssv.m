function [b, dbg] = asp_postcorr_ssv(x, cfg, hyp, nEpoch)
%ASP_POSTCORR_SSV Refine the spoofing signature from despread correlator outputs.
%
%   [B, DBG] = ASP_POSTCORR_SSV(X, CFG, HYP, NEPOCH)
%
%   X       nAnt x M array samples
%   HYP     struct array of despreading hypotheses, one per spoofed PRN:
%             .prn, .tau (s), .fd (Hz), .phi0 (rad, optional)
%           In a product these come from the receiver's tracking loops - the
%           same Doppler feedback the paper's power-maximisation unit
%           already assumes is available.
%   NEPOCH  number of 1 ms epochs to accumulate coherently
%
%   B       nAnt x 1 refined spoofing spatial signature
%   DBG.R   post-correlation spatial covariance
%   DBG.lam its eigenvalues (a rank indicator: rank > 1 means the spoofer
%           has a resolvable multipath component)
%   DBG.snapshotSnrDB  per-snapshot array SNR
%
%   WHY THIS IS THE SINGLE LARGEST AVAILABLE IMPROVEMENT
%   ---------------------------------------------------
%   The pre-correlation estimator saturates.  Measured on the full waveform
%   model, the eigen-based covariance estimator reaches -23.7 dB at 1 ms and
%   only -26.8 dB at 20 ms, where a purely noise-limited estimator would
%   have gained 13 dB over that range.  The shortfall is a bias, not noise:
%   the covariance contains sum_m p_m a_m a_m' from the authentic
%   constellation, which tilts the principal eigenvector by an amount set by
%   p_auth/P_spoof and not by the dwell length.  No amount of integration
%   removes it.
%
%   Despreading raises the wanted term by the full despreading gain,
%   10*log10(K) ~ 36-42 dB, and replaces the authentic-constellation bias
%   with a much smaller one: at the spoofed PRN's code phase, every other
%   signal contributes only a Gold cross-correlation sidelobe.
%
%   IT DOES NOT MAKE THE ESTIMATOR NOISE LIMITED.  It substitutes a
%   different, lower bias floor, set by the C/A code itself.  Measured over
%   a 20 dB sweep of spoofing power (10 scenes each, 20 ms coherent):
%
%       SAPR    0.0 dB  ->  null -26.96 dB
%       SAPR    5.5 dB  ->  null -32.47 dB
%       SAPR   12.0 dB  ->  null -38.96 dB
%       SAPR   20.0 dB  ->  null -46.92 dB
%
%   which is
%
%       null depth  =  -(SAPR + 27.0) dB                                (*)
%
%   to better than 0.1 dB across the whole sweep.  The 27.0 is the C/A
%   cross-correlation bound of 23.9 dB (65/1023) plus about 3.1 dB, because
%   the 65/1023 figure is a worst-case peak and the RMS over random code
%   phases sits below it.  The 1 dB-per-dB slope is simply the authentic
%   contribution falling relative to the spoofer.
%
%   Two things follow from (*).  First, this stage gets BETTER as the attack
%   gets stronger, which is the behaviour you want and the opposite of what
%   a threshold-based detector does.  Second, if you ever need to go deeper
%   than -(SAPR+27) dB, the limiting resource is the CODE, not the array or
%   the dwell - which points at the modernised signals (L1C and L5 pilots,
%   10230 chips, roughly 10 dB better cross-correlation) or at successive
%   interference cancellation.
%
%   ARCHITECTURAL CONSEQUENCE
%   -------------------------
%   The system should be two-stage, and the stages have different jobs and
%   different failure modes:
%
%     Stage 1, pre-correlation, open loop, ~1 ms latency.  Covariance +
%     eigen null.  Reaches -24 dB at 1 ms and saturates near -27 dB, which
%     is enough to drop the spoofer below the authentic signals so the
%     receiver can acquire the REAL peaks.  It needs no receiver state, so
%     it works from a cold start and during reacquisition - exactly when an
%     attack is most effective.  Bias limited by the authentic
%     constellation.
%
%     Stage 2, post-correlation, closed loop, ~20-100 ms latency.  Reaches
%     -(SAPR+27) dB.  Bias limited by the C/A cross-correlation.  It depends
%     on the receiver having tracked something, so it cannot stand alone.
%
%   Be clear about the size of the prize: at a 5.5 dB SAPR stage 2 is worth
%   about 6 dB of extra null depth over stage 1's saturation point, not the
%   20 dB one might hope for.  The deeper null is a bonus.  The real value
%   of stage 2 is what it can do that stage 1 CANNOT at any dwell length:
%
%     * PER-PRN ATTRIBUTION.  Comparing each PRN's despread spatial
%       snapshot against the estimated spoofing signature tells you WHICH
%       PRNs are counterfeit.  Stage 1 sees only an aggregate.  This is what
%       lets the receiver exclude specific measurements rather than
%       discarding the whole solution, and it is the basis of any credible
%       integrity claim.
%     * PER-SATELLITE OPTIMAL COMBINING with no Doppler ambiguity, because
%       the residual Doppler after the tracking loop is near zero by
%       construction.  The paper's equation (20) samples the Doppler phasor
%       once per code period, i.e. at 1 kHz, against Dopplers spanning
%       +/-5 kHz - it aliases.
%     * RANK AND MULTIPATH STRUCTURE of the spoofing source.
%
%   The paper implements neither stage cleanly: its SSV estimator is stage 1
%   with a weaker statistic, and its power-maximisation unit is a partial
%   stage 2 that estimates the AUTHENTIC signatures but never revisits the
%   spoofing signature it nulled with.

nAnt = size(x,1);
K = cfg.K;
Racc = zeros(nAnt);
nSnap = 0;
sigPow = 0;

for hIdx = 1:numel(hyp)
    hh = hyp(hIdx);
    code = asp_ca_code(hh.prn);

    z = zeros(nAnt,1);
    for e = 0:nEpoch-1
        n0 = e*K;
        idx = n0 + (1:K);
        if idx(end) > size(x,2), break; end
        t = (n0 + (0:K-1))/cfg.fs;

        chips = (t - hh.tau) * cfg.chipRate * (1 + hh.fd/cfg.fc);
        ci = mod(floor(chips), 1023) + 1;
        ref = code(ci) .* exp(1i*2*pi*hh.fd*t);

        % Coherent accumulation across epochs: the wanted term adds
        % coherently, the cross-correlation terms rotate at their own
        % Doppler offset and average away.
        z = z + x(:, idx) * conj(ref(:));
    end

    Racc = Racc + z*z';
    sigPow = sigPow + real(z'*z);
    nSnap = nSnap + 1;
end

Racc = Racc / max(nSnap,1);
Racc = (Racc + Racc')/2;

[U, lam] = asp_evd_herm(Racc, cfg.est.jacobiSweeps);
b = U(:,1);

dbg.R = Racc;
dbg.lam = lam;
dbg.U = U;
dbg.nSnap = nSnap;
dbg.rankRatio = lam(1)/max(lam(2), eps);
dbg.snapshotSnrDB = 10*log10(max(lam(1)/max(mean(lam(2:end)),eps) - 1, eps));

end
