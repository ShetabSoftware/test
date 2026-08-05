function [y, dbg] = asp_ssv_from_cov(R, mode, opt)
%ASP_SSV_FROM_COV Spoofing spatial signature from the sample spatial covariance.
%
%   [Y, DBG] = ASP_SSV_FROM_COV(R, MODE, OPT)
%
%   R     N x N sample spatial covariance, R = (1/K) sum_n r(n) r(n)'
%   MODE  'evd'    principal eigenvector of the whitened covariance
%                  (recommended; supports rank > 1)
%         'gamma'  magnitude from the covariance DIAGONAL, phase from the
%                  reference column: y_i = sqrt(R_ii)*exp(j*angle(R_i1)).
%                  This is the paper's equation (7) with beta replaced by
%                  the diagonal - same structure, strictly better statistic,
%                  no epoch-delay memory.  The cheapest useful variant.
%         'column' raw reference column of R, y_i = R_i1.  Included as a
%                  CAUTIONARY baseline: it does not work, and the reason is
%                  instructive.  R_11 = sum|r_1|^2 is dominated by NOISE
%                  power, while R_i1 for i != 1 contains only the spoofer
%                  cross-term, so the raw column is not proportional to b -
%                  its first entry is larger than the rest by sigma^2/P_s,
%                  about 12 dB at 5.5 dB SAPR.  Measured accuracy is
%                  rho = 0.74 and null depth -6.5 dB.  Splitting magnitude
%                  from phase, as the paper does, is not a stylistic choice;
%                  it is what makes the estimator work at all.
%   OPT   .refIdx        reference element for the 'gamma'/'column' modes
%         .rank          number of spoofing subspace dimensions to return
%         .jacobiSweeps  sweeps for the eigen-decomposition
%         .fmt           optional fx_fmt for a finite-precision EVD
%
%   Y     N x rank spatial signature basis (columns need not be normalised;
%         everything downstream is scale invariant)
%
%   ------------------------------------------------------------------
%   THE UNIFYING OBSERVATION
%   ------------------------------------------------------------------
%   Write D = diag(sqrt(R_11), ..., sqrt(R_NN)) and let Rw = D^-1 R D^-1 be
%   the whitened (unit-diagonal) covariance.  Then the paper's estimator is
%   exactly
%
%       y_paper  ~  D * (first column of Rw)
%
%   because R_i1 carries the phase of C_i b_i (its noise term vanishes for
%   i != 1) and sqrt(R_ii) carries the magnitude |C_i| (the mismatch is
%   post-LNA, so it scales signal and noise together and the common factor
%   sqrt(P_total + sigma0^2) cancels out of the projector).
%
%   Taking the first COLUMN of Rw is precisely one step of power iteration
%   started from the unit vector e_ref.  The estimator this module
%   recommends replaces that single step with the converged answer:
%
%       y_evd    =  D * (principal eigenvector of Rw)
%
%   Everything the paper achieves is retained, and the following is gained
%   for the cost of 6 extra complex multiply-accumulators on an N = 4 array
%   (N(N+1)/2 = 10 correlators instead of 2N-1 = 7):
%
%     * All N eigenvalues, hence an actual spoofing DETECTOR.  The paper has
%       none: its pipeline nulls unconditionally, so in a clean environment
%       it steers a null into the strongest authentic satellite.  Both
%       provided scripts contain comments acknowledging this and neither
%       fixes it.  This is the single most serious defect in the design as
%       it stands, because the failure is silent and it degrades the system
%       exactly when there is no threat.
%     * Rank > 1 nulling, needed for the spoofer's ground bounce and for a
%       simultaneous jammer.
%     * R^-1 for free from the eigen-decomposition, hence MVDR/LCMV weights
%       and diagonal loading without a matrix inversion block.
%     * No epoch-delay memory (about 58 BRAM36 saved at fs = 16.368 MHz).
%     * Immunity to navigation-bit transitions.
%     * Immunity to the random-walk fading of the paper's constant d.
%
%   ------------------------------------------------------------------
%   WHY WHITEN AT ALL
%   ------------------------------------------------------------------
%   With post-LNA channel mismatch, R = C*(R0)*C' where C = diag(C_i) and R0
%   is the covariance seen at the antenna terminals, including the noise
%   term sigma0^2*I.  The principal eigenvector of R is NOT C*b unless C is
%   a scalar multiple of a unitary.  Whitening by the measured diagonal
%   turns C into a pure phase diagonal, which IS unitary, so the principal
%   eigenvector of Rw is the phase-only steering vector and multiplying back
%   by D recovers C*b in the measured domain.  This restores the paper's
%   calibration-free property to the eigen-based estimator, which is not a
%   property it has if you simply run eig(R).

if nargin < 2 || isempty(mode)
    mode = 'evd';
end
if nargin < 3
    opt = struct();
end
if ~isfield(opt,'refIdx') || isempty(opt.refIdx),       opt.refIdx = 1;       end
if ~isfield(opt,'rank') || isempty(opt.rank),           opt.rank = 1;         end
if ~isfield(opt,'jacobiSweeps') || isempty(opt.jacobiSweeps), opt.jacobiSweeps = 6; end
if ~isfield(opt,'fmt'),                                 opt.fmt = [];         end

n = size(R,1);
R = (R + R')/2;

d = real(diag(R));
d(d <= 0) = eps;
dsqrt = sqrt(d);

% Whitening.  In hardware this is one reciprocal-square-root per channel per
% millisecond (4 CORDIC operations at 1 kHz), not a per-sample cost.
Dinv = 1 ./ dsqrt;
Rw = bsxfun(@times, Dinv, bsxfun(@times, R, Dinv.'));
Rw = (Rw + Rw')/2;

switch lower(mode)
    case 'column'
        y = R(:, opt.refIdx);
        dbg.lam = [];
        dbg.U = [];

    case 'gamma'
        % Magnitude from the diagonal, phase from the reference column -
        % the paper's equation (7), with sqrt(|beta_i|) replaced by
        % sqrt(R_ii).  Both estimate |C_i| up to a factor common to all
        % elements, which the projector removes; the diagonal does it with
        % far lower variance and with no code-period delay line.
        col = R(:, opt.refIdx);
        y = dsqrt .* exp(1i*angle(col));
        dbg.lam = [];
        dbg.U = [];

    case 'evd'
        [U, lam] = asp_evd_herm(Rw, opt.jacobiSweeps, opt.fmt);
        p = min(max(opt.rank,1), n);
        y = bsxfun(@times, dsqrt, U(:,1:p));
        dbg.lam = lam;
        dbg.U   = U;

    otherwise
        error('asp_ssv_from_cov:mode', 'Unknown mode "%s".', mode);
end

dbg.Rw    = Rw;
dbg.dsqrt = dsqrt;
dbg.mode  = mode;

end
