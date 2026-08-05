function [f, dbg] = asp_weights(method, varargin)
%ASP_WEIGHTS Beamformer weight computation.
%
%   f = ASP_WEIGHTS('project', Y, H)
%   f = ASP_WEIGHTS('mvdr',    U, LAM, A, LOADDB)
%   f = ASP_WEIGHTS('lcmv',    U, LAM, C, G, LOADDB)
%   f = ASP_WEIGHTS('mrc',     Q)
%
%   All weights are returned UNNORMALISED except for a power-of-two block
%   floating point scaling.  See FX_BFP_SCALE for why every Euclidean-norm
%   normalisation in the original scripts is removable.
%
%   'project'  Orthogonal projection away from the columns of Y:
%                  f = h - Y*(Y\h)     (least squares, no explicit inverse)
%              For rank 1 this is exactly the paper's equation (14),
%                  f = (I - y y'/(y'y)) h,
%              but computed without forming the N x N projector.  Forming P
%              explicitly costs N^2 storage and N^2 multiplies per beam;
%              the factored form costs 2*N*rank.  At N = 4, rank = 1 that is
%              8 multiplies instead of 16, and the saving grows with N.
%
%   'mvdr'     Minimum variance distortionless response toward steering
%              vector A, using the eigen-decomposition already computed:
%                  Rinv = U*diag(1./(lam + delta))*U'
%                  f    = Rinv*a / (a'*Rinv*a)
%              Requires a calibrated manifold (A must be known), which the
%              projection method does not.  Given a calibrated array, MVDR
%              strictly dominates projection because it also suppresses
%              residual interference and multipath rather than only the one
%              direction that was explicitly nulled.
%
%   'lcmv'     Linearly constrained minimum variance, f = Rinv*C*(C'*Rinv*C)^-1*g.
%              The natural formulation for this product: constrain unit gain
%              toward the satellite AND a hard null toward the spoofer, and
%              let the remaining degrees of freedom minimise output power.
%              With N elements the degrees of freedom ledger is
%                  1 (desired) + rank_spoof (nulls) + (N - 1 - rank_spoof) adaptive
%              which is the sharpest single argument for N = 4 over N = 3:
%              against a spoofer with one ground bounce, N = 3 has ZERO
%              adaptive degrees of freedom left and N = 4 has one.
%
%   'mrc'     Maximal-ratio combining from post-correlation snapshots Q.
%              f = Q (already in the projected subspace if Q came from
%              projected data), which is the exact optimum for white noise.
%
%   DIAGONAL LOADING
%   ----------------
%   LOADDB is relative to trace(R)/N.  Loading is not optional in a shipping
%   product: it bounds the white-noise gain, it desensitises the weights to
%   steering-vector error and to the eigenvector error of a short dwell, and
%   it is the difference between a beamformer that survives a calibration
%   drift and one that does not.  -20 dB is a reasonable default; the
%   classical rule is to load to roughly the level of the steering vector
%   uncertainty (Carlson, IEEE T-AES 1988).

switch lower(method)

    case 'project'
        Y = varargin{1};
        h = varargin{2};
        % Least-squares projection, no explicit inverse, no normalisation.
        if isempty(Y) || size(Y,2) == 0
            f = h;
            dbg.coef = [];
        else
            % Orthonormalise Y first (modified Gram-Schmidt: N*rank^2 ops,
            % numerically far better than forming Y'*Y and inverting it).
            Q = mgs(Y);
            coef = Q' * h;
            f = h - Q*coef;
            dbg.coef = coef;
            dbg.Q = Q;
        end

    case 'mvdr'
        U = varargin{1}; lam = real(varargin{2}(:)); a = varargin{3};
        loadDB = getarg(varargin, 4, -20);
        delta = 10^(loadDB/10) * mean(lam);
        Rinv = U * diag(1./(lam + delta)) * U';
        num = Rinv * a;
        den = real(a' * num);
        f = num / max(den, eps);
        dbg.delta = delta;

    case 'lcmv'
        U = varargin{1}; lam = real(varargin{2}(:));
        C = varargin{3}; g = varargin{4};
        loadDB = getarg(varargin, 5, -20);
        delta = 10^(loadDB/10) * mean(lam);
        Rinv = U * diag(1./(lam + delta)) * U';
        Mmat = C' * Rinv * C;
        % Solve rather than invert; C is N x nConstraints with
        % nConstraints <= N, so this is a 2x2 or 3x3 solve.
        f = Rinv * C * (Mmat \ g(:));
        dbg.delta = delta;
        dbg.constraintCond = cond(Mmat);

    case 'mrc'
        f = varargin{1};
        dbg = struct();

    otherwise
        error('asp_weights:method', 'Unknown method "%s".', method);
end

[f, sh] = fx_bfp_scale(f);
dbg.bfpShift = sh;

end

% -------------------------------------------------------------------------
function Q = mgs(Y)
%MGS Modified Gram-Schmidt orthonormalisation.
[n, p] = size(Y);
Q = zeros(n, p);
kept = 0;
for j = 1:p
    v = Y(:,j);
    for i = 1:kept
        v = v - Q(:,i)*(Q(:,i)'*v);
    end
    nv = norm(v);
    if nv > 1e-12 * max(norm(Y(:,j)), eps)
        kept = kept + 1;
        Q(:,kept) = v / nv;
    end
end
Q = Q(:,1:kept);
end

function v = getarg(args, k, dflt)
if numel(args) >= k && ~isempty(args{k})
    v = args{k};
else
    v = dflt;
end
end
