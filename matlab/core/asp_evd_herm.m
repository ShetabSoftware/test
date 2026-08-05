function [U, lam, diag_] = asp_evd_herm(R, nSweeps, fmt)
%ASP_EVD_HERM Hermitian eigen-decomposition by cyclic Jacobi rotations.
%
%   [U, LAM] = ASP_EVD_HERM(R, NSWEEPS)
%   [U, LAM, DIAG] = ASP_EVD_HERM(R, NSWEEPS, FMT)
%
%   R        N x N Hermitian matrix
%   NSWEEPS  FIXED number of cyclic sweeps (default 6).  Fixed, not
%            tolerance-based: a hardware block must have deterministic
%            latency, and a convergence test that can run long is a timing
%            closure hazard and a real-time scheduling hazard.
%   FMT      optional fx_fmt applied to the working matrix and the
%            accumulated eigenvector matrix after every rotation, modelling
%            a finite-precision datapath.
%
%   U    N x N unitary, columns are eigenvectors, ordered by DESCENDING
%        eigenvalue
%   LAM  N x 1 real eigenvalues, descending
%   DIAG.offNorm       off-diagonal Frobenius norm after each sweep
%        .frobDrift    relative change in Frobenius norm (should be ~eps;
%                      a growing value means the fixed-point rotation is
%                      losing orthogonality)
%        .rotations    number of 2x2 rotations applied
%
%   WHY JACOBI, AND NOT ANYTHING ELSE, FOR THIS BLOCK
%   -------------------------------------------------
%   Every operation is a unitary similarity transform, so
%       ||A||_F is invariant, exactly, at every step.
%   That single property gives:
%     * zero dynamic-range growth, hence a fixed word length with no
%       intermediate rescaling and no overflow analysis beyond the input;
%     * unconditional numerical stability, independent of the conditioning
%       of R - which matters because R IS ill-conditioned in this
%       application (eigenvalue spread 1.05 with no interferer, 10^6 with a
%       60 dB jammer);
%     * a direct CORDIC mapping: each rotation is one vectoring operation
%       (to find the phase that makes the off-diagonal real) and two
%       rotation operations (to apply it), with no multipliers at all if
%       desired;
%     * natural parallelism: for N = 4 the three disjoint index pairs
%       {(1,2),(3,4)}, {(1,3),(2,4)}, {(1,4),(2,3)} can each be applied
%       simultaneously, so one sweep is three parallel double-rotation
%       steps.
%
%   The alternatives are worse here.  Explicit inversion or Cholesky of R
%   needs a divider and has dynamic range proportional to the condition
%   number.  Power iteration converges at a rate set by the eigenvalue gap,
%   which for a 0 dB SAPR spoofer is about 1.07, so it would need ~200
%   iterations - the very case where the answer matters most.
%
%   ROTATION MATHEMATICS
%   --------------------
%   For the (p,q) pair with a = R(p,p), b = R(q,q) real and c = R(p,q)
%   complex, apply V = D*G where
%       D = I except D(q,q) = exp(-j*angle(c))   makes the entry real, |c|
%       G = I except G(p,p)=G(q,q)=cos(th), G(p,q)=sin(th), G(q,p)=-sin(th)
%       th = 0.5*atan2(2*|c|, a-b)
%   Then (V'*R*V)(p,q) = 0 and the larger eigenvalue lands at index p.

if nargin < 2 || isempty(nSweeps)
    nSweeps = 6;
end
doFx = (nargin >= 3) && ~isempty(fmt);

n = size(R,1);
if size(R,2) ~= n
    error('asp_evd_herm:square', 'R must be square.');
end

A = (R + R')/2;                    % enforce Hermitian
U = eye(n);

frob0 = norm(A,'fro');
diag_.offNorm  = zeros(1,nSweeps);
diag_.frobDrift = zeros(1,nSweeps);
diag_.rotations = 0;

for sweep = 1:nSweeps
    for p = 1:n-1
        for q = p+1:n
            c = A(p,q);
            magc = abs(c);

            % Deterministic: the rotation is applied unconditionally.  When
            % magc is already zero the rotation is the identity, which costs
            % the same number of cycles as any other and keeps the schedule
            % fixed.
            if magc > 0
                alpha = angle(c);
            else
                alpha = 0;
            end

            a = real(A(p,p));
            b = real(A(q,q));
            th = 0.5*atan2(2*magc, a-b);

            cs = cos(th);
            sn = sin(th);
            if doFx
                cs = fx_quant(cs, fmt);
                sn = fx_quant(sn, fmt);
            end

            % V acts on columns p and q.  Build it implicitly.
            eja = exp(-1i*alpha);

            % Column update of A: A <- A*V
            Ap = A(:,p);
            Aq = A(:,q) * eja;
            A(:,p) =  Ap*cs + Aq*sn;
            A(:,q) = -Ap*sn + Aq*cs;

            % Row update of A: A <- V'*A
            Rp = A(p,:);
            Rq = A(q,:) * conj(eja);
            A(p,:) =  Rp*cs + Rq*sn;
            A(q,:) = -Rp*sn + Rq*cs;

            % Accumulate eigenvectors: U <- U*V
            Up = U(:,p);
            Uq = U(:,q) * eja;
            U(:,p) =  Up*cs + Uq*sn;
            U(:,q) = -Up*sn + Uq*cs;

            if doFx
                A = fx_quant(A, fmt);
                U = fx_quant(U, fmt);
            end

            diag_.rotations = diag_.rotations + 1;
        end
    end

    off = A - diag(diag(A));
    diag_.offNorm(sweep)   = norm(off,'fro');
    diag_.frobDrift(sweep) = abs(norm(A,'fro') - frob0)/max(frob0, eps);
end

lam = real(diag(A));
[lam, ord] = sort(lam, 'descend');
U = U(:,ord);

% Fix the phase of each eigenvector so the largest-magnitude component is
% real positive.  Purely cosmetic for the projector (which is phase
% invariant) but it makes golden-vector comparison against RTL possible.
for k = 1:n
    [~, imax] = max(abs(U(:,k)));
    ph = angle(U(imax,k));
    U(:,k) = U(:,k) * exp(-1i*ph);
end

end
