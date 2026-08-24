function ok = test_evd_herm()
%TEST_EVD_HERM Validate the fixed-sweep Jacobi Hermitian EVD.
%
%   Checks:
%     1. eigenvalues match eig() to high accuracy for random Hermitian
%        matrices at N = 3, 4, 7, 8
%     2. eigenvectors are unitary to machine precision
%     3. the Frobenius norm is invariant across sweeps (the property that
%        guarantees no dynamic-range growth in fixed point)
%     4. convergence with sweep count, for the ill-conditioned case that
%        actually occurs here: a rank-one source only 0.3 dB above the noise
%     5. accuracy is preserved with an 18-bit rotation datapath
%     6. the eigenvector of the dominant source is recovered to better than
%        the estimator's own statistical accuracy, i.e. the EVD is not the
%        bottleneck

ok = true;
fprintf('test_evd_herm:\n');

rng(7);

% --- 1,2,3: general random Hermitian
maxEigErr = 0; maxUnitErr = 0; maxDrift = 0;
for n = [3 4 7 8]
    for trial = 1:20
        A = randn(n) + 1i*randn(n);
        R = A*A'/n;
        [U, lam, dg] = asp_evd_herm(R, 8);
        lamRef = sort(real(eig(R)), 'descend');
        maxEigErr  = max(maxEigErr,  max(abs(lam - lamRef))/max(lamRef));
        maxUnitErr = max(maxUnitErr, norm(U'*U - eye(n), 'fro'));
        maxDrift   = max(maxDrift,   dg.frobDrift(end));
    end
end
ok = report(ok, maxEigErr < 1e-10, sprintf('eigenvalues match eig() (max rel err %.2e)', maxEigErr));
ok = report(ok, maxUnitErr < 1e-12, sprintf('eigenvectors unitary (max ||U''U-I|| %.2e)', maxUnitErr));
ok = report(ok, maxDrift  < 1e-13, sprintf('Frobenius norm invariant (max drift %.2e)', maxDrift));

% --- 4: sweep count on the operating-point matrix
% One rank-one source 0.3 dB above a white floor: eigen gap 1.07, the
% regime where power iteration is useless and Jacobi is not.
n = 4;
b = exp(1i*2*pi*rand(n,1));
R = eye(n) + 0.07/n * (b*b');
lamRef = sort(real(eig(R)),'descend');
errBySweep = zeros(1,6);
for s = 1:6
    [~, lam] = asp_evd_herm(R, s);
    errBySweep(s) = max(abs(lam - lamRef))/lamRef(1);
end
fprintf('        sweep-count convergence (eigen gap %.3f): ', lamRef(1)/lamRef(end));
fprintf('%.1e ', errBySweep); fprintf('\n');
ok = report(ok, errBySweep(4) < 1e-12, ...
    sprintf('4 sweeps suffice at N=4 (err %.2e); 6 chosen for margin', errBySweep(4)));

% --- 5: finite-precision rotation datapath
fmt = fx_fmt(18,16);
worstVec = 0;
for trial = 1:50
    b = exp(1i*2*pi*rand(n,1));
    R = eye(n) + 0.5/n*(b*b');
    [Uq, ~] = asp_evd_herm(R, 6, fmt);
    [Uf, ~] = asp_evd_herm(R, 6);
    c = abs(Uq(:,1)'*Uf(:,1))/(norm(Uq(:,1))*norm(Uf(:,1)));
    worstVec = max(worstVec, sqrt(max(1-c^2,0)));
end
ok = report(ok, worstVec < 1e-3, ...
    sprintf('18-bit rotations track float to sin(theta) < %.1e (= %.0f dB null floor)', ...
    worstVec, 20*log10(max(worstVec,realmin))));

% --- 6: dominant eigenvector recovery vs statistical accuracy
% Statistical accuracy of a K-snapshot estimate of a rank-one subspace is
% sin^2 = (N-1)/K * lam1*lamN/(lam1-lamN)^2.  The EVD must be far better
% than that or it becomes the bottleneck.
K = 16368; Ps = 0.25;
lam1 = 1 + n*Ps/n; lamN = 1;
statSin = sqrt((n-1)/K * lam1*lamN/(lam1-lamN)^2);
ok = report(ok, worstVec < statSin/20, ...
    sprintf('EVD error (%.1e) is >20x below the statistical floor (%.1e)', worstVec, statSin));

end

function ok = report(ok, cond, msg)
if cond, fprintf('  PASS  %s\n', msg);
else,    fprintf('  FAIL  %s\n', msg); ok = false; end
end
