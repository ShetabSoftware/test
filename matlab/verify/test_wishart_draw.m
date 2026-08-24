function ok = test_wishart_draw()
%TEST_WISHART_DRAW Check the Bartlett draw against brute-force sampling.
%
%   Compares the mean and the eigenvalue distribution of asp_wishart_draw
%   against explicitly generating and correlating K snapshots.

ok = true;
fprintf('test_wishart_draw:\n');

rng(3);
n = 4; K = 4000; nTrial = 400;

b = exp(1i*2*pi*rand(n,1));
Rtrue = eye(n) + 0.25*(b*b');
L = chol(Rtrue,'lower');

lamB = zeros(nTrial,n);
lamD = zeros(nTrial,n);
accB = zeros(n); accD = zeros(n);

for t = 1:nTrial
    Rb = asp_wishart_draw(Rtrue, K, L);
    lamB(t,:) = sort(real(eig(Rb)),'descend').';
    accB = accB + Rb;

    X = L*(randn(n,K)+1i*randn(n,K))/sqrt(2);
    Rd = (X*X')/K;
    lamD(t,:) = sort(real(eig(Rd)),'descend').';
    accD = accD + Rd;
end

meanErr = norm(accB/nTrial - accD/nTrial,'fro')/norm(Rtrue,'fro');
ok = report(ok, meanErr < 0.01, sprintf('sample means agree to %.4f relative', meanErr));

mB = mean(lamB,1); mD = mean(lamD,1);
sB = std(lamB,0,1); sD = std(lamD,0,1);
fprintf('        E[lambda] Bartlett: %s\n', sprintf('%.4f ', mB));
fprintf('        E[lambda] direct  : %s\n', sprintf('%.4f ', mD));
fprintf('        sd[lambda] Bartlett: %s\n', sprintf('%.4f ', sB));
fprintf('        sd[lambda] direct  : %s\n', sprintf('%.4f ', sD));

ok = report(ok, max(abs(mB-mD)./mD) < 0.01, ...
    sprintf('eigenvalue means agree to %.4f relative', max(abs(mB-mD)./mD)));
ok = report(ok, max(abs(sB-sD)./sD) < 0.15, ...
    sprintf('eigenvalue spreads agree to %.3f relative', max(abs(sB-sD)./sD)));

end

function ok = report(ok, cond, msg)
if cond, fprintf('  PASS  %s\n', msg);
else,    fprintf('  FAIL  %s\n', msg); ok = false; end
end
