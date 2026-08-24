function ok = test_fixedpoint()
%TEST_FIXEDPOINT Numerical hazards of the fixed-point covariance datapath.
%
%   HEADLINE RESULT
%   ---------------
%   Rounding the sample products BEFORE they enter the covariance
%   accumulator injects an error into R that is
%     (a) DETERMINISTIC - a fixed matrix, not noise;
%     (b) low rank, with a fixed phantom eigenvector set purely by the
%         arithmetic conventions, not by anything physical;
%     (c) INDEPENDENT OF DWELL LENGTH, so it cannot be integrated away.
%
%   (c) is what makes it dangerous.  Every other error in this estimator
%   falls as 1/sqrt(K), so a designer who validates at 1 ms reasonably
%   expects improvement at 100 ms.  With product rounding the floor is flat:
%   the system silently stops improving, and nothing in a floating-point
%   simulation predicts it.
%
%   The phantom sits at BORESIGHT, which is the worst possible place: it is
%   where the satellites are.  A truncating implementation therefore
%   synthesises a fake source at zenith and steers the null into the sky.
%   Measured alignment with the all-ones (boresight) vector is 1.000 when
%   the RTL computes both triangles and Hermitises, and 0.906 when it
%   computes the upper triangle and mirrors the conjugate; a random
%   direction would give 0.5.  The reason is immediate once stated: a
%   constant bias c on every product puts the SAME c in every entry of R,
%   and the all-ones matrix is rank one with the boresight steering vector
%   as its eigenvector.
%
%   THE FIX IS FREE.  A DSP48 holds the exact 32-bit product of two 16-bit
%   words and accumulates it in the 48-bit P register with no rounding at
%   all.  Rounding inside the accumulation loop is a design choice, and it
%   is the wrong one.  If a narrower accumulator is unavoidable, convergent
%   rounding converts the flat floor back into a 1/sqrt(K) error.

ok = true;
fprintf('test_fixedpoint:\n');

cfg = asp_config();
rng(11);

n = 4;
Kmax = 80000;
K = 20000;

compositeRms = 10^(-cfg.fx.agcBackoffDB/20);
sigma2 = 2*compositeRms^2;
PsOverSigma2 = cfg.nSpoof * cfg.spoofSnrSample;
Ps = PsOverSigma2 * sigma2;

b = exp(1i*2*pi*rand(n,1));
Rtrue = sigma2*eye(n) + Ps*(b*b');
L = chol(Rtrue,'lower');
xAll = L * (randn(n,Kmax)+1i*randn(n,Kmax))/sqrt(2);

sampleFmt = cfg.fx.sample;
xqAll = fx_quant(xAll, sampleFmt);
xq = xqAll(:,1:K);
Rexact = accumExact(xq, sampleFmt);
[Uc,~] = asp_evd_herm(Rexact, 8);
onesVec = ones(n,1)/sqrt(n);

% ---------------------------------------------------------------- test 1
fprintf('        product-retention sweep, K = %d, spoofer covariance term %.2e\n', K, Ps);
fprintf('        %-8s %-11s %-12s %-13s %-13s %s\n', ...
    'prodFL', 'meanBias', 'bias/Ps dB', 'rho(herm,e)', 'rho(mirr,e)', 'nullFloor dB');

prodFLs = [10 14 18 22 26 30];
nullFloorMirror = zeros(size(prodFLs));
rhoHerm  = zeros(size(prodFLs));
rhoMirr  = zeros(size(prodFLs));

for kk = 1:numel(prodFLs)
    Rherm = accumRounded(xq, sampleFmt, prodFLs(kk), 'floor', 'hermitise');
    Rmirr = accumRounded(xq, sampleFmt, prodFLs(kk), 'floor', 'mirror');

    Eh = Rherm - Rexact;
    Em = Rmirr - Rexact;

    % The dominant mode of the error is the one with the largest ABSOLUTE
    % eigenvalue.  For the Hermitised convention that eigenvalue is
    % NEGATIVE (E = -c*ones(N), eigenvalue -c*N), so sorting by signed value
    % would pick a null-space direction instead.
    rhoHerm(kk) = dominantAlignment(Eh, onesVec);
    rhoMirr(kk) = dominantAlignment(Em, onesVec);

    bias = mean(real(Eh(:)));

    [Ut,~] = asp_evd_herm(Rmirr, 8);
    c = min(abs(Uc(:,1)'*Ut(:,1)), 1);
    nullFloorMirror(kk) = 20*log10(max(sqrt(max(1-c^2,0)), 1e-16));

    fprintf('        %-8d %-11.3e %-12.1f %-13.4f %-13.4f %.1f\n', ...
        prodFLs(kk), bias, 20*log10(max(abs(bias),realmin)/Ps), ...
        rhoHerm(kk), rhoMirr(kk), nullFloorMirror(kk));
end

% For n = 4 a RANDOM unit vector has E|u'e|^2 = 1/n, i.e. rho = 0.5.
% For n = 4 a RANDOM unit vector has E|u'e|^2 = 1/n, i.e. rho = 0.5.  Both
% conventions land far above that, so the phantom is a boresight source in
% both cases: exactly (Hermitised) or predominantly (mirrored, where the
% surviving antisymmetric imaginary bias tilts it slightly).
ok = report(ok, rhoHerm(1) > 0.99, ...
    sprintf('Hermitised convention: bias IS exactly the boresight vector (rho = %.4f; random gives 0.5)', rhoHerm(1)));
ok = report(ok, rhoMirr(1) > 0.85, ...
    sprintf('mirrored convention: phantom is still boresight-dominated (rho = %.4f)', rhoMirr(1)));
ok = report(ok, nullFloorMirror(1) > -30, ...
    sprintf('10-bit product retention caps the null at %.1f dB', nullFloorMirror(1)));
ok = report(ok, nullFloorMirror(end) < -100, ...
    sprintf('full-width products remove the floor entirely (%.1f dB)', nullFloorMirror(end)));

slope = (nullFloorMirror(1) - nullFloorMirror(5)) / (prodFLs(5) - prodFLs(1));
ok = report(ok, abs(slope - 6.02) < 1.5, ...
    sprintf('null floor improves %.2f dB per retained product bit (theory 6.02)', slope));

% ---------------------------------------------------------------- test 2
Ks = [5000 20000 80000];
biasVsK = zeros(size(Ks));
for kk = 1:numel(Ks)
    xk = xqAll(:,1:Ks(kk));
    Ee = accumRounded(xk, sampleFmt, 14, 'floor', 'hermitise') - accumExact(xk, sampleFmt);
    biasVsK(kk) = mean(real(Ee(:)));
end
spread = max(biasVsK)/min(biasVsK);
ok = report(ok, spread < 1.05, ...
    sprintf('bias is dwell-independent: %.3e / %.3e / %.3e over 16x in K', ...
    biasVsK(1), biasVsK(2), biasVsK(3)));

% ---------------------------------------------------------------- test 3
% Convergent rounding: removes the deterministic component, so the residual
% error is zero mean and DOES average down with K.
biasFloorMode = zeros(1,numel(Ks));
biasConvMode  = zeros(1,numel(Ks));
for kk = 1:numel(Ks)
    xk = xqAll(:,1:Ks(kk));
    Rx = accumExact(xk, sampleFmt);
    biasFloorMode(kk) = norm(accumRounded(xk, sampleFmt, 14, 'floor',      'hermitise') - Rx, 'fro');
    biasConvMode(kk)  = norm(accumRounded(xk, sampleFmt, 14, 'convergent', 'hermitise') - Rx, 'fro');
end
fprintf('        K                  : %8d %8d %8d\n', Ks);
fprintf('        ||E||_F, truncation: %.2e %.2e %.2e\n', biasFloorMode);
fprintf('        ||E||_F, convergent: %.2e %.2e %.2e\n', biasConvMode);
ok = report(ok, biasFloorMode(3)/biasFloorMode(1) > 0.8, ...
    sprintf('truncation error does NOT fall with dwell (%.2fx over 16x in K)', ...
    biasFloorMode(3)/biasFloorMode(1)));
ok = report(ok, biasConvMode(3)/biasConvMode(1) < 0.45, ...
    sprintf('convergent-rounding error DOES fall with dwell (%.2fx over 16x in K; 1/sqrt(K) predicts 0.25)', ...
    biasConvMode(3)/biasConvMode(1)));

% ---------------------------------------------------------------- test 4
[acc, info] = fx_cov_accum(xq, sampleFmt, cfg.fx.covAccWl);
err = max(abs(acc.R(:) - Rexact(:)));
ok = report(ok, err < 1e-15, sprintf('fx_cov_accum is bit-exact (max err %.1e)', err));
ok = report(ok, ~info.overflow, ...
    sprintf('accumulator uses %d of %d budgeted bits (%d spare)', ...
    info.bitsUsed, cfg.fx.covAccWl, info.headroom));

% ---------------------------------------------------------------- test 5
trials = 300;
achieved = zeros(1,trials);
for t = 1:trials
    bb = exp(1i*2*pi*rand(n,1));
    h  = ones(n,1)/sqrt(n);
    f  = asp_weights('project', bb, h);
    fq = fx_quant(f, cfg.fx.weight);
    achieved(t) = 10*log10(abs(fq'*bb)^2 / real(fq'*fq));
end
measured = 10*log10(mean(10.^(achieved/10)));
ok = report(ok, abs(measured - cfg.fx.predictedWeightNullFloorDB) < 4, ...
    sprintf('weight-quantisation null floor: measured %.1f dB, predicted %.1f dB', ...
    measured, cfg.fx.predictedWeightNullFloorDB));

% ---------------------------------------------------------------- test 6
fprintf('        converter word length vs usable J/N (14 dB backoff, 1 dB loss):\n');
for wl = [8 10 12 14 16]
    fprintf('          %2d bits -> %5.1f dB\n', wl, asp_adc_jn_limit(wl, 14, 1.0));
end
ok = report(ok, asp_adc_jn_limit(16,14,1.0) - asp_adc_jn_limit(12,14,1.0) > 20, ...
    '16-bit converters buy >20 dB of J/N over the AD9361''s 12-bit');

end

% -------------------------------------------------------------------------
function R = accumExact(xq, fmt)
acc = fx_cov_accum(xq, fmt, 64);
R = acc.R;
end

function rho = dominantAlignment(E, v)
%DOMINANTALIGNMENT Alignment of E's largest-|eigenvalue| eigenvector with v.
if norm(E,'fro') < 1e-18
    rho = NaN;
    return;
end
[U, lam] = asp_evd_herm(E, 8);
[~, idx] = max(abs(lam));
rho = abs(U(:,idx)' * v);
end

function R = accumRounded(xq, fmt, prodFL, rmode, sym)
%ACCUMROUNDED Naive implementation: round each product before accumulating.
%   sym = 'hermitise' : compute all N^2 entries, then (R+R')/2
%   sym = 'mirror'    : compute the upper triangle only, mirror the conjugate
pf = fx_fmt(prodFL + 10, prodFL, 'round', rmode, 'overflow', 'saturate');
[n, K] = size(xq);
R = zeros(n);
chunk = 5000;
for s = 1:chunk:K
    e = min(s+chunk-1, K);
    blk = xq(:, s:e);
    for i = 1:n
        if strcmp(sym,'mirror'), jlist = i:n; else, jlist = 1:n; end
        for j = jlist
            p  = blk(i,:) .* conj(blk(j,:));
            R(i,j) = R(i,j) + sum(fx_quant(p, pf));
        end
    end
end
if strcmp(sym,'mirror')
    for i = 1:n
        for j = i+1:n
            R(j,i) = conj(R(i,j));
        end
    end
else
    R = (R + R')/2;
end
R = R / K;
end

function ok = report(ok, cond, msg)
if cond, fprintf('  PASS  %s\n', msg);
else,    fprintf('  FAIL  %s\n', msg); ok = false; end
end
