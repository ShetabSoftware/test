function asp_demo_single_file()
%ASP_DEMO_SINGLE_FILE Self-contained demo of the anti-spoofing array processor.
%
%   Copy this ONE file anywhere and run it.  No other files, no toolboxes,
%   no path setup.  MATLAB R2018b+ or Octave 7+.
%
%       >> asp_demo_single_file
%
%   It is a condensed version of the full reference model in this repo,
%   included so the algorithm can be read and run in one sitting.  It
%   demonstrates the four changes that matter:
%
%     1. the spatial covariance estimator that replaces the paper's
%        gamma/beta construction  (asp_ssv_from_cov in the full model)
%     2. the eigenvalue DETECTOR the paper does not have (asp_detect)
%     3. the fixed-sweep Jacobi eigen-decomposition that maps onto CORDIC
%        with no dynamic-range growth (asp_evd_herm)
%     4. exact covariance accumulation, and what happens without it
%        (fx_cov_accum / verify/test_fixedpoint)
%
%   For the full model - impairment modelling, bit-exact fixed point,
%   Monte Carlo studies, RTL vector export - use the matlab/ tree and read
%   docs/07-using-the-code.md.

clc;
fprintf('\n=====================================================================\n');
fprintf(' GNSS ANTI-SPOOFING ARRAY PROCESSOR - self-contained demo\n');
fprintf('=====================================================================\n\n');

% ------------------------------------------------------------------ setup
fc     = 1575.42e6;
c      = 299792458;
lambda = c/fc;
fs     = 4*1.023e6;          % 4 samples/chip; 4092 samples per 1 ms epoch
K      = round(fs*1e-3);
nAnt   = 4;
nAuth  = 9;
nSpoof = 9;
cn0    = 45;                 % dB-Hz per authentic PRN
saprDB = 5.5;                % spoof-to-authentic power ratio, per PRN

% CRITICAL: powers are derived from C/N0, never set as arbitrary
% amplitudes.  Per-sample SNR = (C/N0)/fs.  At 45 dB-Hz and 4.092 MHz this
% is -26.2 dB: every GPS signal sits 26 dB BELOW the thermal noise in each
% sample.  Getting this wrong by even 20 dB makes the whole problem look
% trivially easy and every reported null depth meaningless.
authSnr  = 10^(cn0/10)/fs;
spoofSnr = authSnr*10^(saprDB/10);
fprintf('Operating point: per-sample SNR %.1f dB (authentic), %.1f dB (spoof)\n', ...
    10*log10(authSnr), 10*log10(spoofSnr));

% Array: 3 elements on a ring at 120 deg + 1 at the centre, radius 0.45*lambda.
% A triangular baseline lattice aliases only at 1.155*(lambda/2) where a
% square lattice aliases at lambda/2, so this fits 0.78*lambda of baseline
% inside the same radome that limits a 2x2 array to 0.53*lambda.
R0  = 0.45*lambda;
ang = [90 210 330]*pi/180;
pos = [ [R0*cos(ang); R0*sin(ang); zeros(1,3)], [0;0;0] ];

rng(20120926);

% Per-channel gain/phase mismatch (cables, filters, transceiver channels).
% The algorithm never uses the array manifold, so it absorbs this - and
% mutual coupling too - without any calibration.
Cdiag = (10.^(0.5/20*(2*rand(1,nAnt)-1))) .* exp(1i*deg2rad(5*(2*rand(1,nAnt)-1)));
Cdiag = Cdiag(:);

% ------------------------------------------------- scene: 9 auth + 9 spoof
spoofAz = 45; spoofEl = 15;          % terrestrial spoofer => LOW elevation
b = Cdiag .* steer(spoofAz, spoofEl, pos, lambda);

az = rand(1,nAuth)*360;
el = asind(sind(5) + (1-sind(5))*rand(1,nAuth));   % uniform in SOLID ANGLE
A  = zeros(nAnt,nAuth);
for m = 1:nAuth
    A(:,m) = Cdiag .* steer(az(m), el(m), pos, lambda);
end

prns = [2 5 10 12 21 25 29 30 31];
nMs  = 5;
fprintf('Scene: %d authentic PRNs, %d spoofing PRNs from az %d el %d, %d ms\n\n', ...
    nAuth, nSpoof, spoofAz, spoofEl, nMs);

% ------------------------------------------------------- generate signals
x      = genRx(true,  prns, A, b, authSnr, spoofSnr, K, nMs, fs, nAnt);
xClean = genRx(false, prns, A, b, authSnr, spoofSnr, K, nMs, fs, nAnt);

h = ones(nAnt,1)/sqrt(nAnt);

% ============================================================== PART 1
fprintf('--- 1. SPOOFER PRESENT -----------------------------------------\n');
Rhat = (x*x')/size(x,2);
[y, lam] = ssvFromCov(Rhat);
det1 = detect(lam, size(x,2));
f1 = project(y, h);

fprintf('  eigenvalues (whitened) : %s\n', sprintf('%.4f ', lam));
fprintf('  detector statistic     : %.4f  (threshold %.4f) -> %s\n', ...
    det1.stat, det1.thr, tern(det1.detected,'SPOOFING DETECTED','no threat'));
fprintf('  SSV accuracy rho       : %.5f\n', abs(y'*b)/(norm(y)*norm(b)));
fprintf('  null toward spoofer    : %+7.2f dB   (re one antenna element)\n', gainDb(f1,b));
fprintf('  quiescent beam toward it: %+6.2f dB   -> suppression %.1f dB\n', ...
    gainDb(h,b), gainDb(h,b)-gainDb(f1,b));
fprintf('  mean authentic gain    : %+7.2f dB   (quiescent beam: %+.2f dB)\n', ...
    meanGainDb(f1,A), meanGainDb(h,A));
fprintf('  Note h = ones/sqrt(N) is a ZENITH-pointing beam, so its mean gain\n');
fprintf('  over satellites spread across the sky is already below 0 dB.  What\n');
fprintf('  matters is the CHANGE: %+.2f dB.  All of the array gain in this\n', ...
    meanGainDb(f1,A)-meanGainDb(h,A));
fprintf('  method comes from per-satellite power maximisation, not from nulling.\n');

% ============================================================== PART 2
fprintf('\n--- 2. NO SPOOFER (the case the original design gets wrong) ----\n');
RClean = (xClean*xClean')/size(xClean,2);
[yc, lamc] = ssvFromCov(RClean);
det2 = detect(lamc, size(xClean,2));
fprintf('  detector statistic     : %.4f  (threshold %.4f) -> %s\n', ...
    det2.stat, det2.thr, tern(det2.detected,'SPOOFING DETECTED','no threat, BYPASS'));
if det2.detected
    fc2 = project(yc, h);
else
    fc2 = h;                       % pass the quiescent beam unchanged
end
fprintf('  gated weights, mean authentic gain : %+7.2f dB\n', meanGainDb(fc2,A));
fprintf('  if we had nulled anyway (no gate)  : %+7.2f dB   <-- the damage\n', ...
    meanGainDb(project(yc,h),A));
fprintf('  An anti-spoofing box is NOT under attack ~100%% of its operating\n');
fprintf('  hours.  Without the gate it degrades the receiver the entire time.\n');

% ============================================================== PART 3
fprintf('\n--- 3. ESTIMATOR COMPARISON ------------------------------------\n');
fprintf('  %-28s %10s %12s\n','method','rho','null dB');
yp = ssvPaper(x, K, nAnt);
fprintf('  %-28s %10.5f %12.2f\n','paper gamma/beta (2012)', ...
    abs(yp'*b)/(norm(yp)*norm(b)), gainDb(project(yp,h),b));
yr = Rhat(:,1);
fprintf('  %-28s %10.5f %12.2f\n','raw covariance column', ...
    abs(yr'*b)/(norm(yr)*norm(b)), gainDb(project(yr,h),b));
fprintf('  %-28s %10.5f %12.2f\n','whitened EVD (recommended)', ...
    abs(y'*b)/(norm(y)*norm(b)), gainDb(f1,b));
fprintf('  The raw column fails because R(1,1) is dominated by NOISE power\n');
fprintf('  while R(i,1) holds only the spoofer term.  The paper''s split of\n');
fprintf('  magnitude from phase is essential, not stylistic.\n');

% ============================================================== PART 4
fprintf('\n--- 4. THE FIXED-POINT TRAP ------------------------------------\n');
fprintf('  Rounding sample products BEFORE the covariance accumulator adds a\n');
fprintf('  constant to every entry of R.  ones(N) is rank one and its\n');
fprintf('  eigenvector is the BORESIGHT steering vector - so a truncating\n');
fprintf('  implementation invents a source at zenith, where the satellites are.\n');
fprintf('  The error is DWELL-INDEPENDENT, so integrating longer never helps.\n\n');
% Work in ADC units, not in noise-normalised units.  An AGC holds the
% COMPOSITE signal at a fixed backoff below full scale (14 dB gives 5 sigma
% of crest headroom), which scales the signal down while the quantiser LSB
% stays put.  Doing this analysis in noise-normalised units understates the
% problem by about 30 dB - a mistake that is easy to make and expensive.
backoffDB = 14;
g    = 10^(-backoffDB/20) / sqrt(mean(abs(x(:)).^2)/2);
Radc = (g^2) * Rhat;
Ps   = (g^2) * nSpoof * spoofSnr;
fprintf('  %-22s %14s %14s\n','product bits kept','bias/Ps dB','null floor dB');
[Uc,~] = jacobiEig(Radc, 8);
for fl = [10 14 18 22 26]
    bias = 0.5*2^-fl;                       % mean truncation error, 1 product
    [Ue,~] = jacobiEig(Radc - bias*ones(nAnt), 8);
    s = min(abs(Uc(:,1)'*Ue(:,1)),1);
    fprintf('  %-22d %14.1f %14.1f\n', fl, 20*log10(bias/Ps), ...
        20*log10(max(sqrt(max(1-s^2,0)),1e-16)));
end
fprintf('\n  6.02 dB per retained bit, and INDEPENDENT of dwell length: the bias\n');
fprintf('  is a fixed matrix, so averaging longer does not shrink it.  Every\n');
fprintf('  other error here falls as 1/sqrt(K), which is exactly why this one\n');
fprintf('  is missed - the system quietly stops improving with integration.\n');
fprintf('  (verify/test_fixedpoint.m measures this bit-exactly, not analytically.)\n');
fprintf('\n  Fix: a DSP48 holds the exact 32-bit product of two 16-bit words and\n');
fprintf('  accumulates in its 48-bit P register with NO rounding.  Free.\n');

fprintf('\n=====================================================================\n');
fprintf(' Full model: matlab/  |  Usage guide: docs/07-using-the-code.md\n');
fprintf('=====================================================================\n\n');

end

% =========================================================================
%  CORE ALGORITHM
% =========================================================================

function [y, lam] = ssvFromCov(R)
%SSVFROMCOV Spoofing signature from the whitened covariance.
%
%   The paper's estimator is exactly ONE power-iteration step applied to the
%   whitened covariance, started from the reference element:
%       y_paper ~ D * (first column of Rw)
%   Running the iteration to convergence is the whole improvement:
%       y_evd   = D * (principal eigenvector of Rw)
%
%   Whitening by the measured diagonal is what preserves the paper's
%   calibration-free property.  With post-LNA channel mismatch, R = C*R0*C',
%   and the principal eigenvector of R is NOT C*b unless C is a scalar times
%   a unitary.  Whitening turns C into a pure-phase diagonal, which IS
%   unitary.  Plain eig(R) does not have this property.
R = (R+R')/2;
d = real(diag(R)); d(d<=0) = eps;
ds = sqrt(d);
Rw = (R ./ (ds*ds.'));
Rw = (Rw+Rw')/2;
[U, lam] = jacobiEig(Rw, 6);
y = ds .* U(:,1);
end

function [U, lam] = jacobiEig(R, nSweeps)
%JACOBIEIG Hermitian EVD by cyclic Jacobi rotations, FIXED sweep count.
%
%   Every operation is a unitary similarity transform, so ||A||_F is
%   invariant EXACTLY.  Three consequences that matter for hardware:
%     - zero dynamic-range growth, so one word length for the whole block;
%     - unconditional stability regardless of conditioning (and R here IS
%       badly conditioned: the signal is a ~25% perturbation of the identity);
%     - every step is a CORDIC rotation, so no multipliers are needed at all.
%   Fixed sweeps, not a tolerance test: hardware needs deterministic latency.
n = size(R,1);
A = (R+R')/2;
U = eye(n);
for s = 1:nSweeps
    for p = 1:n-1
        for q = p+1:n
            cpq = A(p,q);
            m   = abs(cpq);
            if m > 0, alpha = angle(cpq); else, alpha = 0; end
            th = 0.5*atan2(2*m, real(A(p,p))-real(A(q,q)));
            cs = cos(th); sn = sin(th);
            ej = exp(-1i*alpha);
            Ap = A(:,p); Aq = A(:,q)*ej;
            A(:,p) =  Ap*cs + Aq*sn;  A(:,q) = -Ap*sn + Aq*cs;
            Rp = A(p,:); Rq = A(q,:)*conj(ej);
            A(p,:) =  Rp*cs + Rq*sn;  A(q,:) = -Rp*sn + Rq*cs;
            Up = U(:,p); Uq = U(:,q)*ej;
            U(:,p) =  Up*cs + Uq*sn;  U(:,q) = -Up*sn + Uq*cs;
        end
    end
end
lam = real(diag(A));
[lam, ord] = sort(lam,'descend');
U = U(:,ord);
end

function d = detect(lam, K)
%DETECT Eigenvalue spoofing detector - the block the paper does not have.
%
%   The statistic separates spoofing from authentic signals by spatial RANK,
%   not by power:
%     spoofer:   sum_k p_k*b*b'   -> rank 1, all power in one eigenvalue
%     authentic: sum_m p_m*a*a'  ~= p_tot*I, isotropic, raises the FLOOR
%   which is why it still works at 0 dB SAPR, where the two have equal power.
%
%   The threshold below is a simplified analytic stand-in.  The full model
%   calibrates it by Monte Carlo against a realistic H0 that INCLUDES the
%   authentic constellation, because that contribution does not fall as
%   1/sqrt(K) and a threshold set from the white-noise formula false-alarms
%   on a clean sky at long dwells.  See analysis/asp_detect_threshold.m.
n = numel(lam);
floorLam = mean(lam(2:end));
d.stat = lam(1)/max(floorLam,eps);
d.thr  = (1 + sqrt(n/K))^2 + 6*sqrt(n/K);
d.detected = d.stat > d.thr;
end

function f = project(y, h)
%PROJECT Orthogonal projection away from y, WITHOUT forming the projector.
%
%   f = h - Q*(Q'*h) costs 2*N*rank multiplies; forming P = I - yy'/(y'y)
%   and applying it costs N^2.  No normalisation: the projector is scale
%   invariant, the beamformer output feeds a correlator that is invariant to
%   a constant complex gain, and in hardware a power-of-two shift replaces
%   the inverse square root entirely.
if isempty(y), f = h; return; end
Q = zeros(size(y,1),0);
for j = 1:size(y,2)
    v = y(:,j);
    for i = 1:size(Q,2), v = v - Q(:,i)*(Q(:,i)'*v); end
    if norm(v) > 1e-12, Q(:,end+1) = v/norm(v); end %#ok<AGROW>
end
f = h - Q*(Q'*h);
if norm(f) < eps, f = h; end          % graceful degradation, never error()
end

function y = ssvPaper(x, K, nAnt)
%SSVPAPER The 2012 estimator, equations (7)-(12), as the comparison baseline.
%   Phase from the same-epoch cross-correlation (noise is spatially
%   uncorrelated), magnitude from the one-epoch-lagged self-correlation
%   (noise is temporally uncorrelated).  Correct, but beta is both weaker
%   than the covariance diagonal and needs a full code period of delay
%   memory per channel: ~58 BRAM36 at fs = 16.368 MHz.
g = zeros(nAnt,1); bta = zeros(nAnt,1);
nEp = floor(size(x,2)/K);
for e = 2:nEp
    cur  = x(:,(e-1)*K+(1:K));
    prev = x(:,(e-2)*K+(1:K));
    g   = g   + sum(cur .* conj(repmat(cur(1,:),nAnt,1)), 2);
    bta = bta + sum(cur .* conj(prev), 2);
end
y = sqrt(abs(bta)) .* exp(1i*angle(g));
end

% =========================================================================
%  SIGNAL MODEL (simulation only - none of this goes in the FPGA)
% =========================================================================

function x = genRx(withSpoof, prns, A, b, authSnr, spoofSnr, K, nMs, fs, nAnt)
N = K*nMs;
t = (0:N-1)/fs;
x = zeros(nAnt, N);
rng(7);
for m = 1:numel(prns)
    x = x + emit(prns(m), A(:,m), sqrt(authSnr), t, fs, rand*1e-3, (rand*2-1)*5000, rand*2*pi);
end
if withSpoof
    % A real spoofer synthesises its constellation from ONE clock, so every
    % PRN shares a carrier phase reference.  That is what makes the spoofing
    % energy add coherently in space, and it is the property the whole
    % method depends on.
    phi = rand*2*pi;
    for m = 1:numel(prns)
        x = x + emit(prns(m), b, sqrt(spoofSnr), t, fs, rand*1e-3, (rand*2-1)*5000, phi);
    end
end
x = x + (randn(nAnt,N) + 1i*randn(nAnt,N))/sqrt(2);   % unit-variance noise
end

function s = emit(prn, a, amp, t, fs, tau, fd, phi)
code  = caCode(prn);
chips = (t - tau)*1.023e6*(1 + fd/1575.42e6);          % includes code Doppler
cs    = code(mod(floor(chips),1023) + 1);
s     = a * (amp * cs .* exp(1i*(2*pi*fd*t + phi)));
end

function code = caCode(prn)
%CACODE GPS L1 C/A Gold code, IS-GPS-200 section 3.3.2.3.  Logic 1 -> -1.
%   Verified against the ICD "first 10 chips (octal)" table for all 32 PRNs
%   in verify/test_ca_code.m, together with code balance, the three-valued
%   autocorrelation and the -23.9 dB cross-correlation bound.
tapTab = [2 6;3 7;4 8;5 9;1 9;2 10;1 8;2 9;3 10;2 3;3 4;5 6;6 7;7 8;8 9;9 10; ...
          1 4;2 5;3 6;4 7;5 8;6 9;1 3;4 6;5 7;6 8;7 9;8 10;1 6;2 7;3 8;4 9];
g1 = -ones(1,10); g2 = -ones(1,10); code = zeros(1,1023); tp = tapTab(prn,:);
for k = 1:1023
    code(k) = g1(10)*g2(tp(1))*g2(tp(2));
    g1 = [g1(3)*g1(10), g1(1:9)];
    g2 = [g2(2)*g2(3)*g2(6)*g2(8)*g2(9)*g2(10), g2(1:9)];
end
end

function a = steer(azDeg, elDeg, pos, lambda)
u = [cosd(elDeg)*cosd(azDeg); cosd(elDeg)*sind(azDeg); sind(elDeg)];
a = exp(1i*2*pi*(pos.'*u)/lambda);
end

% =========================================================================
%  METRICS
% =========================================================================

function g = gainDb(f, a)
%GAINDB Gain toward a, in dB relative to ONE antenna element.
%   |f'a|^2/||f||^2 equals 1 for a single element (|a_i| = 1) and has unit
%   mean over the sky, so dB values read directly as "gain over one antenna"
%   and a beamformer can only redistribute, never add, on average.
g = 10*log10(max(abs(f'*a)^2,realmin)/real(f'*f));
end

function g = meanGainDb(f, A)
v = zeros(1,size(A,2));
for k = 1:size(A,2), v(k) = abs(f'*A(:,k))^2/real(f'*f); end
g = 10*log10(mean(v));
end

function s = tern(c,a,b)
if c, s = a; else, s = b; end
end
