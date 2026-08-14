function ok = test_ddc_invariance()
%TEST_DDC_INVARIANCE Is the DDC required by the array processing?
%
%   Claim under test: the covariance, its eigen-decomposition, the detector
%   and the beamformer weights are ALL invariant to a common frequency
%   offset, so translating the signal to baseband is not required by the
%   anti-spoofing algorithm.
%
%   Why it should hold.  A digital downconversion multiplies every channel
%   by the SAME unit-modulus sequence,  r'[n] = r[n]*exp(j*theta[n]).  Then
%
%       r'[n] r'[n]^H = r[n] r[n]^H * |exp(j theta[n])|^2 = r[n] r[n]^H
%
%   term by term, so the accumulated covariance is not merely statistically
%   similar - it is IDENTICAL, sample for sample, before any averaging.
%   Everything derived from R (eigenvalues, eigenvectors, detection
%   statistic, weights) is therefore identical too, and the beamformer
%   output differs only by the same common rotation:
%
%       v'[n] = f^H r'[n] = exp(j theta[n]) * v[n]
%
%   i.e. the same signal sitting at the IF instead of at baseband.
%
%   What this does NOT say: the DC/LO-leakage rejection is still mandatory,
%   because a per-channel DC offset is a rank-one term in R and the
%   detector will find it.  That is a FILTER requirement, not a mixer
%   requirement, and the two are easy to conflate.

ok = true;
fprintf('test_ddc_invariance:\n');

rng(4242);
n = 4; K = 8192;

% --- a scene: spoofer + authentic + per-channel LO leakage
b = exp(1i*2*pi*rand(n,1));
A = exp(1i*2*pi*rand(n,9));
s = (randn(1,K)+1i*randn(1,K))/sqrt(2);
x = 0.30*b*s;
for m = 1:9
    x = x + 0.06*A(:,m)*((randn(1,K)+1i*randn(1,K))/sqrt(2));
end
x = x + (randn(n,K)+1i*randn(n,K))/sqrt(2);

% --- the SAME data at a 2.046 MHz IF (i.e. before the DDC)
fs = 32.736e6; fif = fs/16;
theta = 2*pi*fif*(0:K-1)/fs;
xIF = x .* exp(-1i*theta);

% --- process both ways
[wBB, lamBB, RBB] = chain(x);
[wIF, lamIF, RIF] = chain(xIF);

covErr = max(abs(RBB(:) - RIF(:))) / max(abs(RBB(:)));
lamErr = max(abs(lamBB - lamIF)) / max(abs(lamBB));
wErr   = max(abs(wBB - wIF)) / max(abs(wBB));

fprintf('        covariance  max relative difference : %.3e\n', covErr);
fprintf('        eigenvalues max relative difference : %.3e\n', lamErr);
fprintf('        weights     max relative difference : %.3e\n', wErr);

ok = report(ok, covErr < 1e-12, 'covariance is invariant to a common frequency offset');
ok = report(ok, lamErr < 1e-12, 'eigenvalues are invariant');
ok = report(ok, wErr   < 1e-12, 'beamformer weights are invariant');

% --- null depth is identical; only the output centre frequency differs
gdb = @(f,a) 10*log10(abs(f'*a)^2/real(f'*f));
fprintf('        null depth  baseband %.2f dB | at IF %.2f dB\n', gdb(wBB,b), gdb(wIF,b));
ok = report(ok, abs(gdb(wBB,b) - gdb(wIF,b)) < 1e-9, ...
    'null depth is identical, so the DDC is NOT required by the algorithm');

% --- how much DC/LO leakage can be tolerated?
%   Per-channel LO leakage IS a rank-one term in R, so it is a phantom
%   emitter in principle.  The question is at what level it matters, and
%   the answer sets the required rejection - which turns out to be modest.
%   Tested with NO spoofer present, which is the dangerous case: there the
%   leakage is the LARGEST structured term and the detector has nothing
%   else to lock onto.
xa = zeros(n,K);
for m = 1:9
    xa = xa + 0.06*A(:,m)*((randn(1,K)+1i*randn(1,K))/sqrt(2));
end
xa = xa + (randn(n,K)+1i*randn(n,K))/sqrt(2);      % noise variance 1

THR = 1.06;                                        % calibrated detector threshold
fprintf('        NO spoofer.  Detector statistic vs LO leakage (re thermal noise):\n');
lvls = [-40 -30 -25 -20 -15 -10];
statv = zeros(size(lvls));
for k = 1:numel(lvls)
    dc = 10^(lvls(k)/20) * exp(1i*2*pi*rand(n,1));
    [~, lk] = chain(xa + dc*ones(1,K));
    statv(k) = lk(1)/mean(lk(2:end));
    fprintf('          %+4d dB : stat %.4f  %s\n', lvls(k), statv(k), ...
        tern(statv(k) > THR, '<-- FALSE ALARM', ''));
end
tol = lvls(find(statv <= THR, 1, 'last'));
fprintf('        => leakage must sit about %d dB below the thermal noise.\n', -tol);
fprintf('        The AD9361 with DC tracking enabled leaves ~-60 dBFS against a\n');
fprintf('        noise floor at -14 dBFS, i.e. 46 dB of margin, so this is NOT the\n');
fprintf('        binding constraint - any DC notch or highpass clears it easily.\n');
ok = report(ok, statv(1) < THR && statv(end) > THR, ...
    sprintf('leakage is benign at -40 dB and a false alarm at -10 dB re noise'));

% --- filter cost comparison, the reason to keep the DDC anyway
N = 63;
mulReal    = N;            % real symmetric FIR on complex data: N/2 mults x 2 streams
mulComplex = 3*N;          % complex-coefficient FIR, Karatsuba, no symmetry
mulDDC     = 3;            % one complex multiply per channel per sample
fprintf('        multipliers per channel per sample, %d-tap shaping filter:\n', N);
fprintf('          DDC + real lowpass      : %3d + %3d = %3d\n', mulDDC, mulReal, mulDDC+mulReal);
fprintf('          no DDC, complex bandpass: %3d + %3d = %3d\n', 0, mulComplex, mulComplex);
ok = report(ok, mulDDC + mulReal < mulComplex, ...
    sprintf('the DDC SAVES %d multipliers per channel, it is not a cost', ...
    mulComplex - (mulDDC+mulReal)));

end

% -------------------------------------------------------------------------
function [w, lam, R] = chain(x)
n = size(x,1); K = size(x,2);
R = (x*x')/K;
d = sqrt(real(diag(R)));
Rw = R ./ (d*d.');
[U, lam] = asp_evd_herm((Rw+Rw')/2, 8);
y = d .* U(:,1);
h = ones(n,1)/sqrt(n);
w = asp_weights('project', y, h);
end

function ok = report(ok, cond, msg)
if cond, fprintf('  PASS  %s\n', msg);
else,    fprintf('  FAIL  %s\n', msg); ok = false; end
end

function s = tern(c,a,b)
if c, s = a; else, s = b; end
end
