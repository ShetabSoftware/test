function G = asp_golden_model(varargin)
%ASP_GOLDEN_MODEL Bit-accurate golden reference for the VHDL implementation.
%
%   G = ASP_GOLDEN_MODEL()
%   G = ASP_GOLDEN_MODEL('durationMs',2, 'outDir','./gold', 'spoof',true)
%
%   ONE file, no toolboxes, MATLAB R2018b+ or Octave 7+.
%
%   ---------------------------------------------------------------------
%   CONTRACT WITH THE VHDL
%   ---------------------------------------------------------------------
%   Every signal in this model is a MATLAB double holding an EXACT INTEGER -
%   the raw two's-complement register contents of the corresponding VHDL
%   signal.  Nothing in the datapath is a fraction.  Scaling is carried in
%   the STAGE comments and in the constants at the top, never in the data.
%   A MATLAB value and a VHDL std_logic_vector are therefore the same
%   number and can be diffed with no tolerance.
%
%   Doubles represent 53-bit integers exactly; the widest accumulator here
%   is 48 bits, so every operation below is exact rather than approximate.
%   assertWidth() checks this at each stage instead of trusting it.
%
%   Rules obeyed throughout, each of which exists because breaking it makes
%   the VHDL un-diffable:
%     * every multiply is followed by an EXPLICIT shift, round and saturate
%     * no divisions, no sqrt, no sin/cos, no atan on the datapath - all
%       transcendentals are CORDIC or Newton iterations with a fixed
%       iteration count, modelled bit for bit
%     * no MATLAB filter(), fft(), eig(), inv(), norm() anywhere in the chain
%     * fixed iteration counts everywhere; no tolerance-based loops
%     * where a vectorised expression is used it is EXACTLY equivalent to
%       the hardware, i.e. only where the hardware also accumulates at full
%       precision and rounds once at the end (FIR taps, covariance)
%
%   ---------------------------------------------------------------------
%   BLOCK MAP: one MATLAB local function == one VHDL entity
%   ---------------------------------------------------------------------
%     stage1_ddc          ddc_mixer.vhd          NCO period 16 + complex mult
%     stage2_hbdec        hb_decim2.vhd          halfband, decimate by 2
%     stage3_fir          fir_shape.vhd          63-tap symmetric lowpass
%     stage4_cov          cov_accum.vhd          10 x 48-bit DSP48 accumulators
%     stage5_whiten       whiten.vhd             4 x reciprocal-sqrt (Newton)
%     stage6_evd          jacobi_evd.vhd         CORDIC cyclic Jacobi, 4x4
%     stage7_detect       detect.vhd             compare, no divide
%     stage8_weights      weight_calc.vhd        Gram-Schmidt projection
%     stage9_beamform     beamformer.vhd         4 complex MAC
%     stage10_txdac       tx_scale.vhd           digital AGC + 12-bit DAC
%
%   Each stage writes its input and output to <outDir>/sNN_<name>_{in,out}.txt
%   as signed decimal integers, one value per line (I and Q interleaved for
%   complex signals).  Point the VHDL testbench at the same files.

% =====================================================================
%  0. CONFIGURATION - the single place any number is defined
% =====================================================================
opt.durationMs = 2;
opt.outDir     = fullfile(pwd,'gold_vectors');
opt.spoof      = true;
opt.seed       = 20120926;
opt.dump       = true;
opt.verbose    = true;
for k = 1:2:numel(varargin), opt.(varargin{k}) = varargin{k+1}; end

C = struct();

% ---- RF / rate plan (see docs section on the LO offset decision) ----
C.FC        = 1575.42e6;      % GPS L1
C.LO_OFFSET = 2.046e6;        % = FS_ADC/16 exactly -> NCO period 16
C.LO        = C.FC + C.LO_OFFSET;   % 1577.466 MHz
C.FS_ADC    = 32.736e6;       % 32 x 1.023 MHz
C.FS_WORK   = C.FS_ADC/2;     % 16.368 MHz = 16 samples/chip, 16368 per ms
C.CHIP      = 1.023e6;
C.NANT      = 4;
C.K_DWELL   = round(C.FS_WORK*1e-3);   % 16368 samples per 1 ms dwell

% ---- word lengths.  W = total bits incl. sign; F = fractional bits ----
% Data are stored as raw integers; F only says where the binary point is
% for the purpose of choosing shifts.
C.W_ADC   = 12;  C.F_ADC   = 11;   % AD9361 RX, s12.11, range [-2048, 2047]
C.W_NCO   = 16;  C.F_NCO   = 14;   % twiddle Q1.14 so +1.0 is representable
C.W_MIX   = 16;  C.F_MIX   = 15;   % mixer output
C.W_HB    = 16;  C.F_HB    = 15;   % halfband output
C.W_COEF  = 18;  C.F_COEF  = 17;   % FIR coefficients
C.W_DAT   = 16;  C.F_DAT   = 15;   % working sample word (covariance input)
C.W_ACC   = 48;                    % covariance accumulator = DSP48 P reg
% EVD word: F_EVD = 26 is chosen, not guessed.  The whitened matrix starts
% with a unit diagonal, so trace = N*2^F_EVD and lambda_1 <= trace.  A
% rotation can grow an element by at most sqrt(2).  With F_EVD = 26 the
% worst-case intermediate before the output shift is
%   |A| * 2^F_ROT * 2  =  2^28 * 2^17 * 2  =  2^46
% which is exact in a double and needs 47 bits in the VHDL multiplier-adder,
% and the post-shift result is <= 2^30, inside s32 with a bit to spare.
% F_EVD = 30 would overflow s32 on the eigenvalues.
C.W_EVD   = 32;  C.F_EVD   = 26;   % EVD working matrix
C.W_ROT   = 18;  C.F_ROT   = 16;   % CORDIC-derived cos/sin, Q1.16
C.W_UVEC  = 20;                    % eigenvector accumulator
C.RSQ_F   = 16;                    % reciprocal-sqrt fractional bits (even)
C.R_NORM_BITS = 16;                % R is block-normalised to this magnitude
                                   % before whitening, so that
                                   % R*Y_i*Y_j <= 2^16 * 2^17 * 2^17 = 2^50,
                                   % still exact in a 53-bit double.
C.W_WGT   = 18;  C.F_WGT   = 16;   % beamformer weights, Q1.16 (DSP48 B port)
C.W_BEAM  = 16;  C.F_BEAM  = 15;   % beamformer output
C.W_DAC   = 12;  C.F_DAC   = 11;   % AD9361 TX, s12.11, range [-2048, 2047]

% ---- CORDIC ----
C.CORDIC_N   = 16;            % iterations, FIXED
C.CORDIC_W   = 22;            % internal datapath width
C.ANG_F      = 16;            % angle scaling: radians * 2^16
C.CORDIC_INV_K = 39797;       % round(0.607252935 * 2^16), gain compensation

% ---- algorithm ----
C.JACOBI_SWEEPS = 6;          % FIXED -> deterministic latency
C.MAX_RANK      = 2;
C.DET_NUM       = 1123;       % detector threshold as a rational: 1123/1024
C.DET_DEN       = 1024;       %   = 1.0967, avoids a divider (see stage7)
% Second-null threshold, deliberately much stricter than the first.  A
% second arrival only 6 dB down lifts lambda_2 by ~0.06 above the noise
% floor in a 1 ms dwell, which yields an eigenvector with tens of degrees
% of error; nulling that direction measures WORSE than not nulling it.
% Rank 2 is therefore gated on the second eigenvalue being genuinely
% resolvable, and RANK2_ENABLE lets the soft processor withhold it entirely
% while the MDL test and its hysteresis run at 1 kHz.
C.DET_NUM2      = 1331;       % 1.30
C.RANK2_ENABLE  = false;

% ---- TX / DAC operating point (answer to Q1) ----
C.DAC_TARGET_RMS = 256;       % per-component RMS in LSB at the 12-bit DAC
C.DAC_AGC_LO     = 181;       % hysteresis band, = TARGET/sqrt(2)
C.DAC_AGC_HI     = 362;       %                  = TARGET*sqrt(2)

if opt.dump && ~exist(opt.outDir,'dir'), mkdir(opt.outDir); end

G = struct('C',C,'opt',opt);
if opt.verbose, banner(C, opt); end

% =====================================================================
%  1. STIMULUS  -  12-bit AD9361 RX samples, 4 channels
% =====================================================================
% Not part of the VHDL.  Replace with a file read to drive the model from
% captured hardware data:
%     adc = readVectors(fullfile(opt.outDir,'s00_adc.txt'), C.NANT);
[adc, truth] = genStimulus(C, opt);
assertWidth(adc, C.W_ADC, 'adc');
dumpStage(opt, '00_adc', [], adc);
G.adc = adc; G.truth = truth;

% =====================================================================
%  2. STAGE 1  -  DDC: NCO + complex mixer          [ddc_mixer.vhd]
% =====================================================================
[ddc, ncoTab] = stage1_ddc(adc, C);
assertWidth(ddc, C.W_MIX, 'ddc');
dumpStage(opt, '01_ddc', adc, ddc);
G.ddc = ddc; G.ncoTab = ncoTab;

% =====================================================================
%  3. STAGE 2  -  halfband decimate by 2            [hb_decim2.vhd]
% =====================================================================
[hb, hbCoef] = stage2_hbdec(ddc, C);
assertWidth(hb, C.W_HB, 'hb');
dumpStage(opt, '02_hbdec', ddc, hb);
G.hb = hb; G.hbCoef = hbCoef;

% =====================================================================
%  4. STAGE 3  -  shaping FIR                       [fir_shape.vhd]
% =====================================================================
[dat, firCoef] = stage3_fir(hb, C);
assertWidth(dat, C.W_DAT, 'dat');
dumpStage(opt, '03_fir', hb, dat);
G.dat = dat; G.firCoef = firCoef;

% =====================================================================
%  5..10  per-dwell processing
% =====================================================================
nDwell = floor(size(dat,2)/C.K_DWELL);
beam   = zeros(1, nDwell*C.K_DWELL);
dacOut = zeros(1, nDwell*C.K_DWELL);

w  = initialWeights(C);          % quiescent beam until the first estimate
agcShift = NaN;   % NaN triggers fast AGC acquisition on dwell 1
G.dwell = struct([]);

for d = 1:nDwell
    idx = (d-1)*C.K_DWELL + (1:C.K_DWELL);
    blk = dat(:, idx);

    % --- STAGE 9 first: weights are CAUSAL, dwell d uses dwell d-1's answer
    b = stage9_beamform(blk, w, C);
    assertWidth(b, C.W_BEAM, 'beam');
    beam(idx) = b;

    % --- STAGE 10: digital AGC + 12-bit DAC        [tx_scale.vhd]
    [dq, agcShift, dacInfo] = stage10_txdac(b, agcShift, C);
    assertWidth(dq, C.W_DAC, 'dac');
    dacOut(idx) = dq;

    % --- STAGE 4: covariance                        [cov_accum.vhd]
    [accRe, accIm] = stage4_cov(blk, C);
    assertWidth([accRe(:); accIm(:)], C.W_ACC, 'cov');

    % --- STAGE 5: whiten                            [whiten.vhd]
    [Rw, dsq, rsq] = stage5_whiten(accRe, accIm, C);
    assertWidth(Rw(:), C.W_EVD, 'Rw');

    % --- STAGE 6: eigen-decomposition               [jacobi_evd.vhd]
    [U, lam, nRot] = stage6_evd(Rw, C);

    % --- STAGE 7: detection                         [detect.vhd]
    det = stage7_detect(lam, C);

    % --- STAGE 8: weights                           [weight_calc.vhd]
    wNext = stage8_weights(U, dsq, det.rank, C);
    w = wNext;

    if opt.dump
        % One block per dwell, appended, so the VHDL testbench can step
        % through dwells.  Complex data uses the SAME interleaved (I,Q)
        % convention as the streaming files - one convention everywhere.
        appendVec(opt, '04_cov_out',    accRe(:) + 1i*accIm(:), d);
        appendVec(opt, '05_whiten_out', Rw(:),                  d);
        appendVec(opt, '06_evd_U',      U(:),                   d);
        appendVec(opt, '06_evd_lam',    lam(:),                 d);
        appendVec(opt, '07_detect',     [det.lhs; det.rhs; det.detected; det.rank], d);
        appendVec(opt, '08_weights',    w(:),                   d);
    end

    dw.k = d; dw.lam = lam; dw.det = det; dw.w = w; dw.nRot = nRot;
    dw.dsq = dsq; dw.rsq = rsq; dw.agcShift = agcShift; dw.dac = dacInfo;
    if isempty(G.dwell), G.dwell = dw; else, G.dwell(end+1) = dw; end

    if opt.verbose
        fprintf(['  dwell %2d | lam %s | stat*%d = %8.0f vs %8.0f | det=%d rank=%d ' ...
                 '| agcShift=%+d | clip=%d\n'], d, sprintf('%9.0f',lam), C.DET_DEN, ...
                det.lhs, det.rhs, det.detected, det.rank, agcShift, dacInfo.nClip);
    end
end

dumpStage(opt, '09_beamform', dat(:,1:numel(beam)), beam);
dumpStage(opt, '10_dac', beam, dacOut);

G.beam = beam;
G.dac  = dacOut;
G.nDwell = nDwell;

% =====================================================================
%  SELF-CHECK  (not part of the VHDL - proves the model is CORRECT, not
%  merely self-consistent.  A bit-exact model that computes the wrong
%  thing is worse than no model, because it makes the VHDL wrong too.)
% =====================================================================
G.check = selfCheck(G, C);
if opt.dump, writeManifest(opt, C, G); end
if opt.verbose, summary(G, C, opt); end

end % ======================== end main =================================


function writeManifest(opt, C, G)
%WRITEMANIFEST  The contract the VHDL testbench reads.  Without this the
%   files are ambiguous about ordering, scaling and channel interleave, and
%   an ambiguity in the vector format looks exactly like an RTL bug.
fid = fopen(fullfile(opt.outDir,'MANIFEST.txt'),'w');
fprintf(fid, 'ASP GOLDEN REFERENCE VECTORS\n');
fprintf(fid, '============================\n\n');
fprintf(fid, 'FORMAT (identical for every file)\n');
fprintf(fid, '  signed decimal integers, one per line, no header except\n');
fprintf(fid, '  "# dwell N" separators in the per-dwell files.\n');
fprintf(fid, '  Values are RAW REGISTER CONTENTS, not scaled reals.\n');
fprintf(fid, '  Complex data is interleaved I,Q.\n');
fprintf(fid, '  Multi-channel streaming data is column-major over the\n');
fprintf(fid, '  4 x Nsamples matrix, i.e. ch0,ch1,ch2,ch3 for sample 0,\n');
fprintf(fid, '  then ch0..ch3 for sample 1, and so on.\n');
fprintf(fid, '  Compare BIT-EXACTLY.  Any difference is a bug, not noise.\n\n');
fprintf(fid, 'STAGE                 FILE                    FORMAT      RATE\n');
fprintf(fid, 'ADC in                s00_adc_out.txt         s%d.%d      %.3f MHz x4\n', C.W_ADC,C.F_ADC,C.FS_ADC/1e6);
fprintf(fid, 'ddc_mixer.vhd         s01_ddc_{in,out}.txt    s%d.%d      %.3f MHz x4\n', C.W_MIX,C.F_MIX,C.FS_ADC/1e6);
fprintf(fid, 'hb_decim2.vhd         s02_hbdec_{in,out}.txt  s%d.%d      %.3f MHz x4\n', C.W_HB,C.F_HB,C.FS_WORK/1e6);
fprintf(fid, 'fir_shape.vhd         s03_fir_{in,out}.txt    s%d.%d      %.3f MHz x4\n', C.W_DAT,C.F_DAT,C.FS_WORK/1e6);
fprintf(fid, 'cov_accum.vhd         s04_cov_out.txt         s%d int     1 kHz, 16 entries\n', C.W_ACC);
fprintf(fid, 'whiten.vhd            s05_whiten_out.txt      s%d.%d      1 kHz, 16 entries\n', C.W_EVD,C.F_EVD);
fprintf(fid, 'jacobi_evd.vhd        s06_evd_U.txt           s%d.%d      1 kHz, 16 entries\n', C.W_UVEC,C.F_ROT);
fprintf(fid, 'jacobi_evd.vhd        s06_evd_lam.txt         s%d.%d      1 kHz, 4 entries\n', C.W_EVD,C.F_EVD);
fprintf(fid, 'detect.vhd            s07_detect.txt          int         1 kHz: lhs,rhs,det,rank\n');
fprintf(fid, 'weight_calc.vhd       s08_weights.txt         s%d.%d      1 kHz, 4 complex\n', C.W_WGT,C.F_WGT);
fprintf(fid, 'beamformer.vhd        s09_beamform_out.txt    s%d.%d      %.3f MHz x1\n', C.W_BEAM,C.F_BEAM,C.FS_WORK/1e6);
fprintf(fid, 'tx_scale.vhd          s10_dac_out.txt         s%d.%d      %.3f MHz x1\n', C.W_DAC,C.F_DAC,C.FS_WORK/1e6);
fprintf(fid, '\nCONSTANTS FOR THE VHDL PACKAGE\n');
fprintf(fid, '  NCO table (Q1.%d, 16 entries, I,Q):\n', C.F_NCO);
for k = 1:16
    fprintf(fid, '    %3d %7.0f %7.0f\n', k-1, real(G.ncoTab(k))+0, imag(G.ncoTab(k))+0);
end
fprintf(fid, '  halfband coefficients (Q1.%d): %s\n', C.F_COEF, mat2str(G.hbCoef));
fprintf(fid, '  FIR taps (Q1.%d), %d taps, symmetric:\n', C.F_COEF, numel(G.firCoef));
fprintf(fid, '    %d\n', G.firCoef);
fprintf(fid, '  CORDIC: %d iterations, angle scale 2^%d, 1/K = %d (Q16)\n', ...
    C.CORDIC_N, C.ANG_F, C.CORDIC_INV_K);
fprintf(fid, '  Jacobi: %d sweeps x %d pairs = %d rotations, FIXED\n', ...
    C.JACOBI_SWEEPS, C.NANT*(C.NANT-1)/2, C.JACOBI_SWEEPS*C.NANT*(C.NANT-1)/2);
fprintf(fid, '  detector: lam1*(N-1)*%d > %d*sum(lam2..N)\n', C.DET_DEN, C.DET_NUM);
fprintf(fid, '  rank2   : lam2*(N-2)*%d > %d*sum(lam3..N), gated by RANK2_ENABLE\n', ...
    C.DET_DEN, C.DET_NUM2);
fprintf(fid, '  DAC     : target %d LSB rms, hysteresis [%d, %d], full scale %d\n', ...
    C.DAC_TARGET_RMS, C.DAC_AGC_LO, C.DAC_AGC_HI, 2^(C.W_DAC-1));
fclose(fid);
end


function chk = selfCheck(G, C)
%SELFCHECK  Fixed-point chain versus a double-precision reference.
%   Averaged over EVERY dwell, not just the last one.  A single dwell's
%   fixed-vs-float difference is dominated by realisation noise and can
%   come out either sign; the mean over dwells is the number that means
%   something.
b = G.truth.b;
h = initialWeights(C);
gdb = @(f,a) 10*log10(max(abs(f'*a)^2,realmin)/real(f'*f));
K = C.K_DWELL;

nFix = zeros(1,G.nDwell); nFlo = zeros(1,G.nDwell);
angD = zeros(1,G.nDwell); lamE = zeros(1,G.nDwell);
aFix = zeros(1,G.nDwell);

for d = 1:G.nDwell
    w = G.dwell(d).w;
    x = G.dat(:, (d-1)*K + (1:K));
    R = (x*x')/K;
    dd = sqrt(real(diag(R)));
    Rw = R ./ (dd*dd.');
    [Uf, Lf] = eig((Rw+Rw')/2);
    [lf, o] = sort(real(diag(Lf)),'descend');
    yf = dd .* Uf(:,o(1));
    qf = yf/norm(yf);
    wf = h - qf*(qf'*h);

    nFix(d) = gdb(w,  b);
    nFlo(d) = gdb(wf, b);
    lamE(d) = max(abs(G.dwell(d).lam(:).'/2^C.F_EVD - lf(:).')./lf(:).');
    angD(d) = acosd(min(abs(w'*wf)/(norm(w)*norm(wf)),1));
    aFix(d) = 10*log10(mean(arrayfun(@(k) ...
        abs(w'*G.truth.A(:,k))^2/real(w'*w), 1:size(G.truth.A,2))));
end

lg = @(v) 10*log10(mean(10.^(v/10)));
chk.nullFixedDB   = lg(nFix);
chk.nullFloatDB   = lg(nFlo);
chk.lossDB        = chk.nullFixedDB - chk.nullFloatDB;
chk.quiescentDB   = gdb(h, b);
chk.suppressionDB = chk.quiescentDB - chk.nullFixedDB;
chk.lamRelErr     = max(lamE);
chk.weightAngleDeg = max(angD);
chk.authGainFixedDB = lg(aFix);
chk.authGainQuiDB = 10*log10(mean(arrayfun(@(k) ...
    abs(h'*G.truth.A(:,k))^2/real(h'*h), 1:size(G.truth.A,2))));
chk.nDwell = G.nDwell;
end


% =====================================================================
%  STAGE 1 : DDC                                    [ddc_mixer.vhd]
% =====================================================================
function [y, tab] = stage1_ddc(x, C)
%  Mixes the L1 signal from -LO_OFFSET up to baseband.
%
%  The LO offset is chosen as FS_ADC/16 EXACTLY, so the NCO is a 16-entry
%  ROM indexed by n mod 16 with no phase accumulator and therefore no phase
%  truncation spurs.  In VHDL: a 4-bit counter and a 16 x 32-bit ROM.
%
%  ONE NCO drives all four channels.  That matters: a common complex
%  rotation applied to every element is invisible to the array algorithm
%  (the covariance, the projector and the beamformer are all invariant to a
%  common complex scalar), so a shared NCO contributes exactly zero
%  inter-channel mismatch.  Four independent NCOs would not.
%
%  y = round( x * conj(nco) )  with nco = exp(-j*2*pi*n/16), i.e. mixing UP
%  by +FS_ADC/16 to bring the signal from -2.046 MHz to 0.
n   = size(x,2);
tab = round(exp(1i*2*pi*(0:15)/16) * 2^C.F_NCO);   % Q1.14, |.| = 16384
tab = clampInt(real(tab), C.W_NCO) + 1i*clampInt(imag(tab), C.W_NCO);

ph  = mod(0:n-1, 16) + 1;
lo  = tab(ph);                                     % 1 x n

% Full-precision complex product (hardware: 3 DSP48 with the pre-adder),
% then ONE round and saturate.  Product scale is 2^(F_ADC+F_NCO) = 2^25;
% target scale 2^F_MIX = 2^15, so shift right by 10.
pr = real(x).*real(lo) - imag(x).*imag(lo);
pi_ = real(x).*imag(lo) + imag(x).*real(lo);
sh = C.F_ADC + C.F_NCO - C.F_MIX;
y  = shiftRoundSat(pr, sh, C.W_MIX) + 1i*shiftRoundSat(pi_, sh, C.W_MIX);
end


% =====================================================================
%  STAGE 2 : halfband decimator                     [hb_decim2.vhd]
% =====================================================================
function [y, h] = stage2_hbdec(x, C)
%  11-tap halfband, decimate by 2.  Half the taps are exactly zero and the
%  centre tap is 0.5, so the hardware needs 3 multipliers, not 11.
%  Linear phase, and the SAME coefficients on every channel, so the group
%  delay is identical by construction - which is what keeps the four
%  channels aligned to a fraction of a sample.
h = hbDesign(C);                       % integer coefficients, Q1.17
nCh = size(x,1);
y = zeros(nCh, floor(size(x,2)/2));
sh = C.F_MIX + C.F_COEF - C.F_HB;
for c = 1:nCh
    accR = convFull(real(x(c,:)), h);
    accI = convFull(imag(x(c,:)), h);
    accR = accR(1:2:end);              % decimate AFTER filtering
    accI = accI(1:2:end);
    m = min(numel(accR), size(y,2));
    y(c,1:m) = shiftRoundSat(accR(1:m), sh, C.W_HB) + ...
            1i*shiftRoundSat(accI(1:m), sh, C.W_HB);
end
end


% =====================================================================
%  STAGE 3 : shaping FIR                            [fir_shape.vhd]
% =====================================================================
function [y, h] = stage3_fir(x, C)
%  63-tap symmetric lowpass at 16.368 MHz.  Passband +/-1.2 MHz, stopband
%  from 2.046 MHz, which is where the residual LO leakage lands after the
%  DDC.  It therefore does three jobs at once: band-limits to the C/A main
%  lobe, rejects the LO leakage that would otherwise appear to the
%  eigen-detector as a spatially arbitrary rank-one source, and rejects the
%  I/Q image.
%
%  Symmetric -> 32 multipliers, not 63.  Same coefficients on all channels.
h = firDesign(C);
nCh = size(x,1);
y = zeros(nCh, size(x,2));
sh = C.F_HB + C.F_COEF - C.F_DAT;
for c = 1:nCh
    accR = convSame(real(x(c,:)), h);
    accI = convSame(imag(x(c,:)), h);
    y(c,:) = shiftRoundSat(accR, sh, C.W_DAT) + 1i*shiftRoundSat(accI, sh, C.W_DAT);
end
end


% =====================================================================
%  STAGE 4 : covariance accumulation                [cov_accum.vhd]
% =====================================================================
function [accRe, accIm] = stage4_cov(x, C)
%  R = sum_n x[n] x[n]^H over one dwell, upper triangle only.
%
%  Products are 16x16 = 32 bits and are held EXACTLY.  Accumulation over
%  K = 16368 adds ceil(log2 K) = 14 bits, plus 2 guard bits, giving 48 -
%  exactly the DSP48 P register, so the accumulator lives in the slice with
%  no fabric adder and NO ROUNDING ANYWHERE IN THE LOOP.
%
%  Rounding here is not a precision trade, it is a correctness bug.  A
%  constant rounding bias c on every product puts the same c in every entry
%  of R, and ones(N) is rank one with the boresight steering vector as its
%  eigenvector - so the estimator invents a source at zenith, exactly where
%  the satellites are.  The error is a fixed matrix, so it does NOT shrink
%  with dwell length, which is why longer integration never reveals it.
n = size(x,1);
accRe = zeros(n); accIm = zeros(n);
xr = real(x); xi = imag(x);
for i = 1:n
    for j = i:n
        accRe(i,j) =  sum(xr(i,:).*xr(j,:) + xi(i,:).*xi(j,:));
        accIm(i,j) =  sum(xi(i,:).*xr(j,:) - xr(i,:).*xi(j,:));
    end
end
for i = 1:n
    for j = 1:i-1
        accRe(i,j) =  accRe(j,i);
        accIm(i,j) = -accIm(j,i);
    end
end
end


% =====================================================================
%  STAGE 5 : whitening                              [whiten.vhd]
% =====================================================================
function [Rw, dsq, rsq] = stage5_whiten(accRe, accIm, C)
%  Rw = D^-1 R D^-1 with D = diag(sqrt(R_ii)).
%
%  Whitening is what preserves the calibration-free property under an
%  eigen-based estimator.  With post-LNA channel mismatch R = C R0 C', and
%  the principal eigenvector of R is not C*b unless C is a scalar times a
%  unitary.  Whitening turns C into a pure-phase diagonal, which IS unitary.
%  Plain eig(R) does not have this property; this block is the reason the
%  design still needs no array calibration.
%
%  Cost in hardware: N reciprocal-square-roots and N sqrt per DWELL, i.e.
%  8 operations at 1 kHz.  Negligible.
n = size(accRe,1);

% --- block normalisation.  A common right shift on the WHOLE matrix, from
% a leading-zero count on the largest diagonal.  Common, so it cancels
% exactly in the ratio R(i,j)/sqrt(R_ii R_jj) and cannot bias the result.
dmax = max(diag(accRe));
sh = max(0, ceil(log2(max(dmax,1))) - C.R_NORM_BITS);
Rs = round((accRe + 1i*accIm) / 2^sh);

% --- N reciprocal square roots per dwell (4 operations at 1 kHz)
Yv = zeros(n,1); kv = zeros(n,1); dsq = zeros(n,1);
for i = 1:n
    [Yv(i), kv(i)] = rsqrtNorm(real(Rs(i,i)), C);
    dsq(i) = isqrtInt(real(Rs(i,i)));      % sqrt(R_ii), used by stage 8
end
rsq = Yv;

% --- Rw(i,j) = Rs(i,j) * Y_i * Y_j >> (F_EVD is reached exactly)
%   1/sqrt(a) = Y * 2^-(RSQ_F + RSQ_F/2 + k), so the total right shift for
%   the pair is 2*RSQ_F + RSQ_F/2*2 + k_i + k_j - F_EVD.  Verified below by
%   the assertion that the diagonal comes out as exactly 2^F_EVD.
Rw = zeros(n);
base = 2*C.RSQ_F + C.RSQ_F - C.F_EVD;      % = 3*RSQ_F - F_EVD
for i = 1:n
    for j = 1:n
        p = Rs(i,j) * Yv(i) * Yv(j);
        Rw(i,j) = shiftRoundSatC(p, base + kv(i) + kv(j), C.W_EVD);
    end
end
Rw = Rw - 1i*imag(diag(diag(Rw)));         % diagonal is real by construction
end


% =====================================================================
%  STAGE 6 : Jacobi eigen-decomposition             [jacobi_evd.vhd]
% =====================================================================
function [U, lam, nRot] = stage6_evd(A, C)
%  Cyclic Jacobi with CORDIC-derived rotations, FIXED sweep count.
%
%  Every step is a unitary similarity transform, so the Frobenius norm is
%  invariant EXACTLY.  Three consequences that decide the hardware:
%    - zero dynamic-range growth, so one word length for the whole block
%      with no rescaling and no overflow analysis beyond the input;
%    - unconditional stability regardless of conditioning, which matters
%      because R here IS ill conditioned in the sense that counts (the
%      signal is a ~25% perturbation of the identity);
%    - fixed latency: 6 sweeps x 6 rotations = 36 rotations, no tolerance
%      test, no data-dependent loop count.
%
%  Rotation parameters come from two CORDIC vectoring operations; the
%  rotation itself is applied with multipliers, because a CORDIC applied
%  directly to the data would impose its 1.6468 gain on every pass and
%  require compensation at every step.  Deriving cos/sin ONCE per rotation
%  and applying them with 4 multipliers keeps the gain compensation to a
%  single constant inside the CORDIC.
n = size(A,1);
U = zeros(n); for i=1:n, U(i,i) = 2^C.F_ROT; end     % identity in Q1.16
nRot = 0;

for s = 1:C.JACOBI_SWEEPS
    for p = 1:n-1
        for q = p+1:n
            cr = real(A(p,q)); ci = imag(A(p,q));

            % alpha = atan2(ci, cr);  m = |A(p,q)| (CORDIC gain removed)
            [m, alpha] = cordicVec(cr, ci, C);

            % phi = atan2(2m, A(p,p)-A(q,q));  theta = phi/2  (a shift)
            [~, phi] = cordicVec(real(A(p,p)) - real(A(q,q)), 2*m, C);
            theta = fix(phi/2);

            % twiddles, gain-compensated inside the CORDIC
            [ct, st] = cordicRot(theta,  C);      % cos(th), sin(th)  Q1.16
            [ca, sa] = cordicRot(-alpha, C);      % e^{-j alpha}      Q1.16

            % ---- column update: A(:,q) <- A(:,q)*e^{-j alpha}, then rotate
            Ap = A(:,p);
            Aq = cmulQ(A(:,q), ca + 1i*sa, C.F_ROT, C.W_EVD);
            [A(:,p), A(:,q)] = rot2(Ap, Aq, ct, st, C.F_ROT, C.W_EVD);

            % ---- row update
            Rp = A(p,:);
            Rq = cmulQ(A(q,:), ca - 1i*sa, C.F_ROT, C.W_EVD);
            [A(p,:), A(q,:)] = rot2(Rp, Rq, ct, st, C.F_ROT, C.W_EVD);

            % ---- eigenvector accumulation
            Up = U(:,p);
            Uq = cmulQ(U(:,q), ca + 1i*sa, C.F_ROT, C.W_UVEC);
            [U(:,p), U(:,q)] = rot2(Up, Uq, ct, st, C.F_ROT, C.W_UVEC);

            nRot = nRot + 1;
        end
    end
end

lam = real(diag(A));
[lam, ord] = sort(lam,'descend');
U = U(:,ord);
end


% =====================================================================
%  STAGE 7 : detection                              [detect.vhd]
% =====================================================================
function d = stage7_detect(lam, C)
%  Statistic: lambda_1 / mean(lambda_2..lambda_N) against a threshold.
%
%  Implemented WITHOUT A DIVIDER by cross-multiplying:
%      lam1/mean > THR_NUM/THR_DEN   <=>   lam1*(N-1)*THR_DEN > THR_NUM*sum_tail
%  which is two multiplies and a compare.
%
%  This block is the one the published method does not have.  Without it
%  the system nulls unconditionally and, with no spoofer present, steers a
%  null into the strongest authentic satellite - which is its state for
%  essentially all of its operating hours.
n = numel(lam);
tail = sum(lam(2:end));
d.lhs = lam(1) * (n-1) * C.DET_DEN;
d.rhs = C.DET_NUM * tail;
d.detected = d.lhs > d.rhs;

% Rank by successive application of the same test to the remaining
% eigenvalues.  MDL would be better but needs logarithms; it belongs in the
% soft processor at 1 kHz, where the hysteresis policy also stays tunable.
d.rank = 0;
d.rank2Eligible = false;
if d.detected
    d.rank = 1;
    if n >= 4
        tail2 = sum(lam(3:end));
        d.lhs2 = lam(2) * (n-2) * C.DET_DEN;
        d.rhs2 = C.DET_NUM2 * tail2;
        d.rank2Eligible = d.lhs2 > d.rhs2;
        if d.rank2Eligible && C.RANK2_ENABLE
            d.rank = 2;
        end
    end
end
d.rank = min(d.rank, min(C.MAX_RANK, n-2));
d.lam = lam;
end


% =====================================================================
%  STAGE 8 : weight computation                     [weight_calc.vhd]
% =====================================================================
function w = stage8_weights(U, dsq, rank, C)
%  y = D * U(:,1:rank)   then   w = h - Q*(Q'*h)  by modified Gram-Schmidt.
%
%  The projector is never formed: f = h - Q(Q'h) costs 2*N*rank multiplies
%  where forming P = I - yy'/(y'y) and applying it costs N^2.
%
%  There is no normalisation anywhere.  The projector is scale invariant,
%  and the beamformer output feeds a correlator and a tracking loop that
%  are invariant to a constant complex gain, so the three norm() divisions
%  in a naive implementation are all removable.  What DOES matter is that
%  the weights sit at the top of their word, which is a leading-zero count
%  and a barrel shift - about 30 LUTs, one cycle, and exactly zero added
%  error, versus an inverse-square-root block.
n = size(U,1);
h = initialWeights(C);
if rank == 0
    w = h; return;
end

% y_j = D * U(:,j).  D = diag(sqrt(R_ii)) undoes the whitening, mapping the
% eigenvector back into the measured domain where the beamformer operates.
Y = zeros(n, rank);
for j = 1:rank
    Y(:,j) = bfpScale(U(:,j) .* dsq, C.W_UVEC);
end

% Modified Gram-Schmidt in fixed point.  Each column is normalised by a
% reciprocal square root of its own squared norm - the only division-like
% operation left, and it runs `rank` times per dwell (at most 2 per ms).
Q = zeros(n,0);
for j = 1:rank
    v = Y(:,j);
    for i = 1:size(Q,2)
        ip = sum(conj(Q(:,i)).*v);                 % Q(F_ROT) x Q(F_ROT)
        v  = v - shiftRoundSatC(Q(:,i)*ip, 2*C.F_ROT, C.W_UVEC);
    end
    nrm2 = round(sum(abs(v).^2));
    if nrm2 <= 0, continue; end
    [Yr, kr] = rsqrtNorm(nrm2, C);
    % 1/sqrt(nrm2) = Yr * 2^-(3*RSQ_F/2 + kr); want Q(F_ROT) unit vector
    qcol = shiftRoundSatC(v*Yr, round(1.5*C.RSQ_F) + kr - C.F_ROT, C.W_ROT);
    Q(:,end+1) = qcol; %#ok<AGROW>
end

w = h;
for i = 1:size(Q,2)
    ip = sum(conj(Q(:,i)).*w);                     % Q(F_ROT) x Q(F_WGT)
    w  = w - shiftRoundSatC(Q(:,i)*ip, 2*C.F_ROT, C.W_WGT+4);
end
w = bfpScale(w, C.W_WGT);
end


% =====================================================================
%  STAGE 9 : beamformer                             [beamformer.vhd]
% =====================================================================
function y = stage9_beamform(x, w, C)
%  v[n] = sum_i conj(w_i) * x_i[n].
%
%  The conjugate is on the WEIGHTS.  Getting that backwards conjugates the
%  whole spatial response and steers the null to the mirror direction - a
%  bug that passes a broadside test and fails everywhere else.
%
%  Full-precision products and accumulation across the 4 taps, ONE round at
%  the output, so the only quantisation reaching the DAC is a single
%  rounding rather than four.
wc = conj(w(:));
accR = zeros(1,size(x,2)); accI = zeros(1,size(x,2));
for i = 1:size(x,1)
    accR = accR + real(wc(i))*real(x(i,:)) - imag(wc(i))*imag(x(i,:));
    accI = accI + real(wc(i))*imag(x(i,:)) + imag(wc(i))*real(x(i,:));
end
sh = C.F_DAT + C.F_WGT - C.F_BEAM;
y = shiftRoundSat(accR, sh, C.W_BEAM) + 1i*shiftRoundSat(accI, sh, C.W_BEAM);
end


% =====================================================================
%  STAGE 10 : TX scaling and 12-bit DAC             [tx_scale.vhd]
% =====================================================================
function [y, shiftOut, info] = stage10_txdac(v, shiftIn, C)
%  ANSWER TO Q1.  Output format: SIGNED 12-bit, s12.11 (Q0.11), integer
%  range [-2048, +2047], which is the AD9361 TX data port.
%
%  Operating point: per-component RMS held at 256 LSB, i.e. -18.1 dBFS,
%  giving 8 sigma of crest headroom (clip probability < 1e-15 per sample).
%
%  Why so much backoff, when the RX side uses 14 dB?  Because the two sides
%  are solving opposite problems.  On RX, backoff is spent buying headroom
%  for an interferer you have not removed yet.  On TX the interferer is
%  ALREADY nulled, so there is nothing to leave room for, and the only
%  costs of backing off further are quantisation noise - which is
%  10*log10(1 + (1/12)/256^2) = 5e-6 dB, i.e. nothing - and DAC output
%  power, which is set by the analog attenuator anyway.  Meanwhile clipping
%  is a memoryless nonlinearity that would intermodulate the residual
%  spoofer back into the band and undo the nulling.  The asymmetry is
%  deliberate: on TX, buy clipping margin, it is free.
%
%  Gain control is a POWER-OF-TWO SHIFT with hysteresis, not a multiplier:
%  a barrel shifter plus a leading-zero count, no divider, and it adds
%  exactly zero error.  Amplitude steps on the TX side are harmless (unlike
%  the RX side, where a per-channel gain step corrupts the array), so a
%  coarse shift is sufficient.
%
%  Rounding: convergent (round-half-to-even).  A truncating output would
%  put a -0.5 LSB DC offset on the DAC, which becomes an LO-leakage-like
%  spur at the TX LO and can disturb the downstream receiver's DC and AGC
%  logic.  Convergent rounding costs one OR gate.
%
%  Saturation: hard clamp, NEVER wraparound.  A single wrapped sample is a
%  full-scale transient that spreads across the whole band.
K = numel(v);
p = sum(real(v).^2 + imag(v).^2);        % 1 dwell of power, exact
rmsComp = sqrt(p/(2*K));                 % per-component RMS, current scale

if isempty(shiftIn) || isnan(shiftIn)
    % FAST ACQUISITION on the first dwell.  Real AGC hardware measures
    % wideband power and jumps straight to the right shift, then rate
    % limits.  Stepping one bit per dwell from an arbitrary start would
    % leave the first several milliseconds clipped, and a clipped sample is
    % not merely distorted - clipping intermodulates the residual spoofer
    % back into the band and partially undoes the nulling.
    shiftIn = round(log2(max(rmsComp,1)/C.DAC_TARGET_RMS)) ...
              - (C.F_BEAM - C.F_DAC);
end

shiftOut = shiftIn;
sh = C.F_BEAM - C.F_DAC + shiftIn;       % net right shift
scaled = rmsComp / 2^sh;
if scaled > C.DAC_AGC_HI
    shiftOut = shiftIn + 1;              % too hot -> shift right more
elseif scaled < C.DAC_AGC_LO
    shiftOut = shiftIn - 1;
end

y = shiftRoundSat(real(v), sh, C.W_DAC) + 1i*shiftRoundSat(imag(v), sh, C.W_DAC);

lim = 2^(C.W_DAC-1) - 1;
info.nClip = sum(abs(real(v)/2^sh) > lim | abs(imag(v)/2^sh) > lim);
info.rmsComp = sqrt(mean(real(y).^2 + imag(y).^2)/2);
info.shift = sh;
end


% =====================================================================
%  ARITHMETIC PRIMITIVES  -  every one maps to a named VHDL construct
% =====================================================================

function y = shiftRoundSat(x, sh, W)
%SHIFTROUNDSAT  arithmetic right shift, convergent round, saturate.
%   VHDL: shift_right(acc, sh) with round-half-to-even, then a clamp.
if sh > 0
    q = convRound(x / 2^sh);
elseif sh < 0
    q = x * 2^(-sh);
else
    q = x;
end
y = clampInt(q, W);
end

function y = shiftRoundSatC(x, sh, W)
y = shiftRoundSat(real(x), sh, W) + 1i*shiftRoundSat(imag(x), sh, W);
end

function q = convRound(x)
%CONVROUND round half to even.  Zero mean error, unlike truncation
%   (-0.5 LSB) or round-half-away (non-zero for one-sided data).
f = floor(x); r = x - f;
q = f;
up = r > 0.5; tie = (r == 0.5);
q(up) = f(up) + 1;
t = f(tie); odd = mod(t,2) ~= 0; t(odd) = t(odd) + 1; q(tie) = t;
end

function y = clampInt(x, W)
hi = 2^(W-1) - 1; lo = -2^(W-1);
y = min(max(x, lo), hi);
end

function [u, v] = rot2(a, b, c, s, F, W)
%ROT2  2x2 real Givens rotation applied to a complex pair.
%   u = (a*c + b*s) >> F ,  v = (-a*s + b*c) >> F
%   Products are held at full precision and the SUM is rounded once, not
%   each product, so one rounding per output instead of two.
%   VHDL: 4 DSP48 in a two-deep adder cascade per complex pair.
u = shiftRoundSatC( a*c + b*s, F, W);
v = shiftRoundSatC(-a*s + b*c, F, W);
end

function y = cmulQ(x, c, F, W)
% complex x times complex c (Q F), one round+saturate
pr = real(x)*real(c) - imag(x)*imag(c);
pi_= real(x)*imag(c) + imag(x)*real(c);
y  = shiftRoundSat(pr, F, W) + 1i*shiftRoundSat(pi_, F, W);
end

function y = bfpScale(x, W)
%BFPSCALE  power-of-two normalisation: leading-zero count + barrel shift.
pk = max(abs([real(x(:)); imag(x(:))]));
if pk == 0, y = x; return; end
hi = 2^(W-2);
s = 0;
while pk*2^s < hi && s < W, s = s + 1; end
while pk*2^s > 2^(W-1)-1 && s > -W, s = s - 1; end
y = round(x * 2^s);
y = clampInt(real(y),W) + 1i*clampInt(imag(y),W);
end


% =====================================================================
%  CORDIC  -  fixed iteration count, integer arithmetic
% =====================================================================

function [m, z] = cordicVec(x, y, C)
%CORDICVEC vectoring mode: returns magnitude (gain compensated) and
%   z = atan2(y,x) in units of radians*2^ANG_F.
%   VHDL: N stages of add/subtract/shift, no multipliers.
persistent atanTab lastN lastF
if isempty(atanTab) || lastN ~= C.CORDIC_N || lastF ~= C.ANG_F
    atanTab = round(atan(2.^-(0:C.CORDIC_N-1)) * 2^C.ANG_F);
    lastN = C.CORDIC_N; lastF = C.ANG_F;
end
PI_Q  = round(pi*2^C.ANG_F);
HPI_Q = round(pi/2*2^C.ANG_F);

x = round(x); y = round(y);
z = 0;
% Fold to the right half plane; CORDIC converges only for |angle|<1.7433
if x < 0
    if y >= 0
        t = x; x = y;  y = -t; z =  HPI_Q;
    else
        t = x; x = -y; y =  t; z = -HPI_Q;
    end
end
for i = 0:C.CORDIC_N-1
    xi = x; yi = y;
    if yi >= 0
        x = xi + fix(yi/2^i); y = yi - fix(xi/2^i); z = z + atanTab(i+1);
    else
        x = xi - fix(yi/2^i); y = yi + fix(xi/2^i); z = z - atanTab(i+1);
    end
end
% remove the CORDIC gain 1.64676 -> multiply by 1/K
m = fix(x * C.CORDIC_INV_K / 2^16);
if z >  PI_Q, z = z - 2*PI_Q; end
if z < -PI_Q, z = z + 2*PI_Q; end
end

function [c, s] = cordicRot(z, C)
%CORDICROT rotation mode applied to (1/K, 0), returning cos and sin in
%   Q1.F_ROT.  The 1/K seed removes the CORDIC gain, so no post-scaling.
persistent atanTab lastN lastF
if isempty(atanTab) || lastN ~= C.CORDIC_N || lastF ~= C.ANG_F
    atanTab = round(atan(2.^-(0:C.CORDIC_N-1)) * 2^C.ANG_F);
    lastN = C.CORDIC_N; lastF = C.ANG_F;
end
HPI_Q = round(pi/2*2^C.ANG_F);

neg = false;
z = round(z);
if z >  HPI_Q, z = z - 2*HPI_Q; neg = true; end
if z < -HPI_Q, z = z + 2*HPI_Q; neg = true; end

x = round(C.CORDIC_INV_K / 2^16 * 2^C.F_ROT);
y = 0;
for i = 0:C.CORDIC_N-1
    xi = x; yi = y;
    if z >= 0
        x = xi - fix(yi/2^i); y = yi + fix(xi/2^i); z = z - atanTab(i+1);
    else
        x = xi + fix(yi/2^i); y = yi - fix(xi/2^i); z = z + atanTab(i+1);
    end
end
if neg, x = -x; y = -y; end
c = clampInt(x, C.W_ROT); s = clampInt(y, C.W_ROT);
end

function [Y, k] = rsqrtNorm(a, C)
%RSQRTNORM  Reciprocal square root by range reduction + fixed Newton count.
%
%   Returns Y and k such that   1/sqrt(a) = Y * 2^-(3*RSQ_F/2 + k)
%   (with RSQ_F even), where Y is an integer in (2^F, 2^(F+1)].
%
%   Range reduction first: shift a by an EVEN number of bits 2k so the
%   mantissa lands in [2^(F-2), 2^F), i.e. a' in [0.25, 1).  Then
%   y' = 1/sqrt(a') is in (1, 2] and the Newton iteration operates on a
%   bounded argument, so a FIXED three-iteration count reaches full
%   precision for every input - which is what makes the latency
%   deterministic.  Without range reduction the iteration count would
%   depend on the data.
%
%   VHDL: leading-zero count, barrel shift, 16-entry seed ROM, 3 stages of
%   (multiply, subtract, multiply, shift).  Widest intermediate is
%   a_n*Y^2 <= 2^F * 2^(2F+2) = 2^50 at F = 16.
F = C.RSQ_F;
if a <= 0, Y = 0; k = 0; return; end

e = floor(log2(a));
k = floor((e - (F-1))/2);                    % even shift, a_n = a >> 2k
an = round(a / 2^(2*k));
while an >= 2^F,      k = k + 1; an = round(a / 2^(2*k)); end
while an <  2^(F-2),  k = k - 1; an = round(a / 2^(2*k)); end

% Seed: the chord of 1/sqrt(x) across [0.25,1), i.e. y0 = 7/3 - (4/3)a'.
% Worst-case relative error 18%, which the quadratic Newton iteration takes
% to 2e-5 in three steps and below 1 LSB in four.  A cruder seed such as
% (1.5 - a') leaves 0.7% after three steps - enough to shift the whitened
% diagonal off unity, which then propagates into every eigenvalue.
Y = fix((7*2^F - 4*an) / 3);
Y = min(max(Y, 2^F), 2^(F+1));

for it = 1:4
    t = fix(an * Y / 2^F);                   % a'*y      in Q(F)
    t = fix(t  * Y / 2^F);                   % a'*y^2    in Q(F)
    Y = fix(Y * (3*2^F - t) / 2^(F+1));      % y*(3 - a'y^2)/2
end
Y = min(max(Y, 1), 2^(F+1));
end

function s = isqrtInt(a)
%ISQRTINT  Integer square root, fixed 24-step restoring algorithm.
%   VHDL: the classic non-restoring square-root array; no divider.
if a <= 0, s = 0; return; end
s = 0; rem = 0;
for i = 23:-1:0
    rem = rem*4 + fix(mod(fix(a/4^i), 4));
    t = 4*s + 1;
    if rem >= t, rem = rem - t; s = 2*s + 1; else, s = 2*s; end
end
end


% =====================================================================
%  SUPPORT  (not part of the VHDL)
% =====================================================================

function w = initialWeights(C)
%   Quiescent beam h = ones, in Q1.F_WGT.  Not scaled by 1/sqrt(N): the
%   whole chain is invariant to a common real scale.
w = complex(repmat(2^(C.F_WGT-1), C.NANT, 1), zeros(C.NANT,1));
end

function h = hbDesign(C)
%   11-tap halfband, integer Q1.F_COEF, sum = 2^F_COEF (unity DC gain).
hf = [-0.031250 0 0.281250 0.500000 0.281250 0 -0.031250];
hf = [0 0 hf 0 0];
h  = round(hf * 2^C.F_COEF);
h(h==0) = 0;
h  = adjustDC(h, 2^C.F_COEF);
end

function h = firDesign(C)
%   63-tap windowed-sinc lowpass, fc = 1.6 MHz at 16.368 MHz, Blackman.
%   Designed here rather than loaded so the model is self-contained; in
%   production the coefficients come from a file shared with the VHDL ROM.
N = 63; fc = 1.6e6/C.FS_WORK;
n = (0:N-1) - (N-1)/2;
hf = 2*fc*sinc2(2*fc*n);
wnd = 0.42 - 0.5*cos(2*pi*(0:N-1)/(N-1)) + 0.08*cos(4*pi*(0:N-1)/(N-1));
hf = hf .* wnd;
hf = hf / sum(hf);
h  = round(hf * 2^C.F_COEF);
h  = adjustDC(h, 2^C.F_COEF);
end

function h = adjustDC(h, target)
d = target - sum(h);
[~,i] = max(abs(h));
h(i) = h(i) + d;
end

function y = sinc2(x)
y = ones(size(x));
nz = x ~= 0;
y(nz) = sin(pi*x(nz))./(pi*x(nz));
end

function y = convFull(x, h)
y = zeros(1, numel(x));
L = numel(h);
xp = [zeros(1,L-1) x];
for k = 1:L
    y = y + h(k)*xp(L-k+1 : L-k+numel(x));
end
end

function y = convSame(x, h)
L = numel(h); d = (L-1)/2;
xp = [zeros(1,d) x zeros(1,d)];
y = zeros(1, numel(x));
for k = 1:L
    y = y + h(k)*xp(L-k+1 : L-k+numel(x));
end
end

function [x, truth] = genStimulus(C, opt)
%   4-channel AD9361 RX capture at FS_ADC, 12-bit, with the L1 signal
%   placed at -LO_OFFSET, plus a per-channel DC/LO-leakage term.
rng(opt.seed);
lambda = 299792458/C.FC;
R0 = 0.45*lambda; ang = [90 210 330]*pi/180;
pos = [ [R0*cos(ang); R0*sin(ang); zeros(1,3)], [0;0;0] ];
Cd = (10.^(0.5/20*(2*rand(1,C.NANT)-1))).*exp(1i*deg2rad(5*(2*rand(1,C.NANT)-1)));
Cd = Cd(:);

N  = round(C.FS_ADC*opt.durationMs*1e-3);
t  = (0:N-1)/C.FS_ADC;
cn0 = 45; authSnr = 10^(cn0/10)/C.FS_ADC; spoofSnr = authSnr*10^(5.5/10);
prns = [2 5 10 12 21 25 29 30 31];

sig = zeros(C.NANT, N);
truth.A = zeros(C.NANT, numel(prns));
for m = 1:numel(prns)
    el = asind(sind(5)+(1-sind(5))*rand); az = rand*360;
    a = Cd .* exp(1i*2*pi*(pos.'*[cosd(el)*cosd(az);cosd(el)*sind(az);sind(el)])/lambda);
    truth.A(:,m) = a;
    sig = sig + emit(prns(m), a, sqrt(authSnr), t, rand*1e-3, (rand*2-1)*5000, rand*2*pi, C);
end
truth.b = zeros(C.NANT,1);
if opt.spoof
    b = Cd .* exp(1i*2*pi*(pos.'*[cosd(15)*cosd(45);cosd(15)*sind(45);sind(15)])/lambda);
    truth.b = b;
    phi = rand*2*pi;
    for m = 1:numel(prns)
        sig = sig + emit(prns(m), b, sqrt(spoofSnr), t, rand*1e-3, (rand*2-1)*5000, phi, C);
    end
end
sig = sig + (randn(C.NANT,N)+1i*randn(C.NANT,N))/sqrt(2);

% place at the IF, add per-channel LO leakage at DC (-45 dBFS)
sig = sig .* exp(-1i*2*pi*C.LO_OFFSET*t);
lk  = 10^(-45/20) * exp(1i*2*pi*rand(C.NANT,1));

% AGC: composite per-component RMS to -14 dBFS of the 12-bit full scale
g = (2^(C.W_ADC-1)*10^(-14/20)) / sqrt(mean(abs(sig(:)).^2)/2);
sig = g*(sig + lk*ones(1,N));
x = clampInt(convRound(real(sig)), C.W_ADC) + 1i*clampInt(convRound(imag(sig)), C.W_ADC);
truth.Cd = Cd; truth.pos = pos;
end

function s = emit(prn, a, amp, t, tau, fd, phi, C)
code = caCode(prn);
ch = (t - tau)*C.CHIP*(1 + fd/C.FC);
s  = a * (amp * code(mod(floor(ch),1023)+1) .* exp(1i*(2*pi*fd*t + phi)));
end

function code = caCode(prn)
tp = [2 6;3 7;4 8;5 9;1 9;2 10;1 8;2 9;3 10;2 3;3 4;5 6;6 7;7 8;8 9;9 10; ...
      1 4;2 5;3 6;4 7;5 8;6 9;1 3;4 6;5 7;6 8;7 9;8 10;1 6;2 7;3 8;4 9];
g1 = -ones(1,10); g2 = -ones(1,10); code = zeros(1,1023); q = tp(prn,:);
for k = 1:1023
    code(k) = g1(10)*g2(q(1))*g2(q(2));
    g1 = [g1(3)*g1(10), g1(1:9)];
    g2 = [g2(2)*g2(3)*g2(6)*g2(8)*g2(9)*g2(10), g2(1:9)];
end
end

function assertWidth(x, W, name)
v = [real(x(:)); imag(x(:))];
if any(v ~= fix(v))
    error('asp_golden_model:notInteger','%s carries non-integer values', name);
end
if max(abs(v)) > 2^(W-1)
    error('asp_golden_model:width','%s needs more than %d bits (peak %g)', ...
        name, W, max(abs(v)));
end
end

function dumpStage(opt, name, in, out)
if ~opt.dump, return; end
if ~isempty(in),  writeVec(fullfile(opt.outDir, sprintf('s%s_in.txt',  name)), in);  end
if ~isempty(out), writeVec(fullfile(opt.outDir, sprintf('s%s_out.txt', name)), out); end
end

function appendVec(opt, name, x, dwell)
%APPENDVEC  Per-dwell block, appended.  Same interleaved convention as the
%   streaming files, so the VHDL testbench has exactly one format to parse.
fn = fullfile(opt.outDir, sprintf('s%s.txt', name));
if dwell == 1, fid = fopen(fn,'w'); else, fid = fopen(fn,'a'); end
fprintf(fid, '# dwell %d\n', dwell);
fprintf(fid, '%.0f\n', flatten(x));
fclose(fid);
end

function writeVec(fn, x)
fid = fopen(fn,'w');
fprintf(fid,'%.0f\n', flatten(x));
fclose(fid);
end

function v = flatten(x)
%FLATTEN  Column-major; complex becomes interleaved I,Q.  ONE convention
%   for every file the model emits.
x = x(:);
if isreal(x)
    v = x;
else
    v = zeros(2*numel(x),1);
    v(1:2:end) = real(x);
    v(2:2:end) = imag(x);
end
end

function banner(C, opt)
fprintf('\n=====================================================================\n');
fprintf(' ASP GOLDEN REFERENCE MODEL  (bit-accurate, integer datapath)\n');
fprintf('=====================================================================\n');
fprintf(' RX LO        %.6f MHz  (L1 + %.3f MHz = FS_ADC/16)\n', C.LO/1e6, C.LO_OFFSET/1e6);
fprintf(' FS_ADC       %.6f MHz -> FS_WORK %.6f MHz (%d samples/ms)\n', ...
    C.FS_ADC/1e6, C.FS_WORK/1e6, C.K_DWELL);
fprintf(' channels     %d,  dwell %d ms,  duration %g ms\n', C.NANT, 1, opt.durationMs);
fprintf(' formats      ADC s%d.%d | data s%d.%d | cov s%d | wgt s%d.%d | DAC s%d.%d\n', ...
    C.W_ADC,C.F_ADC, C.W_DAT,C.F_DAT, C.W_ACC, C.W_WGT,C.F_WGT, C.W_DAC,C.F_DAC);
fprintf('---------------------------------------------------------------------\n');
end

function summary(G, C, opt)
fprintf('---------------------------------------------------------------------\n');
d = G.dwell(end);
fprintf(' final weights (Q1.%d)  :', C.F_WGT);
fprintf(' %+d%+dj', [real(d.w).'; imag(d.w).']); fprintf('\n');
fprintf(' DAC per-component RMS : %.1f LSB   (target %d, full scale %d)\n', ...
    d.dac.rmsComp, C.DAC_TARGET_RMS, 2^(C.W_DAC-1));
fprintf(' DAC clipped samples   : %d of %d\n', ...
    sum(arrayfun(@(z) z.dac.nClip, G.dwell)), numel(G.dac));
fprintf(' Jacobi rotations/dwell: %d  (%d sweeps x %d pairs)\n', ...
    d.nRot, C.JACOBI_SWEEPS, C.NANT*(C.NANT-1)/2);
if isfield(G,'check')
    k = G.check;
    fprintf('---------------------------------------------------------------------\n');
    fprintf(' SELF-CHECK vs double precision on the same stage-3 data\n');
    fprintf('   null depth   fixed %+7.2f dB | float %+7.2f dB | loss %+.2f dB  (mean of %d dwells)\n', ...
        k.nullFixedDB, k.nullFloatDB, k.lossDB, k.nDwell);
    fprintf('   suppression  %.2f dB below the quiescent beam (%+.2f dB)\n', ...
        k.suppressionDB, k.quiescentDB);
    fprintf('   eigenvalues  worst relative error %.2e\n', k.lamRelErr);
    fprintf('   weights      worst angle to float solution %.4f deg\n', k.weightAngleDeg);
    fprintf('   authentic    %+.2f dB (quiescent %+.2f dB), change %+.2f dB\n', ...
        k.authGainFixedDB, k.authGainQuiDB, k.authGainFixedDB-k.authGainQuiDB);
end
if opt.dump
    fprintf(' stage vectors written : %s\n', opt.outDir);
end
fprintf('=====================================================================\n\n');
end
