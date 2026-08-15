function asp_export_unit_vectors(outDir, nRand)
%ASP_EXPORT_UNIT_VECTORS Reference vectors for the RTL arithmetic units.
%
%   ASP_EXPORT_UNIT_VECTORS('build/gold')
%
%   FORMAT: one integer per line, same as every other vector file in the
%   project, with '#' comment lines carrying the field order.  Packing
%   several fields onto one line would need a second parser in the VHDL
%   testbench package for no benefit.
%
%   The stage testbenches compare whole blocks against the golden model's
%   dumped stage vectors, which is the acceptance criterion.  But when a
%   stage fails, the question is WHICH primitive is wrong, and a 4x4
%   whitened matrix does not answer it.  These files pin the arithmetic
%   units down individually:
%
%       unit_rsqrt.txt   a, Y, k          (reciprocal square root)
%       unit_isqrt.txt   a, s             (integer square root)
%       unit_cordic_v.txt x, y, m, z      (CORDIC vectoring)
%       unit_cordic_r.txt z, c, s         (CORDIC rotation)
%
%   The MATLAB side here is a faithful re-implementation of the local
%   functions in asp_golden_model.m.  That duplication is deliberate and
%   is checked, not assumed: the whiten and EVD testbenches compare
%   against the model's OWN dumps, so if this re-implementation ever
%   drifts from the model the stage tests fail even when the unit tests
%   pass.  The unit vectors are a debugging aid, not the contract.

if nargin < 1 || isempty(outDir), outDir = 'build/gold'; end
if nargin < 2 || isempty(nRand),  nRand  = 400; end
if ~exist(outDir,'dir'), mkdir(outDir); end

C.RSQ_F      = 16;
C.CORDIC_N   = 16;
C.ANG_F      = 16;
C.CORDIC_INV_K = 39797;
C.W_ROT      = 18;
C.F_ROT      = 16;

rng(7);

% ---------------------------------------------------------------- rsqrt
% Cover the whole dynamic range the block actually sees: the whitening
% feeds it R_ii after block normalisation (~2^16), the weight stage feeds
% it a squared norm that can reach ~2^40.
a = [1 2 3 4 5 7 8 9 15 16 17 255 256 257 (2^15-1) 2^15 (2^15+1) ...
     (2^16-1) 2^16 (2^16+1) (2^20) (2^24) (2^31) (2^32) (2^39) (2^40)];
a = [a, round(2.^(rand(1,nRand)*40) )];
a = unique(max(a,1));
fid = fopen(fullfile(outDir,'unit_rsqrt.txt'),'w');
fprintf(fid, '# one value per line, repeating: a, Y, k\n');
fprintf(fid, '#   1/sqrt(a) = Y * 2^-(3*RSQ_F/2 + k)\n');
for v = a
    [Y,k] = rsqrtNorm(v, C);
    fprintf(fid, '%.0f\n%d\n%d\n', v, Y, k);
end
fclose(fid);

% ---------------------------------------------------------------- isqrt
b = [0 1 2 3 4 8 15 16 17 99 100 (2^24-1) 2^24 (2^31) (2^40) (2^47)];
b = [b, round(2.^(rand(1,nRand)*47))];
b = unique(max(b,0));
fid = fopen(fullfile(outDir,'unit_isqrt.txt'),'w');
fprintf(fid, '# one value per line, repeating: a, s   (s = floor(sqrt(a)))\n');
for v = b
    fprintf(fid, '%.0f\n%.0f\n', v, isqrtInt(v));
end
fclose(fid);

% ------------------------------------------------------- CORDIC vectoring
% Includes every quadrant and the axes, which is where the pre-rotation
% folding is easiest to get wrong.
% Magnitudes up to 2^29, which is what the Jacobi engine actually
% presents: |A| <= 2^28 for the whitened matrix, and the second
% vectoring is called with 2*m.
xs = [1 -1 0 1000 -1000 32767 -32768 12345 -12345 2^29 -(2^29) 2^28];
ys = [0 0 1 1000 -1000 -32768 32767 -54321 54321 2^29 -(2^29) -(2^28)];
pts = [];
for ix = 1:numel(xs)
    for iy = 1:numel(ys)
        pts(end+1,:) = [xs(ix) ys(iy)]; %#ok<AGROW>
    end
end
r = round((rand(nRand,2)*2-1) .* 2.^(rand(nRand,2)*29));
pts = [pts; r];
fid = fopen(fullfile(outDir,'unit_cordic_v.txt'),'w');
fprintf(fid, '# one value per line, repeating: x, y, m, z\n');
fprintf(fid, '#   m = |.| gain compensated, z = atan2(y,x)*2^ANG_F\n');
for i = 1:size(pts,1)
    [m,z] = cordicVec(pts(i,1), pts(i,2), C);
    fprintf(fid, '%d\n%d\n%.0f\n%.0f\n', pts(i,1), pts(i,2), m, z);
end
fclose(fid);

% ------------------------------------------------------- CORDIC rotation
PI_Q = round(pi*2^C.ANG_F);
zs = [0 1 -1 PI_Q -PI_Q round(PI_Q/2) -round(PI_Q/2) round(PI_Q/4)];
zs = [zs, round((rand(1,nRand)*2-1)*PI_Q)];
fid = fopen(fullfile(outDir,'unit_cordic_r.txt'),'w');
fprintf(fid, '# one value per line, repeating: z, c, s\n');
fprintf(fid, '#   c = cos(z)*2^F_ROT, s = sin(z)*2^F_ROT\n');
for z = zs
    [c,s] = cordicRot(z, C);
    fprintf(fid, '%.0f\n%d\n%d\n', z, c, s);
end
fclose(fid);

fprintf('asp_export_unit_vectors: wrote %d rsqrt, %d isqrt, %d cordicVec, %d cordicRot\n', ...
    numel(a), numel(b), size(pts,1), numel(zs));
end


% =====================================================================
%  Re-implementations, identical to the local functions in the model.
% =====================================================================

function [Y, k] = rsqrtNorm(a, C)
F = C.RSQ_F;
if a <= 0, Y = 0; k = 0; return; end
e = floor(log2(a));
k = floor((e - (F-1))/2);
an = round(a / 2^(2*k));
while an >= 2^F,      k = k + 1; an = round(a / 2^(2*k)); end
while an <  2^(F-2),  k = k - 1; an = round(a / 2^(2*k)); end
Y = fix((7*2^F - 4*an) / 3);
Y = min(max(Y, 2^F), 2^(F+1));
for it = 1:4
    t = fix(an * Y / 2^F);
    t = fix(t  * Y / 2^F);
    Y = fix(Y * (3*2^F - t) / 2^(F+1));
end
Y = min(max(Y, 1), 2^(F+1));
end

function s = isqrtInt(a)
if a <= 0, s = 0; return; end
s = 0; rem = 0;
for i = 23:-1:0
    rem = rem*4 + fix(mod(fix(a/4^i), 4));
    t = 4*s + 1;
    if rem >= t, rem = rem - t; s = 2*s + 1; else, s = 2*s; end
end
end

function [m, z] = cordicVec(x, y, C)
atanTab = round(atan(2.^-(0:C.CORDIC_N-1)) * 2^C.ANG_F);
PI_Q  = round(pi*2^C.ANG_F);
HPI_Q = round(pi/2*2^C.ANG_F);
x = round(x); y = round(y);
z = 0;
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
m = fix(x * C.CORDIC_INV_K / 2^16);
if z >  PI_Q, z = z - 2*PI_Q; end
if z < -PI_Q, z = z + 2*PI_Q; end
end

function [c, s] = cordicRot(z, C)
atanTab = round(atan(2.^-(0:C.CORDIC_N-1)) * 2^C.ANG_F);
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

function y = clampInt(x, W)
hi = 2^(W-1) - 1; lo = -2^(W-1);
y = min(max(x, lo), hi);
end
