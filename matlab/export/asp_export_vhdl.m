function asp_export_vhdl(outFile)
%ASP_EXPORT_VHDL Generate the VHDL constant package from the golden model.
%
%   ASP_EXPORT_VHDL()                 % writes rtl/pkg/asp_coef_pkg.vhd
%   ASP_EXPORT_VHDL('path/pkg.vhd')
%
%   Every word length, coefficient and table in the RTL comes from HERE, and
%   here reads them out of asp_golden_model.  Hand-transcribing them into
%   VHDL is the obvious way to do it and it is wrong: the model and the RTL
%   then drift silently, and the symptom is an RTL "bug" that is actually a
%   constant mismatch.  Regenerate this file whenever the model changes and
%   the divergence becomes a compile-time fact instead of a simulation
%   mystery.

if nargin < 1 || isempty(outFile)
    here = fileparts(fileparts(mfilename('fullpath')));   % <repo>/matlab
    outFile = fullfile(fileparts(here), 'rtl', 'pkg', 'asp_coef_pkg.vhd');
end
d = fileparts(outFile);
if ~isempty(d) && ~exist(d,'dir'), mkdir(d); end

G = asp_golden_model('durationMs',2, 'dump',false, 'verbose',false);
C = G.C;

nco = G.ncoTab(:);
hb  = G.hbCoef(:);
fir = G.firCoef(:);
atanTab = round(atan(2.^-(0:C.CORDIC_N-1)) * 2^C.ANG_F);

fid = fopen(outFile,'w');
w = @(varargin) fprintf(fid, varargin{:});

w('-- =====================================================================\n');
w('--  asp_coef_pkg  -  constants and coefficient tables for the GNSS\n');
w('--                  anti-spoofing array processor.\n');
w('--\n');
w('--  GENERATED FILE - DO NOT EDIT BY HAND.\n');
w('--  Produced by matlab/export/asp_export_vhdl.m from the golden model\n');
w('--  matlab/golden/asp_golden_model.m, which is the sole functional\n');
w('--  reference for this design.  Regenerate after any model change:\n');
w('--      octave-cli --eval "asp_startup; asp_export_vhdl"\n');
w('--\n');
w('--  Every value below is an EXACT INTEGER, identical to the raw register\n');
w('--  contents the model computes, so RTL and model outputs are compared\n');
w('--  with no tolerance.\n');
w('-- =====================================================================\n');
w('library ieee;\n');
w('use ieee.std_logic_1164.all;\n');
w('use ieee.numeric_std.all;\n\n');
w('package asp_coef_pkg is\n\n');

w('  -- ---------------------------------------------------------------\n');
w('  -- Rate plan.  All rates are integer ratios of the ADC clock, so the\n');
w('  -- entire datapath is ONE synchronous clock domain with clock enables\n');
w('  -- and there is no true CDC anywhere between stage 1 and stage 10.\n');
w('  -- ---------------------------------------------------------------\n');
w('  constant FS_ADC_HZ    : natural := %d;   -- %.3f MHz, 32 x 1.023 MHz\n', round(C.FS_ADC), C.FS_ADC/1e6);
w('  constant FS_WORK_HZ   : natural := %d;   -- %.3f MHz after the halfband\n', round(C.FS_WORK), C.FS_WORK/1e6);
w('  constant CLK_DSP_HZ   : natural := %d;  -- 4 x FS_ADC = 8 x FS_WORK\n', round(4*C.FS_ADC));
w('  constant N_ANT        : natural := %d;\n', C.NANT);
w('  constant K_DWELL      : natural := %d;  -- samples per 1 ms dwell at FS_WORK\n', C.K_DWELL);
w('  constant LOG2_K_DWELL : natural := %d;  -- ceil(log2(K_DWELL))\n', ceil(log2(C.K_DWELL)));
w('\n');

w('  -- ---------------------------------------------------------------\n');
w('  -- Word lengths.  W = total bits including sign, F = fractional bits.\n');
w('  -- Data are carried as raw two''s-complement integers; F only records\n');
w('  -- where the binary point sits so the shifts can be derived.\n');
w('  -- ---------------------------------------------------------------\n');
flds = {'W_ADC','F_ADC','W_NCO','F_NCO','W_MIX','F_MIX','W_HB','F_HB', ...
        'W_COEF','F_COEF','W_DAT','F_DAT','W_ACC','W_EVD','F_EVD', ...
        'W_ROT','F_ROT','W_UVEC','RSQ_F','R_NORM_BITS','W_WGT','F_WGT', ...
        'W_BEAM','F_BEAM','W_DAC','F_DAC'};
for k = 1:numel(flds)
    w('  constant %-12s : natural := %d;\n', flds{k}, C.(flds{k}));
end
w('\n');

w('  -- Derived shift amounts.  Each is (input F) + (coefficient F) - (output F).\n');
w('  constant SH_MIX  : natural := %d;  -- F_ADC + F_NCO  - F_MIX\n',  C.F_ADC + C.F_NCO  - C.F_MIX);
w('  constant SH_HB   : natural := %d;  -- F_MIX + F_COEF - F_HB\n',   C.F_MIX + C.F_COEF - C.F_HB);
w('  constant SH_FIR  : natural := %d;  -- F_HB  + F_COEF - F_DAT\n',  C.F_HB  + C.F_COEF - C.F_DAT);
w('  constant SH_BEAM : natural := %d;  -- F_DAT + F_WGT  - F_BEAM\n', C.F_DAT + C.F_WGT  - C.F_BEAM);
w('  constant SH_EVD_BASE : integer := %d;  -- 3*RSQ_F - F_EVD, stage 5\n', 3*C.RSQ_F - C.F_EVD);
w('\n');

w('  -- ---------------------------------------------------------------\n');
w('  -- CORDIC.  Fixed iteration count, so latency is deterministic and\n');
w('  -- the block needs no tolerance test and no data-dependent loop.\n');
w('  -- ---------------------------------------------------------------\n');
w('  constant CORDIC_N     : natural := %d;\n', C.CORDIC_N);
w('  constant CORDIC_W     : natural := %d;\n', C.CORDIC_W);
w('  constant ANG_F        : natural := %d;   -- angle scaling: radians * 2^ANG_F\n', C.ANG_F);
w('  constant CORDIC_INV_K : natural := %d;   -- round(0.607252935 * 2^16)\n', C.CORDIC_INV_K);
w('  constant ANG_PI       : integer := %d;\n', round(pi*2^C.ANG_F));
w('  constant ANG_HPI      : integer := %d;\n', round(pi/2*2^C.ANG_F));
w('\n');

w('  -- ---------------------------------------------------------------\n');
w('  -- Algorithm constants.\n');
w('  -- ---------------------------------------------------------------\n');
w('  constant JACOBI_SWEEPS : natural := %d;\n', C.JACOBI_SWEEPS);
w('  constant JACOBI_PAIRS  : natural := %d;   -- N*(N-1)/2\n', C.NANT*(C.NANT-1)/2);
w('  constant JACOBI_ROTS   : natural := %d;  -- SWEEPS * PAIRS, FIXED\n', C.JACOBI_SWEEPS*C.NANT*(C.NANT-1)/2);
w('  constant MAX_RANK      : natural := %d;\n', C.MAX_RANK);
w('  -- Detector threshold as a rational so no divider is needed.  1280/1024\n');
w('  -- = 5/4 exactly, so stage 7 needs no multiplier either: the test is\n');
w('  -- lam1*3*4 > 5*sum(tail), i.e. (x sll 2) and (x sll 2) + x.\n');
w('  constant DET_NUM       : natural := %d;\n', C.DET_NUM);
w('  constant DET_DEN       : natural := %d;\n', C.DET_DEN);
w('  constant DET_NUM2      : natural := %d;   -- second-null threshold, 1.30\n', C.DET_NUM2);
w('\n');
w('  constant DAC_TARGET_RMS : natural := %d;\n', C.DAC_TARGET_RMS);
w('  constant DAC_AGC_LO     : natural := %d;   -- TARGET / sqrt(2)\n', C.DAC_AGC_LO);
w('  constant DAC_AGC_HI     : natural := %d;   -- TARGET * sqrt(2)\n', C.DAC_AGC_HI);
w('\n');

w('  -- ---------------------------------------------------------------\n');
w('  -- Types.\n');
w('  -- ---------------------------------------------------------------\n');
w('  type int_array_t  is array (natural range <>) of integer;\n');
w('\n');

w('  -- ---------------------------------------------------------------\n');
w('  -- NCO table, Q1.%d.  16 entries indexed by n mod 16, so the mixer\n', C.F_NCO);
w('  -- needs a 4-bit counter and NO phase accumulator - hence no phase\n');
w('  -- truncation spurs at all, which a DDS Compiler could not promise.\n');
w('  -- exp(+j*2*pi*n/16); the mixer multiplies by this directly.\n');
w('  -- ---------------------------------------------------------------\n');
w('  constant NCO_RE : int_array_t(0 to 15) := (\n    ');
w('%s);\n', joinInts(real(nco), 8));
w('  constant NCO_IM : int_array_t(0 to 15) := (\n    ');
w('%s);\n\n', joinInts(imag(nco), 8));

w('  -- ---------------------------------------------------------------\n');
w('  -- Halfband decimator, %d taps, Q1.%d, DC gain exactly 2^%d.\n', numel(hb), C.F_COEF, C.F_COEF);
w('  -- EVERY non-zero coefficient is a sum of at most two powers of two:\n');
for k = 1:numel(hb)
    if hb(k) ~= 0
        w('  --   h(%d) = %+7d = %s\n', k-1, hb(k), powersOfTwo(hb(k)));
    end
end
w('  -- so the halfband is implemented with shifts and adds and uses\n');
w('  -- ZERO DSP48 slices.  That is why it is custom RTL and not FIR\n');
w('  -- Compiler, which would spend multipliers on it.\n');
w('  -- ---------------------------------------------------------------\n');
w('  constant HB_NTAPS : natural := %d;\n', numel(hb));
w('  constant HB_COEF  : int_array_t(0 to %d) := (\n    ', numel(hb)-1);
w('%s);\n\n', joinInts(hb, 8));

w('  -- ---------------------------------------------------------------\n');
w('  -- Shaping FIR, %d taps, symmetric, Q1.%d, DC gain exactly 2^%d.\n', numel(fir), C.F_COEF, C.F_COEF);
w('  -- Passband +/-1.2 MHz, stopband from %.3f MHz (where the residual LO\n', C.LO_OFFSET/1e6);
w('  -- leakage lands after the DDC), >= 70 dB.\n');
w('  -- Symmetric, so %d multipliers rather than %d.\n', (numel(fir)+1)/2, numel(fir));
w('  -- ---------------------------------------------------------------\n');
w('  constant FIR_NTAPS : natural := %d;\n', numel(fir));
w('  constant FIR_NMULT : natural := %d;   -- (NTAPS+1)/2 with the centre tap\n', (numel(fir)+1)/2);
w('  constant FIR_COEF  : int_array_t(0 to %d) := (\n    ', numel(fir)-1);
w('%s);\n\n', joinInts(fir, 8));

w('  -- ---------------------------------------------------------------\n');
w('  -- CORDIC arctangent table: round(atan(2^-i) * 2^%d).\n', C.ANG_F);
w('  -- ---------------------------------------------------------------\n');
w('  constant CORDIC_ATAN : int_array_t(0 to %d) := (\n    ', numel(atanTab)-1);
w('%s);\n\n', joinInts(atanTab, 8));

w('  -- ---------------------------------------------------------------\n');
w('  -- Quiescent beam h = (2^(F_WGT-1), 0) on every element.  NOT scaled\n');
w('  -- by 1/sqrt(N): the whole chain is invariant to a common real scale,\n');
w('  -- so the division would cost hardware and buy nothing.\n');
w('  -- ---------------------------------------------------------------\n');
w('  constant WGT_QUIESCENT : integer := %d;\n\n', 2^(C.F_WGT-1));

w('end package asp_coef_pkg;\n');
fclose(fid);

fprintf('asp_export_vhdl: wrote %s\n', outFile);
fprintf('  %d NCO entries, %d halfband taps, %d FIR taps, %d CORDIC angles\n', ...
    numel(nco), numel(hb), numel(fir), numel(atanTab));
end


function s = joinInts(v, perLine)
v = round(real(v(:)));
parts = cell(1,numel(v));
for k = 1:numel(v)
    parts{k} = sprintf('%d', v(k));
end
s = '';
for k = 1:numel(parts)
    s = [s parts{k}]; %#ok<AGROW>
    if k < numel(parts)
        s = [s ', ']; %#ok<AGROW>
        if mod(k, perLine) == 0, s = [s sprintf('\n    ')]; end %#ok<AGROW>
    end
end
end


function s = powersOfTwo(v)
%POWERSOFTWO  Render an integer as a signed sum of powers of two, which is
%   what the shift-add implementation actually builds.
neg = v < 0; a = abs(v); terms = {};
b = 0;
while a > 0
    if bitand(a,1) == 1, terms{end+1} = sprintf('2^%d', b); end %#ok<AGROW>
    a = floor(a/2); b = b + 1;
end
terms = fliplr(terms);
s = strjoin(terms, ' + ');
if neg, s = ['-(' s ')']; end
end
