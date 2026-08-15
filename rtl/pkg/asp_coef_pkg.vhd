-- =====================================================================
--  asp_coef_pkg  -  constants and coefficient tables for the GNSS
--                  anti-spoofing array processor.
--
--  GENERATED FILE - DO NOT EDIT BY HAND.
--  Produced by matlab/export/asp_export_vhdl.m from the golden model
--  matlab/golden/asp_golden_model.m, which is the sole functional
--  reference for this design.  Regenerate after any model change:
--      octave-cli --eval "asp_startup; asp_export_vhdl"
--
--  Every value below is an EXACT INTEGER, identical to the raw register
--  contents the model computes, so RTL and model outputs are compared
--  with no tolerance.
-- =====================================================================
library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

package asp_coef_pkg is

  -- ---------------------------------------------------------------
  -- Rate plan.  All rates are integer ratios of the ADC clock, so the
  -- entire datapath is ONE synchronous clock domain with clock enables
  -- and there is no true CDC anywhere between stage 1 and stage 10.
  -- ---------------------------------------------------------------
  constant FS_ADC_HZ    : natural := 32736000;   -- 32.736 MHz, 32 x 1.023 MHz
  constant FS_WORK_HZ   : natural := 16368000;   -- 16.368 MHz after the halfband
  constant CLK_DSP_HZ   : natural := 130944000;  -- 4 x FS_ADC = 8 x FS_WORK
  constant N_ANT        : natural := 4;
  constant K_DWELL      : natural := 16368;  -- samples per 1 ms dwell at FS_WORK
  constant LOG2_K_DWELL : natural := 14;  -- ceil(log2(K_DWELL))

  -- ---------------------------------------------------------------
  -- Word lengths.  W = total bits including sign, F = fractional bits.
  -- Data are carried as raw two's-complement integers; F only records
  -- where the binary point sits so the shifts can be derived.
  -- ---------------------------------------------------------------
  constant W_ADC        : natural := 12;
  constant F_ADC        : natural := 11;
  constant W_NCO        : natural := 16;
  constant F_NCO        : natural := 14;
  constant W_MIX        : natural := 16;
  constant F_MIX        : natural := 15;
  constant W_HB         : natural := 16;
  constant F_HB         : natural := 15;
  constant W_COEF       : natural := 18;
  constant F_COEF       : natural := 17;
  constant W_DAT        : natural := 16;
  constant F_DAT        : natural := 15;
  constant W_ACC        : natural := 48;
  constant W_EVD        : natural := 32;
  constant F_EVD        : natural := 26;
  constant W_ROT        : natural := 18;
  constant F_ROT        : natural := 16;
  constant W_UVEC       : natural := 20;
  constant RSQ_F        : natural := 16;
  constant R_NORM_BITS  : natural := 16;
  constant W_WGT        : natural := 18;
  constant F_WGT        : natural := 16;
  constant W_BEAM       : natural := 16;
  constant F_BEAM       : natural := 15;
  constant W_DAC        : natural := 12;
  constant F_DAC        : natural := 11;

  -- Derived shift amounts.  Each is (input F) + (coefficient F) - (output F).
  constant SH_MIX  : natural := 10;  -- F_ADC + F_NCO  - F_MIX
  constant SH_HB   : natural := 17;  -- F_MIX + F_COEF - F_HB
  constant SH_FIR  : natural := 17;  -- F_HB  + F_COEF - F_DAT
  constant SH_BEAM : natural := 16;  -- F_DAT + F_WGT  - F_BEAM
  constant SH_EVD_BASE : integer := 22;  -- 3*RSQ_F - F_EVD, stage 5

  -- ---------------------------------------------------------------
  -- CORDIC.  Fixed iteration count, so latency is deterministic and
  -- the block needs no tolerance test and no data-dependent loop.
  -- ---------------------------------------------------------------
  constant CORDIC_N     : natural := 16;
  constant CORDIC_W     : natural := 22;
  constant ANG_F        : natural := 16;   -- angle scaling: radians * 2^ANG_F
  constant CORDIC_INV_K : natural := 39797;   -- round(0.607252935 * 2^16)
  constant ANG_PI       : integer := 205887;
  constant ANG_HPI      : integer := 102944;

  -- ---------------------------------------------------------------
  -- Algorithm constants.
  -- ---------------------------------------------------------------
  constant JACOBI_SWEEPS : natural := 6;
  constant JACOBI_PAIRS  : natural := 6;   -- N*(N-1)/2
  constant JACOBI_ROTS   : natural := 36;  -- SWEEPS * PAIRS, FIXED
  constant MAX_RANK      : natural := 2;
  -- Detector threshold as a rational so no divider is needed.  1280/1024
  -- = 5/4 exactly, so stage 7 needs no multiplier either: the test is
  -- lam1*3*4 > 5*sum(tail), i.e. (x sll 2) and (x sll 2) + x.
  constant DET_NUM       : natural := 1280;
  constant DET_DEN       : natural := 1024;
  constant DET_NUM2      : natural := 1331;   -- second-null threshold, 1.30

  constant DAC_TARGET_RMS : natural := 256;
  constant DAC_AGC_LO     : natural := 181;   -- TARGET / sqrt(2)
  constant DAC_AGC_HI     : natural := 362;   -- TARGET * sqrt(2)

  -- ---------------------------------------------------------------
  -- Types.
  -- ---------------------------------------------------------------
  type int_array_t  is array (natural range <>) of integer;

  -- ---------------------------------------------------------------
  -- NCO table, Q1.14.  16 entries indexed by n mod 16, so the mixer
  -- needs a 4-bit counter and NO phase accumulator - hence no phase
  -- truncation spurs at all, which a DDS Compiler could not promise.
  -- exp(+j*2*pi*n/16); the mixer multiplies by this directly.
  -- ---------------------------------------------------------------
  constant NCO_RE : int_array_t(0 to 15) := (
    16384, 15137, 11585, 6270, 0, -6270, -11585, -15137, 
    -16384, -15137, -11585, -6270, 0, 6270, 11585, 15137);
  constant NCO_IM : int_array_t(0 to 15) := (
    0, 6270, 11585, 15137, 16384, 15137, 11585, 6270, 
    0, -6270, -11585, -15137, -16384, -15137, -11585, -6270);

  -- ---------------------------------------------------------------
  -- Halfband decimator, 11 taps, Q1.17, DC gain exactly 2^17.
  -- EVERY non-zero coefficient is a sum of at most two powers of two:
  --   h(2) =   -4096 = -(2^12)
  --   h(4) =  +36864 = 2^15 + 2^12
  --   h(5) =  +65536 = 2^16
  --   h(6) =  +36864 = 2^15 + 2^12
  --   h(8) =   -4096 = -(2^12)
  -- so the halfband is implemented with shifts and adds and uses
  -- ZERO DSP48 slices.  That is why it is custom RTL and not FIR
  -- Compiler, which would spend multipliers on it.
  -- ---------------------------------------------------------------
  constant HB_NTAPS : natural := 11;
  constant HB_COEF  : int_array_t(0 to 10) := (
    0, 0, -4096, 0, 36864, 65536, 36864, 0, 
    -4096, 0, 0);

  -- ---------------------------------------------------------------
  -- Shaping FIR, 63 taps, symmetric, Q1.17, DC gain exactly 2^17.
  -- Passband +/-1.2 MHz, stopband from 2.046 MHz (where the residual LO
  -- leakage lands after the DDC), >= 70 dB.
  -- Symmetric, so 32 multipliers rather than 63.
  -- ---------------------------------------------------------------
  constant FIR_NTAPS : natural := 63;
  constant FIR_NMULT : natural := 32;   -- (NTAPS+1)/2 with the centre tap
  constant FIR_COEF  : int_array_t(0 to 62) := (
    0, -1, -5, -13, -19, -10, 22, 75, 
    130, 144, 78, -88, -318, -520, -559, -322, 
    214, 919, 1512, 1640, 1021, -381, -2244, -3885, 
    -4429, -3086, 530, 6165, 12896, 19319, 23940, 25622, 
    23940, 19319, 12896, 6165, 530, -3086, -4429, -3885, 
    -2244, -381, 1021, 1640, 1512, 919, 214, -322, 
    -559, -520, -318, -88, 78, 144, 130, 75, 
    22, -10, -19, -13, -5, -1, 0);

  -- ---------------------------------------------------------------
  -- CORDIC arctangent table: round(atan(2^-i) * 2^16).
  -- ---------------------------------------------------------------
  constant CORDIC_ATAN : int_array_t(0 to 15) := (
    51472, 30386, 16055, 8150, 4091, 2047, 1024, 512, 
    256, 128, 64, 32, 16, 8, 4, 2);

  -- ---------------------------------------------------------------
  -- Quiescent beam h = (2^(F_WGT-1), 0) on every element.  NOT scaled
  -- by 1/sqrt(N): the whole chain is invariant to a common real scale,
  -- so the division would cost hardware and buy nothing.
  -- ---------------------------------------------------------------
  constant WGT_QUIESCENT : integer := 32768;

end package asp_coef_pkg;
