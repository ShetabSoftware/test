-- =====================================================================
--  asp_tx_scale  -  STAGE 10.  Digital AGC and 12-bit DAC formatting.
--
--  MODEL REFERENCE : stage10_txdac() in matlab/golden/asp_golden_model.m
--
--  OUTPUT FORMAT
--    SIGNED 12-bit, s12.11, integer range [-2048, +2047] - the AD9361 TX
--    data port.  Operating point: per-component RMS held at 256 LSB, i.e.
--    -18.1 dBFS, giving 8 sigma of crest headroom (clip probability
--    below 1e-15 per sample).
--
--  WHY 18 dB OF BACKOFF ON TX WHEN RX USES 14
--    The two sides solve opposite problems.  On RX, backoff buys
--    headroom for an interferer that has not been removed yet.  On TX
--    the interferer is ALREADY nulled, so there is nothing to leave room
--    for, and the only costs of backing off further are quantisation
--    noise - 10*log10(1 + (1/12)/256^2) = 5e-6 dB, i.e. nothing - and
--    output power, which the analogue attenuator sets anyway.  Meanwhile
--    clipping is a memoryless nonlinearity that would intermodulate the
--    residual spoofer back into the band and undo the nulling.  On TX,
--    clipping margin is free: buy it.
--
--  GAIN CONTROL IS A SHIFT, NOT A MULTIPLIER
--    A barrel shifter and a comparison, with hysteresis at TARGET*sqrt2
--    and TARGET/sqrt2.  It adds EXACTLY zero error, unlike a fractional
--    gain.  Amplitude steps are harmless here (unlike the RX side, where
--    a per-channel gain step corrupts the array), so a coarse
--    power-of-two step is sufficient.
--
--  ROUNDING AND SATURATION
--    Convergent rounding, because a truncating output would put a
--    -0.5 LSB DC offset on the DAC, which becomes an LO-leakage-like
--    spur at the TX LO and can disturb the downstream receiver's own DC
--    and AGC logic.  Saturation is a hard clamp, NEVER wraparound: a
--    single wrapped sample is a full-scale transient spread across the
--    whole band.
--
--  ================= DEVIATION FROM THE MODEL, DWELL 0 =================
--  The model's first dwell uses "fast acquisition": it measures the RMS
--  of a dwell and then applies the resulting shift TO THAT SAME DWELL.
--  That is not causal and no hardware can do it without buffering the
--  entire dwell (16368 complex samples, ~15 BRAM36) purely to improve
--  the first millisecond after reset.
--
--  This block therefore scales dwell 0 with i_shift_init (an AXI
--  register, default 0) and performs the same acquisition jump at the
--  END of dwell 0.  From dwell 1 onward the shift sequence is identical
--  to the model's, because the acquisition arithmetic and the hysteresis
--  are the same - only the dwell they are first applied to differs.  The
--  testbench skips dwell 0 for exactly this reason and says so.
--
--  Set i_shift_init from a power measurement made during AD9361 RX
--  calibration and even dwell 0 is close; leave it at 0 and dwell 0 may
--  clip, which the o_clip_count status register reports.
--  =====================================================================
--
--  NO SQUARE ROOT ANYWHERE
--    The model works with rmsComp = sqrt(p/(2K)).  Every use of it is a
--    comparison, so this block squares the thresholds instead and stays
--    in exact integer arithmetic:
--        rms > 362*2^sh   <=>   p > 2K*131044*2^(2sh)
--        rms < 181*2^sh   <=>   p < 2K*32761 *2^(2sh)
--    and round(log2(max(rms,1))) = e where 2K*2^(2e-1) <= p < 2K*2^(2e+1),
--    which is a leading-zero comparison rather than a logarithm.  Exact,
--    and it removes both the sqrt and the log from the hardware.
--
--  LATENCY   3 clk per sample; the shift update lands at the dwell tick.
--  RESOURCES 2 DSP48E1 (the power accumulator), ~400 FF, ~700 LUT
--            (dominated by the wide threshold comparators).
-- =====================================================================
library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

library work;
use work.asp_pkg.all;
use work.asp_coef_pkg.all;

entity asp_tx_scale is
  generic (
    G_KDWELL : natural := K_DWELL
  );
  port (
    clk     : in  std_logic;
    rst     : in  std_logic;

    i_shift_init : in  std_logic_vector(7 downto 0);   -- signed, dwell 0

    s_valid : in  std_logic;
    s_re    : in  std_logic_vector(W_BEAM-1 downto 0);
    s_im    : in  std_logic_vector(W_BEAM-1 downto 0);

    m_valid : out std_logic;
    m_re    : out std_logic_vector(W_DAC-1 downto 0);
    m_im    : out std_logic_vector(W_DAC-1 downto 0);

    o_shift : out std_logic_vector(7 downto 0);        -- signed, in use
    o_clip  : out std_logic_vector(31 downto 0);       -- clipped samples
    o_dwell_tick : out std_logic                       -- for observability
  );
end entity asp_tx_scale;


architecture rtl of asp_tx_scale is

  constant PACC_W : natural := 48;
  constant CMPW   : natural := 96;      -- threshold comparison width

  -- 2*K*TARGET^2 scaled constants, all exact integers
  constant HI_C : natural := DAC_AGC_HI * DAC_AGC_HI;   -- 362^2 = 131044
  constant LO_C : natural := DAC_AGC_LO * DAC_AGC_LO;   -- 181^2 =  32761

  -- The dwell boundary is counted HERE rather than taken from stage 4.
  -- Both count K_DWELL from reset off the same sample stream, so they
  -- stay aligned by construction; taking it as an input would instead
  -- make the power accumulation depend on the exact cycle the tick
  -- arrives relative to the sample pipeline, which is a hazard with no
  -- upside.
  signal samp_cnt : unsigned(LOG2_K_DWELL downto 0) := (others => '0');
  signal tick     : std_logic := '0';

  signal pacc  : unsigned(PACC_W-1 downto 0) := (others => '0');
  -- Completed dwell sum, frozen at the boundary.  The accumulator itself
  -- must restart on the very next sample, so the value the AGC reads has
  -- to be captured separately - otherwise the last sample of the dwell
  -- overwrites the sum it was supposed to complete.
  signal pacc_hold : unsigned(PACC_W-1 downto 0) := (others => '0');
  signal shift_cur : integer range -64 to 64 := 0;
  signal acquired  : std_logic := '0';
  signal clip_cnt  : unsigned(31 downto 0) := (others => '0');

  signal d1_val : std_logic := '0';
  signal d1_re, d1_im : signed(W_BEAM-1 downto 0) := (others => '0');
  signal d2_val : std_logic := '0';
  signal d2_re, d2_im : signed(W_DAC-1 downto 0) := (others => '0');
  signal d2_clip : std_logic := '0';

begin

  -- -------------------------------------------------------------------
  -- Per-sample: scale with the CURRENT shift and accumulate power for
  -- the NEXT shift decision.
  -- -------------------------------------------------------------------
  p_sample : process (clk)
    variable sh   : integer;
    variable qr, qi : signed(W_DAC-1 downto 0);
    variable clipped : boolean;
    variable ar, ai  : unsigned(CMPW-1 downto 0);
    variable fs      : unsigned(CMPW-1 downto 0);
  begin
    if rising_edge(clk) then
      if rst = '1' then
        d1_val   <= '0';
        d2_val   <= '0';
        pacc      <= (others => '0');
        pacc_hold <= (others => '0');
        samp_cnt  <= (others => '0');
        tick      <= '0';
      else
        d1_val <= s_valid;
        d1_re  <= signed(s_re);
        d1_im  <= signed(s_im);

        -- Power accumulator, exact.  The dwell tick is raised on the
        -- sample AFTER the last one of a dwell, so the accumulator is
        -- complete when the AGC reads it and restarts on that same
        -- sample: no dead sample and no lost sample at the boundary.
        tick <= '0';
        if d1_val = '1' then
          if samp_cnt = to_unsigned(G_KDWELL-1, samp_cnt'length) then
            -- last sample of the dwell: complete the sum into the hold
            -- register and restart the accumulator from zero, so no
            -- sample is counted twice and none is dropped
            pacc_hold <= pacc + resize(unsigned(d1_re * d1_re) +
                                       unsigned(d1_im * d1_im), PACC_W);
            pacc      <= (others => '0');
            samp_cnt  <= (others => '0');
            tick      <= '1';
          else
            samp_cnt <= samp_cnt + 1;
            pacc     <= pacc + resize(unsigned(d1_re * d1_re) +
                                      unsigned(d1_im * d1_im), PACC_W);
          end if;
        end if;

        -- scale, round, saturate
        sh := F_BEAM - F_DAC + shift_cur;
        d2_val <= d1_val;
        if d1_val = '1' then
          qr := shift_round_sat_v(d1_re, sh, W_DAC);
          qi := shift_round_sat_v(d1_im, sh, W_DAC);
          d2_re <= qr;
          d2_im <= qi;

          -- Clip flag.  The model tests |v| / 2^sh against full scale in
          -- exact real arithmetic, so the comparison is moved to the
          -- other side rather than dividing: |v| > FS * 2^sh.  Dividing
          -- first would truncate and miscount samples that sit exactly
          -- on the boundary.
          ar := resize(abs_ext(d1_re), CMPW);
          ai := resize(abs_ext(d1_im), CMPW);
          fs := to_unsigned(2**(W_DAC-1)-1, CMPW);
          if sh >= 0 then
            fs := shift_left(fs, sh);
          else
            ar := shift_left(ar, -sh);
            ai := shift_left(ai, -sh);
          end if;
          clipped := (ar > fs) or (ai > fs);
          if clipped then
            d2_clip <= '1';
          else
            d2_clip <= '0';
          end if;
        else
          d2_clip <= '0';
        end if;
      end if;
    end if;
  end process p_sample;

  p_clip : process (clk)
  begin
    if rising_edge(clk) then
      if rst = '1' then
        clip_cnt <= (others => '0');
      elsif d2_val = '1' and d2_clip = '1' then
        clip_cnt <= clip_cnt + 1;
      end if;
    end if;
  end process p_clip;

  -- -------------------------------------------------------------------
  -- Dwell boundary: acquire on the first tick, then apply hysteresis.
  -- All comparisons are exact integer, with the RMS squared away.
  -- -------------------------------------------------------------------
  p_agc : process (clk)
    variable p      : unsigned(CMPW-1 downto 0);
    variable pl     : unsigned(CMPW-1 downto 0);
    variable thr    : unsigned(CMPW-1 downto 0);
    variable thr2   : unsigned(CMPW-1 downto 0);
    variable base   : unsigned(CMPW-1 downto 0);
    variable e      : integer;
    variable sh     : integer;
    variable twok   : natural;
    variable newsh  : integer;
  begin
    if rising_edge(clk) then
      if rst = '1' then
        shift_cur <= to_integer(signed(i_shift_init));
        acquired  <= '0';
      elsif tick = '1' then
        twok := 2*G_KDWELL;
        p    := resize(pacc_hold, CMPW);

        if acquired = '0' then
          -- FAST ACQUISITION: e = round(log2(max(rms,1))), found by
          -- locating p between 2K*2^(2e-1) and 2K*2^(2e+1).  No log, no
          -- sqrt - just a search over at most ~24 exact comparisons.
          e := 0;
          for t in 1 to 30 loop
            base := to_unsigned(twok, CMPW);
            thr  := shift_left(base, 2*t-1);
            if p >= thr then
              e := t;
            end if;
          end loop;
          newsh := e - 8 - (F_BEAM - F_DAC);
          acquired <= '1';
        else
          newsh := shift_cur;
        end if;

        -- Hysteresis, evaluated at the shift that WILL be in force.
        -- if/elsif, NOT two independent tests: the model's branches are
        -- mutually exclusive and although HI > LO makes both firing
        -- impossible today, mirroring the structure keeps it that way if
        -- the band is ever re-tuned.
        sh   := F_BEAM - F_DAC + newsh;
        base := to_unsigned(twok, CMPW);
        pl   := p;
        if sh < 0 then
          pl := shift_left(p, -2*sh);
        end if;

        thr := resize(base * to_unsigned(HI_C, 20), CMPW);
        if sh >= 0 then
          thr := shift_left(thr, 2*sh);
        end if;
        thr2 := resize(base * to_unsigned(LO_C, 20), CMPW);
        if sh >= 0 then
          thr2 := shift_left(thr2, 2*sh);
        end if;

        if pl > thr then
          newsh := newsh + 1;          -- too hot: shift right more
        elsif pl < thr2 then
          newsh := newsh - 1;
        end if;

        shift_cur <= newsh;
      end if;
    end if;
  end process p_agc;

  m_valid <= d2_val;
  m_re    <= std_logic_vector(d2_re);
  m_im    <= std_logic_vector(d2_im);
  o_shift <= std_logic_vector(to_signed(shift_cur, 8));
  o_clip  <= std_logic_vector(clip_cnt);
  o_dwell_tick <= tick;

end architecture rtl;
