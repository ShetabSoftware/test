-- =====================================================================
--  asp_beamformer  -  STAGE 9.  v[n] = sum_i conj(w_i) * x_i[n].
--
--  MODEL REFERENCE : stage9_beamform() in matlab/golden/asp_golden_model.m
--
--  THE CONJUGATE IS ON THE WEIGHTS
--    Getting that backwards conjugates the whole spatial response and
--    steers the null to the MIRROR direction.  It passes a broadside
--    test - where the mirror and the original coincide - and fails
--    everywhere else, which makes it one of the most expensive one-
--    character mistakes available in this design.  Expanded:
--        acc_re = sum ( re(w)*re(x) + im(w)*im(x) )
--        acc_im = sum ( re(w)*im(x) - im(w)*re(x) )
--
--  FULL PRECISION ACROSS THE FOUR TAPS, ONE ROUNDING AT THE OUTPUT
--    The four antenna contributions are accumulated at full width and
--    rounded once, so the only quantisation reaching the DAC is a single
--    rounding rather than four.  Rounding per tap would also correlate
--    the error with the weights, which is exactly the direction the null
--    is in.
--
--  CAUSALITY
--    This block simply applies whatever weights are presented on i_w.
--    The rule that dwell k's estimate is applied to dwell k+1 is
--    enforced where it belongs - in asp_dwell_engine, which updates the
--    weight register only on a dwell boundary.  Putting the causality in
--    the datapath would make it invisible; putting it in the sequencer
--    makes it a single, reviewable register update.
--
--  ARCHITECTURE
--    Consumes the same TDM stream every other block uses: one antenna
--    per valid slot, channels cycling 0..3.  Four real multipliers cover
--    one antenna per slot, and the accumulator is cleared on channel 0
--    and emitted on channel 3, so no separate framing signal is needed.
--
--  FIXED POINT
--    data   s16.15
--    weight s18.16
--    product s34, accumulator s40 (8 terms, no rounding)
--    output s16.15 after >> (F_DAT + F_WGT - F_BEAM) = >> 16
--
--  LATENCY   4 clk from the channel-3 slot to o_valid.
--  THROUGHPUT one output sample per 8 clk (16.368 MS/s).
--  RESOURCES 4 DSP48E1, ~250 FF, ~200 LUT.
-- =====================================================================
library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

library work;
use work.asp_pkg.all;
use work.asp_coef_pkg.all;

entity asp_beamformer is
  generic (
    G_NCH : natural := 4
  );
  port (
    clk     : in  std_logic;
    rst     : in  std_logic;

    -- beamformer weights, flattened: {im,re} pairs, element k at
    -- bits ((2k+2)*W_WGT-1 downto 2k*W_WGT) = re then im
    i_w     : in  std_logic_vector(2*G_NCH*W_WGT-1 downto 0);

    s_valid : in  std_logic;
    s_chan  : in  std_logic_vector(1 downto 0);
    s_re    : in  std_logic_vector(W_DAT-1 downto 0);
    s_im    : in  std_logic_vector(W_DAT-1 downto 0);

    m_valid : out std_logic;
    m_re    : out std_logic_vector(W_BEAM-1 downto 0);
    m_im    : out std_logic_vector(W_BEAM-1 downto 0)
  );
end entity asp_beamformer;


architecture rtl of asp_beamformer is

  constant PW  : natural := W_DAT + W_WGT;      -- 34
  constant ACW : natural := 40;

  signal wr, wi : signed(W_WGT-1 downto 0) := (others => '0');
  signal xr, xi : signed(W_DAT-1 downto 0) := (others => '0');
  signal s1_val : std_logic := '0';
  signal s1_first, s1_last : std_logic := '0';

  signal p0, p1, p2, p3 : signed(PW-1 downto 0) := (others => '0');
  signal s2_val : std_logic := '0';
  signal s2_first, s2_last : std_logic := '0';

  signal acc_re, acc_im : signed(ACW-1 downto 0) := (others => '0');
  signal s3_val, s3_last : std_logic := '0';

  signal out_v : std_logic := '0';
  signal out_re, out_im : signed(W_BEAM-1 downto 0) := (others => '0');

begin

  assert SH_BEAM = F_DAT + F_WGT - F_BEAM
    report "asp_beamformer: output shift inconsistent" severity failure;

  -- Stage 1: select this slot's weight and register the operands.
  p_s1 : process (clk)
    variable k : integer range 0 to G_NCH-1;
  begin
    if rising_edge(clk) then
      if rst = '1' then
        s1_val <= '0';
      else
        s1_val <= s_valid;
      end if;
      k  := to_integer(unsigned(s_chan));
      wr <= signed(i_w((2*k+1)*W_WGT-1 downto (2*k)*W_WGT));
      wi <= signed(i_w((2*k+2)*W_WGT-1 downto (2*k+1)*W_WGT));
      xr <= signed(s_re);
      xi <= signed(s_im);
      if k = 0 then
        s1_first <= '1';
      else
        s1_first <= '0';
      end if;
      if k = G_NCH-1 then
        s1_last <= '1';
      else
        s1_last <= '0';
      end if;
    end if;
  end process p_s1;

  -- Stage 2: the four real products, full precision.
  p_s2 : process (clk)
  begin
    if rising_edge(clk) then
      if rst = '1' then
        s2_val <= '0';
      else
        s2_val <= s1_val;
      end if;
      p0 <= wr * xr;
      p1 <= wi * xi;
      p2 <= wr * xi;
      p3 <= wi * xr;
      s2_first <= s1_first;
      s2_last  <= s1_last;
    end if;
  end process p_s2;

  -- Stage 3: accumulate across the four antennas.  Cleared on channel 0
  -- by LOADING rather than adding, so there is no clear pass and no dead
  -- slot between output samples.
  p_s3 : process (clk)
  begin
    if rising_edge(clk) then
      if rst = '1' then
        s3_val  <= '0';
        s3_last <= '0';
        acc_re  <= (others => '0');
        acc_im  <= (others => '0');
      else
        s3_val  <= s2_val;
        s3_last <= s2_last;
        if s2_val = '1' then
          if s2_first = '1' then
            acc_re <= resize(p0, ACW) + resize(p1, ACW);
            acc_im <= resize(p2, ACW) - resize(p3, ACW);
          else
            acc_re <= acc_re + resize(p0, ACW) + resize(p1, ACW);
            acc_im <= acc_im + resize(p2, ACW) - resize(p3, ACW);
          end if;
        end if;
      end if;
    end if;
  end process p_s3;

  -- Stage 4: the single rounding point.
  p_s4 : process (clk)
  begin
    if rising_edge(clk) then
      if rst = '1' then
        out_v <= '0';
      else
        out_v <= s3_val and s3_last;
      end if;
      out_re <= shift_round_sat(acc_re, SH_BEAM, W_BEAM);
      out_im <= shift_round_sat(acc_im, SH_BEAM, W_BEAM);
    end if;
  end process p_s4;

  m_valid <= out_v;
  m_re    <= std_logic_vector(out_re);
  m_im    <= std_logic_vector(out_im);

end architecture rtl;
