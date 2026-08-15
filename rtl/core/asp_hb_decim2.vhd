-- =====================================================================
--  asp_hb_decim2  -  STAGE 2.  Halfband decimate-by-two.
--
--  MODEL REFERENCE : stage2_hbdec() in matlab/golden/asp_golden_model.m
--
--  FUNCTION
--    y[p] = round_even( sum_j h[j] * x[2p - j] >> 17 ),  saturated to s16.
--    11 taps, Q1.17, DC gain exactly 2^17 so the filter is unity gain and
--    no scaling drift accumulates down the chain.
--
--  WHY THIS IS CUSTOM RTL AND NOT FIR COMPILER
--    Every non-zero coefficient is a sum of at most two powers of two:
--        h[2] = h[8] = -4096  = -(2^12)
--        h[4] = h[6] = +36864 = 2^15 + 2^12
--        h[5]        = +65536 = 2^16
--    so the entire filter is shifts and adds and costs ZERO DSP48.  FIR
--    Compiler would spend multipliers on it, and multipliers are the one
--    resource this design is actually tight on once the shaping FIR and
--    the covariance array are placed.  The elaboration assertions below
--    tie that claim to the generated coefficient table: if the model is
--    ever re-tuned so a coefficient stops being a two-term power-of-two
--    sum, this file fails to elaborate instead of silently computing the
--    wrong filter.
--
--  WHY A HALFBAND AND WHY NO CIC
--    Linear phase with IDENTICAL group delay on all four channels by
--    construction, which is what keeps the array aligned.  A CIC would be
--    cheaper still but has passband droop and non-linear phase that would
--    have to be equalised PER CHANNEL - and any per-channel difference in
--    the equaliser is exactly the mismatch the array cannot tolerate.
--
--  DECIMATION PHASE
--    The model keeps the EVEN input indices (accR(1:2:end) in 1-based
--    MATLAB is n = 0, 2, 4, ...).  Getting this off by one shifts every
--    channel by half a sample; because it shifts them all equally it does
--    NOT show up as a broken covariance, only as a slightly wrong filter
--    response - which is why it is worth stating rather than assuming.
--
--  ARCHITECTURE
--    Operates on the dense TDM stream from stage 1: one channel per clk,
--    channels cycling 0..3.  A sample of the same channel is therefore
--    4 slots back, so tap j lives at shift-register position 4j and the
--    four antennas share one delay line and one arithmetic unit.  The
--    output is a TDM stream at half the input rate: valid on four of
--    every eight clocks.
--
--  FIXED POINT
--    input   s16.15
--    coeffs  s18.17 (exact powers of two)
--    accum   s36     full precision, no rounding inside the sum
--    output  s16.15  after >> 17, convergent round, saturate
--
--  LATENCY   5 clk from the input slot to the corresponding output slot.
--  THROUGHPUT one output channel-sample per two clks (65.472 MS/s).
--  RESOURCES 0 DSP48, ~1100 FF (the 33-deep x 32-bit delay line),
--            ~250 LUT, 0 BRAM.
-- =====================================================================
library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

library work;
use work.asp_pkg.all;
use work.asp_coef_pkg.all;

entity asp_hb_decim2 is
  generic (
    G_NCH : natural := 4
  );
  port (
    clk     : in  std_logic;
    rst     : in  std_logic;

    s_valid : in  std_logic;
    s_chan  : in  std_logic_vector(1 downto 0);
    s_re    : in  std_logic_vector(W_MIX-1 downto 0);
    s_im    : in  std_logic_vector(W_MIX-1 downto 0);

    m_valid : out std_logic;
    m_chan  : out std_logic_vector(1 downto 0);
    m_re    : out std_logic_vector(W_HB-1 downto 0);
    m_im    : out std_logic_vector(W_HB-1 downto 0)
  );
end entity asp_hb_decim2;


architecture rtl of asp_hb_decim2 is

  -- Tap j of the filter sits G_NCH*j slots back in the TDM stream.
  constant TAP_A1 : natural := 2 * 4;      -- h[2]
  constant TAP_B1 : natural := 4 * 4;      -- h[4]
  constant TAP_C  : natural := 5 * 4;      -- h[5], centre
  constant TAP_B2 : natural := 6 * 4;      -- h[6] = h[4]
  constant TAP_A2 : natural := 8 * 4;      -- h[8] = h[2]
  constant SR_LEN : natural := TAP_A2 + 1; -- 33

  constant ACC_W  : natural := 36;

  type sr_t is array (0 to SR_LEN-1) of signed(W_MIX-1 downto 0);
  signal sr_re : sr_t := (others => (others => '0'));
  signal sr_im : sr_t := (others => (others => '0'));

  -- Sample parity: toggles once per group of G_NCH slots.  The output is
  -- produced only on even sample index, which is the model's decimation
  -- phase.
  signal parity   : std_logic := '0';
  signal sr_valid : std_logic := '0';
  signal sr_chan  : unsigned(1 downto 0) := (others => '0');
  signal sr_par   : std_logic := '0';

  -- pipeline
  signal p1_a_re, p1_a_im : signed(W_MIX downto 0) := (others => '0');
  signal p1_b_re, p1_b_im : signed(W_MIX downto 0) := (others => '0');
  signal p1_c_re, p1_c_im : signed(W_MIX-1 downto 0) := (others => '0');
  signal p1_val  : std_logic := '0';
  signal p1_chan : unsigned(1 downto 0) := (others => '0');

  signal p2_t1_re, p2_t1_im : signed(ACC_W-1 downto 0) := (others => '0');
  signal p2_t2_re, p2_t2_im : signed(ACC_W-1 downto 0) := (others => '0');
  signal p2_val  : std_logic := '0';
  signal p2_chan : unsigned(1 downto 0) := (others => '0');

  signal p3_re, p3_im : signed(ACC_W-1 downto 0) := (others => '0');
  signal p3_val  : std_logic := '0';
  signal p3_chan : unsigned(1 downto 0) := (others => '0');

  signal p4_re, p4_im : signed(W_HB-1 downto 0) := (others => '0');
  signal p4_val  : std_logic := '0';
  signal p4_chan : unsigned(1 downto 0) := (others => '0');

begin

  -- -------------------------------------------------------------------
  -- The shift-add decomposition above is only valid for these exact
  -- coefficients.  Tie it to the generated table so a model change that
  -- invalidates it stops the build rather than corrupting the filter.
  -- -------------------------------------------------------------------
  assert HB_NTAPS = 11
    report "asp_hb_decim2: halfband length changed; rework the tap map"
    severity failure;
  assert HB_COEF(0) = 0 and HB_COEF(1) = 0 and HB_COEF(3) = 0 and
         HB_COEF(7) = 0 and HB_COEF(9) = 0 and HB_COEF(10) = 0
    report "asp_hb_decim2: a coefficient assumed zero is not zero"
    severity failure;
  assert HB_COEF(2) = -4096 and HB_COEF(8) = -4096
    report "asp_hb_decim2: outer tap is no longer -2^12"
    severity failure;
  assert HB_COEF(4) = 36864 and HB_COEF(6) = 36864
    report "asp_hb_decim2: inner tap is no longer 2^15 + 2^12"
    severity failure;
  assert HB_COEF(5) = 65536
    report "asp_hb_decim2: centre tap is no longer 2^16"
    severity failure;
  assert SH_HB = 17
    report "asp_hb_decim2: output shift changed"
    severity failure;

  -- -------------------------------------------------------------------
  -- Delay line.  One line shared by all four antennas: they traverse the
  -- same registers, so they cannot acquire different group delays.
  -- -------------------------------------------------------------------
  p_sr : process (clk)
  begin
    if rising_edge(clk) then
      if rst = '1' then
        sr_re    <= (others => (others => '0'));
        sr_im    <= (others => (others => '0'));
        parity   <= '0';
        sr_valid <= '0';
        sr_par   <= '0';
      else
        sr_valid <= s_valid;
        if s_valid = '1' then
          sr_re(0) <= signed(s_re);
          sr_im(0) <= signed(s_im);
          for i in 1 to SR_LEN-1 loop
            sr_re(i) <= sr_re(i-1);
            sr_im(i) <= sr_im(i-1);
          end loop;
          sr_chan <= unsigned(s_chan);
          sr_par  <= parity;
          -- advance the sample index after the last channel of a group
          if unsigned(s_chan) = to_unsigned(G_NCH-1, 2) then
            parity <= not parity;
          end if;
        end if;
      end if;
    end if;
  end process p_sr;

  -- -------------------------------------------------------------------
  -- P1 : symmetric pre-adds.  Pairing h[2] with h[8] and h[4] with h[6]
  -- halves the arithmetic and is exact - the model sums the same terms.
  -- -------------------------------------------------------------------
  p_p1 : process (clk)
  begin
    if rising_edge(clk) then
      if rst = '1' then
        -- The decimation decision travels WITH the data rather than
        -- being recomputed downstream, so a stalled or missing input
        -- slot can never slip the output phase by half a sample.
        p1_val <= '0';
      else
        p1_val <= sr_valid and (not sr_par);
      end if;
      p1_a_re <= resize(sr_re(TAP_A1), W_MIX+1) + resize(sr_re(TAP_A2), W_MIX+1);
      p1_a_im <= resize(sr_im(TAP_A1), W_MIX+1) + resize(sr_im(TAP_A2), W_MIX+1);
      p1_b_re <= resize(sr_re(TAP_B1), W_MIX+1) + resize(sr_re(TAP_B2), W_MIX+1);
      p1_b_im <= resize(sr_im(TAP_B1), W_MIX+1) + resize(sr_im(TAP_B2), W_MIX+1);
      p1_c_re <= sr_re(TAP_C);
      p1_c_im <= sr_im(TAP_C);
      p1_chan <= sr_chan;
    end if;
  end process p_p1;

  -- -------------------------------------------------------------------
  -- P2 : the shift-add products.  Written as explicit shifts, not as
  -- multiplications by a constant, so that no DSP48 can be inferred.
  --    t1 = (b << 15) + (b << 12)       =  36864 * b
  --    t2 = (c << 16) - (a << 12)       =  65536 * c - 4096 * a
  -- -------------------------------------------------------------------
  p_p2 : process (clk)
  begin
    if rising_edge(clk) then
      if rst = '1' then
        p2_val <= '0';
      else
        p2_val <= p1_val;
      end if;
      p2_t1_re <= shift_left(resize(p1_b_re, ACC_W), 15) +
                  shift_left(resize(p1_b_re, ACC_W), 12);
      p2_t1_im <= shift_left(resize(p1_b_im, ACC_W), 15) +
                  shift_left(resize(p1_b_im, ACC_W), 12);
      p2_t2_re <= shift_left(resize(p1_c_re, ACC_W), 16) -
                  shift_left(resize(p1_a_re, ACC_W), 12);
      p2_t2_im <= shift_left(resize(p1_c_im, ACC_W), 16) -
                  shift_left(resize(p1_a_im, ACC_W), 12);
      p2_chan  <= p1_chan;
    end if;
  end process p_p2;

  -- -------------------------------------------------------------------
  -- P3 : final accumulate, still full precision.
  -- -------------------------------------------------------------------
  p_p3 : process (clk)
  begin
    if rising_edge(clk) then
      if rst = '1' then
        p3_val <= '0';
      else
        p3_val <= p2_val;
      end if;
      p3_re   <= p2_t1_re + p2_t2_re;
      p3_im   <= p2_t1_im + p2_t2_im;
      p3_chan <= p2_chan;
    end if;
  end process p_p3;

  -- -------------------------------------------------------------------
  -- P4 : the single rounding point of this block.
  -- -------------------------------------------------------------------
  p_p4 : process (clk)
  begin
    if rising_edge(clk) then
      if rst = '1' then
        p4_val <= '0';
      else
        p4_val <= p3_val;
      end if;
      p4_re   <= shift_round_sat(p3_re, SH_HB, W_HB);
      p4_im   <= shift_round_sat(p3_im, SH_HB, W_HB);
      p4_chan <= p3_chan;
    end if;
  end process p_p4;

  m_valid <= p4_val;
  m_chan  <= std_logic_vector(p4_chan);
  m_re    <= std_logic_vector(p4_re);
  m_im    <= std_logic_vector(p4_im);

end architecture rtl;
