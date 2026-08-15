-- =====================================================================
--  asp_ddc_mixer  -  STAGE 1.  Digital down-converter.
--
--  MODEL REFERENCE : stage1_ddc() in matlab/golden/asp_golden_model.m
--
--  FUNCTION
--    y[n] = x[n] * exp(+j*2*pi*n/16), one round and saturate.
--    Mixes the L1 signal from -2.046 MHz (where the offset analogue LO
--    puts it) up to baseband.
--
--  WHY THERE IS NO PHASE ACCUMULATOR
--    The LO offset is FS_ADC/16 EXACTLY, so the NCO is a 16-entry ROM
--    addressed by a 4-bit counter.  No accumulator means no phase
--    truncation, which means no truncation spurs - not "spurs below some
--    level", none at all.  A DDS Compiler instance would be larger AND
--    worse here, which is the whole reason this block is custom RTL.
--
--  WHY ONE NCO FOR ALL FOUR CHANNELS
--    A common complex rotation applied to every element is invisible to
--    the array algorithm: the covariance, the projector and the
--    beamformer are all invariant to a common complex scalar.  A shared
--    NCO therefore contributes EXACTLY ZERO inter-channel mismatch.
--    Four independent NCOs would not, and the mismatch would appear to
--    the eigen-detector as a spurious source.
--
--  ARCHITECTURE
--    Time-division multiplexed over the four antennas.  clk runs at
--    4 x FS_ADC, so each ADC sample period contains exactly four slots
--    and one slot carries one antenna.  That folds four complex
--    multipliers into one:  4 DSP48E1 instead of 16, with no loss of
--    throughput and no FIFO, because the ratio is an exact integer.
--
--    The output is a dense TDM stream - m_valid is high on every clock -
--    and every downstream block consumes that same format.  Choosing one
--    stream format for the whole datapath is what keeps the channel
--    alignment provably identical: the four antennas traverse the SAME
--    logic, not four copies of it, so there is no path for them to
--    acquire different group delays.
--
--  FIXED POINT
--    input   s12.11  (AD9361 RX, raw two's complement)
--    NCO     s16.14  (Q1.14 so +1.0 is representable; |.| = 16384)
--    product s28     held at full precision
--    output  s16.15  after >> (F_ADC + F_NCO - F_MIX) = >> 10,
--                    convergent round, saturate
--
--  LATENCY   4 clk after the slot in which a channel is selected.
--            Deterministic; no data-dependent paths.
--  THROUGHPUT one channel-sample per clk (130.944 MS/s aggregate).
--  RESOURCES 4 DSP48E1, ~120 FF, ~90 LUT, 0 BRAM.
--
--  TIMING
--    The multiplier operands are registered (stage A) and the products
--    are registered (stage B) so Vivado can use the DSP48E1 internal
--    A/B and M registers.  The critical path is therefore inside the
--    DSP slice, which closes above 400 MHz on -1 silicon; at 130.944 MHz
--    this block has ~3x margin and is never the limiting path.
-- =====================================================================
library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

library work;
use work.asp_pkg.all;
use work.asp_coef_pkg.all;

entity asp_ddc_mixer is
  generic (
    G_NCH : natural := 4                       -- antennas; TDM slots per sample
  );
  port (
    clk      : in  std_logic;
    rst      : in  std_logic;                  -- synchronous, active high

    -- Parallel input, one strobe per ADC sample (FS_ADC = clk/G_NCH).
    -- Channel k occupies bits ((k+1)*W_ADC-1 downto k*W_ADC).
    s_valid  : in  std_logic;
    s_re     : in  std_logic_vector(G_NCH*W_ADC-1 downto 0);
    s_im     : in  std_logic_vector(G_NCH*W_ADC-1 downto 0);

    -- Dense TDM output, one channel per clk.
    m_valid  : out std_logic;
    m_chan   : out std_logic_vector(1 downto 0);
    m_re     : out std_logic_vector(W_MIX-1 downto 0);
    m_im     : out std_logic_vector(W_MIX-1 downto 0)
  );
end entity asp_ddc_mixer;


architecture rtl of asp_ddc_mixer is

  -- ---- input holding registers, loaded once per ADC sample ----------
  type adc_arr_t is array (0 to G_NCH-1) of signed(W_ADC-1 downto 0);
  signal hold_re : adc_arr_t := (others => (others => '0'));
  signal hold_im : adc_arr_t := (others => (others => '0'));

  -- ---- NCO -----------------------------------------------------------
  -- 4-bit phase counter, advanced once per ADC sample.  The SAME phase
  -- is applied to all G_NCH channels of that sample, which is what makes
  -- the rotation common mode.
  signal nco_phase : unsigned(3 downto 0) := (others => '0');
  signal lo_re_r   : signed(W_NCO-1 downto 0) := (others => '0');
  signal lo_im_r   : signed(W_NCO-1 downto 0) := (others => '0');

  -- ---- TDM sequencing -------------------------------------------------
  signal slot      : unsigned(1 downto 0) := (others => '0');
  signal run       : std_logic := '0';       -- set by the first s_valid

  -- ---- pipeline -------------------------------------------------------
  -- A : operands selected and registered      (DSP48 A/B registers)
  -- B : four products registered              (DSP48 M register)
  -- C : two sums registered                   (DSP48 P register)
  -- D : shift / convergent round / saturate
  constant PROD_W : natural := W_ADC + W_NCO;             -- 28
  constant SUM_W  : natural := PROD_W + 1;                -- 29

  signal a_re, a_im : signed(W_ADC-1 downto 0) := (others => '0');
  signal a_lr, a_li : signed(W_NCO-1 downto 0) := (others => '0');
  signal a_chan     : unsigned(1 downto 0) := (others => '0');
  signal a_val      : std_logic := '0';

  signal b_p1, b_p2, b_p3, b_p4 : signed(PROD_W-1 downto 0) := (others => '0');
  signal b_chan     : unsigned(1 downto 0) := (others => '0');
  signal b_val      : std_logic := '0';

  signal c_pr, c_pi : signed(SUM_W-1 downto 0) := (others => '0');
  signal c_chan     : unsigned(1 downto 0) := (others => '0');
  signal c_val      : std_logic := '0';

  signal d_re, d_im : signed(W_MIX-1 downto 0) := (others => '0');
  signal d_chan     : unsigned(1 downto 0) := (others => '0');
  signal d_val      : std_logic := '0';

begin

  -- -------------------------------------------------------------------
  -- Input capture and TDM slot sequencing.
  --
  -- s_valid arrives once every G_NCH clocks by construction (the ADC
  -- clock and clk come from the same MMCM with an exact integer ratio).
  -- The slot counter free-runs once started, so a missing strobe cannot
  -- desynchronise the channel order; it would only repeat the previous
  -- sample, which is a benign failure compared with silently permuting
  -- the antennas.
  -- -------------------------------------------------------------------
  p_input : process (clk)
  begin
    if rising_edge(clk) then
      if rst = '1' then
        slot      <= (others => '0');
        nco_phase <= (others => '0');
        run       <= '0';
      else
        if s_valid = '1' then
          for k in 0 to G_NCH-1 loop
            hold_re(k) <= signed(s_re((k+1)*W_ADC-1 downto k*W_ADC));
            hold_im(k) <= signed(s_im((k+1)*W_ADC-1 downto k*W_ADC));
          end loop;
          -- Latch the twiddle for this sample, then advance the phase.
          lo_re_r   <= to_signed(NCO_RE(to_integer(nco_phase)), W_NCO);
          lo_im_r   <= to_signed(NCO_IM(to_integer(nco_phase)), W_NCO);
          nco_phase <= nco_phase + 1;
          slot      <= (others => '0');
          run       <= '1';
        elsif run = '1' then
          slot <= slot + 1;
        end if;
      end if;
    end if;
  end process p_input;

  -- -------------------------------------------------------------------
  -- Stage A : select the slot's channel and register the operands.
  -- -------------------------------------------------------------------
  p_stage_a : process (clk)
  begin
    if rising_edge(clk) then
      if rst = '1' then
        a_val <= '0';
      else
        a_val <= run;
      end if;
      a_re   <= hold_re(to_integer(slot));
      a_im   <= hold_im(to_integer(slot));
      a_lr   <= lo_re_r;
      a_li   <= lo_im_r;
      a_chan <= slot;
    end if;
  end process p_stage_a;

  -- -------------------------------------------------------------------
  -- Stage B : the four real products, full precision, no rounding.
  -- Maps to 4 DSP48E1; the M register is inferred from these registers.
  -- -------------------------------------------------------------------
  p_stage_b : process (clk)
  begin
    if rising_edge(clk) then
      if rst = '1' then
        b_val <= '0';
      else
        b_val <= a_val;
      end if;
      b_p1   <= a_re * a_lr;
      b_p2   <= a_im * a_li;
      b_p3   <= a_re * a_li;
      b_p4   <= a_im * a_lr;
      b_chan <= a_chan;
    end if;
  end process p_stage_b;

  -- -------------------------------------------------------------------
  -- Stage C : accumulate the complex product.  Full precision still -
  -- the model sums the products BEFORE rounding, so rounding here would
  -- change the result by up to 1 LSB and break the bit-exact contract.
  -- -------------------------------------------------------------------
  p_stage_c : process (clk)
  begin
    if rising_edge(clk) then
      if rst = '1' then
        c_val <= '0';
      else
        c_val <= b_val;
      end if;
      c_pr   <= resize(b_p1, SUM_W) - resize(b_p2, SUM_W);
      c_pi   <= resize(b_p3, SUM_W) + resize(b_p4, SUM_W);
      c_chan <= b_chan;
    end if;
  end process p_stage_c;

  -- -------------------------------------------------------------------
  -- Stage D : the single rounding point of this block.
  -- -------------------------------------------------------------------
  p_stage_d : process (clk)
  begin
    if rising_edge(clk) then
      if rst = '1' then
        d_val <= '0';
      else
        d_val <= c_val;
      end if;
      d_re   <= shift_round_sat(c_pr, SH_MIX, W_MIX);
      d_im   <= shift_round_sat(c_pi, SH_MIX, W_MIX);
      d_chan <= c_chan;
    end if;
  end process p_stage_d;

  m_valid <= d_val;
  m_chan  <= std_logic_vector(d_chan);
  m_re    <= std_logic_vector(d_re);
  m_im    <= std_logic_vector(d_im);

end architecture rtl;
