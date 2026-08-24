-- =====================================================================
--  asp_datapath  -  the complete PL signal chain, stages 1 to 10.
--
--  Wires the ten co-simulated blocks into one streaming pipeline and
--  owns the two pieces of sequencing that do not belong inside any of
--  them: the dwell handshake between the estimator blocks, and the
--  causal weight update.
--
--  ------------------------------------------------------------------
--  DATA FLOW
--  ------------------------------------------------------------------
--    4 x s12.11 @ 32.736 MHz   (from the AD9361 pair)
--        |
--    [1] asp_ddc_mixer      -> TDM, dense, 130.944 MS/s
--    [2] asp_hb_decim2      -> TDM, half rate
--    [3] asp_fir_shape      -> TDM, s16.15, 65.472 MS/s
--        |                                    |
--        |                                    +--> [9] asp_beamformer
--    [4] asp_cov_accum  (1 kHz)                        |
--    [5] asp_whiten                              [10] asp_tx_scale
--    [6] asp_jacobi_evd                                |
--    [7] asp_detect                              s12.11 @ 16.368 MHz
--    [8] asp_weight_calc  --(weights)------------------+
--
--  ------------------------------------------------------------------
--  WEIGHT LATENCY: TWO DWELLS, NOT ONE.  READ THIS BEFORE COMPARING
--  AGAINST THE MODEL.
--  ------------------------------------------------------------------
--  The model applies dwell k's estimate to dwell k+1.  Hardware cannot:
--  the covariance for dwell k is only complete AT THE END of dwell k,
--  and the estimator chain then needs
--
--      whiten ~270 + EVD ~4700 + detect 2 + weights ~350  =  ~5400 clk
--
--  which is about 4% of a dwell.  The new weights are therefore ready
--  partway INTO dwell k+1, not at its start.
--
--  Switching weights mid-dwell would put a discontinuity inside a
--  coherent integration period, so this design does not: the weight
--  register is loaded only at a dwell boundary, and dwell k's estimate
--  is applied from the start of dwell k+2.  The cost is that the weights
--  are one dwell (1 ms) staler than the model's.
--
--  That is a real difference and it is worth being explicit about its
--  size: the spatial signature of a spoofer moves on the timescale of
--  platform and satellite motion - seconds - so 1 ms of extra staleness
--  changes the null depth by an amount far below the dwell-to-dwell
--  estimation noise.  What it is NOT is a free choice; it is what the
--  estimator's own latency costs.
--
--  ------------------------------------------------------------------
--  BACKPRESSURE
--  ------------------------------------------------------------------
--  There is none, deliberately.  Every rate in the chain is an exact
--  integer ratio of the ADC clock, so the pipeline is rate-matched by
--  construction and a tready would never deassert.  The one place
--  elasticity IS needed - the burst out of the decimator - has a FIFO
--  inside asp_fir_shape with an overflow flag wired to a status
--  register, so a rate-plan error shows up as a sticky bit rather than
--  as silently corrupted data.
--
--  ------------------------------------------------------------------
--  RESOURCES (estimated, XC7Z020)
--  ------------------------------------------------------------------
--    DSP48E1   ~85 of 220
--    BRAM36      0 of 140  (the delay lines are registers/SRL)
--    FF        ~19000 of 106400
--    LUT       ~13000 of 53200
-- =====================================================================
library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

library work;
use work.asp_pkg.all;
use work.asp_coef_pkg.all;

entity asp_datapath is
  generic (
    G_NCH    : natural := 4;
    G_KDWELL : natural := K_DWELL;
    G_DSW    : natural := 24
  );
  port (
    clk : in std_logic;
    rst : in std_logic;

    -- ---- RX, parallel, one strobe per ADC sample --------------------
    s_valid : in  std_logic;
    s_re    : in  std_logic_vector(G_NCH*W_ADC-1 downto 0);
    s_im    : in  std_logic_vector(G_NCH*W_ADC-1 downto 0);

    -- ---- TX, one complex sample per FS_WORK period ------------------
    m_valid : out std_logic;
    m_re    : out std_logic_vector(W_DAC-1 downto 0);
    m_im    : out std_logic_vector(W_DAC-1 downto 0);

    -- ---- control (AXI4-Lite register file) --------------------------
    i_bypass     : in  std_logic;                     -- force quiescent beam
    i_rank2_en   : in  std_logic;
    i_shift_init : in  std_logic_vector(7 downto 0);

    -- ---- status / telemetry, latched per dwell ----------------------
    o_dwell_tick : out std_logic;
    o_lam        : out std_logic_vector(G_NCH*W_EVD-1 downto 0);
    o_det        : out std_logic;
    o_rank       : out std_logic_vector(1 downto 0);
    o_rank2_el   : out std_logic;
    o_lhs        : out std_logic_vector(47 downto 0);
    o_rhs        : out std_logic_vector(47 downto 0);
    o_w          : out std_logic_vector(2*G_NCH*W_WGT-1 downto 0);
    o_tx_shift   : out std_logic_vector(7 downto 0);
    o_clip       : out std_logic_vector(31 downto 0);
    o_fir_ovf    : out std_logic;
    o_dwell_cnt  : out std_logic_vector(31 downto 0)
  );
end entity asp_datapath;


architecture rtl of asp_datapath is

  -- stage 1 -> 2
  signal x1_v : std_logic;
  signal x1_c : std_logic_vector(1 downto 0);
  signal x1_r, x1_i : std_logic_vector(W_MIX-1 downto 0);

  -- stage 2 -> 3
  signal x2_v : std_logic;
  signal x2_c : std_logic_vector(1 downto 0);
  signal x2_r, x2_i : std_logic_vector(W_HB-1 downto 0);

  -- stage 3 -> 4 and 9
  signal x3_v : std_logic;
  signal x3_c : std_logic_vector(1 downto 0);
  signal x3_r, x3_i : std_logic_vector(W_DAT-1 downto 0);

  -- stage 4 -> 5
  signal c4_v, c4_l, c4_tick : std_logic;
  signal c4_x : std_logic_vector(3 downto 0);
  signal c4_r, c4_i : std_logic_vector(W_ACC-1 downto 0);

  -- stage 5 -> 6
  signal w5_v, w5_done : std_logic;
  signal w5_x : std_logic_vector(3 downto 0);
  signal w5_r, w5_i : std_logic_vector(W_EVD-1 downto 0);
  signal w5_dsq : std_logic_vector(G_NCH*G_DSW-1 downto 0);

  -- stage 6 -> 7, 8
  signal e6_v, e6_done : std_logic;
  signal e6_x : std_logic_vector(3 downto 0);
  signal e6_r, e6_i : std_logic_vector(W_UVEC-1 downto 0);
  signal e6_lam : std_logic_vector(G_NCH*W_EVD-1 downto 0);

  -- stage 7
  signal d7_v, d7_det, d7_el2 : std_logic;
  signal d7_rank : std_logic_vector(1 downto 0);
  signal d7_lhs, d7_rhs : std_logic_vector(47 downto 0);

  -- stage 8
  signal w8_done : std_logic;
  signal w8_w    : std_logic_vector(2*G_NCH*W_WGT-1 downto 0);

  -- applied weights
  signal w_applied : std_logic_vector(2*G_NCH*W_WGT-1 downto 0);
  signal w_pending : std_logic_vector(2*G_NCH*W_WGT-1 downto 0);
  signal w_have    : std_logic := '0';

  -- stage 9 -> 10
  signal b9_v : std_logic;
  signal b9_r, b9_i : std_logic_vector(W_BEAM-1 downto 0);

  signal dwell_cnt : unsigned(31 downto 0) := (others => '0');
  signal quiescent : std_logic_vector(2*G_NCH*W_WGT-1 downto 0);

begin

  -- quiescent beam: h = 2^(F_WGT-1) on every element, imaginary part zero
  g_quies : for k in 0 to G_NCH-1 generate
    quiescent((2*k+1)*W_WGT-1 downto (2*k)*W_WGT) <=
      std_logic_vector(to_signed(WGT_QUIESCENT, W_WGT));
    quiescent((2*k+2)*W_WGT-1 downto (2*k+1)*W_WGT) <=
      (others => '0');
  end generate;

  -- ------------------------------------------------------ stage 1
  u_ddc : entity work.asp_ddc_mixer
    generic map (G_NCH => G_NCH)
    port map (clk => clk, rst => rst,
              s_valid => s_valid, s_re => s_re, s_im => s_im,
              m_valid => x1_v, m_chan => x1_c, m_re => x1_r, m_im => x1_i);

  -- ------------------------------------------------------ stage 2
  u_hb : entity work.asp_hb_decim2
    generic map (G_NCH => G_NCH)
    port map (clk => clk, rst => rst,
              s_valid => x1_v, s_chan => x1_c, s_re => x1_r, s_im => x1_i,
              m_valid => x2_v, m_chan => x2_c, m_re => x2_r, m_im => x2_i);

  -- ------------------------------------------------------ stage 3
  u_fir : entity work.asp_fir_shape
    generic map (G_NCH => G_NCH)
    port map (clk => clk, rst => rst,
              s_valid => x2_v, s_chan => x2_c, s_re => x2_r, s_im => x2_i,
              m_valid => x3_v, m_chan => x3_c, m_re => x3_r, m_im => x3_i,
              o_overflow => o_fir_ovf);

  -- ------------------------------------------------------ stage 4
  u_cov : entity work.asp_cov_accum
    generic map (G_NCH => G_NCH, G_KDWELL => G_KDWELL)
    port map (clk => clk, rst => rst,
              s_valid => x3_v, s_chan => x3_c, s_re => x3_r, s_im => x3_i,
              o_valid => c4_v, o_idx => c4_x, o_re => c4_r, o_im => c4_i,
              o_last => c4_l, o_dwell_tick => c4_tick);

  -- ------------------------------------------------------ stage 5
  u_whiten : entity work.asp_whiten
    generic map (G_N => G_NCH, G_DSW => G_DSW)
    port map (clk => clk, rst => rst,
              s_valid => c4_v, s_idx => c4_x, s_re => c4_r, s_im => c4_i,
              s_last => c4_l,
              o_valid => w5_v, o_idx => w5_x, o_re => w5_r, o_im => w5_i,
              o_done => w5_done, o_dsq => w5_dsq);

  -- ------------------------------------------------------ stage 6
  u_evd : entity work.asp_jacobi_evd
    generic map (G_N => G_NCH)
    port map (clk => clk, rst => rst,
              s_valid => w5_v, s_idx => w5_x, s_re => w5_r, s_im => w5_i,
              s_start => w5_done,
              o_valid => e6_v, o_idx => e6_x, o_re => e6_r, o_im => e6_i,
              o_lam => e6_lam, o_done => e6_done);

  -- ------------------------------------------------------ stage 7
  u_det : entity work.asp_detect
    generic map (G_N => G_NCH)
    port map (clk => clk, rst => rst,
              i_start => e6_done, i_lam => e6_lam, i_rank2_en => i_rank2_en,
              o_valid => d7_v, o_lhs => d7_lhs, o_rhs => d7_rhs,
              o_det => d7_det, o_rank => d7_rank, o_rank2_eligible => d7_el2);

  -- ------------------------------------------------------ stage 8
  -- The eigenvector matrix is captured inside asp_weight_calc as it
  -- streams past, so stage 8 starts as soon as the detector has produced
  -- a rank - no separate buffer between them.
  u_wgt : entity work.asp_weight_calc
    generic map (G_N => G_NCH, G_DSW => G_DSW)
    port map (clk => clk, rst => rst,
              s_valid => e6_v, s_idx => e6_x, s_re => e6_r, s_im => e6_i,
              i_start => d7_v, i_rank => d7_rank, i_dsq => w5_dsq,
              o_w => w8_w, o_done => w8_done);

  -- ------------------------------------------------------ weight update
  -- The ONLY place causality is enforced.  New weights are parked in
  -- w_pending when the estimator finishes and are promoted to
  -- w_applied only on a dwell boundary, so a weight vector is never
  -- changed inside a coherent integration period.
  p_weights : process (clk)
  begin
    if rising_edge(clk) then
      if rst = '1' then
        w_applied <= quiescent;
        w_pending <= quiescent;
        w_have    <= '0';
        dwell_cnt <= (others => '0');
      else
        if w8_done = '1' then
          w_pending <= w8_w;
          w_have    <= '1';
        end if;
        if c4_tick = '1' then
          dwell_cnt <= dwell_cnt + 1;
          if i_bypass = '1' then
            w_applied <= quiescent;
          elsif w_have = '1' then
            w_applied <= w_pending;
            w_have    <= '0';
          end if;
        end if;
      end if;
    end if;
  end process p_weights;

  -- ------------------------------------------------------ stage 9
  u_beam : entity work.asp_beamformer
    generic map (G_NCH => G_NCH)
    port map (clk => clk, rst => rst, i_w => w_applied,
              s_valid => x3_v, s_chan => x3_c, s_re => x3_r, s_im => x3_i,
              m_valid => b9_v, m_re => b9_r, m_im => b9_i);

  -- ------------------------------------------------------ stage 10
  u_tx : entity work.asp_tx_scale
    generic map (G_KDWELL => G_KDWELL)
    port map (clk => clk, rst => rst,
              i_shift_init => i_shift_init,
              s_valid => b9_v, s_re => b9_r, s_im => b9_i,
              m_valid => m_valid, m_re => m_re, m_im => m_im,
              o_shift => o_tx_shift, o_clip => o_clip, o_dwell_tick => open);

  -- ------------------------------------------------------ telemetry
  p_status : process (clk)
  begin
    if rising_edge(clk) then
      if rst = '1' then
        o_det      <= '0';
        o_rank     <= (others => '0');
        o_rank2_el <= '0';
      elsif d7_v = '1' then
        o_det      <= d7_det;
        o_rank     <= d7_rank;
        o_rank2_el <= d7_el2;
        o_lhs      <= d7_lhs;
        o_rhs      <= d7_rhs;
        o_lam      <= e6_lam;
      end if;
    end if;
  end process p_status;

  o_w          <= w_applied;
  o_dwell_tick <= c4_tick;
  o_dwell_cnt  <= std_logic_vector(dwell_cnt);

end architecture rtl;
