-- =====================================================================
--  asp_top  -  PL top level for the GNSS anti-spoofing array processor.
--
--  TARGET: Xilinx Zynq-7000 XC7Z020 (-1 speed grade assumed).
--
--  ------------------------------------------------------------------
--  CLOCKING - THE DECISION EVERYTHING ELSE FOLLOWS FROM
--  ------------------------------------------------------------------
--  clk_dsp = 130.944 MHz = 4 x FS_ADC = 8 x FS_WORK, and it MUST be
--  derived from the AD9361 sample clock, not from the PS PLL.
--
--  Every rate in the datapath is an exact integer ratio of the ADC
--  clock, so every rate change in the design is a clock enable rather
--  than a FIFO, and there is NO true clock-domain crossing anywhere
--  between the ADC pins and the DAC pins.  That removes an entire class
--  of metastability, reconvergence and gray-code bugs rather than
--  managing them.  The price is that the PS AXI clock is asynchronous to
--  the datapath, which is dealt with ONCE by an AXI Clock Converter -
--  a reviewed IP instance instead of a synchroniser per register.
--
--    -- REQUIRED VIVADO IP CORE:
--    -- IP NAME: Clocking Wizard (clk_wiz)
--    -- CONFIGURATION:
--    --   PRIMITIVE=MMCM
--    --   PRIM_IN_FREQ=245.760          (AD9361 DATA_CLK, LVDS DDR mode)
--    --   CLKOUT1_REQUESTED_OUT_FREQ=130.944
--    --   CLKOUT2_REQUESTED_OUT_FREQ=32.736     (ADC-rate enable domain)
--    --   USE_RESET=true, RESET_TYPE=ACTIVE_LOW
--    --   USE_LOCKED=true
--    --   Jitter filter: minimise output jitter
--    -- CONNECTIONS:
--    --   clk_in1  <- AD9361 DATA_CLK through an IBUFDS + BUFG
--    --   clk_out1 -> clk_dsp (this module)
--    --   locked   -> part of the reset tree, see below
--    -- NOTE: 245.760 / 130.944 is not an integer ratio; the MMCM
--    -- realises it as M/D.  Verify the exact M/D the wizard picks
--    -- reports zero frequency error, and if the AD9361 is clocked from
--    -- a 40 MHz reference derive 130.944 = 40 * 3.2736 instead.  The
--    -- one thing that must NOT happen is clk_dsp being asynchronous to
--    -- the ADC data, because the whole enable-based rate plan assumes
--    -- it is not.
--
--  ------------------------------------------------------------------
--  RESET STRATEGY
--  ------------------------------------------------------------------
--  ONE synchronous, active-high reset for the whole PL datapath,
--  asserted asynchronously and released synchronously (a standard reset
--  synchroniser, below).  Synchronous because Xilinx SRL primitives have
--  no reset input and an asynchronous reset would force the FIR and
--  halfband delay lines into 8000+ flip-flops instead.  Active high
--  because that is what the fabric registers want natively.
--
--  Datapath pipeline registers that self-flush are NOT reset - only
--  control state and valid pipelines are.  That keeps the reset net
--  small, which matters on a device this size: a global high-fanout
--  reset is a routing and timing problem, not a safety feature.
--
--  ------------------------------------------------------------------
--  INTERFACES - NOTHING IS ASSUMED TO EXIST
--  ------------------------------------------------------------------
--    s_axi_*        AXI4-Lite slave, 4 kB, control and telemetry.
--                   Reached from the PS M_AXI_GP0 through an AXI
--                   Interconnect and an AXI Clock Converter.
--    rx_*           Parallel 4-antenna ADC input at FS_ADC.  Supplied by
--                   two AD9361 interfaces; see the IP note below.
--    tx_*           Parallel complex DAC output at FS_WORK.
--    m_axis_*       AXI4-Stream master carrying the beamformed samples,
--                   for optional capture to DDR through an AXI DMA.
--                   NOT required for the anti-spoofing function itself.
--    irq            Level interrupt to the PS, one per dwell (1 kHz).
--
--    -- REQUIRED VIVADO IP CORE:
--    -- IP NAME: AXI AD9361 (Analog Devices axi_ad9361), TWO instances
--    -- CONFIGURATION:
--    --   MODE_1R1T=0 (2R2T), so each instance supplies TWO antennas
--    --   DAC_DDS_DISABLE=1, ADC_DATAPATH_DISABLE=0
--    --   DELAY_REFCLK_FREQUENCY=200
--    --   Both devices MUST share one LO and be MCS synchronised;
--    --   an independent LO per pair puts an uncalibrated phase offset
--    --   between antenna pairs that the array cannot detect.
--    -- CONNECTIONS:
--    --   device 0 -> rx_re[1:0], rx_im[1:0]
--    --   device 1 -> rx_re[3:2], rx_im[3:2]
--    --   both     -> a common rx_valid strobe at FS_ADC
--    -- ALTERNATIVE: if the AD9361 reference design is not used, a
--    -- SelectIO ISERDESE2 front end plus a small deserialiser produces
--    -- the same parallel interface; the ports below do not care which.
--
--    -- REQUIRED VIVADO IP CORE (optional capture path only):
--    -- IP NAME: AXI Direct Memory Access (axi_dma)
--    -- CONFIGURATION:
--    --   Enable Scatter Gather=0
--    --   Write Channel only (S2MM), Read Channel=0
--    --   Width of Buffer Length Register=26
--    --   Memory Map Data Width=64, Stream Data Width=32
--    --   Max Burst Size=256
--    -- CONNECTIONS:
--    --   S_AXIS_S2MM <- m_axis_* of this module
--    --   M_AXI_S2MM  -> Zynq PS S_AXI_HP0
--    --   s_axi_lite  <- PS M_AXI_GP0
--    --   mm2s/s2mm introut -> PS IRQF2P
--
--  ------------------------------------------------------------------
--  PS / PL PARTITIONING
--  ------------------------------------------------------------------
--  In the PL, because it is sample rate and must be deterministic:
--    stages 1..10 - mixer, decimator, FIR, covariance, whitening, EVD,
--    detector, weights, beamformer, DAC formatting.
--
--  In the PS, because it is 1 kHz or slower, policy rather than
--  arithmetic, or needs floating point and logarithms:
--    * MDL rank estimation and the detection hysteresis (M-of-N voting).
--      The PL exports both sides of the detector comparison so the PS
--      can apply any policy without needing a divider.
--    * AD9361 configuration over SPI, calibration, and freezing the DC
--      tracking cal - a tracking cal that runs mid-dwell changes the
--      channel response and corrupts the covariance.
--    * threshold re-calibration when fs, the dwell length or the
--      satellite count changes.
--    * logging, health monitoring, and the operator interface.
--
--  Not in either yet: the post-correlation stage (per-PRN spatial
--  clustering).  It needs correlators in the PL and clustering in the
--  PS; the hooks are the AXI4-Stream port and the dwell interrupt.
-- =====================================================================
library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

library work;
use work.asp_pkg.all;
use work.asp_coef_pkg.all;

entity asp_top is
  generic (
    G_NCH       : natural := 4;
    G_KDWELL    : natural := K_DWELL;
    G_AXI_AW    : natural := 12
  );
  port (
    -- ---- datapath clock and asynchronous reset ---------------------
    clk_dsp     : in  std_logic;
    aresetn     : in  std_logic;          -- async assert, from PS + MMCM lock

    -- ---- AXI4-Lite slave (already converted to clk_dsp) ------------
    s_axi_awaddr  : in  std_logic_vector(G_AXI_AW-1 downto 0);
    s_axi_awvalid : in  std_logic;
    s_axi_awready : out std_logic;
    s_axi_wdata   : in  std_logic_vector(31 downto 0);
    s_axi_wstrb   : in  std_logic_vector(3 downto 0);
    s_axi_wvalid  : in  std_logic;
    s_axi_wready  : out std_logic;
    s_axi_bresp   : out std_logic_vector(1 downto 0);
    s_axi_bvalid  : out std_logic;
    s_axi_bready  : in  std_logic;
    s_axi_araddr  : in  std_logic_vector(G_AXI_AW-1 downto 0);
    s_axi_arvalid : in  std_logic;
    s_axi_arready : out std_logic;
    s_axi_rdata   : out std_logic_vector(31 downto 0);
    s_axi_rresp   : out std_logic_vector(1 downto 0);
    s_axi_rvalid  : out std_logic;
    s_axi_rready  : in  std_logic;

    -- ---- RX from the AD9361 pair, parallel, FS_ADC -----------------
    rx_valid    : in  std_logic;
    rx_re       : in  std_logic_vector(G_NCH*W_ADC-1 downto 0);
    rx_im       : in  std_logic_vector(G_NCH*W_ADC-1 downto 0);

    -- ---- TX to the AD9361, FS_WORK ---------------------------------
    tx_valid    : out std_logic;
    tx_re       : out std_logic_vector(W_DAC-1 downto 0);
    tx_im       : out std_logic_vector(W_DAC-1 downto 0);

    -- ---- optional capture stream -----------------------------------
    m_axis_tvalid : out std_logic;
    m_axis_tdata  : out std_logic_vector(31 downto 0);
    m_axis_tlast  : out std_logic;
    m_axis_tready : in  std_logic;

    -- ---- interrupt to the PS ---------------------------------------
    irq         : out std_logic
  );
end entity asp_top;


architecture rtl of asp_top is

  -- reset synchroniser
  signal rst_meta : std_logic := '1';
  signal rst_sync : std_logic := '1';
  signal rst_syncn : std_logic;
  signal soft_rst : std_logic;
  signal rst      : std_logic;

  signal enable, bypass, rank2_en : std_logic;
  signal shift_init : std_logic_vector(7 downto 0);

  signal dwell_tick : std_logic;
  signal lam        : std_logic_vector(G_NCH*W_EVD-1 downto 0);
  signal det, r2el, fir_ovf : std_logic;
  signal rank       : std_logic_vector(1 downto 0);
  signal lhs, rhs   : std_logic_vector(47 downto 0);
  signal wts        : std_logic_vector(2*G_NCH*W_WGT-1 downto 0);
  signal tx_shift   : std_logic_vector(7 downto 0);
  signal clip       : std_logic_vector(31 downto 0);
  signal dwell_cnt  : std_logic_vector(31 downto 0);

  signal dp_valid : std_logic;
  signal dp_re, dp_im : std_logic_vector(W_DAC-1 downto 0);
  signal rx_valid_g : std_logic;

  signal tick_cnt : unsigned(15 downto 0) := (others => '0');

begin

  -- ------------------------------------------------------------------
  -- Reset synchroniser: asynchronous assert, synchronous release.
  -- Two flops is the minimum; the datapath reset then has a clean
  -- release edge and cannot put half the pipeline in reset for one
  -- clock while the other half is running.
  -- ------------------------------------------------------------------
  p_rst : process (clk_dsp, aresetn)
  begin
    if aresetn = '0' then
      rst_meta <= '1';
      rst_sync <= '1';
    elsif rising_edge(clk_dsp) then
      rst_meta <= '0';
      rst_sync <= rst_meta;
    end if;
  end process p_rst;

  rst       <= rst_sync or soft_rst;
  rst_syncn <= not rst_sync;

  -- Gating the input strobe rather than holding the datapath in reset
  -- means a disabled system still tracks its AGC and still reports
  -- telemetry, instead of coming back with stale state when re-enabled.
  rx_valid_g <= rx_valid and enable;

  -- ------------------------------------------------------------------
  u_dp : entity work.asp_datapath
    generic map (G_NCH => G_NCH, G_KDWELL => G_KDWELL, G_DSW => 24)
    port map (
      clk => clk_dsp, rst => rst,
      s_valid => rx_valid_g, s_re => rx_re, s_im => rx_im,
      m_valid => dp_valid, m_re => dp_re, m_im => dp_im,
      i_bypass => bypass, i_rank2_en => rank2_en,
      i_shift_init => shift_init,
      o_dwell_tick => dwell_tick, o_lam => lam, o_det => det,
      o_rank => rank, o_rank2_el => r2el,
      o_lhs => lhs, o_rhs => rhs, o_w => wts,
      o_tx_shift => tx_shift, o_clip => clip, o_fir_ovf => fir_ovf,
      o_dwell_cnt => dwell_cnt);

  -- ------------------------------------------------------------------
  u_regs : entity work.asp_axi_lite_regs
    generic map (G_NCH => G_NCH, G_AW => G_AXI_AW)
    port map (
      s_axi_aclk => clk_dsp, s_axi_aresetn => rst_syncn,
      s_axi_awaddr => s_axi_awaddr, s_axi_awvalid => s_axi_awvalid,
      s_axi_awready => s_axi_awready,
      s_axi_wdata => s_axi_wdata, s_axi_wstrb => s_axi_wstrb,
      s_axi_wvalid => s_axi_wvalid, s_axi_wready => s_axi_wready,
      s_axi_bresp => s_axi_bresp, s_axi_bvalid => s_axi_bvalid,
      s_axi_bready => s_axi_bready,
      s_axi_araddr => s_axi_araddr, s_axi_arvalid => s_axi_arvalid,
      s_axi_arready => s_axi_arready,
      s_axi_rdata => s_axi_rdata, s_axi_rresp => s_axi_rresp,
      s_axi_rvalid => s_axi_rvalid, s_axi_rready => s_axi_rready,
      o_enable => enable, o_bypass => bypass, o_rank2_en => rank2_en,
      o_soft_rst => soft_rst, o_shift_init => shift_init,
      i_dwell_tick => dwell_tick, i_lam => lam, i_det => det,
      i_rank => rank, i_rank2_el => r2el, i_lhs => lhs, i_rhs => rhs,
      i_w => wts, i_tx_shift => tx_shift, i_clip => clip,
      i_fir_ovf => fir_ovf, i_dwell_cnt => dwell_cnt,
      o_irq => irq);

  -- ------------------------------------------------------------------
  -- DAC output.
  -- ------------------------------------------------------------------
  tx_valid <= dp_valid;
  tx_re    <= dp_re;
  tx_im    <= dp_im;

  -- ------------------------------------------------------------------
  -- Capture stream.  No backpressure handling beyond dropping: the DMA
  -- must keep up or the capture is not representative anyway, and
  -- stalling here would back-pressure a datapath that has no elasticity
  -- and would corrupt the covariance.  Dropping a capture sample is
  -- harmless; stalling the array is not.
  -- ------------------------------------------------------------------
  m_axis_tvalid <= dp_valid;
  m_axis_tdata  <= std_logic_vector(resize(signed(dp_im), 16)) &
                   std_logic_vector(resize(signed(dp_re), 16));

  p_tlast : process (clk_dsp)
  begin
    if rising_edge(clk_dsp) then
      if rst = '1' then
        tick_cnt <= (others => '0');
      elsif dp_valid = '1' then
        tick_cnt <= tick_cnt + 1;
      end if;
    end if;
  end process p_tlast;

  -- one AXI4-Stream packet per 4096 samples, so the DMA can use a fixed
  -- descriptor size without needing to know the dwell length
  m_axis_tlast <= '1' when (dp_valid = '1' and tick_cnt = x"0FFF") else '0';

end architecture rtl;
