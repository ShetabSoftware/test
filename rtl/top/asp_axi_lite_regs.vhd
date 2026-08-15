-- =====================================================================
--  asp_axi_lite_regs  -  AXI4-Lite control and status register file.
--
--  ------------------------------------------------------------------
--  CLOCK DOMAIN
--  ------------------------------------------------------------------
--  This slave runs entirely on clk_dsp (130.944 MHz), NOT on the PS AXI
--  clock.  There is no manual clock-domain crossing anywhere in this
--  file, and that is deliberate: the crossing is done once, by an AXI
--  Clock Converter in the block design, where it is a reviewed IP
--  instance rather than a hand-rolled synchroniser per register.
--
--  The reason clk_dsp cannot simply be the PS clock is sample coherence:
--  every rate in the datapath is an exact integer ratio of the AD9361
--  sample clock, so clk_dsp must be derived from that clock and not from
--  the PS PLL.  See the header of asp_top.
--
--    -- REQUIRED VIVADO IP CORE:
--    -- IP NAME: AXI Clock Converter (axi_clock_converter)
--    -- CONFIGURATION:
--    --   PROTOCOL=AXI4LITE
--    --   ADDR_WIDTH=12
--    --   DATA_WIDTH=32
--    --   ID_WIDTH=0
--    --   SI_CLK: FCLK_CLK0 from the Zynq PS (100 MHz)
--    --   MI_CLK: clk_dsp (130.944 MHz)
--    --   ASYNC_CLK=1
--    -- CONNECTIONS:
--    --   S_AXI  <- Zynq PS M_AXI_GP0 (through the AXI Interconnect)
--    --   M_AXI  -> this module's s_axi_* ports
--    --   s_axi_aclk/aresetn from the PS, m_axi_aclk/aresetn from clk_dsp
--
--  ------------------------------------------------------------------
--  REGISTER MAP (byte offsets, 32-bit registers)
--  ------------------------------------------------------------------
--   0x00  ID           RO   0x41535002  ("AS" + revision 2)
--   0x04  CONTROL      RW   [0] enable      [1] bypass (force quiescent)
--                           [2] rank2_en    [3] soft reset (self clearing)
--   0x08  STATUS       RO   [0] detected    [2:1] rank
--                           [3] rank2 eligible
--                           [4] FIR FIFO overflow (sticky)
--   0x0C  DWELL_CNT    RO   dwells since reset
--   0x10  LAMBDA0..3   RO   eigenvalues, s32.26, descending  (0x10..0x1C)
--   0x20  DET_LHS_LO   RO   detector left  side, low  32 bits
--   0x24  DET_LHS_HI   RO   detector left  side, high 16 bits
--   0x28  DET_RHS_LO   RO   detector right side, low  32 bits
--   0x2C  DET_RHS_HI   RO   detector right side, high 16 bits
--   0x30  W0_RE .. 0x4C W3_IM  RO  applied weights, s18.16 sign extended
--   0x50  TX_SHIFT     RO   signed 8-bit AGC shift currently in force
--   0x54  CLIP_CNT     RO   DAC samples clipped since reset
--   0x58  SHIFT_INIT   RW   signed 8-bit shift used for dwell 0
--   0x5C  IRQ_STATUS   RW1C [0] dwell complete
--   0x60  IRQ_ENABLE   RW   [0] dwell complete
--
--  WHY THE DETECTOR OPERANDS ARE EXPORTED AS TWO WORDS EACH
--  The statistic is lambda_1 / mean(lambda_2..N), and the hardware
--  deliberately never forms that ratio - it cross-multiplies to avoid a
--  divider.  Exporting both sides lets the PS apply its own MDL test and
--  its own hysteresis policy at 1 kHz without needing a divider either,
--  and without the PL having to guess what policy the PS wants.
--
--  RESOURCES  ~450 FF, ~600 LUT, 0 DSP, 0 BRAM.
-- =====================================================================
library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

library work;
use work.asp_pkg.all;
use work.asp_coef_pkg.all;

entity asp_axi_lite_regs is
  generic (
    G_NCH  : natural := 4;
    G_AW   : natural := 12                     -- byte address width
  );
  port (
    -- AXI4-Lite slave, all on clk
    s_axi_aclk    : in  std_logic;
    s_axi_aresetn : in  std_logic;

    s_axi_awaddr  : in  std_logic_vector(G_AW-1 downto 0);
    s_axi_awvalid : in  std_logic;
    s_axi_awready : out std_logic;
    s_axi_wdata   : in  std_logic_vector(31 downto 0);
    s_axi_wstrb   : in  std_logic_vector(3 downto 0);
    s_axi_wvalid  : in  std_logic;
    s_axi_wready  : out std_logic;
    s_axi_bresp   : out std_logic_vector(1 downto 0);
    s_axi_bvalid  : out std_logic;
    s_axi_bready  : in  std_logic;
    s_axi_araddr  : in  std_logic_vector(G_AW-1 downto 0);
    s_axi_arvalid : in  std_logic;
    s_axi_arready : out std_logic;
    s_axi_rdata   : out std_logic_vector(31 downto 0);
    s_axi_rresp   : out std_logic_vector(1 downto 0);
    s_axi_rvalid  : out std_logic;
    s_axi_rready  : in  std_logic;

    -- to the datapath
    o_enable      : out std_logic;
    o_bypass      : out std_logic;
    o_rank2_en    : out std_logic;
    o_soft_rst    : out std_logic;
    o_shift_init  : out std_logic_vector(7 downto 0);

    -- from the datapath
    i_dwell_tick  : in  std_logic;
    i_lam         : in  std_logic_vector(G_NCH*W_EVD-1 downto 0);
    i_det         : in  std_logic;
    i_rank        : in  std_logic_vector(1 downto 0);
    i_rank2_el    : in  std_logic;
    i_lhs         : in  std_logic_vector(47 downto 0);
    i_rhs         : in  std_logic_vector(47 downto 0);
    i_w           : in  std_logic_vector(2*G_NCH*W_WGT-1 downto 0);
    i_tx_shift    : in  std_logic_vector(7 downto 0);
    i_clip        : in  std_logic_vector(31 downto 0);
    i_fir_ovf     : in  std_logic;
    i_dwell_cnt   : in  std_logic_vector(31 downto 0);

    -- level interrupt to the PS
    o_irq         : out std_logic
  );
end entity asp_axi_lite_regs;


architecture rtl of asp_axi_lite_regs is

  constant ID_VALUE : std_logic_vector(31 downto 0) := x"41535002";

  signal awready_r, wready_r, bvalid_r : std_logic := '0';
  signal arready_r, rvalid_r : std_logic := '0';
  signal awaddr_r : unsigned(G_AW-1 downto 0) := (others => '0');
  signal araddr_r : unsigned(G_AW-1 downto 0) := (others => '0');
  signal rdata_r  : std_logic_vector(31 downto 0) := (others => '0');

  signal ctrl_r   : std_logic_vector(31 downto 0) := (others => '0');
  signal shinit_r : std_logic_vector(7 downto 0)  := (others => '0');
  signal irq_st   : std_logic_vector(31 downto 0) := (others => '0');
  signal irq_en   : std_logic_vector(31 downto 0) := (others => '0');
  signal soft_rst : std_logic := '0';

  signal ovf_sticky : std_logic := '0';

  -- word index = byte address / 4
  function widx (a : unsigned) return integer is
  begin
    return to_integer(a(G_AW-1 downto 2));
  end function;

begin

  -- ------------------------------------------------------------------
  -- Write channel.  Address and data are accepted independently and the
  -- write commits when both have arrived, which is what AXI4-Lite
  -- requires; a slave that demands them in the same cycle works with
  -- most masters and deadlocks with some.
  -- ------------------------------------------------------------------
  p_write : process (s_axi_aclk)
  begin
    if rising_edge(s_axi_aclk) then
      if s_axi_aresetn = '0' then
        awready_r <= '0';
        wready_r  <= '0';
        bvalid_r  <= '0';
        ctrl_r    <= (others => '0');
        shinit_r  <= (others => '0');
        irq_en    <= (others => '0');
        irq_st    <= (others => '0');
        soft_rst  <= '0';
      else
        soft_rst <= '0';                    -- self-clearing, one cycle

        if awready_r = '0' and s_axi_awvalid = '1' then
          awready_r <= '1';
          awaddr_r  <= unsigned(s_axi_awaddr);
        elsif awready_r = '1' then
          awready_r <= '0';
        end if;

        if wready_r = '0' and s_axi_wvalid = '1' then
          wready_r <= '1';
        elsif wready_r = '1' then
          wready_r <= '0';
        end if;

        if awready_r = '1' and wready_r = '1' then
          case widx(awaddr_r) is
            when 1 =>                                     -- 0x04 CONTROL
              ctrl_r <= s_axi_wdata;
              if s_axi_wdata(3) = '1' then
                soft_rst <= '1';
              end if;
            when 22 =>                                    -- 0x58 SHIFT_INIT
              shinit_r <= s_axi_wdata(7 downto 0);
            when 23 =>                                    -- 0x5C IRQ_STATUS
              -- write-1-to-clear
              irq_st <= irq_st and (not s_axi_wdata);
            when 24 =>                                    -- 0x60 IRQ_ENABLE
              irq_en <= s_axi_wdata;
            when others =>
              null;                                       -- read-only
          end case;
          bvalid_r <= '1';
        elsif bvalid_r = '1' and s_axi_bready = '1' then
          bvalid_r <= '0';
        end if;

        -- the dwell interrupt sets after the write-clear so a tick that
        -- lands in the same cycle as the clear is not lost
        if i_dwell_tick = '1' then
          irq_st(0) <= '1';
        end if;

        if i_fir_ovf = '1' then
          ovf_sticky <= '1';
        end if;
      end if;
    end if;
  end process p_write;

  -- ------------------------------------------------------------------
  -- Read channel.  Telemetry is read combinationally from the datapath's
  -- per-dwell latched registers: they only change at a dwell boundary,
  -- so a read can never catch a half-updated eigenvalue set.
  -- ------------------------------------------------------------------
  p_read : process (s_axi_aclk)
    variable i : integer;
  begin
    if rising_edge(s_axi_aclk) then
      if s_axi_aresetn = '0' then
        arready_r <= '0';
        rvalid_r  <= '0';
        rdata_r   <= (others => '0');
      else
        if arready_r = '0' and s_axi_arvalid = '1' then
          arready_r <= '1';
          araddr_r  <= unsigned(s_axi_araddr);
        elsif arready_r = '1' then
          arready_r <= '0';
        end if;

        if arready_r = '1' and rvalid_r = '0' then
          i := widx(araddr_r);
          rdata_r <= (others => '0');
          case i is
            when 0 => rdata_r <= ID_VALUE;
            when 1 => rdata_r <= ctrl_r;
            when 2 =>
              rdata_r(0) <= i_det;
              rdata_r(2 downto 1) <= i_rank;
              rdata_r(3) <= i_rank2_el;
              rdata_r(4) <= ovf_sticky;
            when 3 => rdata_r <= i_dwell_cnt;
            when 4 | 5 | 6 | 7 =>
              rdata_r <= i_lam((i-3)*W_EVD-1 downto (i-4)*W_EVD);
            when 8  => rdata_r <= i_lhs(31 downto 0);
            when 9  => rdata_r <= std_logic_vector(resize(signed(i_lhs(47 downto 32)), 32));
            when 10 => rdata_r <= i_rhs(31 downto 0);
            when 11 => rdata_r <= std_logic_vector(resize(signed(i_rhs(47 downto 32)), 32));
            when 12 to 19 =>
              -- weights, sign extended from s18.16 to 32 bits so the PS
              -- can read them as plain signed integers
              rdata_r <= std_logic_vector(resize(
                signed(i_w((i-11)*W_WGT-1 downto (i-12)*W_WGT)), 32));
            when 20 => rdata_r <= std_logic_vector(resize(signed(i_tx_shift), 32));
            when 21 => rdata_r <= i_clip;
            when 22 => rdata_r <= std_logic_vector(resize(signed(shinit_r), 32));
            when 23 => rdata_r <= irq_st;
            when 24 => rdata_r <= irq_en;
            when others => rdata_r <= (others => '0');
          end case;
          rvalid_r <= '1';
        elsif rvalid_r = '1' and s_axi_rready = '1' then
          rvalid_r <= '0';
        end if;
      end if;
    end if;
  end process p_read;

  s_axi_awready <= awready_r;
  s_axi_wready  <= wready_r;
  s_axi_bvalid  <= bvalid_r;
  s_axi_bresp   <= "00";                     -- OKAY, always
  s_axi_arready <= arready_r;
  s_axi_rvalid  <= rvalid_r;
  s_axi_rresp   <= "00";
  s_axi_rdata   <= rdata_r;

  o_enable     <= ctrl_r(0);
  o_bypass     <= ctrl_r(1);
  o_rank2_en   <= ctrl_r(2);
  o_soft_rst   <= soft_rst;
  o_shift_init <= shinit_r;

  -- Level interrupt.  The PS IRQF2P inputs are level sensitive and are
  -- synchronised inside the PS, so no extra pulse stretching is needed
  -- provided the level is held until the ISR clears IRQ_STATUS - which
  -- it is.
  o_irq <= '1' when (irq_st and irq_en) /= x"00000000" else '0';

end architecture rtl;
