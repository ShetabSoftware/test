-- =====================================================================
--  tb_asp_datapath  -  the whole PL chain, ADC in to DAC out.
--
--  Drives the raw AD9361 samples the model consumed (s01_ddc_in.txt) and
--  checks the system-level behaviour that no single-stage test can:
--
--    * the WEIGHT SEQUENCE.  Every weight vector the chain produces must
--      equal the model's, bit for bit.  This is the strongest available
--      end-to-end check, because the weights are a function of every
--      stage from the mixer through the EVD - if any of them drifted by
--      one LSB the eigenvector would rotate and the weights would not
--      match.
--
--    * the causal SCHEDULE, including its one unavoidable difference
--      from the model.  The model applies dwell k's estimate to dwell
--      k+1; hardware cannot, because the estimator needs ~5400 clocks
--      after the covariance closes and the covariance only closes at the
--      dwell boundary.  The chain therefore applies dwell k's estimate
--      from dwell k+2.  The check below follows that schedule
--      explicitly rather than papering over it.
--
--    * the detector fires on every dwell (the reference scenario has a
--      spoofer at SAPR +5.5 dB), and the FIR input FIFO never overflows.
--
--  This testbench is slow - it simulates several million clocks - and
--  that is the point: it is the only test that exercises the dwell
--  handshake between blocks, which is where integration bugs live.
-- =====================================================================
library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

library std;
use std.textio.all;

library work;
use work.asp_pkg.all;
use work.asp_coef_pkg.all;
use work.asp_tb_pkg.all;

entity tb_asp_datapath is
  generic (
    G_IN_FILE : string := "gold/s01_ddc_in.txt";
    G_W_FILE  : string := "gold/s08_weights.txt";
    G_MAX     : integer := 400000
  );
end entity tb_asp_datapath;


architecture sim of tb_asp_datapath is

  constant CLK_PER : time := 7.637 ns;
  constant NCH     : natural := 4;

  signal clk     : std_logic := '0';
  signal rst     : std_logic := '1';
  signal running : boolean   := true;

  signal s_valid : std_logic := '0';
  signal s_re    : std_logic_vector(NCH*W_ADC-1 downto 0) := (others => '0');
  signal s_im    : std_logic_vector(NCH*W_ADC-1 downto 0) := (others => '0');

  signal m_valid : std_logic;
  signal m_re    : std_logic_vector(W_DAC-1 downto 0);
  signal m_im    : std_logic_vector(W_DAC-1 downto 0);

  signal o_tick  : std_logic;
  signal o_lam   : std_logic_vector(NCH*W_EVD-1 downto 0);
  signal o_det   : std_logic;
  signal o_rank  : std_logic_vector(1 downto 0);
  signal o_el2   : std_logic;
  signal o_lhs, o_rhs : std_logic_vector(47 downto 0);
  signal o_w     : std_logic_vector(2*NCH*W_WGT-1 downto 0);
  signal o_shift : std_logic_vector(7 downto 0);
  signal o_clip  : std_logic_vector(31 downto 0);
  signal o_ovf   : std_logic;
  signal o_dwell : std_logic_vector(31 downto 0);

  signal stim_done : boolean := false;
  signal n_dac     : integer := 0;

begin

  clk <= not clk after CLK_PER/2 when running else '0';

  dut : entity work.asp_datapath
    generic map (G_NCH => NCH, G_KDWELL => K_DWELL, G_DSW => 24)
    port map (
      clk => clk, rst => rst,
      s_valid => s_valid, s_re => s_re, s_im => s_im,
      m_valid => m_valid, m_re => m_re, m_im => m_im,
      i_bypass => '0', i_rank2_en => '0', i_shift_init => x"00",
      o_dwell_tick => o_tick, o_lam => o_lam, o_det => o_det,
      o_rank => o_rank, o_rank2_el => o_el2,
      o_lhs => o_lhs, o_rhs => o_rhs, o_w => o_w,
      o_tx_shift => o_shift, o_clip => o_clip, o_fir_ovf => o_ovf,
      o_dwell_cnt => o_dwell);

  p_stim : process
    file     fin : text open read_mode is G_IN_FILE;
    variable v   : integer;
    variable ok  : boolean;
    variable cnt : integer := 0;
  begin
    rst <= '1';
    for i in 0 to 15 loop
      wait until rising_edge(clk);
    end loop;
    rst <= '0';
    wait until rising_edge(clk);

    loop
      exit when cnt >= G_MAX;
      for k in 0 to NCH-1 loop
        read_int(fin, v, ok);
        exit when not ok;
        s_re((k+1)*W_ADC-1 downto k*W_ADC) <=
          std_logic_vector(to_signed(v, W_ADC));
        read_int(fin, v, ok);
        exit when not ok;
        s_im((k+1)*W_ADC-1 downto k*W_ADC) <=
          std_logic_vector(to_signed(v, W_ADC));
      end loop;
      exit when not ok;

      s_valid <= '1';
      wait until rising_edge(clk);
      s_valid <= '0';
      for i in 1 to NCH-1 loop
        wait until rising_edge(clk);
      end loop;
      cnt := cnt + 1;
    end loop;
    s_valid <= '0';

    for i in 0 to 4095 loop
      wait until rising_edge(clk);
    end loop;
    stim_done <= true;
    wait;
  end process p_stim;

  -- count DAC samples so a chain that silently stops is visible
  p_dac : process (clk)
  begin
    if rising_edge(clk) then
      if m_valid = '1' then
        n_dac <= n_dac + 1;
      end if;
    end if;
  end process p_dac;

  p_check : process
    file     fw    : text open read_mode is G_W_FILE;
    variable want  : integer;
    variable ok    : boolean;
    variable err   : integer := 0;
    variable shown : integer := 0;
    variable cnt   : integer := 0;
    variable ntick : integer := 0;
    variable ndet  : integer := 0;
    variable nchk  : integer := 0;
  begin
    loop
      wait until rising_edge(clk);
      exit when stim_done;

      if o_tick = '1' then
        ntick := ntick + 1;

        -- Wait for the detector and the weight update that this tick
        -- triggers to settle, then look at what the chain is now using.
        for i in 0 to 63 loop
          wait until rising_edge(clk);
        end loop;

        if o_det = '1' then
          ndet := ndet + 1;
        end if;

        -- Schedule: at the end of dwell k the chain promotes the weights
        -- estimated from dwell k-1.  So after tick 2 the applied weights
        -- are the model's dwell-1 weights, after tick 3 the dwell-2
        -- weights, and so on.
        if ntick >= 2 then
          for k in 0 to NCH-1 loop
            read_int(fw, want, ok);
            exit when not ok;
            check_int("tb_asp_datapath.w.re", cnt,
              to_integer(signed(o_w((2*k+1)*W_WGT-1 downto (2*k)*W_WGT))),
              want, err, shown);
            cnt := cnt + 1;
            read_int(fw, want, ok);
            exit when not ok;
            check_int("tb_asp_datapath.w.im", cnt,
              to_integer(signed(o_w((2*k+2)*W_WGT-1 downto (2*k+1)*W_WGT))),
              want, err, shown);
            cnt := cnt + 1;
          end loop;
          exit when not ok;
          nchk := nchk + 1;
        end if;
      end if;
    end loop;

    assert ntick >= 3
      report "tb_asp_datapath: only " & integer'image(ntick) &
             " dwells completed" severity failure;
    assert nchk >= 2
      report "tb_asp_datapath: only " & integer'image(nchk) &
             " weight vectors checked" severity failure;
    assert ndet = ntick
      report "tb_asp_datapath: detector fired on " & integer'image(ndet) &
             " of " & integer'image(ntick) &
             " dwells; the reference scenario has a spoofer on every one"
      severity failure;
    assert o_ovf = '0'
      report "tb_asp_datapath: the shaping FIR input FIFO overflowed"
      severity failure;
    assert n_dac > 40000
      report "tb_asp_datapath: only " & integer'image(n_dac) &
             " DAC samples produced; the chain stalled"
      severity failure;

    summarise("tb_asp_datapath", cnt, err);
    report "tb_asp_datapath: " & integer'image(ntick) & " dwells, " &
           integer'image(nchk) & " weight vectors bit-identical, detector " &
           "fired on every dwell, " & integer'image(n_dac) &
           " DAC samples" severity note;
    running <= false;
    wait;
  end process p_check;

end architecture sim;
