-- =====================================================================
--  tb_asp_hb_decim2  -  co-simulation of STAGE 2 against the golden model.
--
--  Feeds s02_hbdec_in.txt (which is stage 1's output) and requires the
--  RTL to reproduce s02_hbdec_out.txt exactly.
--
--  The specific failure this is here to catch is the DECIMATION PHASE.
--  Keeping the odd input indices instead of the even ones produces a
--  filter output that still looks correct in a spectrum plot, still has
--  the right group delay difference between channels (zero), and still
--  yields a plausible covariance - it is simply the wrong samples.  Only
--  a bit-exact comparison against the model finds it.
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

entity tb_asp_hb_decim2 is
  generic (
    G_IN_FILE  : string  := "gold/s02_hbdec_in.txt";
    G_OUT_FILE : string  := "gold/s02_hbdec_out.txt";
    G_MAX_SAMP : integer := 400000
  );
end entity tb_asp_hb_decim2;


architecture sim of tb_asp_hb_decim2 is

  constant NCH     : natural := 4;
  constant CLK_PER : time    := 7.637 ns;

  signal clk     : std_logic := '0';
  signal rst     : std_logic := '1';
  signal running : boolean   := true;

  signal s_valid : std_logic := '0';
  signal s_chan  : std_logic_vector(1 downto 0) := (others => '0');
  signal s_re    : std_logic_vector(W_MIX-1 downto 0) := (others => '0');
  signal s_im    : std_logic_vector(W_MIX-1 downto 0) := (others => '0');

  signal m_valid : std_logic;
  signal m_chan  : std_logic_vector(1 downto 0);
  signal m_re    : std_logic_vector(W_HB-1 downto 0);
  signal m_im    : std_logic_vector(W_HB-1 downto 0);

  signal stim_done : boolean := false;

begin

  clk <= not clk after CLK_PER/2 when running else '0';

  dut : entity work.asp_hb_decim2
    generic map (G_NCH => NCH)
    port map (
      clk => clk, rst => rst,
      s_valid => s_valid, s_chan => s_chan, s_re => s_re, s_im => s_im,
      m_valid => m_valid, m_chan => m_chan, m_re => m_re, m_im => m_im
    );

  -- Dense TDM stimulus: one channel every clk, channels cycling 0..3,
  -- exactly as stage 1 emits it.
  p_stim : process
    file     fin : text open read_mode is G_IN_FILE;
    variable v   : integer;
    variable ok  : boolean;
    variable cnt : integer := 0;
    variable ch  : integer := 0;
  begin
    rst <= '1';
    for i in 0 to 7 loop
      wait until rising_edge(clk);
    end loop;
    rst <= '0';

    loop
      exit when cnt >= G_MAX_SAMP;
      read_int(fin, v, ok);
      exit when not ok;
      s_re <= std_logic_vector(to_signed(v, W_MIX));
      read_int(fin, v, ok);
      exit when not ok;
      s_im <= std_logic_vector(to_signed(v, W_MIX));
      s_chan  <= std_logic_vector(to_unsigned(ch, 2));
      s_valid <= '1';
      wait until rising_edge(clk);
      ch  := (ch + 1) mod NCH;
      cnt := cnt + 1;
    end loop;
    s_valid <= '0';

    for i in 0 to 63 loop
      wait until rising_edge(clk);
    end loop;
    stim_done <= true;
    wait;
  end process p_stim;

  p_check : process
    file     fexp  : text open read_mode is G_OUT_FILE;
    variable want  : integer;
    variable ok    : boolean;
    variable err   : integer := 0;
    variable shown : integer := 0;
    variable idx   : integer := 0;
    variable expch : integer := 0;
  begin
    loop
      wait until rising_edge(clk);
      exit when stim_done;
      if m_valid = '1' then
        assert to_integer(unsigned(m_chan)) = expch
          report "tb_asp_hb_decim2: TDM channel order broken at index " &
                 integer'image(idx)
          severity failure;
        expch := (expch + 1) mod NCH;

        read_int(fexp, want, ok);
        exit when not ok;
        check_int("tb_asp_hb_decim2.re", idx,
                  to_integer(signed(m_re)), want, err, shown);
        idx := idx + 1;

        read_int(fexp, want, ok);
        exit when not ok;
        check_int("tb_asp_hb_decim2.im", idx,
                  to_integer(signed(m_im)), want, err, shown);
        idx := idx + 1;
      end if;
    end loop;

    summarise("tb_asp_hb_decim2", idx, err);
    running <= false;
    wait;
  end process p_check;

end architecture sim;
