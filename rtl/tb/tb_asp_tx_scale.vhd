-- =====================================================================
--  tb_asp_tx_scale  -  co-simulation of STAGE 10 against the golden model.
--
--  DWELL 1 IS DELIBERATELY EXCLUDED, and that exclusion is the most
--  important thing this testbench says.
--
--  The model's first dwell uses "fast acquisition": it measures the RMS
--  of a dwell and applies the resulting shift TO THAT SAME DWELL.  No
--  hardware can do that without buffering the whole dwell (16368 complex
--  samples, ~15 BRAM36) purely to improve the first millisecond after
--  reset.  The RTL therefore scales dwell 1 with i_shift_init and
--  performs the identical acquisition jump at the END of dwell 1, so
--  from dwell 2 onward the shift sequence is exactly the model's.
--
--  This test verifies that claim rather than assuming it: it skips the
--  first dwell's samples and then requires every remaining sample to be
--  bit-identical.  If the acquisition arithmetic or the hysteresis
--  differed by even one step, dwell 2 would be off by a factor of two
--  and every sample in it would mismatch.
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

entity tb_asp_tx_scale is
  generic (
    G_IN_FILE  : string := "gold/s10_dac_in.txt";
    G_OUT_FILE : string := "gold/s10_dac_out.txt"
  );
end entity tb_asp_tx_scale;


architecture sim of tb_asp_tx_scale is

  constant CLK_PER : time := 7.637 ns;

  signal clk     : std_logic := '0';
  signal rst     : std_logic := '1';
  signal running : boolean   := true;

  signal s_valid : std_logic := '0';
  signal s_re    : std_logic_vector(W_BEAM-1 downto 0) := (others => '0');
  signal s_im    : std_logic_vector(W_BEAM-1 downto 0) := (others => '0');
  signal m_valid : std_logic;
  signal m_re    : std_logic_vector(W_DAC-1 downto 0);
  signal m_im    : std_logic_vector(W_DAC-1 downto 0);
  signal o_shift : std_logic_vector(7 downto 0);
  signal o_clip  : std_logic_vector(31 downto 0);
  signal o_tick  : std_logic;

  signal stim_done : boolean := false;

begin

  clk <= not clk after CLK_PER/2 when running else '0';

  dut : entity work.asp_tx_scale
    generic map (G_KDWELL => K_DWELL)
    port map (clk => clk, rst => rst,
              i_shift_init => x"00",
              s_valid => s_valid, s_re => s_re, s_im => s_im,
              m_valid => m_valid, m_re => m_re, m_im => m_im,
              o_shift => o_shift, o_clip => o_clip, o_dwell_tick => o_tick);

  p_stim : process
    file     fin : text open read_mode is G_IN_FILE;
    variable v   : integer;
    variable ok  : boolean;
  begin
    rst <= '1';
    for i in 0 to 7 loop
      wait until rising_edge(clk);
    end loop;
    rst <= '0';

    loop
      read_int(fin, v, ok);
      exit when not ok;
      s_re <= std_logic_vector(to_signed(v, W_BEAM));
      read_int(fin, v, ok);
      exit when not ok;
      s_im    <= std_logic_vector(to_signed(v, W_BEAM));
      s_valid <= '1';
      wait until rising_edge(clk);
      s_valid <= '0';
      for i in 1 to 7 loop           -- one sample per 8 clk, as in hardware
        wait until rising_edge(clk);
      end loop;
    end loop;
    s_valid <= '0';

    for i in 0 to 31 loop
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
    variable cmp   : integer := 0;
  begin
    loop
      wait until rising_edge(clk);
      exit when stim_done;
      if m_valid = '1' then
        read_int(fexp, want, ok);
        exit when not ok;
        if idx >= 2*K_DWELL then      -- skip dwell 1 (see the header)
          check_int("tb_asp_tx_scale.re", idx,
                    to_integer(signed(m_re)), want, err, shown);
          cmp := cmp + 1;
        end if;
        idx := idx + 1;

        read_int(fexp, want, ok);
        exit when not ok;
        if idx >= 2*K_DWELL then
          check_int("tb_asp_tx_scale.im", idx,
                    to_integer(signed(m_im)), want, err, shown);
          cmp := cmp + 1;
        end if;
        idx := idx + 1;
      end if;
    end loop;

    assert cmp > 50000
      report "tb_asp_tx_scale: only " & integer'image(cmp) &
             " values compared after the excluded first dwell"
      severity failure;

    summarise("tb_asp_tx_scale", cmp, err);
    report "tb_asp_tx_scale: dwell 1 excluded by design (non-causal fast " &
           "acquisition in the model); " & integer'image(cmp) &
           " values compared from dwell 2 onward" severity note;
    running <= false;
    wait;
  end process p_check;

end architecture sim;
