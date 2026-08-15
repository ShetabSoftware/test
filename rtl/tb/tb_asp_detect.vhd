-- =====================================================================
--  tb_asp_detect  -  co-simulation of STAGE 7 against the golden model.
--
--  Drives the model's eigenvalues (s06_evd_lam.txt) and checks all four
--  outputs against s07_detect.txt: the two sides of the cross-multiplied
--  comparison, the detection flag and the rank.
--
--  Checking lhs and rhs rather than only the flag is deliberate.  The
--  flag is one bit and would pass with a threshold that is wrong by a
--  factor of two on any dwell where the statistic is far from the
--  boundary - which is most of them.  The two operands pin the actual
--  arithmetic down.
--
--  RANK2_ENABLE is driven low, matching the model's default.  Rank 2 is
--  gated on the second eigenvalue being genuinely resolvable and the PS
--  withholds it until its own MDL test agrees.
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

entity tb_asp_detect is
  generic (
    G_LAM_FILE : string := "gold/s06_evd_lam.txt";
    G_DET_FILE : string := "gold/s07_detect.txt"
  );
end entity tb_asp_detect;


architecture sim of tb_asp_detect is

  constant CLK_PER : time := 7.637 ns;
  constant NA      : natural := 4;

  signal clk     : std_logic := '0';
  signal rst     : std_logic := '1';
  signal running : boolean   := true;

  signal i_start : std_logic := '0';
  signal i_lam   : std_logic_vector(NA*W_EVD-1 downto 0) := (others => '0');
  signal o_valid : std_logic;
  signal o_lhs   : std_logic_vector(47 downto 0);
  signal o_rhs   : std_logic_vector(47 downto 0);
  signal o_det   : std_logic;
  signal o_rank  : std_logic_vector(1 downto 0);
  signal o_el2   : std_logic;

begin

  clk <= not clk after CLK_PER/2 when running else '0';

  dut : entity work.asp_detect
    generic map (G_N => NA)
    port map (clk => clk, rst => rst,
              i_start => i_start, i_lam => i_lam, i_rank2_en => '0',
              o_valid => o_valid, o_lhs => o_lhs, o_rhs => o_rhs,
              o_det => o_det, o_rank => o_rank, o_rank2_eligible => o_el2);

  p_main : process
    file     flam : text open read_mode is G_LAM_FILE;
    file     fdet : text open read_mode is G_DET_FILE;
    variable lv   : integer;
    variable want : big_t;
    variable wi   : integer;
    variable ok   : boolean;
    variable err  : integer := 0;
    variable shown: integer := 0;
    variable cnt  : integer := 0;
    variable ndw  : integer := 0;
  begin
    rst <= '1';
    for i in 0 to 7 loop
      wait until rising_edge(clk);
    end loop;
    rst <= '0';
    wait until rising_edge(clk);

    loop
      for i in 0 to NA-1 loop
        read_int(flam, lv, ok);
        exit when not ok;
        i_lam((i+1)*W_EVD-1 downto i*W_EVD) <=
          std_logic_vector(to_signed(lv, W_EVD));
      end loop;
      exit when not ok;

      i_start <= '1';
      wait until rising_edge(clk);
      i_start <= '0';
      loop
        wait until rising_edge(clk);
        exit when o_valid = '1';
      end loop;

      -- s07_detect.txt holds lhs, rhs, detected, rank per dwell.  lhs and
      -- rhs exceed 32 bits, so they are read as wide values.
      read_big(fdet, want, ok);
      exit when not ok;
      check_big("tb_asp_detect.lhs", cnt, resize(signed(o_lhs), 64), want, err, shown);
      cnt := cnt + 1;

      read_big(fdet, want, ok);
      exit when not ok;
      check_big("tb_asp_detect.rhs", cnt, resize(signed(o_rhs), 64), want, err, shown);
      cnt := cnt + 1;

      read_int(fdet, wi, ok);
      exit when not ok;
      if o_det = '1' then
        check_int("tb_asp_detect.detected", cnt, 1, wi, err, shown);
      else
        check_int("tb_asp_detect.detected", cnt, 0, wi, err, shown);
      end if;
      cnt := cnt + 1;

      read_int(fdet, wi, ok);
      exit when not ok;
      check_int("tb_asp_detect.rank", cnt,
                to_integer(unsigned(o_rank)), wi, err, shown);
      cnt := cnt + 1;

      ndw := ndw + 1;
    end loop;

    assert ndw >= 2
      report "tb_asp_detect: only " & integer'image(ndw) & " dwells checked"
      severity failure;

    summarise("tb_asp_detect", cnt, err);
    running <= false;
    wait;
  end process p_main;

end architecture sim;
