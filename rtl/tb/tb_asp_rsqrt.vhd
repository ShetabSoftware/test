-- =====================================================================
--  tb_asp_rsqrt  -  the two square-root units against MATLAB.
--
--  Covers asp_rsqrt (reciprocal square root, range reduction + four
--  Newton steps) and asp_isqrt (restoring integer square root), because
--  they are companions: stage 5 calls both on the same R_ii.
--
--  The vectors deliberately include the awkward inputs rather than only
--  random ones - powers of two and their neighbours, where the range
--  reduction changes k; 1, where the reduction shifts LEFT; and the top
--  of the range, where the weight stage's squared norm lives.  A block
--  that is right on random data and wrong at 2^16 is the normal failure
--  mode for range-reduced iterations, so random data alone would not be
--  a test.
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

entity tb_asp_rsqrt is
  generic (
    G_RSQ_FILE : string := "gold/unit_rsqrt.txt";
    G_ISQ_FILE : string := "gold/unit_isqrt.txt"
  );
end entity tb_asp_rsqrt;


architecture sim of tb_asp_rsqrt is

  constant CLK_PER : time := 7.637 ns;
  constant AW      : natural := 48;

  signal clk     : std_logic := '0';
  signal rst     : std_logic := '1';
  signal running : boolean   := true;

  signal r_start : std_logic := '0';
  signal r_a     : std_logic_vector(AW-1 downto 0) := (others => '0');
  signal r_done  : std_logic;
  signal r_y     : std_logic_vector(RSQ_F+1 downto 0);
  signal r_k     : std_logic_vector(7 downto 0);

  signal i_start : std_logic := '0';
  signal i_a     : std_logic_vector(AW-1 downto 0) := (others => '0');
  signal i_done  : std_logic;
  signal i_s     : std_logic_vector(AW/2-1 downto 0);

begin

  clk <= not clk after CLK_PER/2 when running else '0';

  u_rsqrt : entity work.asp_rsqrt
    generic map (G_AW => AW)
    port map (clk => clk, rst => rst, i_start => r_start, i_a => r_a,
              o_done => r_done, o_y => r_y, o_k => r_k);

  u_isqrt : entity work.asp_isqrt
    generic map (G_AW => AW)
    port map (clk => clk, rst => rst, i_start => i_start, i_a => i_a,
              o_done => i_done, o_s => i_s);

  p_main : process
    file     frs   : text open read_mode is G_RSQ_FILE;
    file     fis   : text open read_mode is G_ISQ_FILE;
    variable a     : big_t;
    variable wy, wk, ws : integer;
    variable ok    : boolean;
    variable err   : integer := 0;
    variable shown : integer := 0;
    variable n     : integer := 0;
  begin
    rst <= '1';
    for i in 0 to 7 loop
      wait until rising_edge(clk);
    end loop;
    rst <= '0';
    wait until rising_edge(clk);

    -- ---------------------------------------------------- rsqrt
    loop
      read_big(frs, a, ok);
      exit when not ok;
      read_int(frs, wy, ok);
      exit when not ok;
      read_int(frs, wk, ok);
      exit when not ok;

      r_a     <= std_logic_vector(a(AW-1 downto 0));
      r_start <= '1';
      wait until rising_edge(clk);
      r_start <= '0';
      loop
        wait until rising_edge(clk);
        exit when r_done = '1';
      end loop;

      check_int("tb_asp_rsqrt.Y", n, to_integer(unsigned(r_y)), wy, err, shown);
      check_int("tb_asp_rsqrt.k", n, to_integer(signed(r_k)),   wk, err, shown);
      n := n + 1;
    end loop;

    report "tb_asp_rsqrt: " & integer'image(n) &
           " reciprocal-square-root cases checked" severity note;

    -- ---------------------------------------------------- isqrt
    n := 0;
    loop
      read_big(fis, a, ok);
      exit when not ok;
      read_int(fis, ws, ok);
      exit when not ok;

      i_a     <= std_logic_vector(a(AW-1 downto 0));
      i_start <= '1';
      wait until rising_edge(clk);
      i_start <= '0';
      loop
        wait until rising_edge(clk);
        exit when i_done = '1';
      end loop;

      check_int("tb_asp_isqrt.s", n, to_integer(unsigned(i_s)), ws, err, shown);
      n := n + 1;
    end loop;

    summarise("tb_asp_rsqrt+isqrt", n, err);
    running <= false;
    wait;
  end process p_main;

end architecture sim;
