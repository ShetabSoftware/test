-- =====================================================================
--  tb_asp_cordic  -  both CORDIC modes against MATLAB.
--
--  Vectoring is checked on every quadrant and on both axes, because the
--  left-half-plane pre-rotation is where a hand-written CORDIC usually
--  goes wrong and a first-quadrant-only test cannot see it.  Magnitudes
--  run up to 2^29, which is what the Jacobi engine actually presents.
--
--  Rotation is checked across the full (-pi, pi] range including the
--  fold points at +/-pi/2, where the sign flip on both outputs is
--  applied.
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

entity tb_asp_cordic is
  generic (
    G_VEC_FILE : string := "gold/unit_cordic_v.txt";
    G_ROT_FILE : string := "gold/unit_cordic_r.txt"
  );
end entity tb_asp_cordic;


architecture sim of tb_asp_cordic is

  constant CLK_PER : time := 7.637 ns;
  constant W  : natural := 40;
  constant ZW : natural := 22;

  signal clk     : std_logic := '0';
  signal rst     : std_logic := '1';
  signal running : boolean   := true;

  signal c_start : std_logic := '0';
  signal c_mode  : std_logic := '0';
  signal c_x, c_y: std_logic_vector(W-1 downto 0) := (others => '0');
  signal c_zi    : std_logic_vector(ZW-1 downto 0) := (others => '0');
  signal c_done  : std_logic;
  signal c_m     : std_logic_vector(W-1 downto 0);
  signal c_zo    : std_logic_vector(ZW-1 downto 0);
  signal c_c     : std_logic_vector(W_ROT-1 downto 0);
  signal c_s     : std_logic_vector(W_ROT-1 downto 0);

begin

  clk <= not clk after CLK_PER/2 when running else '0';

  dut : entity work.asp_cordic
    generic map (G_W => W, G_ZW => ZW)
    port map (clk => clk, rst => rst,
              i_start => c_start, i_mode => c_mode,
              i_x => c_x, i_y => c_y, i_z => c_zi,
              o_done => c_done, o_m => c_m, o_z => c_zo,
              o_c => c_c, o_s => c_s);

  p_main : process
    file     fv    : text open read_mode is G_VEC_FILE;
    file     fr    : text open read_mode is G_ROT_FILE;
    variable bx, by, bm, bz : big_t;
    variable wz, wc, ws : integer;
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

    -- ------------------------------------------------ vectoring
    loop
      read_big(fv, bx, ok);  exit when not ok;
      read_big(fv, by, ok);  exit when not ok;
      read_big(fv, bm, ok);  exit when not ok;
      read_big(fv, bz, ok);  exit when not ok;

      c_mode  <= '0';
      c_x     <= std_logic_vector(resize(bx, W));
      c_y     <= std_logic_vector(resize(by, W));
      c_start <= '1';
      wait until rising_edge(clk);
      c_start <= '0';
      loop
        wait until rising_edge(clk);
        exit when c_done = '1';
      end loop;

      check_big("tb_asp_cordic.vec.m", n, resize(signed(c_m), 64), bm, err, shown);
      check_big("tb_asp_cordic.vec.z", n, resize(signed(c_zo), 64), bz, err, shown);
      n := n + 1;
    end loop;

    report "tb_asp_cordic: " & integer'image(n) &
           " vectoring cases checked" severity note;

    -- ------------------------------------------------ rotation
    n := 0;
    loop
      read_int(fr, wz, ok);  exit when not ok;
      read_int(fr, wc, ok);  exit when not ok;
      read_int(fr, ws, ok);  exit when not ok;

      c_mode  <= '1';
      c_zi    <= std_logic_vector(to_signed(wz, ZW));
      c_start <= '1';
      wait until rising_edge(clk);
      c_start <= '0';
      loop
        wait until rising_edge(clk);
        exit when c_done = '1';
      end loop;

      check_int("tb_asp_cordic.rot.c", n, to_integer(signed(c_c)), wc, err, shown);
      check_int("tb_asp_cordic.rot.s", n, to_integer(signed(c_s)), ws, err, shown);
      n := n + 1;
    end loop;

    summarise("tb_asp_cordic", n, err);
    running <= false;
    wait;
  end process p_main;

end architecture sim;
