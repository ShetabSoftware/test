-- =====================================================================
--  tb_asp_weight_calc  -  co-simulation of STAGE 8 against the golden model.
--
--  Drives the model's eigenvector matrix (s06_evd_U.txt), the square
--  roots of the covariance diagonal (s05_whiten_dsq.txt) and the rank
--  from the detector dump, then checks the weights against
--  s08_weights.txt.
--
--  These are the numbers that actually steer the array, so the check is
--  exact on all four complex weights.  A weight vector that is merely
--  close produces a null that is merely deep-ish, and the difference
--  between -24 dB and -18 dB of suppression is the difference between
--  the system working and not.
--
--  The rank-0 path is exercised too, because it is the NORMAL state: for
--  almost all operating hours there is no spoofer, the detector does not
--  fire, and this block must hand back the quiescent beam completely
--  untouched.  A rank-0 path that quietly returns a scaled or rotated
--  version of h would degrade the receiver every second it is not under
--  attack.
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

entity tb_asp_weight_calc is
  generic (
    G_U_FILE   : string := "gold/s06_evd_U.txt";
    G_DSQ_FILE : string := "gold/s05_whiten_dsq.txt";
    G_DET_FILE : string := "gold/s07_detect.txt";
    G_W_FILE   : string := "gold/s08_weights.txt"
  );
end entity tb_asp_weight_calc;


architecture sim of tb_asp_weight_calc is

  constant CLK_PER : time := 7.637 ns;
  constant NA      : natural := 4;
  constant DSW     : natural := 24;

  signal clk     : std_logic := '0';
  signal rst     : std_logic := '1';
  signal running : boolean   := true;

  signal s_valid : std_logic := '0';
  signal s_idx   : std_logic_vector(3 downto 0) := (others => '0');
  signal s_re    : std_logic_vector(W_UVEC-1 downto 0) := (others => '0');
  signal s_im    : std_logic_vector(W_UVEC-1 downto 0) := (others => '0');

  signal i_start : std_logic := '0';
  signal i_rank  : std_logic_vector(1 downto 0) := (others => '0');
  signal i_dsq   : std_logic_vector(NA*DSW-1 downto 0) := (others => '0');
  signal o_w     : std_logic_vector(2*NA*W_WGT-1 downto 0);
  signal o_done  : std_logic;

begin

  clk <= not clk after CLK_PER/2 when running else '0';

  dut : entity work.asp_weight_calc
    generic map (G_N => NA, G_DSW => DSW)
    port map (clk => clk, rst => rst,
              s_valid => s_valid, s_idx => s_idx, s_re => s_re, s_im => s_im,
              i_start => i_start, i_rank => i_rank, i_dsq => i_dsq,
              o_w => o_w, o_done => o_done);

  p_main : process
    file     fu   : text open read_mode is G_U_FILE;
    file     fd   : text open read_mode is G_DSQ_FILE;
    file     fdet : text open read_mode is G_DET_FILE;
    file     fw   : text open read_mode is G_W_FILE;
    variable v    : integer;
    variable dummy: big_t;
    variable want : integer;
    variable ok   : boolean;
    variable err  : integer := 0;
    variable shown: integer := 0;
    variable cnt  : integer := 0;
    variable ndw  : integer := 0;
    variable rk   : integer;
  begin
    rst <= '1';
    for i in 0 to 7 loop
      wait until rising_edge(clk);
    end loop;
    rst <= '0';
    wait until rising_edge(clk);

    loop
      -- ---- eigenvector matrix
      for e in 0 to 15 loop
        read_int(fu, v, ok);
        exit when not ok;
        s_re <= std_logic_vector(to_signed(v, W_UVEC));
        read_int(fu, v, ok);
        exit when not ok;
        s_im    <= std_logic_vector(to_signed(v, W_UVEC));
        s_idx   <= std_logic_vector(to_unsigned(e, 4));
        s_valid <= '1';
        wait until rising_edge(clk);
        s_valid <= '0';
      end loop;
      exit when not ok;

      -- ---- sqrt of the covariance diagonal
      for k in 0 to NA-1 loop
        read_int(fd, v, ok);
        exit when not ok;
        i_dsq((k+1)*DSW-1 downto k*DSW) <=
          std_logic_vector(to_unsigned(v, DSW));
      end loop;
      exit when not ok;

      -- ---- rank: the 4th field of the detector record
      read_big(fdet, dummy, ok);  exit when not ok;   -- lhs
      read_big(fdet, dummy, ok);  exit when not ok;   -- rhs
      read_int(fdet, v, ok);      exit when not ok;   -- detected
      read_int(fdet, rk, ok);     exit when not ok;   -- rank
      i_rank <= std_logic_vector(to_unsigned(rk, 2));

      wait until rising_edge(clk);
      i_start <= '1';
      wait until rising_edge(clk);
      i_start <= '0';

      loop
        wait until rising_edge(clk);
        exit when o_done = '1';
      end loop;

      for k in 0 to NA-1 loop
        read_int(fw, want, ok);
        exit when not ok;
        check_int("tb_asp_weight_calc.re", cnt,
                  to_integer(signed(o_w((2*k+1)*W_WGT-1 downto (2*k)*W_WGT))),
                  want, err, shown);
        cnt := cnt + 1;
        read_int(fw, want, ok);
        exit when not ok;
        check_int("tb_asp_weight_calc.im", cnt,
                  to_integer(signed(o_w((2*k+2)*W_WGT-1 downto (2*k+1)*W_WGT))),
                  want, err, shown);
        cnt := cnt + 1;
      end loop;
      exit when not ok;

      ndw := ndw + 1;
    end loop;

    -- ---- rank 0 must return the quiescent beam untouched
    i_rank  <= "00";
    wait until rising_edge(clk);
    i_start <= '1';
    wait until rising_edge(clk);
    i_start <= '0';
    loop
      wait until rising_edge(clk);
      exit when o_done = '1';
    end loop;
    for k in 0 to NA-1 loop
      assert to_integer(signed(o_w((2*k+1)*W_WGT-1 downto (2*k)*W_WGT)))
             = WGT_QUIESCENT
        report "tb_asp_weight_calc: rank 0 did not return the quiescent beam "
             & "on element " & integer'image(k)
        severity failure;
      assert to_integer(signed(o_w((2*k+2)*W_WGT-1 downto (2*k+1)*W_WGT))) = 0
        report "tb_asp_weight_calc: rank 0 produced a non-zero imaginary part"
        severity failure;
    end loop;

    assert ndw >= 2
      report "tb_asp_weight_calc: only " & integer'image(ndw) & " dwells checked"
      severity failure;

    summarise("tb_asp_weight_calc", cnt, err);
    report "tb_asp_weight_calc: rank 0 returns the quiescent beam exactly"
      severity note;
    running <= false;
    wait;
  end process p_main;

end architecture sim;
