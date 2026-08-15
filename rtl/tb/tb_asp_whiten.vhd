-- =====================================================================
--  tb_asp_whiten  -  co-simulation of STAGE 5 against the golden model.
--
--  Drives the model's own covariance dump (s04_cov_out.txt) and checks
--  the whitened matrix against s05_whiten_out.txt, dwell by dwell.
--
--  Two structural properties are asserted on top of the bit comparison,
--  because both are relied on by the block downstream and neither is
--  visible from a single entry:
--
--    * the whitened diagonal must sit at 2^F_EVD to within a few parts
--      in 10^5.  It is NOT exactly 2^F_EVD and cannot be: Y carries 17
--      bits, so Y^2 * Rs_ii reproduces unity only to about 1e-6 per
--      factor, and the measured spread over the reference dwells is
--      -1250 .. +1912 LSB, i.e. 2.9e-5 relative.  The bound below is
--      2^12, comfortably above that and far below the percent-level
--      error any real whitening bug would produce.  What the check is
--      really protecting is trace(Rw) ~ N*2^F_EVD, which is the bound
--      the Jacobi engine's word-length proof rests on.
--    * the imaginary diagonal must be exactly zero, or the matrix handed
--      to the Jacobi engine is not Hermitian and its stability argument
--      no longer holds.
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

entity tb_asp_whiten is
  generic (
    G_IN_FILE  : string := "gold/s04_cov_out.txt";
    G_OUT_FILE : string := "gold/s05_whiten_out.txt"
  );
end entity tb_asp_whiten;


architecture sim of tb_asp_whiten is

  constant CLK_PER : time := 7.637 ns;
  constant DSW     : natural := 24;

  signal clk     : std_logic := '0';
  signal rst     : std_logic := '1';
  signal running : boolean   := true;

  signal s_valid : std_logic := '0';
  signal s_idx   : std_logic_vector(3 downto 0) := (others => '0');
  signal s_re    : std_logic_vector(W_ACC-1 downto 0) := (others => '0');
  signal s_im    : std_logic_vector(W_ACC-1 downto 0) := (others => '0');
  signal s_last  : std_logic := '0';

  signal o_valid : std_logic;
  signal o_idx   : std_logic_vector(3 downto 0);
  signal o_re    : std_logic_vector(W_EVD-1 downto 0);
  signal o_im    : std_logic_vector(W_EVD-1 downto 0);
  signal o_done  : std_logic;
  signal o_dsq   : std_logic_vector(4*DSW-1 downto 0);

  signal fed_all : boolean := false;

begin

  clk <= not clk after CLK_PER/2 when running else '0';

  dut : entity work.asp_whiten
    generic map (G_N => 4, G_DSW => DSW)
    port map (clk => clk, rst => rst,
              s_valid => s_valid, s_idx => s_idx, s_re => s_re,
              s_im => s_im, s_last => s_last,
              o_valid => o_valid, o_idx => o_idx, o_re => o_re,
              o_im => o_im, o_done => o_done, o_dsq => o_dsq);

  p_main : process
    file     fin   : text open read_mode is G_IN_FILE;
    file     fexp  : text open read_mode is G_OUT_FILE;
    variable br, bi: big_t;
    variable wr, wi: integer;
    variable ok    : boolean;
    variable err   : integer := 0;
    variable shown : integer := 0;
    variable n     : integer := 0;
    variable ndw   : integer := 0;
    variable got_re, got_im : integer;
  begin
    rst <= '1';
    for i in 0 to 7 loop
      wait until rising_edge(clk);
    end loop;
    rst <= '0';
    wait until rising_edge(clk);

    loop
      -- ---- feed one dwell's covariance matrix (16 entries)
      for e in 0 to 15 loop
        read_big(fin, br, ok);
        exit when not ok;
        read_big(fin, bi, ok);
        exit when not ok;
        s_idx   <= std_logic_vector(to_unsigned(e, 4));
        s_re    <= std_logic_vector(resize(br, W_ACC));
        s_im    <= std_logic_vector(resize(bi, W_ACC));
        s_valid <= '1';
        if e = 15 then
          s_last <= '1';
        end if;
        wait until rising_edge(clk);
        s_valid <= '0';
        s_last  <= '0';
      end loop;
      exit when not ok;

      -- ---- collect the whitened matrix
      for e in 0 to 15 loop
        loop
          wait until rising_edge(clk);
          exit when o_valid = '1';
        end loop;
        assert to_integer(unsigned(o_idx)) = e
          report "tb_asp_whiten: whitened matrix streamed out of order"
          severity failure;

        got_re := to_integer(signed(o_re));
        got_im := to_integer(signed(o_im));

        read_int(fexp, wr, ok);
        exit when not ok;
        check_int("tb_asp_whiten.re", n, got_re, wr, err, shown);
        n := n + 1;
        read_int(fexp, wi, ok);
        exit when not ok;
        check_int("tb_asp_whiten.im", n, got_im, wi, err, shown);
        n := n + 1;

        if (e mod 4) = (e / 4) then
          assert got_im = 0
            report "tb_asp_whiten: diagonal entry " & integer'image(e) &
                   " has non-zero imaginary part"
            severity failure;
          assert abs(got_re - 2**F_EVD) <= 2**12
            report "tb_asp_whiten: diagonal entry " & integer'image(e) &
                   " is " & integer'image(got_re) & ", which is " &
                   integer'image(got_re - 2**F_EVD) &
                   " LSB from 2^" & integer'image(F_EVD) &
                   "; the reciprocal-sqrt cannot be that far out"
            severity failure;
        end if;
      end loop;
      exit when not ok;

      loop
        wait until rising_edge(clk);
        exit when o_done = '1';
      end loop;
      ndw := ndw + 1;
    end loop;

    assert ndw >= 2
      report "tb_asp_whiten: only " & integer'image(ndw) & " dwells checked"
      severity failure;

    summarise("tb_asp_whiten", n, err);
    report "tb_asp_whiten: " & integer'image(ndw) &
           " dwells, whitened diagonal within 2^12 LSB of 2^" &
           integer'image(F_EVD) & " and exactly real on every one"
      severity note;
    running <= false;
    wait;
  end process p_main;

end architecture sim;
