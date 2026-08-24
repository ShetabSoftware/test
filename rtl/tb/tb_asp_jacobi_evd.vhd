-- =====================================================================
--  tb_asp_jacobi_evd  -  co-simulation of STAGE 6 against the golden model.
--
--  Drives the model's whitened matrix (s05_whiten_out.txt) and checks
--  BOTH outputs: the eigenvalues (s06_evd_lam.txt) and the eigenvector
--  matrix (s06_evd_U.txt), bit for bit, for every dwell.
--
--  Checking the eigenvectors matters more than it looks.  The
--  eigenvalues alone would pass even if the eigenvector accumulation
--  were wrong, because they come off the diagonal of A and never touch
--  U - and U is what the beamformer actually steers with.  A design that
--  detects perfectly and nulls in the wrong direction is the worst
--  possible failure, so U is compared element by element.
--
--  Two invariants are asserted alongside:
--    * the eigenvalues must come out in descending order, since every
--      downstream block indexes them positionally;
--    * their sum must be close to trace(Rw) ~ N*2^F_EVD, because Jacobi
--      rotations are unitary similarity transforms and preserve the
--      trace exactly in exact arithmetic.  A drift here is the signature
--      of a rotation applied on one side only, which is a similarity
--      transform no longer.
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

entity tb_asp_jacobi_evd is
  generic (
    G_IN_FILE  : string := "gold/s05_whiten_out.txt";
    G_LAM_FILE : string := "gold/s06_evd_lam.txt";
    G_U_FILE   : string := "gold/s06_evd_U.txt"
  );
end entity tb_asp_jacobi_evd;


architecture sim of tb_asp_jacobi_evd is

  constant CLK_PER : time := 7.637 ns;
  constant N       : natural := 4;

  signal clk     : std_logic := '0';
  signal rst     : std_logic := '1';
  signal running : boolean   := true;

  signal s_valid : std_logic := '0';
  signal s_idx   : std_logic_vector(3 downto 0) := (others => '0');
  signal s_re    : std_logic_vector(W_EVD-1 downto 0) := (others => '0');
  signal s_im    : std_logic_vector(W_EVD-1 downto 0) := (others => '0');
  signal s_start : std_logic := '0';

  signal o_valid : std_logic;
  signal o_idx   : std_logic_vector(3 downto 0);
  signal o_re    : std_logic_vector(W_UVEC-1 downto 0);
  signal o_im    : std_logic_vector(W_UVEC-1 downto 0);
  signal o_lam   : std_logic_vector(N*W_EVD-1 downto 0);
  signal o_done  : std_logic;

begin

  clk <= not clk after CLK_PER/2 when running else '0';

  dut : entity work.asp_jacobi_evd
    generic map (G_N => N)
    port map (clk => clk, rst => rst,
              s_valid => s_valid, s_idx => s_idx, s_re => s_re,
              s_im => s_im, s_start => s_start,
              o_valid => o_valid, o_idx => o_idx, o_re => o_re,
              o_im => o_im, o_lam => o_lam, o_done => o_done);

  p_main : process
    file     fin  : text open read_mode is G_IN_FILE;
    file     flam : text open read_mode is G_LAM_FILE;
    file     fu   : text open read_mode is G_U_FILE;
    variable vr, vi : integer;
    variable want : integer;
    variable ok   : boolean;
    variable err  : integer := 0;
    variable shown: integer := 0;
    variable nval : integer := 0;
    variable ndw  : integer := 0;
    variable lv   : integer;
    variable prev : integer;
    variable tr   : integer;
  begin
    rst <= '1';
    for i in 0 to 7 loop
      wait until rising_edge(clk);
    end loop;
    rst <= '0';
    wait until rising_edge(clk);

    loop
      -- ---- feed one whitened matrix
      for e in 0 to 15 loop
        read_int(fin, vr, ok);
        exit when not ok;
        read_int(fin, vi, ok);
        exit when not ok;
        s_idx   <= std_logic_vector(to_unsigned(e, 4));
        s_re    <= std_logic_vector(to_signed(vr, W_EVD));
        s_im    <= std_logic_vector(to_signed(vi, W_EVD));
        s_valid <= '1';
        wait until rising_edge(clk);
        s_valid <= '0';
      end loop;
      exit when not ok;

      s_start <= '1';
      wait until rising_edge(clk);
      s_start <= '0';

      -- ---- collect the eigenvector matrix
      for e in 0 to 15 loop
        loop
          wait until rising_edge(clk);
          exit when o_valid = '1';
        end loop;
        assert to_integer(unsigned(o_idx)) = e
          report "tb_asp_jacobi_evd: U streamed out of order" severity failure;

        read_int(fu, want, ok);
        exit when not ok;
        check_int("tb_asp_jacobi_evd.U.re", nval,
                  to_integer(signed(o_re)), want, err, shown);
        nval := nval + 1;
        read_int(fu, want, ok);
        exit when not ok;
        check_int("tb_asp_jacobi_evd.U.im", nval,
                  to_integer(signed(o_im)), want, err, shown);
        nval := nval + 1;
      end loop;
      exit when not ok;

      loop
        wait until rising_edge(clk);
        exit when o_done = '1';
      end loop;

      -- ---- eigenvalues
      prev := integer'high;
      tr   := 0;
      for i in 0 to N-1 loop
        lv := to_integer(signed(o_lam((i+1)*W_EVD-1 downto i*W_EVD)));
        read_int(flam, want, ok);
        exit when not ok;
        check_int("tb_asp_jacobi_evd.lam", nval, lv, want, err, shown);
        nval := nval + 1;

        assert lv <= prev
          report "tb_asp_jacobi_evd: eigenvalues are not descending at " &
                 integer'image(i)
          severity failure;
        prev := lv;
        tr   := tr + lv/1024;         -- scaled to stay inside INTEGER
      end loop;
      exit when not ok;

      -- trace is preserved by unitary similarity transforms
      assert abs(tr - N*(2**F_EVD)/1024) < (2**F_EVD)/1024
        report "tb_asp_jacobi_evd: trace drifted to " & integer'image(tr*1024) &
               ", expected about " & integer'image(N*(2**F_EVD)) &
               "; a rotation is being applied on one side only"
        severity failure;

      ndw := ndw + 1;
    end loop;

    assert ndw >= 2
      report "tb_asp_jacobi_evd: only " & integer'image(ndw) & " dwells checked"
      severity failure;

    summarise("tb_asp_jacobi_evd", nval, err);
    report "tb_asp_jacobi_evd: " & integer'image(ndw) &
           " dwells; eigenvalues descending and trace preserved on every one"
      severity note;
    running <= false;
    wait;
  end process p_main;

end architecture sim;
