-- =====================================================================
--  tb_asp_cov_accum  -  co-simulation of STAGE 4 against the golden model.
--
--  Drives the model's own stage-3 output (s03_fir_out.txt) so the
--  covariance is compared in isolation, then checks all 16 entries of
--  every dwell against s04_cov_out.txt.
--
--  This is the testbench that has to be trusted most, because the error
--  it guards against is invisible everywhere else.  A rounding bias in
--  the accumulator adds a constant to every entry of R; ones(N) is
--  rank one with the boresight steering vector as its eigenvector, so
--  the estimator would invent a source at zenith and null the
--  satellites.  The bias does not shrink with dwell length, so no amount
--  of averaging exposes it, and a floating-point model never shows it.
--  Only an exact 48-bit integer comparison does.
--
--  It also checks the Hermitian structure explicitly: the lower triangle
--  must be the exact conjugate of the upper, and the diagonal must be
--  exactly real.  Both are structural properties that a transposed index
--  map would break while still producing plausible magnitudes.
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

entity tb_asp_cov_accum is
  generic (
    G_IN_FILE  : string  := "gold/s03_fir_out.txt";
    G_OUT_FILE : string  := "gold/s04_cov_out.txt";
    G_MAX_SAMP : integer := 400000
  );
end entity tb_asp_cov_accum;


architecture sim of tb_asp_cov_accum is

  constant NCH     : natural := 4;
  constant CLK_PER : time    := 7.637 ns;

  signal clk     : std_logic := '0';
  signal rst     : std_logic := '1';
  signal running : boolean   := true;

  signal s_valid : std_logic := '0';
  signal s_chan  : std_logic_vector(1 downto 0) := (others => '0');
  signal s_re    : std_logic_vector(W_DAT-1 downto 0) := (others => '0');
  signal s_im    : std_logic_vector(W_DAT-1 downto 0) := (others => '0');

  signal o_valid : std_logic;
  signal o_idx   : std_logic_vector(3 downto 0);
  signal o_re    : std_logic_vector(W_ACC-1 downto 0);
  signal o_im    : std_logic_vector(W_ACC-1 downto 0);
  signal o_last  : std_logic;
  signal o_tick  : std_logic;

  signal stim_done : boolean := false;
  signal n_dwell   : integer := 0;

begin

  clk <= not clk after CLK_PER/2 when running else '0';

  dut : entity work.asp_cov_accum
    generic map (G_NCH => NCH, G_KDWELL => K_DWELL)
    port map (
      clk => clk, rst => rst,
      s_valid => s_valid, s_chan => s_chan, s_re => s_re, s_im => s_im,
      o_valid => o_valid, o_idx => o_idx, o_re => o_re, o_im => o_im,
      o_last => o_last, o_dwell_tick => o_tick
    );

  -- One channel every two clocks, channels cycling 0..3: the exact rate
  -- and phasing stage 3 delivers.
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
      s_re <= std_logic_vector(to_signed(v, W_DAT));
      read_int(fin, v, ok);
      exit when not ok;
      s_im <= std_logic_vector(to_signed(v, W_DAT));
      s_chan  <= std_logic_vector(to_unsigned(ch, 2));
      s_valid <= '1';
      wait until rising_edge(clk);
      s_valid <= '0';
      wait until rising_edge(clk);
      ch  := (ch + 1) mod NCH;
      cnt := cnt + 1;
    end loop;
    s_valid <= '0';

    for i in 0 to 127 loop
      wait until rising_edge(clk);
    end loop;
    stim_done <= true;
    wait;
  end process p_stim;

  p_check : process
    file     fexp  : text open read_mode is G_OUT_FILE;
    variable want  : big_t;
    variable ok    : boolean;
    variable err   : integer := 0;
    variable shown : integer := 0;
    variable idx   : integer := 0;
    variable expidx : integer := 0;
    -- kept for the structural checks
    type acc16_t is array (0 to 15) of big_t;
    variable seen_re, seen_im : acc16_t := (others => (others => '0'));
    variable ndw : integer := 0;
  begin
    loop
      wait until rising_edge(clk);
      exit when stim_done;

      if o_valid = '1' then
        assert to_integer(unsigned(o_idx)) = expidx
          report "tb_asp_cov_accum: matrix streamed out of order"
          severity failure;

        -- The golden file holds the FULL 4x4 in column-major order with
        -- I and Q interleaved, which is what R(:) produces in MATLAB.
        read_big(fexp, want, ok);
        exit when not ok;
        check_big("tb_asp_cov_accum.re", idx,
                  resize(signed(o_re), 64), want, err, shown);
        seen_re(expidx) := resize(signed(o_re), 64);
        idx := idx + 1;

        read_big(fexp, want, ok);
        exit when not ok;
        check_big("tb_asp_cov_accum.im", idx,
                  resize(signed(o_im), 64), want, err, shown);
        seen_im(expidx) := resize(signed(o_im), 64);
        idx := idx + 1;

        if expidx = 15 then
          assert o_last = '1'
            report "tb_asp_cov_accum: o_last missing on the final entry"
            severity failure;
          ndw := ndw + 1;
          n_dwell <= ndw;

          -- Structural checks.  A transposed or mis-mapped index table
          -- can still produce the right magnitudes; these cannot pass
          -- unless the Hermitian structure is actually right.
          for d in 0 to 3 loop
            assert seen_im(4*d + d) = to_signed(0, 64)
              report "tb_asp_cov_accum: diagonal entry " & integer'image(d) &
                     " has non-zero imaginary part"
              severity failure;
            assert seen_re(4*d + d) >= to_signed(0, 64)
              report "tb_asp_cov_accum: diagonal entry " & integer'image(d) &
                     " is negative; it is a sum of squares"
              severity failure;
          end loop;
          for r in 0 to 3 loop
            for c in 0 to 3 loop
              assert seen_re(4*c + r) = seen_re(4*r + c)
                report "tb_asp_cov_accum: real part is not symmetric"
                severity failure;
              assert seen_im(4*c + r) = -seen_im(4*r + c)
                report "tb_asp_cov_accum: imaginary part is not antisymmetric"
                severity failure;
            end loop;
          end loop;
        end if;

        expidx := (expidx + 1) mod 16;
      end if;
    end loop;

    assert ndw >= 2
      report "tb_asp_cov_accum: only " & integer'image(ndw) &
             " dwells completed; expected at least 2"
      severity failure;

    summarise("tb_asp_cov_accum", idx, err);
    report "tb_asp_cov_accum: " & integer'image(ndw) & " dwells checked, " &
           "Hermitian structure verified on every one" severity note;
    running <= false;
    wait;
  end process p_check;

end architecture sim;
