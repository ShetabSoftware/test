-- =====================================================================
--  tb_asp_beamformer  -  co-simulation of STAGE 9 against the golden model.
--
--  Drives the model's stage-3 output (s09_beamform_in.txt) with the
--  CAUSAL weight sequence the model uses: dwell 1 runs on the quiescent
--  beam, dwell 2 on the weights estimated from dwell 1, and so on.
--
--  Reproducing that sequence here rather than holding one weight vector
--  is the point of the test.  The causality rule - dwell k's estimate is
--  applied to dwell k+1 - is the difference between a realisable
--  streaming design and the non-causal "estimate from the whole record,
--  then apply it to that same record" that the original MATLAB did, and
--  which flatters the null depth by several dB.
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

entity tb_asp_beamformer is
  generic (
    G_IN_FILE  : string := "gold/s09_beamform_in.txt";
    G_W_FILE   : string := "gold/s08_weights.txt";
    G_OUT_FILE : string := "gold/s09_beamform_out.txt"
  );
end entity tb_asp_beamformer;


architecture sim of tb_asp_beamformer is

  constant CLK_PER : time := 7.637 ns;
  constant NCH     : natural := 4;

  signal clk     : std_logic := '0';
  signal rst     : std_logic := '1';
  signal running : boolean   := true;

  signal i_w     : std_logic_vector(2*NCH*W_WGT-1 downto 0) := (others => '0');
  signal s_valid : std_logic := '0';
  signal s_chan  : std_logic_vector(1 downto 0) := (others => '0');
  signal s_re    : std_logic_vector(W_DAT-1 downto 0) := (others => '0');
  signal s_im    : std_logic_vector(W_DAT-1 downto 0) := (others => '0');
  signal m_valid : std_logic;
  signal m_re    : std_logic_vector(W_BEAM-1 downto 0);
  signal m_im    : std_logic_vector(W_BEAM-1 downto 0);

  signal stim_done : boolean := false;

begin

  clk <= not clk after CLK_PER/2 when running else '0';

  dut : entity work.asp_beamformer
    generic map (G_NCH => NCH)
    port map (clk => clk, rst => rst, i_w => i_w,
              s_valid => s_valid, s_chan => s_chan,
              s_re => s_re, s_im => s_im,
              m_valid => m_valid, m_re => m_re, m_im => m_im);

  p_stim : process
    file     fin : text open read_mode is G_IN_FILE;
    file     fw  : text open read_mode is G_W_FILE;
    variable v   : integer;
    variable ok  : boolean;
    variable n   : integer := 0;
  begin
    rst <= '1';
    -- dwell 1 runs on the quiescent beam
    for k in 0 to NCH-1 loop
      i_w((2*k+1)*W_WGT-1 downto (2*k)*W_WGT) <=
        std_logic_vector(to_signed(WGT_QUIESCENT, W_WGT));
      i_w((2*k+2)*W_WGT-1 downto (2*k+1)*W_WGT) <=
        std_logic_vector(to_signed(0, W_WGT));
    end loop;
    for i in 0 to 7 loop
      wait until rising_edge(clk);
    end loop;
    rst <= '0';

    loop
      read_int(fin, v, ok);
      exit when not ok;
      s_re <= std_logic_vector(to_signed(v, W_DAT));
      read_int(fin, v, ok);
      exit when not ok;
      s_im    <= std_logic_vector(to_signed(v, W_DAT));
      s_chan  <= std_logic_vector(to_unsigned(n mod NCH, 2));
      s_valid <= '1';
      wait until rising_edge(clk);
      s_valid <= '0';
      wait until rising_edge(clk);
      n := n + 1;

      -- at every dwell boundary, load the weights the model estimated
      -- from the dwell that just finished
      if (n mod (NCH*K_DWELL)) = 0 then
        for k in 0 to NCH-1 loop
          read_int(fw, v, ok);
          exit when not ok;
          i_w((2*k+1)*W_WGT-1 downto (2*k)*W_WGT) <=
            std_logic_vector(to_signed(v, W_WGT));
          read_int(fw, v, ok);
          exit when not ok;
          i_w((2*k+2)*W_WGT-1 downto (2*k+1)*W_WGT) <=
            std_logic_vector(to_signed(v, W_WGT));
        end loop;
      end if;
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
  begin
    loop
      wait until rising_edge(clk);
      exit when stim_done;
      if m_valid = '1' then
        read_int(fexp, want, ok);
        exit when not ok;
        check_int("tb_asp_beamformer.re", idx,
                  to_integer(signed(m_re)), want, err, shown);
        idx := idx + 1;
        read_int(fexp, want, ok);
        exit when not ok;
        check_int("tb_asp_beamformer.im", idx,
                  to_integer(signed(m_im)), want, err, shown);
        idx := idx + 1;
      end if;
    end loop;

    assert idx > 100000
      report "tb_asp_beamformer: only " & integer'image(idx) & " values compared"
      severity failure;
    summarise("tb_asp_beamformer", idx, err);
    running <= false;
    wait;
  end process p_check;

end architecture sim;
