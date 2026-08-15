-- =====================================================================
--  tb_asp_fir_shape  -  co-simulation of STAGE 3 against the golden model.
--
--  Feeds s03_fir_in.txt in the same bursty pattern stage 2 produces
--  (four channels back to back, then four idle clocks) so that the input
--  FIFO and the two-cycle engine schedule are exercised the way they
--  will be in hardware, not with a convenient steady stream.
--
--  The RTL computes the CAUSAL convolution and discards 31 warm-up
--  outputs per channel; the model computes the CENTRED one.  After that
--  discard the two are sample-aligned and must agree bit for bit.  The
--  final 31 model outputs per channel depend on input past the end of
--  the file (the model zero-pads, the RTL has not been given those
--  samples yet), so the comparison stops when the RTL runs out - the
--  checker reports how many were compared so a silent early exit is
--  visible rather than looking like a pass.
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

entity tb_asp_fir_shape is
  generic (
    G_IN_FILE  : string  := "gold/s03_fir_in.txt";
    G_OUT_FILE : string  := "gold/s03_fir_out.txt";
    G_MAX_SAMP : integer := 400000
  );
end entity tb_asp_fir_shape;


architecture sim of tb_asp_fir_shape is

  constant NCH     : natural := 4;
  constant CLK_PER : time    := 7.637 ns;

  signal clk     : std_logic := '0';
  signal rst     : std_logic := '1';
  signal running : boolean   := true;

  signal s_valid : std_logic := '0';
  signal s_chan  : std_logic_vector(1 downto 0) := (others => '0');
  signal s_re    : std_logic_vector(W_HB-1 downto 0) := (others => '0');
  signal s_im    : std_logic_vector(W_HB-1 downto 0) := (others => '0');

  signal m_valid : std_logic;
  signal m_chan  : std_logic_vector(1 downto 0);
  signal m_re    : std_logic_vector(W_DAT-1 downto 0);
  signal m_im    : std_logic_vector(W_DAT-1 downto 0);
  signal m_ovf   : std_logic;

  signal stim_done : boolean := false;

begin

  clk <= not clk after CLK_PER/2 when running else '0';

  dut : entity work.asp_fir_shape
    generic map (G_NCH => NCH)
    port map (
      clk => clk, rst => rst,
      s_valid => s_valid, s_chan => s_chan, s_re => s_re, s_im => s_im,
      m_valid => m_valid, m_chan => m_chan, m_re => m_re, m_im => m_im,
      o_overflow => m_ovf
    );

  -- Bursty stimulus: 4 channels back to back, then 4 idle clocks, which
  -- is exactly what a decimate-by-two stage upstream produces.
  p_stim : process
    file     fin : text open read_mode is G_IN_FILE;
    variable v   : integer;
    variable ok  : boolean;
    variable cnt : integer := 0;
  begin
    rst <= '1';
    for i in 0 to 7 loop
      wait until rising_edge(clk);
    end loop;
    rst <= '0';

    loop
      exit when cnt >= G_MAX_SAMP;
      for ch in 0 to NCH-1 loop
        read_int(fin, v, ok);
        exit when not ok;
        s_re <= std_logic_vector(to_signed(v, W_HB));
        read_int(fin, v, ok);
        exit when not ok;
        s_im <= std_logic_vector(to_signed(v, W_HB));
        s_chan  <= std_logic_vector(to_unsigned(ch, 2));
        s_valid <= '1';
        wait until rising_edge(clk);
        cnt := cnt + 1;
      end loop;
      s_valid <= '0';
      exit when not ok;
      for i in 0 to NCH-1 loop           -- the decimator's idle half
        wait until rising_edge(clk);
      end loop;
    end loop;
    s_valid <= '0';

    for i in 0 to 255 loop
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
          report "tb_asp_fir_shape: TDM channel order broken at index " &
                 integer'image(idx)
          severity failure;
        expch := (expch + 1) mod NCH;

        read_int(fexp, want, ok);
        exit when not ok;
        check_int("tb_asp_fir_shape.re", idx,
                  to_integer(signed(m_re)), want, err, shown);
        idx := idx + 1;

        read_int(fexp, want, ok);
        exit when not ok;
        check_int("tb_asp_fir_shape.im", idx,
                  to_integer(signed(m_im)), want, err, shown);
        idx := idx + 1;
      end if;
    end loop;

    assert m_ovf = '0'
      report "tb_asp_fir_shape: input FIFO overflowed - the rate plan is wrong"
      severity failure;
    assert idx > 100000
      report "tb_asp_fir_shape: only " & integer'image(idx) &
             " values compared; the stimulus ended early"
      severity failure;

    summarise("tb_asp_fir_shape", idx, err);
    running <= false;
    wait;
  end process p_check;

end architecture sim;
