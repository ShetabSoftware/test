-- =====================================================================
--  tb_asp_ddc_mixer  -  co-simulation of STAGE 1 against the golden model.
--
--  Drives the exact ADC samples the MATLAB model consumed
--  (s01_ddc_in.txt) and requires the RTL output to equal the model
--  output (s01_ddc_out.txt) to the LAST BIT.  No tolerance is allowed,
--  because both sides are integers: any difference is a bug.
--
--  What this catches that a "looks about right" test would not:
--    * a wrong rounding mode           - fails on ~50% of ties only
--    * NCO conjugation the wrong way   - passes at phase 0 and 8
--    * a sign slip in the imaginary
--      product                         - passes for real-only input
--    * channel/slot permutation        - passes if all channels equal
--  Each of those is invisible to an eyeball check of a spectrum plot.
--
--  Run:  see rtl/sim/run_sim.sh
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

entity tb_asp_ddc_mixer is
  generic (
    G_IN_FILE  : string  := "gold/s01_ddc_in.txt";
    G_OUT_FILE : string  := "gold/s01_ddc_out.txt";
    G_MAX_SAMP : integer := 200000            -- safety stop
  );
end entity tb_asp_ddc_mixer;


architecture sim of tb_asp_ddc_mixer is

  constant NCH     : natural := 4;
  constant CLK_PER : time    := 7.637 ns;     -- 130.944 MHz

  signal clk     : std_logic := '0';
  signal rst     : std_logic := '1';
  signal s_valid : std_logic := '0';
  signal s_re    : std_logic_vector(NCH*W_ADC-1 downto 0) := (others => '0');
  signal s_im    : std_logic_vector(NCH*W_ADC-1 downto 0) := (others => '0');

  signal m_valid : std_logic;
  signal m_chan  : std_logic_vector(1 downto 0);
  signal m_re    : std_logic_vector(W_MIX-1 downto 0);
  signal m_im    : std_logic_vector(W_MIX-1 downto 0);

  signal stim_done : boolean := false;
  signal running   : boolean := true;
  signal n_in      : integer := 0;

begin

  -- Gated clock generator.  VHDL-93 has no std.env.stop, so the
  -- simulation is ended by stopping the only free-running process.
  clk <= not clk after CLK_PER/2 when running else '0';

  dut : entity work.asp_ddc_mixer
    generic map (G_NCH => NCH)
    port map (
      clk => clk, rst => rst,
      s_valid => s_valid, s_re => s_re, s_im => s_im,
      m_valid => m_valid, m_chan => m_chan, m_re => m_re, m_im => m_im
    );

  -- -------------------------------------------------------------------
  -- Stimulus: one ADC sample every NCH clocks, which is the exact rate
  -- ratio the hardware sees (clk = 4 x FS_ADC from the same MMCM).
  -- -------------------------------------------------------------------
  p_stim : process
    file     fin  : text open read_mode is G_IN_FILE;
    variable v    : integer;
    variable ok   : boolean;
    variable cnt  : integer := 0;
  begin
    rst <= '1';
    for i in 0 to 7 loop
      wait until rising_edge(clk);
    end loop;
    rst <= '0';
    wait until rising_edge(clk);

    loop
      exit when cnt >= G_MAX_SAMP;
      -- assemble one ADC sample: ch0.I ch0.Q ch1.I ch1.Q ... ch3.Q
      for k in 0 to NCH-1 loop
        read_int(fin, v, ok);
        exit when not ok;
        s_re((k+1)*W_ADC-1 downto k*W_ADC) <=
          std_logic_vector(to_signed(v, W_ADC));
        read_int(fin, v, ok);
        exit when not ok;
        s_im((k+1)*W_ADC-1 downto k*W_ADC) <=
          std_logic_vector(to_signed(v, W_ADC));
      end loop;
      exit when not ok;

      s_valid <= '1';
      wait until rising_edge(clk);
      s_valid <= '0';
      for i in 1 to NCH-1 loop
        wait until rising_edge(clk);
      end loop;
      cnt := cnt + 1;
      n_in <= cnt;
    end loop;

    -- let the pipeline drain
    for i in 0 to 63 loop
      wait until rising_edge(clk);
    end loop;
    stim_done <= true;
    wait;
  end process p_stim;

  -- -------------------------------------------------------------------
  -- Checker: every accepted output word, in TDM order, against the model.
  -- -------------------------------------------------------------------
  p_check : process
    file     fexp   : text open read_mode is G_OUT_FILE;
    variable want   : integer;
    variable ok     : boolean;
    variable err    : integer := 0;
    variable shown  : integer := 0;
    variable idx    : integer := 0;
    variable expect_chan : integer := 0;
  begin
    loop
      wait until rising_edge(clk);
      exit when stim_done;
      if m_valid = '1' then
        -- channel order must be 0,1,2,3 repeating: a permutation here
        -- silently transposes the array geometry downstream.
        assert to_integer(unsigned(m_chan)) = expect_chan
          report "tb_asp_ddc_mixer: TDM channel order broken at index " &
                 integer'image(idx)
          severity failure;
        expect_chan := (expect_chan + 1) mod NCH;

        read_int(fexp, want, ok);
        exit when not ok;
        check_int("tb_asp_ddc_mixer.re", idx,
                  to_integer(signed(m_re)), want, err, shown);
        idx := idx + 1;

        read_int(fexp, want, ok);
        exit when not ok;
        check_int("tb_asp_ddc_mixer.im", idx,
                  to_integer(signed(m_im)), want, err, shown);
        idx := idx + 1;
      end if;
    end loop;

    summarise("tb_asp_ddc_mixer", idx, err);
    running <= false;
    wait;
  end process p_check;

end architecture sim;
