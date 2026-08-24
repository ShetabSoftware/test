-- =====================================================================
--  asp_detect  -  STAGE 7.  Eigenvalue spoofing detector.
--
--  MODEL REFERENCE : stage7_detect() in matlab/golden/asp_golden_model.m
--
--  FUNCTION
--    Statistic  lambda_1 / mean(lambda_2..lambda_N)  against a threshold,
--    implemented WITHOUT A DIVIDER by cross-multiplying:
--        lam1/mean > NUM/DEN   <=>   lam1*(N-1)*DEN > NUM*sum_tail
--
--  THIS IS THE BLOCK THE PUBLISHED METHOD DOES NOT HAVE
--    Without it the array nulls unconditionally, and with no spoofer
--    present it steers a null into the strongest authentic satellite -
--    which is its state for essentially all of its operating hours.
--    Everything else in this design improves a number; this block is the
--    difference between a system that helps and one that harms.
--
--  THE THRESHOLD IS 5/4 AND THAT IS NOT AN ACCIDENT
--    NUM/DEN = 1280/1024 = 5/4 exactly, so BOTH multiplies degenerate:
--        lam1*3*4 > 5*sum_tail
--    and 4x is (x sll 2) while 5x is (x sll 2) + x.  No DSP48, no
--    divider, one comparator.  The value itself was calibrated against
--    the REALISTIC null hypothesis - thermal noise plus the sample
--    correlation the shaping FIR imposes plus the authentic
--    constellation - and measured 0/320 false alarms over 80
--    constellations.  Calibrating on white noise instead puts the
--    threshold below the H0 median and gives a ~60% false alarm rate on
--    a clean sky, which is invisible because a false alarm produces no
--    symptom at all.
--
--    The general constant multiply is written below rather than the
--    hard-coded shifts, so the block stays correct if the threshold is
--    ever re-calibrated; with the present constants the synthesiser
--    reduces it to exactly those shifts.  An assertion pins the claim.
--
--  RANK
--    Rank 2 is gated on the second eigenvalue being genuinely
--    RESOLVABLE, never on prior knowledge that a second source exists.
--    A second arrival only 6 dB down lifts lambda_2 by ~0.06 above the
--    noise floor in a 1 ms dwell, which yields an eigenvector with tens
--    of degrees of error - and nulling that direction measures WORSE
--    than not nulling it.  RANK2_ENABLE lets the PS withhold rank 2
--    entirely while its MDL test and hysteresis run at 1 kHz.
--
--  OUTPUTS
--    lhs and rhs are exported at the model's scaling (both multiplied by
--    DEN) so the RTL and the model dumps compare directly.  They are
--    also what the PS reads to run its own hysteresis policy without
--    needing a divider either.
--
--  LATENCY   2 clk from i_start.  Purely combinational arithmetic behind
--            one pipeline register; it is never the critical path.
--  RESOURCES 0 DSP48, ~200 FF, ~250 LUT.
-- =====================================================================
library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

library work;
use work.asp_pkg.all;
use work.asp_coef_pkg.all;

entity asp_detect is
  generic (
    G_N : natural := 4
  );
  port (
    clk     : in  std_logic;
    rst     : in  std_logic;

    i_start : in  std_logic;
    i_lam   : in  std_logic_vector(G_N*W_EVD-1 downto 0);  -- descending

    -- runtime policy, from the AXI4-Lite register file
    i_rank2_en : in  std_logic;

    o_valid : out std_logic;
    o_lhs   : out std_logic_vector(47 downto 0);
    o_rhs   : out std_logic_vector(47 downto 0);
    o_det   : out std_logic;
    o_rank  : out std_logic_vector(1 downto 0);
    o_rank2_eligible : out std_logic
  );
end entity asp_detect;


architecture rtl of asp_detect is

  constant AW : natural := 48;

  signal lhs_r, rhs_r   : signed(AW-1 downto 0) := (others => '0');
  signal lhs2_r, rhs2_r : signed(AW-1 downto 0) := (others => '0');
  signal det_r  : std_logic := '0';
  signal el2_r  : std_logic := '0';
  signal rank_r : unsigned(1 downto 0) := (others => '0');
  signal val_r  : std_logic := '0';
  signal stage1 : std_logic := '0';

begin

  assert DET_NUM = 1280 and DET_DEN = 1024
    report "asp_detect: threshold is no longer 5/4; the multiplier-free "
         & "claim in the header no longer holds (the logic is still correct)"
    severity note;
  assert G_N = 4
    report "asp_detect: the rank-2 test assumes N = 4" severity failure;

  p_calc : process (clk)
    variable lam  : signed(W_EVD-1 downto 0);
    variable tail : signed(W_EVD+2 downto 0);
    variable tal2 : signed(W_EVD+2 downto 0);
  begin
    if rising_edge(clk) then
      if rst = '1' then
        val_r  <= '0';
        stage1 <= '0';
        det_r  <= '0';
        rank_r <= (others => '0');
      else
        stage1 <= i_start;
        val_r  <= stage1;

        if i_start = '1' then
          -- tails
          tail := (others => '0');
          for i in 1 to G_N-1 loop
            tail := tail + resize(signed(i_lam((i+1)*W_EVD-1 downto i*W_EVD)),
                                  tail'length);
          end loop;
          tal2 := (others => '0');
          for i in 2 to G_N-1 loop
            tal2 := tal2 + resize(signed(i_lam((i+1)*W_EVD-1 downto i*W_EVD)),
                                  tal2'length);
          end loop;

          -- Products are formed at their natural width (32x13 and 35x12)
          -- and then resized, not formed at 48x13; the synthesiser is
          -- otherwise entitled to build a 61-bit multiplier and throw
          -- most of it away.
          lam := signed(i_lam(W_EVD-1 downto 0));                 -- lambda_1
          lhs_r <= resize(lam  * to_signed((G_N-1)*DET_DEN, 13), AW);
          rhs_r <= resize(tail * to_signed(DET_NUM, 12), AW);

          lam := signed(i_lam(2*W_EVD-1 downto W_EVD));           -- lambda_2
          lhs2_r <= resize(lam  * to_signed((G_N-2)*DET_DEN, 13), AW);
          rhs2_r <= resize(tal2 * to_signed(DET_NUM2, 12), AW);
        end if;

        if stage1 = '1' then
          if lhs_r > rhs_r then
            det_r <= '1';
            if lhs2_r > rhs2_r then
              el2_r <= '1';
              if i_rank2_en = '1' then
                rank_r <= to_unsigned(2, 2);
              else
                rank_r <= to_unsigned(1, 2);
              end if;
            else
              el2_r  <= '0';
              rank_r <= to_unsigned(1, 2);
            end if;
          else
            det_r  <= '0';
            el2_r  <= '0';
            rank_r <= (others => '0');
          end if;
        end if;
      end if;
    end if;
  end process p_calc;

  o_valid <= val_r;
  o_lhs   <= std_logic_vector(lhs_r);
  o_rhs   <= std_logic_vector(rhs_r);
  o_det   <= det_r;
  o_rank  <= std_logic_vector(rank_r);
  o_rank2_eligible <= el2_r;

end architecture rtl;
