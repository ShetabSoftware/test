-- =====================================================================
--  asp_isqrt  -  integer square root, fixed 24-step restoring algorithm.
--
--  MODEL REFERENCE : isqrtInt() in matlab/golden/asp_golden_model.m
--
--  FUNCTION
--    s = floor(sqrt(a)) for a >= 0, a < 2^48, s < 2^24.
--
--  ALGORITHM
--    The classic restoring (paper-and-pencil) square root, taken two
--    bits of radicand at a time from the top:
--        rem = 4*rem + next two bits
--        t   = 4*s + 1
--        if rem >= t then rem -= t; s = 2s + 1 else s = 2s
--    24 iterations for a 48-bit radicand.  No divider, no multiplier -
--    the "4*s + 1" is a shift and a set bit, and the comparison is a
--    subtract whose borrow is the answer bit.
--
--  WHY NOT A CORDIC OR A NEWTON ITERATION
--    Both would need either a seed table or gain compensation, and both
--    give an approximation that then has to be proved to match the
--    model's floor().  The restoring algorithm IS floor(sqrt(a)) by
--    construction - it is exact for every input by the structure of the
--    recurrence, not by an error bound - so there is nothing to verify
--    beyond the loop count.
--
--  WHERE IT IS USED
--    d_i = sqrt(R_ii) in stage 5, four times per dwell.  D = diag(d)
--    undoes the whitening in stage 8, mapping the eigenvector back into
--    the measured domain where the beamformer actually operates.  Skip
--    it and the null is steered in the whitened domain, which is a
--    different direction.
--
--  LATENCY   26 clk.  o_done is a single-cycle pulse.
--  RESOURCES 0 DSP48, ~120 FF, ~150 LUT.
-- =====================================================================
library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

entity asp_isqrt is
  generic (
    G_AW : natural := 48                 -- must be even
  );
  port (
    clk     : in  std_logic;
    rst     : in  std_logic;
    i_start : in  std_logic;
    i_a     : in  std_logic_vector(G_AW-1 downto 0);      -- unsigned
    o_done  : out std_logic;
    o_s     : out std_logic_vector(G_AW/2-1 downto 0)     -- unsigned
  );
end entity asp_isqrt;


architecture rtl of asp_isqrt is

  constant NIT : natural := G_AW/2;                        -- 24

  signal rad   : unsigned(G_AW-1 downto 0) := (others => '0');
  signal rem_r : unsigned(G_AW/2+1 downto 0) := (others => '0');
  signal root  : unsigned(G_AW/2-1 downto 0) := (others => '0');
  signal it    : integer range 0 to NIT := 0;
  signal busy  : std_logic := '0';
  signal done_r: std_logic := '0';

begin

  assert (G_AW mod 2) = 0
    report "asp_isqrt: radicand width must be even" severity failure;

  p_fsm : process (clk)
    variable r2 : unsigned(G_AW/2+1 downto 0);
    variable t  : unsigned(G_AW/2+1 downto 0);
  begin
    if rising_edge(clk) then
      if rst = '1' then
        busy   <= '0';
        done_r <= '0';
        root   <= (others => '0');
      else
        done_r <= '0';

        if busy = '0' then
          if i_start = '1' then
            rad   <= unsigned(i_a);
            rem_r <= (others => '0');
            root  <= (others => '0');
            it    <= 0;
            busy  <= '1';
          end if;
        else
          -- shift the top two bits of the radicand into the remainder
          r2 := shift_left(rem_r, 2);
          r2(1 downto 0) := rad(G_AW-1 downto G_AW-2);

          t := shift_left(resize(root, t'length), 2) + 1;

          if r2 >= t then
            rem_r <= r2 - t;
            root  <= shift_left(root, 1) or to_unsigned(1, root'length);
          else
            rem_r <= r2;
            root  <= shift_left(root, 1);
          end if;

          rad <= shift_left(rad, 2);

          if it = NIT-1 then
            busy   <= '0';
            done_r <= '1';
          else
            it <= it + 1;
          end if;
        end if;
      end if;
    end if;
  end process p_fsm;

  o_done <= done_r;
  o_s    <= std_logic_vector(root);

end architecture rtl;
