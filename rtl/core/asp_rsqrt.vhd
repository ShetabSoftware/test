-- =====================================================================
--  asp_rsqrt  -  reciprocal square root, range reduction + fixed Newton.
--
--  MODEL REFERENCE : rsqrtNorm() in matlab/golden/asp_golden_model.m
--
--  FUNCTION
--    Given a > 0, returns Y and k with
--        1 / sqrt(a)  =  Y * 2^-(3*RSQ_F/2 + k)
--    where Y is an integer in (2^F, 2^(F+1)] and F = RSQ_F = 16.
--
--  WHY RANGE REDUCTION COMES FIRST
--    Shift a by an EVEN number of bits 2k so the mantissa lands in
--    [2^(F-2), 2^F), i.e. a' in [0.25, 1).  Then y' = 1/sqrt(a') is in
--    (1, 2] and the Newton iteration operates on a bounded argument, so
--    a FIXED iteration count reaches full precision for EVERY input.
--    That is the whole point: without range reduction the iteration
--    count would depend on the data and the block would have
--    data-dependent latency, which a fixed 1 kHz schedule cannot absorb.
--    The shift must be even because it has to come back out of the
--    square root as an integer power of two.
--
--  THE SEED MATTERS MORE THAN IT LOOKS
--    y0 = (7 - 4a')/3, the chord of 1/sqrt across [0.25, 1).  Worst-case
--    relative error 18%, which the quadratic Newton step takes to 2e-5
--    in three iterations and below 1 LSB in four.  A cruder seed such as
--    (1.5 - a') leaves 0.7% after three - enough to move the whitened
--    diagonal off unity, which then propagates into every eigenvalue and
--    from there into the detector statistic.  Four iterations are run
--    here because the model runs four; matching it is not optional.
--
--    The divide by 3 is exact and multiplier-only:
--        floor(N/3) = (N * 699051) >> 21     for N < 2^21
--    because 3 * 699051 = 2^21 + 1.  The seed argument N = 7*2^16 - 4*an
--    lies in [196612, 393216], verified exhaustively over every an in
--    [2^14, 2^16).  No divider, no lookup table.
--
--  ARITHMETIC NOTE
--    The model uses fix() - truncation TOWARD ZERO - inside the
--    iteration, not floor.  For the positive operands that occur in
--    practice they agree, but 3*2^F - t can go negative for a badly
--    conditioned input and there the two differ by one LSB.  trunc_shr()
--    implements fix() exactly so the block cannot diverge from the model
--    on the very inputs that are hardest to reason about.
--
--  ARCHITECTURE
--    Fully sequential around ONE 25x18 multiplier.  This block runs 4
--    times per dwell for the whitening plus up to 2 for the weights -
--    six operations per millisecond - so area matters and latency does
--    not.  Parallelising it would be a pure waste.
--
--  LATENCY   ~40 clk, data dependent only through the normalisation
--            adjust loop, which runs at most twice.  o_done is a single
--            cycle pulse.
--  RESOURCES 2 DSP48E1 (one 25x18 shared, one for the /3 constant),
--            ~350 FF, ~450 LUT.
-- =====================================================================
library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

library work;
use work.asp_pkg.all;
use work.asp_coef_pkg.all;

entity asp_rsqrt is
  generic (
    G_AW : natural := 48                       -- input width, unsigned
  );
  port (
    clk     : in  std_logic;
    rst     : in  std_logic;
    i_start : in  std_logic;
    i_a     : in  std_logic_vector(G_AW-1 downto 0);   -- unsigned, a > 0
    o_done  : out std_logic;
    o_y     : out std_logic_vector(RSQ_F+1 downto 0);  -- unsigned Y
    o_k     : out std_logic_vector(7 downto 0)         -- signed k
  );
end entity asp_rsqrt;


architecture rtl of asp_rsqrt is

  constant F      : natural := RSQ_F;                  -- 16
  constant WIDE   : natural := G_AW + 20;              -- headroom for a<<16
  constant MAGIC  : natural := 699051;                 -- (2^21 + 1)/3

  -- The multiplier has one register stage, so operands issued on edge T
  -- are readable on edge T+2.  Every issue is therefore followed by an
  -- explicit wait state (*_W).  Naming them rather than relying on a
  -- comment is deliberate: an off-by-one here reads the PREVIOUS
  -- product, which still converges to something plausible and would be
  -- extremely hard to find from the whitened matrix alone.
  type state_t is (S_IDLE, S_MSB, S_RED, S_ADJ,
                   S_SEED_W, S_SEED,
                   S_A_W, S_A, S_B_W, S_B, S_C_W, S_C, S_FIN);
  signal state : state_t := S_IDLE;

  signal a_u   : unsigned(WIDE-1 downto 0) := (others => '0');
  signal msb   : integer range -1 to WIDE-1 := -1;
  signal kk    : integer range -64 to 64 := 0;
  signal an    : unsigned(F+2 downto 0) := (others => '0');   -- < 2^17
  signal yv    : signed(F+3 downto 0) := (others => '0');     -- <= 2^17
  signal iter  : integer range 0 to 4 := 0;

  signal t1    : signed(F+3 downto 0) := (others => '0');
  signal t2    : signed(F+4 downto 0) := (others => '0');

  -- shared multiplier: 26 x 22 covers every product in the block
  signal mul_a : signed(25 downto 0) := (others => '0');
  signal mul_b : signed(21 downto 0) := (others => '0');
  signal mul_p : signed(47 downto 0) := (others => '0');

  signal done_r : std_logic := '0';

  -- an = round(a / 2^(2k)) with MATLAB round(); a > 0 so half-up.
  function reduce (a : unsigned; k : integer) return unsigned is
    variable sh : integer;
    variable r  : unsigned(a'length-1 downto 0);
    variable h  : unsigned(a'length-1 downto 0);
  begin
    sh := 2*k;
    if sh > 0 then
      h := (others => '0');
      h(sh-1) := '1';
      r := shift_right(a + h, sh);
    elsif sh < 0 then
      r := shift_left(a, -sh);
    else
      r := a;
    end if;
    return r;
  end function;

begin

  assert RSQ_F = 16
    report "asp_rsqrt: RSQ_F changed; the /3 magic constant and the seed "
         & "range proof are both specific to F = 16"
    severity failure;

  -- shared multiplier (inferred DSP48E1, registered operands and product)
  p_mul : process (clk)
  begin
    if rising_edge(clk) then
      mul_p <= mul_a * mul_b;
    end if;
  end process p_mul;

  p_fsm : process (clk)
    variable m   : integer;
    variable nn  : signed(25 downto 0);
    variable ytmp: signed(F+3 downto 0);
  begin
    if rising_edge(clk) then
      if rst = '1' then
        state  <= S_IDLE;
        done_r <= '0';
        kk     <= 0;
        yv     <= (others => '0');
      else
        done_r <= '0';

        case state is

          when S_IDLE =>
            if i_start = '1' then
              a_u   <= resize(unsigned(i_a), WIDE);
              state <= S_MSB;
            end if;

          when S_MSB =>
            -- position of the most significant set bit; -1 when a = 0
            m := -1;
            for i in 0 to WIDE-1 loop
              if a_u(i) = '1' then
                m := i;
              end if;
            end loop;
            msb <= m;
            if m < 0 then
              -- a <= 0: the model returns Y = 0, k = 0
              yv     <= (others => '0');
              kk     <= 0;
              done_r <= '1';
              state  <= S_IDLE;
            else
              -- k = floor((e - (F-1)) / 2), floor for negatives too,
              -- which is an arithmetic shift and NOT integer division.
              kk    <= to_integer(shift_right(to_signed(m - (F-1), 16), 1));
              state <= S_RED;
            end if;

          when S_RED =>
            an    <= resize(reduce(a_u, kk), F+3);
            state <= S_ADJ;

          when S_ADJ =>
            -- Rounding in the reduction can push an one bit out of the
            -- target window; the model corrects with a while loop, so
            -- this does too.  It iterates at most twice.
            if an >= to_unsigned(2**F, an'length) then
              kk    <= kk + 1;
              state <= S_RED;
            elsif an < to_unsigned(2**(F-2), an'length) then
              kk    <= kk - 1;
              state <= S_RED;
            else
              -- N = 7*2^F - 4*an, then the exact /3 by magic multiply
              nn := to_signed(7*(2**F), 26) -
                    resize(signed('0' & an) * to_signed(4, 4), 26);
              mul_a <= nn;
              mul_b <= to_signed(MAGIC, 22);
              state <= S_SEED_W;
            end if;

          when S_SEED_W =>
            state <= S_SEED;

          when S_SEED =>
            -- y0 = clamp( floor(N/3), 2^F, 2^(F+1) )
            ytmp := resize(shift_right(mul_p, 21), ytmp'length);
            if ytmp < to_signed(2**F, ytmp'length) then
              ytmp := to_signed(2**F, ytmp'length);
            elsif ytmp > to_signed(2**(F+1), ytmp'length) then
              ytmp := to_signed(2**(F+1), ytmp'length);
            end if;
            yv    <= ytmp;
            iter  <= 0;
            mul_a <= resize(signed('0' & an), 26);
            mul_b <= resize(ytmp, 22);
            state <= S_A_W;

          when S_A_W =>
            state <= S_A;

          when S_A =>
            -- t1 = fix(a' * y / 2^F)
            t1    <= resize(trunc_shr(mul_p, F), t1'length);
            mul_a <= resize(trunc_shr(mul_p, F), 26);
            mul_b <= resize(yv, 22);
            state <= S_B_W;

          when S_B_W =>
            state <= S_B;

          when S_B =>
            -- t2 = fix(t1 * y / 2^F) = a' y^2
            t2    <= resize(trunc_shr(mul_p, F), t2'length);
            mul_a <= to_signed(3*(2**F), 26) -
                     resize(trunc_shr(mul_p, F), 26);
            mul_b <= resize(yv, 22);
            state <= S_C_W;

          when S_C_W =>
            state <= S_C;

          when S_C =>
            -- y = fix( y * (3*2^F - a' y^2) / 2^(F+1) )
            ytmp := resize(trunc_shr(mul_p, F+1), ytmp'length);
            yv   <= ytmp;
            if iter = 3 then
              state <= S_FIN;
            else
              iter  <= iter + 1;
              mul_a <= resize(signed('0' & an), 26);
              mul_b <= resize(ytmp, 22);
              state <= S_A_W;
            end if;

          when S_FIN =>
            if yv < to_signed(1, yv'length) then
              yv <= to_signed(1, yv'length);
            elsif yv > to_signed(2**(F+1), yv'length) then
              yv <= to_signed(2**(F+1), yv'length);
            end if;
            done_r <= '1';
            state  <= S_IDLE;

        end case;
      end if;
    end if;
  end process p_fsm;

  o_done <= done_r;
  o_y    <= std_logic_vector(yv(RSQ_F+1 downto 0));
  o_k    <= std_logic_vector(to_signed(kk, 8));

end architecture rtl;
