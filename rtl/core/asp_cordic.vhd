-- =====================================================================
--  asp_cordic  -  CORDIC in vectoring and rotation mode, shared hardware.
--
--  MODEL REFERENCE : cordicVec() and cordicRot() in
--                    matlab/golden/asp_golden_model.m
--
--  VECTORING (i_mode = '0')
--    (x, y)  ->  m = |(x,y)| with the CORDIC gain removed,
--                z = atan2(y, x) in units of radians * 2^ANG_F.
--
--  ROTATION  (i_mode = '1')
--    z       ->  c = cos(z), s = sin(z) in Q1.F_ROT.
--    The iteration is seeded with x = 1/K rather than 1, so the CORDIC
--    gain is cancelled by construction and there is no post-scaling
--    multiply at all.
--
--  WHY BOTH MODES SHARE ONE ENTITY
--    They are the same recurrence with a different decision variable
--    (sign of y versus sign of z) and a different seed.  Sharing costs
--    one mux on the decision and saves a whole adder/shifter datapath.
--    The Jacobi engine issues them strictly one at a time - two
--    vectorings, then two rotations, per plane rotation - so there is
--    nothing to gain from having two instances.
--
--  WHY THE DATA IS NOT ROTATED BY A CORDIC DIRECTLY
--    Applying a CORDIC to the matrix elements would impose its 1.6468
--    gain on every pass and require compensation at every step.
--    Deriving cos/sin ONCE per plane rotation and applying them with
--    four multipliers keeps the gain compensation to a single constant
--    inside this block.  That is why stage 6 has multipliers in it at
--    all.
--
--  ARITHMETIC
--    The shifted addend is fix(y / 2^i) - truncation TOWARD ZERO, not an
--    arithmetic shift.  The two differ by one LSB for negative operands,
--    every iteration, and the error compounds across 16 iterations into
--    a visibly wrong angle.  This is the single most common way a
--    hand-written CORDIC fails to match a reference.
--
--    Internal width is generous on purpose.  The MATLAB reference works
--    in unbounded doubles, so any saturation the RTL introduces that the
--    model does not have is a divergence.  In the Jacobi engine the
--    largest operand is 2*m with |A| <= 2^28, i.e. 2^29, and the CORDIC
--    grows x by the gain to ~2^29.7; G_W = 40 leaves more than 10 bits
--    of headroom above that.
--
--  ANGLE FOLDING
--    Vectoring converges only for |angle| < 1.7433 rad, so a left
--    half-plane input is pre-rotated by +/-90 degrees and the angle is
--    seeded with -/+ pi/2.  Rotation folds |z| > pi/2 by pi and negates
--    both outputs.  Both foldings are exact - they are swaps and sign
--    changes, not approximations.
--
--  LATENCY   vectoring 21 clk, rotation 19 clk.  o_done pulses once.
--  RESOURCES 2 DSP48E1 (the gain-compensation multiply, vectoring only),
--            ~400 FF, ~700 LUT (two barrel shifters dominate).
-- =====================================================================
library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

library work;
use work.asp_pkg.all;
use work.asp_coef_pkg.all;

entity asp_cordic is
  generic (
    G_W  : natural := 40;                 -- internal x/y width
    G_ZW : natural := 22                  -- angle width
  );
  port (
    clk     : in  std_logic;
    rst     : in  std_logic;

    i_start : in  std_logic;
    i_mode  : in  std_logic;              -- '0' vectoring, '1' rotation
    i_x     : in  std_logic_vector(G_W-1 downto 0);
    i_y     : in  std_logic_vector(G_W-1 downto 0);
    i_z     : in  std_logic_vector(G_ZW-1 downto 0);

    o_done  : out std_logic;
    o_m     : out std_logic_vector(G_W-1 downto 0);   -- vectoring
    o_z     : out std_logic_vector(G_ZW-1 downto 0);  -- vectoring
    o_c     : out std_logic_vector(W_ROT-1 downto 0); -- rotation
    o_s     : out std_logic_vector(W_ROT-1 downto 0)  -- rotation
  );
end entity asp_cordic;


architecture rtl of asp_cordic is

  constant PW : natural := G_W + 20;      -- gain-compensation product

  type state_t is (S_IDLE, S_ITER, S_MUL_W, S_MUL, S_FIN);
  signal state : state_t := S_IDLE;

  signal x_r, y_r : signed(G_W-1 downto 0) := (others => '0');
  signal z_r      : signed(G_ZW-1 downto 0) := (others => '0');
  signal it       : integer range 0 to CORDIC_N := 0;
  signal mode_r   : std_logic := '0';
  signal neg_r    : std_logic := '0';

  signal mul_a : signed(G_W-1 downto 0) := (others => '0');
  signal mul_b : signed(19 downto 0) := (others => '0');
  signal mul_p : signed(PW-1 downto 0) := (others => '0');

  signal m_o   : signed(G_W-1 downto 0) := (others => '0');
  signal done_r: std_logic := '0';

begin

  p_mul : process (clk)
  begin
    if rising_edge(clk) then
      mul_p <= mul_a * mul_b;
    end if;
  end process p_mul;

  p_fsm : process (clk)
    variable xv, yv : signed(G_W-1 downto 0);
    variable zv     : signed(G_ZW-1 downto 0);
    variable dx, dy : signed(G_W-1 downto 0);
  begin
    if rising_edge(clk) then
      if rst = '1' then
        state  <= S_IDLE;
        done_r <= '0';
        neg_r  <= '0';
      else
        done_r <= '0';

        case state is

          when S_IDLE =>
            if i_start = '1' then
              mode_r <= i_mode;
              neg_r  <= '0';
              it     <= 0;
              if i_mode = '0' then
                -- vectoring: fold the left half plane into the right
                xv := signed(i_x);
                yv := signed(i_y);
                if xv < 0 then
                  if yv >= 0 then
                    x_r <= yv;
                    y_r <= -xv;
                    z_r <= to_signed(ANG_HPI, G_ZW);
                  else
                    x_r <= -yv;
                    y_r <= xv;
                    z_r <= to_signed(-ANG_HPI, G_ZW);
                  end if;
                else
                  x_r <= xv;
                  y_r <= yv;
                  z_r <= (others => '0');
                end if;
              else
                -- rotation: seed with 1/K so the gain cancels, and fold
                -- |z| > pi/2 by pi with a sign flip on both outputs
                zv := signed(i_z);
                if zv > to_signed(ANG_HPI, G_ZW) then
                  zv := zv - to_signed(2*ANG_HPI, G_ZW);
                  neg_r <= '1';
                end if;
                if zv < to_signed(-ANG_HPI, G_ZW) then
                  zv := zv + to_signed(2*ANG_HPI, G_ZW);
                  neg_r <= '1';
                end if;
                x_r <= to_signed(CORDIC_INV_K, G_W);
                y_r <= (others => '0');
                z_r <= zv;
              end if;
              state <= S_ITER;
            end if;

          when S_ITER =>
            dx := trunc_shr(x_r, it);
            dy := trunc_shr(y_r, it);
            if mode_r = '0' then
              -- vectoring: drive y to zero
              if y_r >= 0 then
                x_r <= x_r + dy;
                y_r <= y_r - dx;
                z_r <= z_r + to_signed(CORDIC_ATAN(it), G_ZW);
              else
                x_r <= x_r - dy;
                y_r <= y_r + dx;
                z_r <= z_r - to_signed(CORDIC_ATAN(it), G_ZW);
              end if;
            else
              -- rotation: drive z to zero
              if z_r >= 0 then
                x_r <= x_r - dy;
                y_r <= y_r + dx;
                z_r <= z_r - to_signed(CORDIC_ATAN(it), G_ZW);
              else
                x_r <= x_r + dy;
                y_r <= y_r - dx;
                z_r <= z_r + to_signed(CORDIC_ATAN(it), G_ZW);
              end if;
            end if;

            if it = CORDIC_N-1 then
              if mode_r = '0' then
                state <= S_MUL_W;
              else
                state <= S_FIN;
              end if;
            else
              it <= it + 1;
            end if;

          when S_MUL_W =>
            -- x_r has just settled; issue the gain-compensation multiply
            mul_a <= x_r;
            mul_b <= to_signed(CORDIC_INV_K, 20);
            state <= S_MUL;

          when S_MUL =>
            state <= S_FIN;

          when S_FIN =>
            if mode_r = '0' then
              m_o <= resize(trunc_shr(mul_p, 16), G_W);
              -- wrap the angle into (-pi, pi]
              if z_r > to_signed(ANG_PI, G_ZW) then
                z_r <= z_r - to_signed(2*ANG_PI, G_ZW);
              elsif z_r < to_signed(-ANG_PI, G_ZW) then
                z_r <= z_r + to_signed(2*ANG_PI, G_ZW);
              end if;
            else
              if neg_r = '1' then
                x_r <= -x_r;
                y_r <= -y_r;
              end if;
            end if;
            done_r <= '1';
            state  <= S_IDLE;

        end case;
      end if;
    end if;
  end process p_fsm;

  o_done <= done_r;
  o_m    <= std_logic_vector(m_o);
  o_z    <= std_logic_vector(z_r);
  o_c    <= std_logic_vector(clamp_s(x_r, W_ROT));
  o_s    <= std_logic_vector(clamp_s(y_r, W_ROT));

end architecture rtl;
