-- =====================================================================
--  asp_whiten  -  STAGE 5.  Rw = D^-1 R D^-1, D = diag(sqrt(R_ii)).
--
--  MODEL REFERENCE : stage5_whiten() in matlab/golden/asp_golden_model.m
--
--  WHY THIS BLOCK IS THE REASON THE SYSTEM NEEDS NO ARRAY CALIBRATION
--    With post-LNA channel mismatch the measured covariance is
--    R = C R0 C^H for some diagonal complex C.  The principal
--    eigenvector of R is NOT C*b unless C is a scalar times a unitary -
--    so a plain eig(R) estimator silently depends on the array being
--    calibrated.  Normalising by the diagonal turns C into a PURE PHASE
--    diagonal, which IS unitary, and the eigenvector comes back right.
--    Delete this block and the design acquires a calibration
--    requirement that nothing else in the chain will reveal.
--
--  BLOCK NORMALISATION
--    One COMMON right shift is applied to the whole matrix, derived from
--    a leading-zero count on the largest diagonal entry.  Because it is
--    common it cancels exactly in the ratio R(i,j)/sqrt(R_ii R_jj) and
--    therefore cannot bias the result - which a per-entry scaling would.
--
--  ARITHMETIC
--    Rs   = round(R / 2^sh)                MATLAB round(), half away
--                                          from zero, so negative
--                                          off-diagonals round the same
--                                          way the model rounds them
--    Y_i  = 1/sqrt(Rs_ii) from asp_rsqrt   (Y, k) pair
--    d_i  = floor(sqrt(Rs_ii)) from asp_isqrt, kept for stage 8
--    Rw   = shift_round_sat( Rs(i,j) * Y_i * Y_j,
--                            SH_EVD_BASE + k_i + k_j, W_EVD )
--
--    The two multiplies are chained at FULL precision (the product
--    reaches 2^50) and rounded ONCE, which is what the model does.
--    Rounding the intermediate would put a bias on the whitened
--    diagonal, and the EVD's word-length proof depends on that diagonal
--    landing at 2^F_EVD.  It lands there to about 3 parts in 10^5 - not
--    exactly, because Y carries only 17 bits - which is what leaves
--    trace(Rw) close enough to N*2^F_EVD for the bound to hold.
--
--  DIAGONAL
--    Forced exactly real on output.  It is real by construction, but
--    "by construction" means "up to the rounding of two independent
--    products", and a non-zero imaginary diagonal makes the matrix
--    non-Hermitian, which breaks the Jacobi invariant that keeps the
--    EVD numerically stable.
--
--  ARCHITECTURE
--    Fully sequential.  This runs ONCE per dwell - 1 kHz - so it shares
--    a single wide multiplier and a single rsqrt/isqrt pair.  Total
--    ~260 clocks out of the 130944 available in a dwell, i.e. 0.2% duty.
--
--  LATENCY   ~270 clk from the last covariance word to o_done.
--  RESOURCES 3 DSP48E1 (shared 40x20 multiply) plus the rsqrt/isqrt
--            instances, ~1900 FF (the 16-entry 96-bit matrix buffer
--            dominates), ~1200 LUT.
-- =====================================================================
library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

library work;
use work.asp_pkg.all;
use work.asp_coef_pkg.all;

entity asp_whiten is
  generic (
    G_N   : natural := 4;
    G_DSW : natural := 24                       -- sqrt(R_ii) width
  );
  port (
    clk     : in  std_logic;
    rst     : in  std_logic;

    -- covariance matrix stream from stage 4 (column-major, 16 entries)
    s_valid : in  std_logic;
    s_idx   : in  std_logic_vector(3 downto 0);
    s_re    : in  std_logic_vector(W_ACC-1 downto 0);
    s_im    : in  std_logic_vector(W_ACC-1 downto 0);
    s_last  : in  std_logic;

    -- whitened matrix stream, same order
    o_valid : out std_logic;
    o_idx   : out std_logic_vector(3 downto 0);
    o_re    : out std_logic_vector(W_EVD-1 downto 0);
    o_im    : out std_logic_vector(W_EVD-1 downto 0);
    o_done  : out std_logic;

    -- d_i = floor(sqrt(Rs_ii)), flattened, stable from o_done until the
    -- next dwell.  Consumed by stage 8 to undo the whitening.
    o_dsq   : out std_logic_vector(G_N*G_DSW-1 downto 0)
  );
end entity asp_whiten;


architecture rtl of asp_whiten is

  constant RSW  : natural := 20;                 -- Rs word (|Rs| <= 2^16)
  constant YW   : natural := RSQ_F + 2;          -- 18, Y <= 2^17
  constant MW   : natural := 60;                 -- chained product

  type acc_arr_t is array (0 to 15) of signed(W_ACC-1 downto 0);
  signal r_re, r_im : acc_arr_t := (others => (others => '0'));

  type rs_arr_t is array (0 to 15) of signed(RSW-1 downto 0);
  signal rs_re, rs_im : rs_arr_t := (others => (others => '0'));

  type y_arr_t is array (0 to G_N-1) of unsigned(YW-1 downto 0);
  signal yv : y_arr_t := (others => (others => '0'));
  type k_arr_t is array (0 to G_N-1) of integer range -64 to 64;
  signal kv : k_arr_t := (others => 0);
  type d_arr_t is array (0 to G_N-1) of unsigned(G_DSW-1 downto 0);
  signal dv : d_arr_t := (others => (others => '0'));

  signal shift_n : natural range 0 to 63 := 0;

  type state_t is (S_IDLE, S_NORM, S_NORM_SH, S_RS, S_SQ_START, S_SQ_WAIT,
                   S_M1_ISS, S_M1_W, S_M1, S_M2_W, S_M2, S_EMIT, S_FIN);
  signal state : state_t := S_IDLE;
  -- Registered peak diagonal so ceil_log2 runs alone on the next cycle
  -- (otherwise the max-of-four tree and the 48-bit priority encoder sit
  -- in one combinational path at clk_dsp).
  signal dmax_r : unsigned(W_ACC-1 downto 0) := (others => '0');

  signal ent   : integer range 0 to 16 := 0;      -- matrix entry index
  signal chn   : integer range 0 to 4 := 0;       -- diagonal index
  signal doing_im : std_logic := '0';

  -- shared multiplier
  signal mul_a : signed(39 downto 0) := (others => '0');
  signal mul_b : signed(19 downto 0) := (others => '0');
  signal mul_p : signed(MW-1 downto 0) := (others => '0');

  signal t_re, t_im : signed(39 downto 0) := (others => '0');
  signal w_re, w_im : signed(W_EVD-1 downto 0) := (others => '0');

  -- sqrt units
  signal sq_start : std_logic := '0';
  signal sq_a     : std_logic_vector(47 downto 0) := (others => '0');
  signal rs_done  : std_logic;
  signal rs_y     : std_logic_vector(RSQ_F+1 downto 0);
  signal rs_k     : std_logic_vector(7 downto 0);
  signal is_done  : std_logic;
  signal is_s     : std_logic_vector(23 downto 0);
  signal rs_seen, is_seen : std_logic := '0';

  signal out_val : std_logic := '0';
  signal out_idx : unsigned(3 downto 0) := (others => '0');
  signal done_r  : std_logic := '0';

begin

  assert SH_EVD_BASE = 3*RSQ_F - F_EVD
    report "asp_whiten: SH_EVD_BASE inconsistent with RSQ_F/F_EVD"
    severity failure;

  u_rsqrt : entity work.asp_rsqrt
    generic map (G_AW => 48)
    port map (clk => clk, rst => rst, i_start => sq_start, i_a => sq_a,
              o_done => rs_done, o_y => rs_y, o_k => rs_k);

  u_isqrt : entity work.asp_isqrt
    generic map (G_AW => 48)
    port map (clk => clk, rst => rst, i_start => sq_start, i_a => sq_a,
              o_done => is_done, o_s => is_s);

  p_mul : process (clk)
  begin
    if rising_edge(clk) then
      mul_p <= mul_a * mul_b;
    end if;
  end process p_mul;

  p_fsm : process (clk)
    variable dmax : unsigned(W_ACC-1 downto 0);
    variable cl   : natural;
    variable sh   : integer;
  begin
    if rising_edge(clk) then
      if rst = '1' then
        state    <= S_IDLE;
        out_val  <= '0';
        done_r   <= '0';
        sq_start <= '0';
      else
        out_val  <= '0';
        done_r   <= '0';
        sq_start <= '0';

        case state is

          when S_IDLE =>
            if s_valid = '1' then
              r_re(to_integer(unsigned(s_idx))) <= signed(s_re);
              r_im(to_integer(unsigned(s_idx))) <= signed(s_im);
              if s_last = '1' then
                state <= S_NORM;
              end if;
            end if;

          when S_NORM =>
            -- Peak of the DIAGONAL entries only (Cauchy-Schwarz: they
            -- dominate the off-diagonal magnitudes; they are also
            -- non-negative by construction).  Registered into dmax_r so
            -- the priority encoder below is not stacked on the same
            -- combinational path.
            dmax := (others => '0');
            for i in 0 to G_N-1 loop
              if r_re(5*i) > signed(dmax) then
                dmax := unsigned(r_re(5*i));
              end if;
            end loop;
            if dmax = 0 then
              dmax := to_unsigned(1, dmax'length);
            end if;
            dmax_r <= dmax;
            state  <= S_NORM_SH;

          when S_NORM_SH =>
            -- sh = max(0, ceil_log2(max(dmax,1)) - R_NORM_BITS)
            cl := ceil_log2_u(dmax_r);
            if cl > R_NORM_BITS then
              shift_n <= cl - R_NORM_BITS;
            else
              shift_n <= 0;
            end if;
            ent   <= 0;
            state <= S_RS;

          when S_RS =>
            -- Rs = round(R / 2^sh), MATLAB round(): half AWAY from zero.
            -- Cauchy-Schwarz bounds |R(i,j)| by the largest diagonal, so
            -- after the common shift every entry fits in 17 bits and the
            -- clamp never fires.  It is here so that a pathological input
            -- (an all-zero dwell, a dead antenna) saturates instead of
            -- wrapping, because a wrapped covariance entry is a
            -- full-scale phantom source.
            rs_re(ent) <= clamp_s(round_away_shr(r_re(ent), shift_n), RSW);
            rs_im(ent) <= clamp_s(round_away_shr(r_im(ent), shift_n), RSW);
            if ent = 15 then
              ent   <= 0;
              chn   <= 0;
              state <= S_SQ_START;
            else
              ent <= ent + 1;
            end if;

          when S_SQ_START =>
            sq_a     <= std_logic_vector(resize(rs_re(5*chn), 48));
            sq_start <= '1';
            rs_seen  <= '0';
            is_seen  <= '0';
            state    <= S_SQ_WAIT;

          when S_SQ_WAIT =>
            if rs_done = '1' then
              yv(chn) <= unsigned(rs_y);
              kv(chn) <= to_integer(signed(rs_k));
              rs_seen <= '1';
            end if;
            if is_done = '1' then
              dv(chn) <= unsigned(is_s);
              is_seen <= '1';
            end if;
            if (rs_done = '1' or rs_seen = '1') and
               (is_done = '1' or is_seen = '1') then
              if chn = G_N-1 then
                ent   <= 0;
                doing_im <= '0';
                state <= S_M1_ISS;
              else
                chn   <= chn + 1;
                state <= S_SQ_START;
              end if;
            end if;

          when S_M1_ISS =>
            -- first multiply: Rs(i,j) * Y_i
            if doing_im = '0' then
              mul_a <= resize(rs_re(ent), 40);
            else
              mul_a <= resize(rs_im(ent), 40);
            end if;
            mul_b <= signed(resize(yv(ent mod 4), 20));
            state <= S_M1_W;

          when S_M1_W =>
            state <= S_M1;

          when S_M1 =>
            -- second multiply: (Rs * Y_i) * Y_j
            mul_a <= resize(mul_p, 40);
            mul_b <= signed(resize(yv(ent / 4), 20));
            state <= S_M2_W;

          when S_M2_W =>
            state <= S_M2;

          when S_M2 =>
            sh := SH_EVD_BASE + kv(ent mod 4) + kv(ent / 4);
            if doing_im = '0' then
              w_re     <= shift_round_sat_v(mul_p, sh, W_EVD);
              doing_im <= '1';
              state    <= S_M1_ISS;
            else
              -- the diagonal is real by construction; force it so that
              -- the matrix handed to the Jacobi engine is exactly
              -- Hermitian, which its stability argument requires
              if (ent mod 4) = (ent / 4) then
                w_im <= (others => '0');
              else
                w_im <= shift_round_sat_v(mul_p, sh, W_EVD);
              end if;
              doing_im <= '0';
              state    <= S_EMIT;
            end if;

          when S_EMIT =>
            out_val <= '1';
            out_idx <= to_unsigned(ent, 4);
            if ent = 15 then
              state <= S_FIN;
            else
              ent   <= ent + 1;
              state <= S_M1_ISS;
            end if;

          when S_FIN =>
            done_r <= '1';
            state  <= S_IDLE;

        end case;
      end if;
    end if;
  end process p_fsm;

  o_valid <= out_val;
  o_idx   <= std_logic_vector(out_idx);
  o_re    <= std_logic_vector(w_re);
  o_im    <= std_logic_vector(w_im);
  o_done  <= done_r;

  g_dsq : for i in 0 to G_N-1 generate
    o_dsq((i+1)*G_DSW-1 downto i*G_DSW) <= std_logic_vector(dv(i));
  end generate;

end architecture rtl;
