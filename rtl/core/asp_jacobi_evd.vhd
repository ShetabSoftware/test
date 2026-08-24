-- =====================================================================
--  asp_jacobi_evd  -  STAGE 6.  Hermitian eigen-decomposition.
--
--  MODEL REFERENCE : stage6_evd() in matlab/golden/asp_golden_model.m
--
--  FUNCTION
--    Cyclic Jacobi with CORDIC-derived rotations and a FIXED sweep
--    count: 6 sweeps x 6 planes = 36 rotations, always, regardless of
--    the data.  Returns the eigenvalues sorted descending and the
--    matching eigenvector matrix.
--
--  WHY JACOBI AND NOT ANYTHING ELSE
--    Every step is a unitary similarity transform, so the Frobenius norm
--    is invariant EXACTLY.  Three consequences, and each one decides a
--    piece of the hardware:
--      * zero dynamic-range growth, so ONE word length covers the whole
--        block with no rescaling and no overflow analysis beyond the
--        input;
--      * unconditional stability regardless of conditioning - which
--        matters here because R genuinely is ill conditioned in the
--        sense that counts: the signal is a ~25% perturbation of the
--        identity;
--      * fixed latency with no tolerance test and no data-dependent
--        loop count, which is what lets the 1 kHz schedule be static.
--    A QR iteration or a power method would give up at least one of
--    those three.
--
--  WHY THE ROTATION IS APPLIED WITH MULTIPLIERS
--    A CORDIC applied directly to the matrix would impose its 1.6468
--    gain on EVERY pass and need compensation at every step.  Deriving
--    cos/sin once per plane rotation and applying them with four
--    multipliers confines the gain compensation to a single constant
--    inside asp_cordic.  That is the whole reason this block contains
--    multipliers at all.
--
--  SCHEDULE PER PLANE ROTATION (p,q)
--    1. vectoring on A(p,q)                 -> m, alpha
--    2. vectoring on (A_pp - A_qq, 2m)      -> phi
--    3. theta = fix(phi/2)                  (a shift, exactly)
--    4. rotation on theta                   -> ct, st
--    5. rotation on -alpha                  -> ca, sa
--    6. column update  A(:,q) *= e^{-j a};  rotate (A(:,p), A(:,q))
--    7. row update     A(q,:) *= e^{+j a};  rotate (A(p,:), A(q,:))
--    8. eigenvector    U(:,q) *= e^{-j a};  rotate (U(:,p), U(:,q))
--    Steps 6 and 7 are strictly ordered: the row update reads elements
--    the column update has already written.  Overlapping them would be
--    a different algorithm.
--
--  ARITHMETIC DATAPATH
--    Four real multipliers, 32 x 18 each.  One issue performs one
--    complex operation:
--      OP_CMUL   b  = x * (c + js)      -> (m0-m1, m2+m3)
--      OP_U      u  =  a*ct + b*st      -> (m0+m1, m2+m3)
--      OP_V      v  = -a*st + b*ct      -> (m1-m0, m3-m2)
--    Products are held at FULL precision and the SUM is rounded once,
--    not each product - one rounding per output instead of two, which
--    is what the model does.
--
--  LATENCY   ~130 clk per plane rotation, 36 rotations, so ~4700 clk
--            per dwell out of the 130944 available: 3.6% duty.  Fully
--            deterministic - the only data dependence in the entire
--            block is the final sort.
--  RESOURCES 8 DSP48E1 (four 32x18 multipliers) plus the shared CORDIC,
--            ~3400 FF (A is 16x64, U is 16x40), ~2500 LUT.
-- =====================================================================
library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

library work;
use work.asp_pkg.all;
use work.asp_coef_pkg.all;

entity asp_jacobi_evd is
  generic (
    G_N : natural := 4
  );
  port (
    clk     : in  std_logic;
    rst     : in  std_logic;

    -- whitened matrix stream from stage 5 (column-major, 16 entries)
    s_valid : in  std_logic;
    s_idx   : in  std_logic_vector(3 downto 0);
    s_re    : in  std_logic_vector(W_EVD-1 downto 0);
    s_im    : in  std_logic_vector(W_EVD-1 downto 0);
    s_start : in  std_logic;                 -- begin after the last word

    -- eigenvector matrix stream, column-major, sorted with the eigenvalues
    o_valid : out std_logic;
    o_idx   : out std_logic_vector(3 downto 0);
    o_re    : out std_logic_vector(W_UVEC-1 downto 0);
    o_im    : out std_logic_vector(W_UVEC-1 downto 0);

    -- eigenvalues, descending, flattened; stable from o_done
    o_lam   : out std_logic_vector(G_N*W_EVD-1 downto 0);
    o_done  : out std_logic
  );
end entity asp_jacobi_evd;


architecture rtl of asp_jacobi_evd is

  constant CW  : natural := 40;               -- CORDIC datapath
  constant ZW  : natural := 22;
  constant PW  : natural := W_EVD + W_ROT;    -- 50, one product
  constant SW  : natural := PW + 2;           -- 52, sum of two

  type a_arr_t is array (0 to 15) of signed(W_EVD-1 downto 0);
  signal a_re, a_im : a_arr_t := (others => (others => '0'));

  type u_arr_t is array (0 to 15) of signed(W_UVEC-1 downto 0);
  signal u_re, u_im : u_arr_t := (others => (others => '0'));

  type lam_t is array (0 to G_N-1) of signed(W_EVD-1 downto 0);
  signal lam : lam_t := (others => (others => '0'));

  -- ---- CORDIC ---------------------------------------------------------
  signal cd_start : std_logic := '0';
  signal cd_mode  : std_logic := '0';
  signal cd_x, cd_y : std_logic_vector(CW-1 downto 0) := (others => '0');
  signal cd_zi    : std_logic_vector(ZW-1 downto 0) := (others => '0');
  signal cd_done  : std_logic;
  signal cd_m     : std_logic_vector(CW-1 downto 0);
  signal cd_zo    : std_logic_vector(ZW-1 downto 0);
  signal cd_c     : std_logic_vector(W_ROT-1 downto 0);
  signal cd_s     : std_logic_vector(W_ROT-1 downto 0);

  signal m_v      : signed(CW-1 downto 0) := (others => '0');
  signal alpha_v  : signed(ZW-1 downto 0) := (others => '0');
  signal theta_v  : signed(ZW-1 downto 0) := (others => '0');
  signal ct, st   : signed(W_ROT-1 downto 0) := (others => '0');
  signal ca, sa   : signed(W_ROT-1 downto 0) := (others => '0');

  -- ---- multiplier bank ------------------------------------------------
  type mx_t is array (0 to 3) of signed(W_EVD-1 downto 0);
  type my_t is array (0 to 3) of signed(W_ROT-1 downto 0);
  signal mx : mx_t := (others => (others => '0'));
  signal my : my_t := (others => (others => '0'));
  type mp_t is array (0 to 3) of signed(PW-1 downto 0);
  signal mp : mp_t := (others => (others => '0'));

  constant OP_CMUL : integer := 0;
  constant OP_U    : integer := 1;
  constant OP_V    : integer := 2;
  signal op_sel   : integer range 0 to 2 := 0;
  signal op_sel_d : integer range 0 to 2 := 0;

  signal sum_a, sum_b : signed(SW-1 downto 0) := (others => '0');

  -- ---- sequencer ------------------------------------------------------
  type state_t is (S_LOAD, S_ROT_INIT,
                   S_V1_ISS, S_V1_WAIT,
                   S_V2_ISS, S_V2_WAIT,
                   S_R1_ISS, S_R1_WAIT,
                   S_R2_ISS, S_R2_WAIT,
                   S_UPD_CM, S_UPD_CM_W, S_UPD_CM_R,
                   S_UPD_U,  S_UPD_U_W,  S_UPD_U_R,
                   S_UPD_V,  S_UPD_V_W,  S_UPD_V_R,
                   S_NEXT_K, S_NEXT_PLANE,
                   S_SORT, S_SORT_INS, S_EMIT, S_DONE);
  signal state : state_t := S_LOAD;
  signal sort_i : integer range 0 to G_N := 1;

  signal sweep : integer range 0 to JACOBI_SWEEPS := 0;
  signal pp    : integer range 0 to 3 := 0;
  signal qq    : integer range 0 to 3 := 1;
  signal kk    : integer range 0 to 4 := 0;
  signal part  : integer range 0 to 3 := 0;   -- 0 col, 1 row, 2 uvec

  -- b = the phase-rotated partner element, held between the three issues
  signal b_re, b_im : signed(W_EVD-1 downto 0) := (others => '0');
  signal a_hold_re, a_hold_im : signed(W_EVD-1 downto 0) := (others => '0');
  signal u_new_re, u_new_im   : signed(W_EVD-1 downto 0) := (others => '0');

  signal emit_i : integer range 0 to 16 := 0;
  signal out_v  : std_logic := '0';
  signal out_i  : unsigned(3 downto 0) := (others => '0');
  signal out_re : signed(W_UVEC-1 downto 0) := (others => '0');
  signal out_im : signed(W_UVEC-1 downto 0) := (others => '0');
  signal done_r : std_logic := '0';

  -- index of element k in column c / row r, column-major storage
  function cidx (k, c : integer) return integer is
  begin
    return 4*c + k;
  end function;
  function ridx (k, r : integer) return integer is
  begin
    return 4*k + r;
  end function;

begin

  assert JACOBI_ROTS = JACOBI_SWEEPS * JACOBI_PAIRS
    report "asp_jacobi_evd: rotation count inconsistent" severity failure;

  u_cordic : entity work.asp_cordic
    generic map (G_W => CW, G_ZW => ZW)
    port map (clk => clk, rst => rst,
              i_start => cd_start, i_mode => cd_mode,
              i_x => cd_x, i_y => cd_y, i_z => cd_zi,
              o_done => cd_done, o_m => cd_m, o_z => cd_zo,
              o_c => cd_c, o_s => cd_s);

  -- ---------------------------------------------------------------------
  -- Multiplier bank: four 32x18 products, registered in and out, then the
  -- two sums.  The sums are full precision; only the caller rounds.
  -- ---------------------------------------------------------------------
  p_mul : process (clk)
  begin
    if rising_edge(clk) then
      for i in 0 to 3 loop
        mp(i) <= mx(i) * my(i);
      end loop;
      op_sel_d <= op_sel;
    end if;
  end process p_mul;

  -- The two sums are COMBINATIONAL on the multiplier outputs, not
  -- registered.  That keeps the issue-to-result distance at two clocks
  -- (operands -> product -> sum readable) instead of three, and a 52-bit
  -- add after a DSP48 output register is a short path at 130.944 MHz.
  -- The alternative, registering them, needs a third wait state in every
  -- one of the three issues and buys nothing.
  p_sum : process (mp, op_sel_d)
  begin
    case op_sel_d is
      when OP_CMUL =>
        sum_a <= resize(mp(0), SW) - resize(mp(1), SW);
        sum_b <= resize(mp(2), SW) + resize(mp(3), SW);
      when OP_U =>
        sum_a <= resize(mp(0), SW) + resize(mp(1), SW);
        sum_b <= resize(mp(2), SW) + resize(mp(3), SW);
      when others =>   -- OP_V
        sum_a <= resize(mp(1), SW) - resize(mp(0), SW);
        sum_b <= resize(mp(3), SW) - resize(mp(2), SW);
    end case;
  end process p_sum;

  p_fsm : process (clk)
    variable ia, ib : integer range 0 to 15;
    variable best   : integer range 0 to 3;
    variable tmp_l  : signed(W_EVD-1 downto 0);
    variable tmp_ur, tmp_ui : signed(W_UVEC-1 downto 0);
  begin
    if rising_edge(clk) then
      if rst = '1' then
        state    <= S_LOAD;
        cd_start <= '0';
        out_v    <= '0';
        done_r   <= '0';
      else
        cd_start <= '0';
        out_v    <= '0';
        done_r   <= '0';

        case state is

          -- ------------------------------------------------- load
          when S_LOAD =>
            if s_valid = '1' then
              a_re(to_integer(unsigned(s_idx))) <= signed(s_re);
              a_im(to_integer(unsigned(s_idx))) <= signed(s_im);
            end if;
            if s_start = '1' then
              -- U starts as the identity in Q1.F_ROT
              for i in 0 to 15 loop
                u_re(i) <= (others => '0');
                u_im(i) <= (others => '0');
              end loop;
              for d in 0 to G_N-1 loop
                u_re(5*d) <= to_signed(2**F_ROT, W_UVEC);
              end loop;
              sweep <= 0;
              pp    <= 0;
              qq    <= 1;
              state <= S_ROT_INIT;
            end if;

          -- ------------------------------------------------- rotation
          when S_ROT_INIT =>
            state <= S_V1_ISS;

          when S_V1_ISS =>
            cd_mode  <= '0';
            cd_x     <= std_logic_vector(resize(a_re(cidx(pp,qq)), CW));
            cd_y     <= std_logic_vector(resize(a_im(cidx(pp,qq)), CW));
            cd_start <= '1';
            state    <= S_V1_WAIT;

          when S_V1_WAIT =>
            if cd_done = '1' then
              m_v     <= signed(cd_m);
              alpha_v <= signed(cd_zo);
              state   <= S_V2_ISS;
            end if;

          when S_V2_ISS =>
            cd_mode  <= '0';
            cd_x     <= std_logic_vector(
                          resize(a_re(cidx(pp,pp)), CW) -
                          resize(a_re(cidx(qq,qq)), CW));
            cd_y     <= std_logic_vector(shift_left(m_v, 1));
            cd_start <= '1';
            state    <= S_V2_WAIT;

          when S_V2_WAIT =>
            if cd_done = '1' then
              -- theta = fix(phi/2): truncation toward zero, not a shift
              theta_v <= trunc_shr(signed(cd_zo), 1);
              state   <= S_R1_ISS;
            end if;

          when S_R1_ISS =>
            cd_mode  <= '1';
            cd_zi    <= std_logic_vector(theta_v);
            cd_start <= '1';
            state    <= S_R1_WAIT;

          when S_R1_WAIT =>
            if cd_done = '1' then
              ct    <= signed(cd_c);
              st    <= signed(cd_s);
              state <= S_R2_ISS;
            end if;

          when S_R2_ISS =>
            cd_mode  <= '1';
            cd_zi    <= std_logic_vector(-alpha_v);
            cd_start <= '1';
            state    <= S_R2_WAIT;

          when S_R2_WAIT =>
            if cd_done = '1' then
              ca    <= signed(cd_c);
              sa    <= signed(cd_s);
              kk    <= 0;
              part  <= 0;
              state <= S_UPD_CM;
            end if;

          -- ------------------------------------------------- updates
          -- Three parts, each four elements, each element three issues:
          --   CM : b = partner * e^{-j alpha}   (row part uses e^{+j})
          --   U  : new value of the p element
          --   V  : new value of the q element
          when S_UPD_CM =>
            case part is
              when 0 =>  ia := cidx(kk, pp);  ib := cidx(kk, qq);
              when 1 =>  ia := ridx(kk, pp);  ib := ridx(kk, qq);
              when others => ia := cidx(kk, pp); ib := cidx(kk, qq);
            end case;

            if part = 2 then
              a_hold_re <= resize(u_re(ia), W_EVD);
              a_hold_im <= resize(u_im(ia), W_EVD);
              mx(0) <= resize(u_re(ib), W_EVD);
              mx(1) <= resize(u_im(ib), W_EVD);
              mx(2) <= resize(u_re(ib), W_EVD);
              mx(3) <= resize(u_im(ib), W_EVD);
            else
              a_hold_re <= a_re(ia);
              a_hold_im <= a_im(ia);
              mx(0) <= a_re(ib);
              mx(1) <= a_im(ib);
              mx(2) <= a_re(ib);
              mx(3) <= a_im(ib);
            end if;

            -- column and eigenvector parts use (ca + j sa); the row part
            -- uses the conjugate (ca - j sa)
            my(0) <= ca;
            my(3) <= ca;
            if part = 1 then
              my(1) <= -sa;
              my(2) <= -sa;
            else
              my(1) <= sa;
              my(2) <= sa;
            end if;
            op_sel <= OP_CMUL;
            state  <= S_UPD_CM_W;

          when S_UPD_CM_W =>
            state <= S_UPD_CM_R;

          when S_UPD_CM_R =>
            if part = 2 then
              b_re <= resize(shift_round_sat(sum_a, F_ROT, W_UVEC), W_EVD);
              b_im <= resize(shift_round_sat(sum_b, F_ROT, W_UVEC), W_EVD);
            else
              b_re <= shift_round_sat(sum_a, F_ROT, W_EVD);
              b_im <= shift_round_sat(sum_b, F_ROT, W_EVD);
            end if;
            state <= S_UPD_U;

          when S_UPD_U =>
            mx(0) <= a_hold_re;  my(0) <= ct;
            mx(1) <= b_re;       my(1) <= st;
            mx(2) <= a_hold_im;  my(2) <= ct;
            mx(3) <= b_im;       my(3) <= st;
            op_sel <= OP_U;
            state  <= S_UPD_U_W;

          when S_UPD_U_W =>
            state <= S_UPD_U_R;

          when S_UPD_U_R =>
            if part = 2 then
              u_new_re <= resize(shift_round_sat(sum_a, F_ROT, W_UVEC), W_EVD);
              u_new_im <= resize(shift_round_sat(sum_b, F_ROT, W_UVEC), W_EVD);
            else
              u_new_re <= shift_round_sat(sum_a, F_ROT, W_EVD);
              u_new_im <= shift_round_sat(sum_b, F_ROT, W_EVD);
            end if;
            state <= S_UPD_V;

          when S_UPD_V =>
            mx(0) <= a_hold_re;  my(0) <= st;
            mx(1) <= b_re;       my(1) <= ct;
            mx(2) <= a_hold_im;  my(2) <= st;
            mx(3) <= b_im;       my(3) <= ct;
            op_sel <= OP_V;
            state  <= S_UPD_V_W;

          when S_UPD_V_W =>
            state <= S_UPD_V_R;

          when S_UPD_V_R =>
            case part is
              when 0 =>  ia := cidx(kk, pp);  ib := cidx(kk, qq);
              when 1 =>  ia := ridx(kk, pp);  ib := ridx(kk, qq);
              when others => ia := cidx(kk, pp); ib := cidx(kk, qq);
            end case;
            if part = 2 then
              u_re(ia) <= resize(u_new_re, W_UVEC);
              u_im(ia) <= resize(u_new_im, W_UVEC);
              u_re(ib) <= shift_round_sat(sum_a, F_ROT, W_UVEC);
              u_im(ib) <= shift_round_sat(sum_b, F_ROT, W_UVEC);
            else
              a_re(ia) <= u_new_re;
              a_im(ia) <= u_new_im;
              a_re(ib) <= shift_round_sat(sum_a, F_ROT, W_EVD);
              a_im(ib) <= shift_round_sat(sum_b, F_ROT, W_EVD);
            end if;
            state <= S_NEXT_K;

          when S_NEXT_K =>
            if kk = G_N-1 then
              kk <= 0;
              if part = 2 then
                state <= S_NEXT_PLANE;
              else
                part  <= part + 1;
                state <= S_UPD_CM;
              end if;
            else
              kk    <= kk + 1;
              state <= S_UPD_CM;
            end if;

          when S_NEXT_PLANE =>
            if qq = G_N-1 then
              if pp = G_N-2 then
                if sweep = JACOBI_SWEEPS-1 then
                  state <= S_SORT;
                else
                  sweep <= sweep + 1;
                  pp <= 0; qq <= 1;
                  state <= S_ROT_INIT;
                end if;
              else
                pp <= pp + 1;
                qq <= pp + 2;
                state <= S_ROT_INIT;
              end if;
            else
              qq <= qq + 1;
              state <= S_ROT_INIT;
            end if;

          -- ------------------------------------------------- sort
          when S_SORT =>
            for d in 0 to G_N-1 loop
              lam(d) <= a_re(5*d);
            end loop;
            sort_i <= 1;
            state  <= S_SORT_INS;

          -- INSERTION sort, descending, one insertion per clock.
          --
          -- Insertion rather than selection because MATLAB's sort() is
          -- STABLE and a selection sort is not: equal eigenvalues would
          -- come back with their eigenvector columns in a different
          -- order than the model produced.  Ties have probability zero
          -- in real data, which is exactly why this would survive every
          -- test and then differ on some captured dwell.
          --
          -- One insertion per clock rather than a single-cycle sorting
          -- network because the network is three levels of compare-and-
          -- mux over 32-bit eigenvalues AND their 4 x 40-bit eigenvector
          -- columns; at 130.944 MHz that is a needless critical path for
          -- an operation with 130000 clocks of slack.
          when S_SORT_INS =>
            -- find the FIRST position holding a value smaller than the
            -- key; inserting there is what makes the sort stable
            best := sort_i;
            for j in G_N-1 downto 0 loop
              if j < sort_i then
                if lam(j) < lam(sort_i) then
                  best := j;
                end if;
              end if;
            end loop;

            if best /= sort_i then
              tmp_l := lam(sort_i);
              for j in 0 to G_N-2 loop
                if j >= best and j < sort_i then
                  lam(j+1) <= lam(j);          -- shift the prefix right
                  for k in 0 to G_N-1 loop
                    u_re(cidx(k,j+1)) <= u_re(cidx(k,j));
                    u_im(cidx(k,j+1)) <= u_im(cidx(k,j));
                  end loop;
                end if;
              end loop;
              lam(best) <= tmp_l;
              for k in 0 to G_N-1 loop
                u_re(cidx(k,best)) <= u_re(cidx(k,sort_i));
                u_im(cidx(k,best)) <= u_im(cidx(k,sort_i));
              end loop;
            end if;

            if sort_i = G_N-1 then
              emit_i <= 0;
              state  <= S_EMIT;
            else
              sort_i <= sort_i + 1;
            end if;

          when S_EMIT =>
            out_v  <= '1';
            out_i  <= to_unsigned(emit_i, 4);
            out_re <= u_re(emit_i);
            out_im <= u_im(emit_i);
            if emit_i = 15 then
              state <= S_DONE;
            else
              emit_i <= emit_i + 1;
            end if;

          when S_DONE =>
            done_r <= '1';
            state  <= S_LOAD;

        end case;
      end if;
    end if;
  end process p_fsm;

  o_valid <= out_v;
  o_idx   <= std_logic_vector(out_i);
  o_re    <= std_logic_vector(out_re);
  o_im    <= std_logic_vector(out_im);
  o_done  <= done_r;

  g_lam : for i in 0 to G_N-1 generate
    o_lam((i+1)*W_EVD-1 downto i*W_EVD) <= std_logic_vector(lam(i));
  end generate;

end architecture rtl;
