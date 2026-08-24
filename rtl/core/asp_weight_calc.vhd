-- =====================================================================
--  asp_weight_calc  -  STAGE 8.  Beamformer weights by projection.
--
--  MODEL REFERENCE : stage8_weights() in matlab/golden/asp_golden_model.m
--
--  FUNCTION
--    y_j = D * U(:,j)                    undo the whitening
--    Q   = orthonormalise(y)             modified Gram-Schmidt
--    w   = h - Q (Q^H h)                 project the quiescent beam
--    w   = block-float w to the top of its word
--
--  THE PROJECTOR IS NEVER FORMED
--    f = h - Q(Q^H h) costs 2*N*rank multiplies.  Forming P = I - yy^H/
--    (y^H y) and applying it costs N^2 plus a division.  For N = 4,
--    rank 1 that is 8 multiplies against 16 plus a divider - and the
--    divider is the expensive part, not the multiplies.
--
--  D = diag(sqrt(R_ii)) IS NOT OPTIONAL
--    It maps the eigenvector back out of the whitened domain into the
--    measured domain where the beamformer actually operates.  Skip it
--    and the null is steered in the whitened domain, which is a
--    different direction - and one that still looks like a null in any
--    whitened-domain plot.
--
--  THERE IS NO NORMALISATION ANYWHERE
--    The projector is scale invariant, and the beamformer output feeds a
--    correlator and a tracking loop that are invariant to a constant
--    complex gain.  So the three norm() divisions a naive implementation
--    would have are all removable.  What DOES matter is that the weights
--    sit at the top of their word, and that is a leading-zero count and
--    a barrel shift - about 30 LUTs, no added error at all, versus an
--    inverse-square-root block.
--
--    The one division-like operation left is the reciprocal square root
--    of each Gram-Schmidt column's own squared norm, which runs at most
--    MAX_RANK times per dwell - twice per millisecond.
--
--  RANK 0
--    Returns the quiescent beam unchanged.  This is the normal state:
--    with no spoofer the detector does not fire, and the array must then
--    behave exactly like a single antenna rather than steering anything.
--
--  ROUNDING
--    round_away_shr (MATLAB round) inside the block-float scaling,
--    conv_round via shift_round_sat everywhere else.  The model uses
--    round() in bfpScale and convRound() in the datapath, and the two
--    differ on exact ties - which occur constantly here because the
--    block-float shift is a power of two.
--
--  ARCHITECTURE
--    One 48x26 multiplier, fully sequential.  Runs once per dwell, about
--    300 clocks out of 130944.
--
--  LATENCY   ~350 clk from i_start to o_done (rank 1).
--  RESOURCES 6 DSP48E1 (the shared wide multiply) plus one asp_rsqrt,
--            ~900 FF, ~1400 LUT.
-- =====================================================================
library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

library work;
use work.asp_pkg.all;
use work.asp_coef_pkg.all;

entity asp_weight_calc is
  generic (
    G_N   : natural := 4;
    G_DSW : natural := 24
  );
  port (
    clk     : in  std_logic;
    rst     : in  std_logic;

    -- eigenvector matrix stream from stage 6 (column-major)
    s_valid : in  std_logic;
    s_idx   : in  std_logic_vector(3 downto 0);
    s_re    : in  std_logic_vector(W_UVEC-1 downto 0);
    s_im    : in  std_logic_vector(W_UVEC-1 downto 0);

    i_start : in  std_logic;
    i_rank  : in  std_logic_vector(1 downto 0);
    i_dsq   : in  std_logic_vector(G_N*G_DSW-1 downto 0);

    -- {re,im} per element, element k at bits ((2k+2)*W_WGT-1 .. 2k*W_WGT)
    o_w     : out std_logic_vector(2*G_N*W_WGT-1 downto 0);
    o_done  : out std_logic
  );
end entity asp_weight_calc;


architecture rtl of asp_weight_calc is

  constant WW  : natural := 48;      -- block-float working width
  constant CMP : natural := 72;      -- block-float comparison width
  constant PW  : natural := 74;      -- multiplier product

  type u_arr_t is array (0 to 15) of signed(W_UVEC-1 downto 0);
  signal u_re, u_im : u_arr_t := (others => (others => '0'));

  type d_arr_t is array (0 to G_N-1) of unsigned(G_DSW-1 downto 0);
  signal dsq : d_arr_t := (others => (others => '0'));

  -- working vector, wide enough for U*dsq before block-floating
  type wv_t is array (0 to G_N-1) of signed(WW-1 downto 0);
  signal vr, vi : wv_t := (others => (others => '0'));

  -- orthonormal basis Q, at most MAX_RANK columns
  type q_t is array (0 to 2*G_N-1) of signed(W_ROT-1 downto 0);
  signal qr, qi : q_t := (others => (others => '0'));

  -- output weights
  type w_t is array (0 to G_N-1) of signed(W_WGT+4-1 downto 0);
  signal wr, wi : w_t := (others => (others => '0'));
  type wo_t is array (0 to G_N-1) of signed(W_WGT-1 downto 0);
  signal wor, woi : wo_t := (others => (others => '0'));

  -- shared multiplier
  signal mul_a : signed(47 downto 0) := (others => '0');
  signal mul_b : signed(25 downto 0) := (others => '0');
  signal mul_p : signed(PW-1 downto 0) := (others => '0');

  signal acc_r, acc_i : signed(63 downto 0) := (others => '0');

  -- rsqrt
  signal rq_start : std_logic := '0';
  signal rq_a     : std_logic_vector(47 downto 0) := (others => '0');
  signal rq_done  : std_logic;
  signal rq_y     : std_logic_vector(RSQ_F+1 downto 0);
  signal rq_k     : std_logic_vector(7 downto 0);
  signal yr_s     : unsigned(RSQ_F+1 downto 0) := (others => '0');
  signal kr_s     : integer range -64 to 64 := 0;

  type state_t is (
    S_IDLE,
    S_Y0, S_Y_W, S_Y_R,                     -- Y = U(:,j) .* dsq
    S_BF_PK, S_BF_PK_CH,                    -- peak mag, one channel/clk
    S_BF_S1, S_BF_S2, S_BF_AP,              -- block-float scale
    S_IP0, S_IP_W, S_IP_D, S_IP_R,          -- ip  = Q_i^H x
    S_SB0, S_SB_W, S_SB_D, S_SB_R,          -- x  -= (Q_i * ip) >> 2*F_ROT
    S_IP_NEXT,
    S_NRM0, S_NRM_W, S_NRM_R,               -- nrm2 = |v|^2
    S_RSQ, S_RSQ_W,
    S_QCOL0, S_QCOL_W, S_QCOL_R,            -- q = v * Yr >> (8+kr)
    S_NEXT_J,
    S_DONE);
  signal state : state_t := S_IDLE;

  -- The Gram-Schmidt deflation and the final projection are the SAME
  -- operation on different vectors:  x -= Q_i (Q_i^H x).  One engine
  -- serves both; op_tgt selects whether x is the working vector v (width
  -- W_UVEC, during orthonormalisation) or the weight vector w (width
  -- W_WGT+4, during projection).  Writing it twice would be two chances
  -- to get the conjugate the wrong way round.
  signal op_tgt : std_logic := '0';         -- '0' = v, '1' = w
  signal ph     : integer range 0 to 3 := 0;
  signal ipr, ipi : signed(47 downto 0) := (others => '0');
  signal tr_s, ti_s : signed(63 downto 0) := (others => '0');

  signal jj, ii, kk : integer range 0 to 4 := 0;
  signal nq         : integer range 0 to 2 := 0;    -- columns of Q so far
  signal rank_r     : integer range 0 to 2 := 0;
  signal bf_target  : integer range 0 to 32 := W_UVEC;
  signal bf_after   : integer range 0 to 3 := 0;    -- 0 = Y stage, 1 = final
  signal bf_s       : integer range -64 to 64 := 0;
  signal bf_pk      : unsigned(CMP-1 downto 0) := (others => '0');
  signal done_r     : std_logic := '0';
  signal do_im      : std_logic := '0';

begin

  u_rsqrt : entity work.asp_rsqrt
    generic map (G_AW => 48)
    port map (clk => clk, rst => rst, i_start => rq_start, i_a => rq_a,
              o_done => rq_done, o_y => rq_y, o_k => rq_k);

  p_mul : process (clk)
  begin
    if rising_edge(clk) then
      mul_p <= mul_a * mul_b;
    end if;
  end process p_mul;

  p_fsm : process (clk)
    variable pk   : unsigned(CMP-1 downto 0);
    variable a    : unsigned(CMP-1 downto 0);
    variable hi   : unsigned(CMP-1 downto 0);
    variable lim  : unsigned(CMP-1 downto 0);
    variable sh   : integer;
    variable src  : signed(WW-1 downto 0);
  begin
    if rising_edge(clk) then
      if rst = '1' then
        state    <= S_IDLE;
        done_r   <= '0';
        rq_start <= '0';
        nq       <= 0;
      else
        done_r   <= '0';
        rq_start <= '0';

        -- capture the eigenvector matrix as it streams past
        if s_valid = '1' then
          u_re(to_integer(unsigned(s_idx))) <= signed(s_re);
          u_im(to_integer(unsigned(s_idx))) <= signed(s_im);
        end if;

        case state is

          when S_IDLE =>
            if i_start = '1' then
              for k in 0 to G_N-1 loop
                dsq(k) <= unsigned(i_dsq((k+1)*G_DSW-1 downto k*G_DSW));
                -- w starts at the quiescent beam h
                wr(k) <= to_signed(WGT_QUIESCENT, W_WGT+4);
                wi(k) <= (others => '0');
              end loop;
              rank_r <= to_integer(unsigned(i_rank));
              nq     <= 0;
              jj     <= 0;
              if unsigned(i_rank) = 0 then
                -- No threat: hand back the quiescent beam UNTOUCHED.
                -- The model returns h and returns - it does NOT
                -- block-float it, and block-floating h would double it
                -- (h = 2^15 sits one bit below the top of an 18-bit
                -- word).  A doubled quiescent beam is not wrong in any
                -- way the array can see, since the whole chain is
                -- invariant to a common real scale - which is precisely
                -- why this would never show up as a symptom.
                -- This is the normal operating state, not an edge case.
                for k in 0 to G_N-1 loop
                  vr(k) <= to_signed(WGT_QUIESCENT, WW);
                  vi(k) <= (others => '0');
                end loop;
                state <= S_DONE;
              else
                kk    <= 0;
                state <= S_Y0;
              end if;
            end if;

          -- ---------------------------------------------- Y = U(:,j).*dsq
          when S_Y0 =>
            mul_a <= resize(u_re(4*jj + kk), 48);
            mul_b <= signed(resize(dsq(kk), 26));
            do_im <= '0';
            state <= S_Y_W;

          when S_Y_W =>
            state <= S_Y_R;

          when S_Y_R =>
            if do_im = '0' then
              vr(kk) <= resize(mul_p, WW);
              mul_a  <= resize(u_im(4*jj + kk), 48);
              mul_b  <= signed(resize(dsq(kk), 26));
              do_im  <= '1';
              state  <= S_Y_W;
            else
              vi(kk) <= resize(mul_p, WW);
              if kk = G_N-1 then
                kk        <= 0;
                bf_target <= W_UVEC;
                bf_after  <= 0;
                state     <= S_BF_PK;
              else
                kk    <= kk + 1;
                state <= S_Y0;
              end if;
            end if;

          -- ---------------------------------------------- block float
          -- Literal transcription of the model's bfpScale: find the peak
          -- magnitude, then step s until the peak sits in the top octave
          -- of the target word.  The two while loops and their s bounds
          -- are reproduced exactly, INCLUDING the -W clamp, because that
          -- clamp is what makes a very large input saturate rather than
          -- shift forever - and saturation is a behaviour the model has.
          --
          -- Peak search is one channel per clock (two abs + two compares)
          -- rather than eight 48-bit magnitudes in one cycle.  At 1 kHz
          -- the extra three clocks are free; at 130.944 MHz they are the
          -- difference between a depth-3 comparator tree and a path that
          -- closes timing without a multicycle exception.
          when S_BF_PK =>
            bf_pk <= (others => '0');
            bf_s  <= 0;
            kk    <= 0;
            state <= S_BF_PK_CH;

          when S_BF_PK_CH =>
            pk := bf_pk;
            a  := resize(abs_ext(vr(kk)), CMP);
            if a > pk then pk := a; end if;
            a := resize(abs_ext(vi(kk)), CMP);
            if a > pk then pk := a; end if;
            bf_pk <= pk;
            if kk = G_N-1 then
              if pk = 0 then
                state <= S_BF_AP;          -- all zero: leave unscaled
              else
                state <= S_BF_S1;
              end if;
            else
              kk <= kk + 1;
            end if;

          when S_BF_S1 =>
            hi := shift_left(to_unsigned(1, CMP), bf_target-2);
            if shift_left(bf_pk, bf_s) < hi and bf_s < bf_target then
              bf_s <= bf_s + 1;
            else
              state <= S_BF_S2;
            end if;

          when S_BF_S2 =>
            lim := shift_left(to_unsigned(1, CMP), bf_target-1) - 1;
            if bf_s >= 0 then
              a := shift_left(bf_pk, bf_s);
            else
              a := bf_pk;
              lim := shift_left(lim, -bf_s);
            end if;
            if a > lim and bf_s > -bf_target then
              bf_s <= bf_s - 1;
            else
              state <= S_BF_AP;
            end if;

          when S_BF_AP =>
            -- y = clamp(round(x * 2^s)), MATLAB round(): half away from
            -- zero, which is NOT the convergent round the datapath uses.
            -- The block-float shift is a power of two, so exact ties
            -- occur constantly here and the two rules differ on every
            -- one of them.
            --
            -- The clamp width is selected between two CONSTANTS rather
            -- than passed as a signal: a saturation whose width depends
            -- on a run-time value is not synthesisable, and writing it
            -- that way would simulate correctly and fail in Vivado.
            for k in 0 to G_N-1 loop
              if bf_s >= 0 then
                src := shift_left(vr(k), bf_s);
              else
                src := round_away_shr(vr(k), -bf_s);
              end if;
              if bf_target = W_UVEC then
                vr(k) <= resize(clamp_s(src, W_UVEC), WW);
              else
                vr(k) <= resize(clamp_s(src, W_WGT), WW);
              end if;

              if bf_s >= 0 then
                src := shift_left(vi(k), bf_s);
              else
                src := round_away_shr(vi(k), -bf_s);
              end if;
              if bf_target = W_UVEC then
                vi(k) <= resize(clamp_s(src, W_UVEC), WW);
              else
                vi(k) <= resize(clamp_s(src, W_WGT), WW);
              end if;
            end loop;
            if bf_after = 0 then
              ii     <= 0;
              op_tgt <= '0';
              state  <= S_IP0;
            else
              state <= S_DONE;
            end if;

          -- ------------------------------ shared:  x -= Q_i (Q_i^H x)
          when S_IP0 =>
            if ii >= nq then
              if op_tgt = '0' then
                kk    <= 0;
                acc_r <= (others => '0');
                state <= S_NRM0;
              else
                -- projection complete: block-float the weights
                bf_target <= W_WGT;
                bf_after  <= 1;
                for k in 0 to G_N-1 loop
                  vr(k) <= resize(wr(k), WW);
                  vi(k) <= resize(wi(k), WW);
                end loop;
                state <= S_BF_PK;
              end if;
            else
              kk    <= 0;
              ph    <= 0;
              acc_r <= (others => '0');
              acc_i <= (others => '0');
              state <= S_IP_W;
            end if;

          when S_IP_W =>
            -- issue the product for (kk, ph)
            --   ph 0: Qr*xr -> ipr   ph 1: Qi*xi -> ipr
            --   ph 2: Qr*xi -> ipi   ph 3: Qi*xr -> ipi (subtracted)
            if ph = 0 or ph = 2 then
              mul_a <= resize(qr(4*ii + kk), 48);
            else
              mul_a <= resize(qi(4*ii + kk), 48);
            end if;
            if op_tgt = '0' then
              if ph = 0 or ph = 3 then
                mul_b <= resize(vr(kk)(25 downto 0), 26);
              else
                mul_b <= resize(vi(kk)(25 downto 0), 26);
              end if;
            else
              if ph = 0 or ph = 3 then
                mul_b <= resize(wr(kk), 26);
              else
                mul_b <= resize(wi(kk), 26);
              end if;
            end if;
            state <= S_IP_D;

          when S_IP_D =>
            -- The multiplier has one register stage, so operands issued
            -- on edge T are readable on edge T+2.  Without this state the
            -- accumulator would sum the PREVIOUS product - which still
            -- produces a plausible weight vector pointing the wrong way.
            state <= S_IP_R;

          when S_IP_R =>
            -- conj(Q) . x = (Qr*xr + Qi*xi) + j (Qr*xi - Qi*xr)
            case ph is
              when 0 => acc_r <= acc_r + resize(mul_p, 64);
              when 1 => acc_r <= acc_r + resize(mul_p, 64);
              when 2 => acc_i <= acc_i + resize(mul_p, 64);
              when others => acc_i <= acc_i - resize(mul_p, 64);
            end case;
            if ph = 3 then
              ph <= 0;
              if kk = G_N-1 then
                kk    <= 0;
                state <= S_SB0;
              else
                kk    <= kk + 1;
                state <= S_IP_W;
              end if;
            else
              ph    <= ph + 1;
              state <= S_IP_W;
            end if;

          when S_SB0 =>
            ipr   <= resize(acc_r, 48);
            ipi   <= resize(acc_i, 48);
            kk    <= 0;
            ph    <= 0;
            tr_s  <= (others => '0');
            ti_s  <= (others => '0');
            state <= S_SB_W;

          when S_SB_W =>
            -- Q_i * ip = (Qr*ipr - Qi*ipi) + j (Qr*ipi + Qi*ipr)
            if ph = 0 or ph = 2 then
              mul_a <= ipr;
            else
              mul_a <= ipi;
            end if;
            if ph = 0 or ph = 3 then
              mul_b <= resize(qr(4*ii + kk), 26);
            else
              mul_b <= resize(qi(4*ii + kk), 26);
            end if;
            state <= S_SB_D;

          when S_SB_D =>
            state <= S_SB_R;

          when S_SB_R =>
            case ph is
              when 0 => tr_s <= resize(mul_p, 64);
              when 1 => tr_s <= tr_s - resize(mul_p, 64);
              when 2 => ti_s <= resize(mul_p, 64);
              when others => ti_s <= ti_s + resize(mul_p, 64);
            end case;
            if ph = 3 then
              ph <= 0;
              state <= S_IP_NEXT;
            else
              ph    <= ph + 1;
              state <= S_SB_W;
            end if;

          when S_IP_NEXT =>
            -- one rounding on the SUM of the two products, as the model
            -- does, then subtract from the target vector
            if op_tgt = '0' then
              vr(kk) <= vr(kk) -
                        resize(shift_round_sat(tr_s, 2*F_ROT, W_UVEC), WW);
              vi(kk) <= vi(kk) -
                        resize(shift_round_sat(ti_s, 2*F_ROT, W_UVEC), WW);
            else
              wr(kk) <= wr(kk) - shift_round_sat(tr_s, 2*F_ROT, W_WGT+4);
              wi(kk) <= wi(kk) - shift_round_sat(ti_s, 2*F_ROT, W_WGT+4);
            end if;
            if kk = G_N-1 then
              ii    <= ii + 1;
              state <= S_IP0;
            else
              kk    <= kk + 1;
              state <= S_SB_W;
            end if;

          -- ---------------------------------------------- norm
          when S_NRM0 =>
            mul_a <= resize(vr(kk), 48);
            mul_b <= resize(vr(kk)(25 downto 0), 26);
            do_im <= '0';
            state <= S_NRM_W;

          when S_NRM_W =>
            state <= S_NRM_R;

          when S_NRM_R =>
            acc_r <= acc_r + resize(mul_p, 64);
            if do_im = '0' then
              mul_a <= resize(vi(kk), 48);
              mul_b <= resize(vi(kk)(25 downto 0), 26);
              do_im <= '1';
              state <= S_NRM_W;
            else
              if kk = G_N-1 then
                state <= S_RSQ;
              else
                kk    <= kk + 1;
                state <= S_NRM0;
              end if;
            end if;

          when S_RSQ =>
            if acc_r <= 0 then
              -- degenerate column: the model skips it (continue), which
              -- leaves Q one column shorter and the beam correspondingly
              -- less constrained.  Not an error - it is what happens
              -- when an eigenvector is numerically zero.
              state <= S_NEXT_J;
            else
              rq_a     <= std_logic_vector(acc_r(47 downto 0));
              rq_start <= '1';
              state    <= S_RSQ_W;
            end if;

          when S_RSQ_W =>
            if rq_done = '1' then
              yr_s  <= unsigned(rq_y);
              kr_s  <= to_integer(signed(rq_k));
              kk    <= 0;
              state <= S_QCOL0;
            end if;

          when S_QCOL0 =>
            mul_a <= resize(vr(kk), 48);
            mul_b <= signed(resize(yr_s, 26));
            do_im <= '0';
            state <= S_QCOL_W;

          when S_QCOL_W =>
            state <= S_QCOL_R;

          when S_QCOL_R =>
            -- 1/sqrt(nrm2) = Yr * 2^-(3*RSQ_F/2 + kr); the model wants a
            -- Q(F_ROT) unit vector, so the shift is 3*RSQ_F/2 + kr - F_ROT
            sh := (3*RSQ_F)/2 + kr_s - F_ROT;
            if do_im = '0' then
              qr(4*nq + kk) <= shift_round_sat_v(mul_p, sh, W_ROT);
              mul_a <= resize(vi(kk), 48);
              mul_b <= signed(resize(yr_s, 26));
              do_im <= '1';
              state <= S_QCOL_W;
            else
              qi(4*nq + kk) <= shift_round_sat_v(mul_p, sh, W_ROT);
              if kk = G_N-1 then
                nq    <= nq + 1;
                state <= S_NEXT_J;
              else
                kk    <= kk + 1;
                state <= S_QCOL0;
              end if;
            end if;

          when S_NEXT_J =>
            if jj = rank_r-1 then
              ii     <= 0;
              op_tgt <= '1';               -- same engine, now on w
              state  <= S_IP0;
            else
              jj    <= jj + 1;
              kk    <= 0;
              state <= S_Y0;
            end if;

          when S_DONE =>
            for k in 0 to G_N-1 loop
              wor(k) <= resize(vr(k), W_WGT);
              woi(k) <= resize(vi(k), W_WGT);
            end loop;
            done_r <= '1';
            state  <= S_IDLE;

        end case;
      end if;
    end if;
  end process p_fsm;

  g_out : for k in 0 to G_N-1 generate
    o_w((2*k+1)*W_WGT-1 downto (2*k)*W_WGT)   <= std_logic_vector(wor(k));
    o_w((2*k+2)*W_WGT-1 downto (2*k+1)*W_WGT) <= std_logic_vector(woi(k));
  end generate;

  o_done <= done_r;

end architecture rtl;
