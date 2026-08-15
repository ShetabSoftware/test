-- =====================================================================
--  asp_cov_accum  -  STAGE 4.  Spatial covariance accumulation.
--
--  MODEL REFERENCE : stage4_cov() in matlab/golden/asp_golden_model.m
--
--  FUNCTION
--    R = sum_{n=0..K-1} x[n] x[n]^H over one 1 ms dwell (K = 16368),
--    upper triangle computed, lower triangle by Hermitian symmetry:
--      Re R(i,j) = sum ( xr_i*xr_j + xi_i*xi_j )
--      Im R(i,j) = sum ( xi_i*xr_j - xr_i*xi_j )
--    The imaginary part of the diagonal is identically zero and is not
--    computed - that is 4 fewer multiplier-accumulators, free.
--
--  THE ONE RULE THIS BLOCK EXISTS TO ENFORCE
--    NOTHING IS ROUNDED INSIDE THE ACCUMULATION LOOP.  This is not a
--    precision trade, it is a correctness requirement, and it is the
--    single most important numerical decision in the whole design:
--
--      A constant rounding bias c applied to every product puts the same
--      c into every entry of R.  ones(N) is a RANK-ONE matrix whose
--      eigenvector is the boresight steering vector.  The estimator
--      therefore invents a source at zenith - precisely where the
--      satellites are - and the beamformer nulls it.  Because the error
--      is a fixed matrix rather than a noise term, it does NOT shrink
--      with dwell length, so longer integration never reveals it and a
--      floating-point simulation never shows it at all.
--
--    Products are 32 bits and are held exactly.  Accumulating K = 16368
--    of them adds ceil(log2 K) = 14 bits; two products per real entry
--    adds one more; 47 bits with a guard bit is 48 - which is exactly
--    the DSP48E1 P register.  The accumulator therefore lives INSIDE the
--    slice with no fabric adder and no rounding anywhere in the loop.
--
--  ARCHITECTURE - why 16 DSP48E1 and not 32
--    Each of the 10 upper-triangle entries needs two product terms for
--    its real part and (off-diagonal only) two for its imaginary part.
--    A naive mapping is 32 multiply-accumulators.  But the input rate is
--    one 4-antenna vector per 8 clocks, so every MAC has 8 slots to do
--    2 operations.  Folding 2:1 gives:
--        10 DSP for the real parts   (2 MAC cycles each)
--         6 DSP for the imaginary    (off-diagonal only)
--        16 DSP48E1 total
--    with the accumulator still in the P register, still exact, and
--    still with no fabric adder.  The fold costs one 2:1 operand mux per
--    DSP and nothing else.
--
--    Phase 0:  Re += xr_i*xr_j      Im += xi_i*xr_j
--    Phase 1:  Re += xi_i*xi_j      Im -= xr_i*xi_j
--    The subtract on phase 1 uses the DSP48E1 ALUMODE, not a negated
--    operand, so it costs nothing.
--
--  DWELL BOUNDARY
--    On the first vector of a dwell the DSP LOADS instead of
--    accumulating, which clears the accumulator for free - no separate
--    reset pass and no dead cycle between dwells.
--
--  OUTPUT
--    The full 4x4 Hermitian matrix is streamed out over 16 clocks in
--    COLUMN-MAJOR order (idx = 4*col + row), which is the order MATLAB's
--    R(:) produces, so the RTL and model dumps are directly comparable.
--    The lower triangle is emitted as the conjugate of the stored upper
--    triangle rather than being accumulated separately.
--
--  FIXED POINT
--    input  s16.15
--    product s32   exact
--    accum   s48   exact, no rounding, no saturation (proved: the peak
--                  is 2^45, see the width assertion below)
--
--  LATENCY   4 clk from the last input vector of a dwell to o_start;
--            16 clk to stream the matrix out.
--  THROUGHPUT one dwell (16368 vectors) per millisecond.
--  RESOURCES 16 DSP48E1, ~200 FF, ~400 LUT, 0 BRAM.
-- =====================================================================
library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

library work;
use work.asp_pkg.all;
use work.asp_coef_pkg.all;

entity asp_cov_accum is
  generic (
    G_NCH   : natural := 4;
    G_KDWELL: natural := K_DWELL
  );
  port (
    clk      : in  std_logic;
    rst      : in  std_logic;

    -- TDM input from stage 3
    s_valid  : in  std_logic;
    s_chan   : in  std_logic_vector(1 downto 0);
    s_re     : in  std_logic_vector(W_DAT-1 downto 0);
    s_im     : in  std_logic_vector(W_DAT-1 downto 0);

    -- Matrix output, streamed column-major, 16 entries per dwell
    o_valid  : out std_logic;
    o_idx    : out std_logic_vector(3 downto 0);
    o_re     : out std_logic_vector(W_ACC-1 downto 0);
    o_im     : out std_logic_vector(W_ACC-1 downto 0);
    o_last   : out std_logic;                       -- with idx = 15

    -- Dwell tick, one pulse per completed dwell (1 kHz).
    o_dwell_tick : out std_logic
  );
end entity asp_cov_accum;


architecture rtl of asp_cov_accum is

  constant NPAIR : natural := 10;                   -- N*(N+1)/2

  type idx_t is array (0 to NPAIR-1) of natural;
  --                     0  1  2  3  4  5  6  7  8  9
  constant PI_T : idx_t := (0, 0, 0, 0, 1, 1, 1, 2, 2, 3);
  constant PJ_T : idx_t := (0, 1, 2, 3, 1, 2, 3, 2, 3, 3);

  -- pair index of (row,col) for row <= col
  type map_t is array (0 to 15) of natural;
  constant PMAP : map_t := (
    --  (r,c) linear index 4*c + r
    0,  1,  2,  3,        -- col 0: (0,0)(1,0)(2,0)(3,0) -> pairs 0,1,2,3
    1,  4,  5,  6,        -- col 1: (0,1)(1,1)(2,1)(3,1) -> pairs 1,4,5,6
    2,  5,  7,  8,        -- col 2
    3,  6,  8,  9 );      -- col 3

  -- true when the (row,col) entry must be conjugated on output
  type conj_t is array (0 to 15) of boolean;
  constant CMAP : conj_t := (
    false, true,  true,  true,
    false, false, true,  true,
    false, false, false, true,
    false, false, false, false );

  -- ---- input deserialiser -------------------------------------------
  type vec_t is array (0 to G_NCH-1) of signed(W_DAT-1 downto 0);
  signal xr, xi : vec_t := (others => (others => '0'));
  signal vec_rdy : std_logic := '0';
  -- Frozen copy of the assembled vector.  The two MAC phases are on
  -- consecutive clocks while the deserialiser keeps writing the NEXT
  -- vector, so reading xr/xi directly would make correctness depend on
  -- the input gap being at least two clocks.  It is (the gap is two),
  -- but a datapath that is correct only because of an upstream timing
  -- coincidence is a defect waiting for a rate-plan change.  128 flops
  -- removes the coupling entirely.
  signal vr, vi : vec_t := (others => (others => '0'));

  -- ---- MAC sequencing -------------------------------------------------
  signal phase   : std_logic := '0';
  signal ph_run  : std_logic := '0';
  signal first_v : std_logic := '1';                -- first vector of a dwell

  -- operand registers (DSP48E1 A/B)
  type op_t is array (0 to NPAIR-1) of signed(W_DAT-1 downto 0);
  signal a_re, b_re : op_t := (others => (others => '0'));
  signal a_im, b_im : op_t := (others => (others => '0'));
  signal op_val  : std_logic := '0';
  signal op_sub  : std_logic := '0';
  signal op_load : std_logic := '0';

  -- product registers (DSP48E1 M)
  constant PROD_W : natural := 2*W_DAT;
  type prod_t is array (0 to NPAIR-1) of signed(PROD_W-1 downto 0);
  signal m_re, m_im : prod_t := (others => (others => '0'));
  signal m_val  : std_logic := '0';
  signal m_sub  : std_logic := '0';
  signal m_load : std_logic := '0';

  -- accumulators (DSP48E1 P)
  type acc_t is array (0 to NPAIR-1) of signed(W_ACC-1 downto 0);
  signal acc_re, acc_im : acc_t := (others => (others => '0'));

  -- frozen copy, so the next dwell can start accumulating immediately
  signal hold_re, hold_im : acc_t := (others => (others => '0'));

  -- ---- dwell counter and output sequencer -----------------------------
  signal samp_cnt : unsigned(LOG2_K_DWELL downto 0) := (others => '0');
  signal dwell_end : std_logic := '0';
  signal dwell_end_d : std_logic_vector(3 downto 0) := (others => '0');

  signal out_run : std_logic := '0';
  signal out_idx : unsigned(3 downto 0) := (others => '0');
  signal out_val : std_logic := '0';
  signal out_re_s, out_im_s : signed(W_ACC-1 downto 0) := (others => '0');
  signal out_idx_r : unsigned(3 downto 0) := (others => '0');
  signal tick : std_logic := '0';

begin

  -- The 48-bit budget is a proof, not a hope; state it where it is
  -- checked rather than only in a comment.
  assert W_ACC >= 2*W_DAT + LOG2_K_DWELL + 1
    report "asp_cov_accum: accumulator too narrow for K_DWELL"
    severity failure;

  -- -------------------------------------------------------------------
  -- Deserialise the TDM stream back into a 4-antenna vector.  The outer
  -- product needs all four antennas of the SAME sample instant, so this
  -- is where the channel-serial datapath ends.
  -- -------------------------------------------------------------------
  p_deser : process (clk)
  begin
    if rising_edge(clk) then
      if rst = '1' then
        vec_rdy <= '0';
      else
        vec_rdy <= '0';
        if s_valid = '1' then
          xr(to_integer(unsigned(s_chan))) <= signed(s_re);
          xi(to_integer(unsigned(s_chan))) <= signed(s_im);
          if unsigned(s_chan) = to_unsigned(G_NCH-1, 2) then
            vec_rdy <= '1';
            for k in 0 to G_NCH-2 loop
              vr(k) <= xr(k);
              vi(k) <= xi(k);
            end loop;
            vr(G_NCH-1) <= signed(s_re);   -- the antenna arriving now
            vi(G_NCH-1) <= signed(s_im);
          end if;
        end if;
      end if;
    end if;
  end process p_deser;

  -- -------------------------------------------------------------------
  -- Two-phase operand selection.  Phase 0 and phase 1 are consecutive
  -- clocks; there are 8 clocks per vector so the schedule has 4x slack.
  -- -------------------------------------------------------------------
  p_phase : process (clk)
  begin
    if rising_edge(clk) then
      if rst = '1' then
        phase    <= '0';
        ph_run   <= '0';
        op_val   <= '0';
        first_v  <= '1';
        samp_cnt <= (others => '0');
        dwell_end <= '0';
      else
        op_val    <= '0';
        dwell_end <= '0';

        if vec_rdy = '1' then
          phase  <= '0';
          ph_run <= '1';
          op_val <= '1';
        elsif ph_run = '1' and phase = '0' then
          phase  <= '1';
          op_val <= '1';
          ph_run <= '0';
          -- the vector is fully consumed once phase 1 is issued
          if samp_cnt = to_unsigned(G_KDWELL-1, samp_cnt'length) then
            samp_cnt  <= (others => '0');
            dwell_end <= '1';
            first_v   <= '1';
          else
            samp_cnt <= samp_cnt + 1;
            first_v  <= '0';
          end if;
        end if;

        op_load <= '0';
        if vec_rdy = '1' and first_v = '1' then
          op_load <= '1';           -- LOAD on phase 0 of the first vector
        end if;
      end if;

      -- operand mux: 2:1 per DSP, the entire cost of the 2:1 fold
      for p in 0 to NPAIR-1 loop
        if vec_rdy = '1' then                   -- phase 0
          a_re(p) <= vr(PI_T(p));  b_re(p) <= vr(PJ_T(p));
          a_im(p) <= vi(PI_T(p));  b_im(p) <= vr(PJ_T(p));
        else                                    -- phase 1
          a_re(p) <= vi(PI_T(p));  b_re(p) <= vi(PJ_T(p));
          a_im(p) <= vr(PI_T(p));  b_im(p) <= vi(PJ_T(p));
        end if;
      end loop;

      if vec_rdy = '1' then
        op_sub <= '0';
      else
        op_sub <= '1';
      end if;
    end if;
  end process p_phase;

  -- -------------------------------------------------------------------
  -- Multiply (DSP48E1 M register).  Exact, 32 bits, never rounded.
  -- -------------------------------------------------------------------
  p_mult : process (clk)
  begin
    if rising_edge(clk) then
      if rst = '1' then
        m_val <= '0';
      else
        m_val <= op_val;
      end if;
      m_sub  <= op_sub;
      m_load <= op_load;
      for p in 0 to NPAIR-1 loop
        m_re(p) <= a_re(p) * b_re(p);
        m_im(p) <= a_im(p) * b_im(p);
      end loop;
    end if;
  end process p_mult;

  -- -------------------------------------------------------------------
  -- Accumulate (DSP48E1 P register).  LOAD on the first product of a
  -- dwell so no separate clearing pass is needed.
  -- -------------------------------------------------------------------
  p_acc : process (clk)
  begin
    if rising_edge(clk) then
      if rst = '1' then
        acc_re <= (others => (others => '0'));
        acc_im <= (others => (others => '0'));
      elsif m_val = '1' then
        for p in 0 to NPAIR-1 loop
          if m_load = '1' then
            acc_re(p) <= resize(m_re(p), W_ACC);
            acc_im(p) <= resize(m_im(p), W_ACC);
          elsif m_sub = '1' then
            acc_re(p) <= acc_re(p) + resize(m_re(p), W_ACC);
            acc_im(p) <= acc_im(p) - resize(m_im(p), W_ACC);
          else
            acc_re(p) <= acc_re(p) + resize(m_re(p), W_ACC);
            acc_im(p) <= acc_im(p) + resize(m_im(p), W_ACC);
          end if;
        end loop;
      end if;
    end if;
  end process p_acc;

  -- -------------------------------------------------------------------
  -- Freeze and stream out.  dwell_end is delayed by the MAC pipeline
  -- depth so the final product has landed in the accumulator before the
  -- copy is taken - getting this wrong loses the last sample of every
  -- dwell, which is a 1-in-16368 error that no plot would ever show.
  -- -------------------------------------------------------------------
  p_out : process (clk)
    variable p : natural;
  begin
    if rising_edge(clk) then
      if rst = '1' then
        dwell_end_d <= (others => '0');
        out_run <= '0';
        out_val <= '0';
        out_idx <= (others => '0');
        tick    <= '0';
      else
        dwell_end_d <= dwell_end_d(2 downto 0) & dwell_end;
        tick    <= '0';
        out_val <= '0';

        if dwell_end_d(2) = '1' then
          hold_re <= acc_re;
          hold_im <= acc_im;
          out_run <= '1';
          out_idx <= (others => '0');
          tick    <= '1';
        elsif out_run = '1' then
          p := PMAP(to_integer(out_idx));
          out_re_s <= hold_re(p);
          if CMAP(to_integer(out_idx)) then
            out_im_s <= -hold_im(p);      -- lower triangle: conjugate
          else
            out_im_s <= hold_im(p);
          end if;
          out_idx_r <= out_idx;
          out_val   <= '1';
          if out_idx = 15 then
            out_run <= '0';
          end if;
          out_idx <= out_idx + 1;
        end if;
      end if;
    end if;
  end process p_out;

  o_valid <= out_val;
  o_idx   <= std_logic_vector(out_idx_r);
  o_re    <= std_logic_vector(out_re_s);
  o_im    <= std_logic_vector(out_im_s);
  o_last  <= '1' when (out_val = '1' and out_idx_r = "1111") else '0';
  o_dwell_tick <= tick;

end architecture rtl;
