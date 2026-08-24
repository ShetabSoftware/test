-- =====================================================================
--  asp_fir_shape  -  STAGE 3.  63-tap symmetric shaping FIR.
--
--  MODEL REFERENCE : stage3_fir() in matlab/golden/asp_golden_model.m
--
--  FUNCTION
--    acc[m] = sum_{j=0..62} h[j] * x[m-j]      (causal, full precision)
--    y      = round_even(acc >> 17), saturated to s16.15
--
--    The model uses a CENTRED convolution, y_model[n] = sum_j h[j]
--    x[n+31-j], so y_model[n] = acc[n+31].  The RTL computes the causal
--    form and discards the first (NTAPS-1)/2 = 31 outputs per channel.
--    That is the whole difference between the two, and it is stated here
--    rather than buried: a reader who assumes the RTL is sample-aligned
--    with the model will chase a 31-sample offset that is not a bug.
--
--  WHAT THE FILTER IS FOR - three jobs, not one
--    1. band-limit to the C/A main lobe;
--    2. reject the residual LO leakage, which after the DDC sits at
--       2.046 MHz.  This is the job that matters most, because a
--       per-channel DC term is a rank-one contribution to the spatial
--       covariance and the eigen-detector will find it and null it -
--       a phantom emitter built out of nothing but analogue mismatch;
--    3. reject the I/Q image.
--
--  ARCHITECTURE AND THE RESOURCE ARGUMENT
--    Rate budget: 4 antennas x 16.368 MHz x (I and Q) x 32 unique
--    coefficients = 4190 M MAC/s.  At clk = 130.944 MHz that is exactly
--    32 multipliers - which is why the clock was chosen as 8 x FS_WORK.
--    A fully parallel implementation would need 256 and does not fit in
--    an XC7Z020 alongside the covariance array.
--
--    The engine therefore runs a strict two-cycle schedule per channel
--    sample: cycle 0 computes the I output, cycle 1 the Q output, both
--    from the same 32 multipliers and the same delay line.  Symmetry
--    (h[j] = h[62-j]) supplies the factor of two that takes 63 taps down
--    to 32 multipliers, using the DSP48E1 PRE-ADDER so the 31 pair sums
--    cost no extra fabric.
--
--  DELAY LINE
--    One line shared by all four antennas.  In the TDM stream a sample
--    of the same channel is 4 slots back, so tap j is at sr(4j-1) and
--    tap 0 is the sample currently being processed.  Sharing the line is
--    not just a saving: it makes it structurally impossible for the four
--    antennas to acquire different group delays, which is the one error
--    the array processing downstream cannot tolerate or detect.
--
--    The line is NOT reset.  It is initialised to zero by the bitstream,
--    which is what lets Vivado map the 4-deep runs between taps into
--    SRL16 primitives; a reset would force 7936 flip-flops instead.  The
--    zero initial state is also exactly the model's zero history, so the
--    warm-up transient matches bit for bit.
--
--  INPUT BURST
--    Stage 2 emits four channels back to back and then idles for four
--    clocks (it decimates by two).  The engine consumes one sample every
--    two clocks steadily, so a small FIFO absorbs the burst.  Depth 8 is
--    twice the worst-case backlog of 4.
--
--  FIXED POINT
--    input   s16.15
--    coeffs  s18.17, DC gain exactly 2^17
--    pre-add s17
--    product s35
--    accum   s40  (worst case 2^15 * sum|h| = 6.38e9 needs 34 bits)
--    output  s16.15 after >> 17, convergent round, saturate
--
--  LATENCY   12 clk from the sample entering the engine to its output.
--  THROUGHPUT one channel-sample per two clk (65.472 MS/s aggregate).
--  RESOURCES 32 DSP48E1 (constant coefficients, pre-adder used),
--            ~8000 FF or ~2000 LUT of SRL for the delay line,
--            ~1500 LUT for the adder tree, 0 BRAM.
--
--  ALTERNATIVE
--    -- OPTIONAL VIVADO IP CORE (not used here, documented for reference):
--    -- IP NAME: FIR Compiler 7.2
--    -- CONFIGURATION:
--    --   Filter Type=Single_Rate, Number of Channels=8,
--    --   Input Sample Frequency=16.368, Clock Frequency=130.944,
--    --   Coefficient Type=Signed, Quantization=Integer_Coefficients,
--    --   Coefficient Width=18, Data Width=16, Output Rounding Mode=
--    --   Convergent_Rounding_to_Even, Output Width=16,
--    --   Coefficient Structure=Inferred (symmetry is detected),
--    --   Coefficient File=fir_shape.coe
--    -- Rejected here because bit-exactness against the golden model is
--    -- the acceptance criterion and the IP's internal accumulation
--    -- order is not user-visible.  It is a legitimate substitute if the
--    -- co-simulation in rtl/tb/tb_asp_fir_shape.vhd passes with it.
-- =====================================================================
library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

library work;
use work.asp_pkg.all;
use work.asp_coef_pkg.all;

entity asp_fir_shape is
  generic (
    G_NCH : natural := 4
  );
  port (
    clk     : in  std_logic;
    rst     : in  std_logic;

    s_valid : in  std_logic;
    s_chan  : in  std_logic_vector(1 downto 0);
    s_re    : in  std_logic_vector(W_HB-1 downto 0);
    s_im    : in  std_logic_vector(W_HB-1 downto 0);

    m_valid : out std_logic;
    m_chan  : out std_logic_vector(1 downto 0);
    m_re    : out std_logic_vector(W_DAT-1 downto 0);
    m_im    : out std_logic_vector(W_DAT-1 downto 0);

    -- Asserted if the input FIFO ever overflows.  Sticky; cleared only
    -- by rst.  Wired to a status register so a rate-plan error shows up
    -- as a flag rather than as silently corrupted data.
    o_overflow : out std_logic
  );
end entity asp_fir_shape;


architecture rtl of asp_fir_shape is

  constant NT      : natural := FIR_NTAPS;              -- 63
  constant NMULT   : natural := FIR_NMULT;              -- 32
  constant NPAIR   : natural := (NT-1)/2;               -- 31
  constant CENTRE  : natural := NPAIR;                  -- 31
  constant WARMUP  : natural := (NT-1)/2;               -- 31 discarded
  constant SR_LEN  : natural := G_NCH*(NT-1);           -- 248
  constant PA_W    : natural := W_HB + 1;               -- 17
  constant ACC_W   : natural := 40;

  -- ---- input FIFO ----------------------------------------------------
  constant FIFO_LOG : natural := 3;
  constant FIFO_N   : natural := 2**FIFO_LOG;
  subtype  fifo_word_t is std_logic_vector(2 + 2*W_HB - 1 downto 0);
  type     fifo_t is array (0 to FIFO_N-1) of fifo_word_t;
  signal   fifo     : fifo_t := (others => (others => '0'));
  signal   wr_ptr   : unsigned(FIFO_LOG-1 downto 0) := (others => '0');
  signal   rd_ptr   : unsigned(FIFO_LOG-1 downto 0) := (others => '0');
  signal   fcount   : unsigned(FIFO_LOG downto 0)   := (others => '0');
  signal   ovf      : std_logic := '0';

  -- ---- engine --------------------------------------------------------
  type state_t is (ST_IDLE, ST_I, ST_Q);
  signal state : state_t := ST_IDLE;

  signal cur_re  : signed(W_HB-1 downto 0) := (others => '0');
  signal cur_im  : signed(W_HB-1 downto 0) := (others => '0');
  signal cur_ch  : unsigned(1 downto 0) := (others => '0');

  -- ---- delay line: past samples only, tap j (j>=1) is sr(4j-1) -------
  type sr_t is array (0 to SR_LEN-1) of signed(W_HB-1 downto 0);
  signal sr_re : sr_t := (others => (others => '0'));
  signal sr_im : sr_t := (others => (others => '0'));

  -- ---- arithmetic pipeline -------------------------------------------
  type pa_t   is array (0 to NMULT-1) of signed(PA_W-1 downto 0);
  type prod_t is array (0 to NMULT-1) of signed(ACC_W-1 downto 0);

  signal pa      : pa_t   := (others => (others => '0'));
  signal pa_val  : std_logic := '0';
  signal pa_isq  : std_logic := '0';
  signal pa_ch   : unsigned(1 downto 0) := (others => '0');

  signal pr      : prod_t := (others => (others => '0'));
  signal pr_val  : std_logic := '0';
  signal pr_isq  : std_logic := '0';
  signal pr_ch   : unsigned(1 downto 0) := (others => '0');

  signal t1, t2, t3, t4 : prod_t := (others => (others => '0'));
  signal t1_val, t2_val, t3_val, t4_val : std_logic := '0';
  signal t1_isq, t2_isq, t3_isq, t4_isq : std_logic := '0';
  signal t1_ch, t2_ch, t3_ch, t4_ch : unsigned(1 downto 0) := (others => '0');

  signal acc     : signed(ACC_W-1 downto 0) := (others => '0');
  signal acc_val : std_logic := '0';
  signal acc_isq : std_logic := '0';
  signal acc_ch  : unsigned(1 downto 0) := (others => '0');

  signal hold_re : signed(W_DAT-1 downto 0) := (others => '0');

  -- ---- warm-up suppression, per channel ------------------------------
  type warm_t is array (0 to G_NCH-1) of unsigned(6 downto 0);
  signal warm : warm_t := (others => (others => '0'));

  signal out_val : std_logic := '0';
  signal out_ch  : unsigned(1 downto 0) := (others => '0');
  signal out_re  : signed(W_DAT-1 downto 0) := (others => '0');
  signal out_im  : signed(W_DAT-1 downto 0) := (others => '0');

begin

  assert FIR_NTAPS = 63
    report "asp_fir_shape: tap count changed; check SR_LEN and WARMUP"
    severity failure;
  assert SH_FIR = 17
    report "asp_fir_shape: output shift changed" severity failure;

  -- -------------------------------------------------------------------
  -- Input FIFO.  Absorbs stage 2's four-on / four-off burst.
  -- -------------------------------------------------------------------
  p_fifo : process (clk)
    variable pop : std_logic;
  begin
    if rising_edge(clk) then
      if rst = '1' then
        wr_ptr <= (others => '0');
        rd_ptr <= (others => '0');
        fcount <= (others => '0');
        ovf    <= '0';
      else
        pop := '0';
        if (state = ST_Q) or (state = ST_IDLE) then
          if fcount /= 0 then
            pop := '1';
          end if;
        end if;

        if s_valid = '1' then
          fifo(to_integer(wr_ptr)) <= s_chan & s_re & s_im;
          wr_ptr <= wr_ptr + 1;
          if fcount = FIFO_N then
            ovf <= '1';        -- sticky: a rate-plan error, not a glitch
          end if;
        end if;

        if pop = '1' then
          rd_ptr <= rd_ptr + 1;
        end if;

        if s_valid = '1' and pop = '0' then
          fcount <= fcount + 1;
        elsif s_valid = '0' and pop = '1' then
          fcount <= fcount - 1;
        end if;
      end if;
    end if;
  end process p_fifo;

  o_overflow <= ovf;

  -- -------------------------------------------------------------------
  -- Engine sequencer and delay line.
  --   ST_I : compute the I output of the held sample
  --   ST_Q : compute the Q output, then shift the held sample into the
  --          delay line and fetch the next one
  -- -------------------------------------------------------------------
  p_engine : process (clk)
    variable w : fifo_word_t;
  begin
    if rising_edge(clk) then
      if rst = '1' then
        state <= ST_IDLE;
      else
        case state is

          when ST_IDLE =>
            if fcount /= 0 then
              w      := fifo(to_integer(rd_ptr));
              cur_ch <= unsigned(w(2*W_HB+1 downto 2*W_HB));
              cur_re <= signed(w(2*W_HB-1 downto W_HB));
              cur_im <= signed(w(W_HB-1 downto 0));
              state  <= ST_I;
            end if;

          when ST_I =>
            state <= ST_Q;

          when ST_Q =>
            -- Shift the sample just processed into the line.  Both the
            -- I and the Q pass read the SAME line state, which is why
            -- the shift happens here and not on entry.
            sr_re(0) <= cur_re;
            sr_im(0) <= cur_im;
            for i in 1 to SR_LEN-1 loop
              sr_re(i) <= sr_re(i-1);
              sr_im(i) <= sr_im(i-1);
            end loop;

            if fcount /= 0 then
              w      := fifo(to_integer(rd_ptr));
              cur_ch <= unsigned(w(2*W_HB+1 downto 2*W_HB));
              cur_re <= signed(w(2*W_HB-1 downto W_HB));
              cur_im <= signed(w(W_HB-1 downto 0));
              state  <= ST_I;
            else
              state <= ST_IDLE;
            end if;

        end case;
      end if;
    end if;
  end process p_engine;

  -- -------------------------------------------------------------------
  -- Symmetric pre-adds.  pa(k) = t[k] + t[62-k] for k < 31,
  -- pa(31) = t[31] (the centre tap has no partner).
  -- Maps to the DSP48E1 pre-adder, so these 31 adds cost no fabric.
  -- -------------------------------------------------------------------
  p_preadd : process (clk)
    variable ta, tb : signed(W_HB-1 downto 0);
    variable isq    : std_logic;
  begin
    if rising_edge(clk) then
      if rst = '1' then
        pa_val <= '0';
      elsif state = ST_I or state = ST_Q then
        pa_val <= '1';
      else
        pa_val <= '0';
      end if;

      if state = ST_Q then
        isq := '1';
      else
        isq := '0';
      end if;
      pa_isq <= isq;
      pa_ch  <= cur_ch;

      for k in 0 to NMULT-1 loop
        -- low-index tap of the pair
        if k = 0 then
          if isq = '0' then ta := cur_re; else ta := cur_im; end if;
        else
          if isq = '0' then ta := sr_re(G_NCH*k - 1);
          else                ta := sr_im(G_NCH*k - 1); end if;
        end if;

        if k = CENTRE then
          pa(k) <= resize(ta, PA_W);        -- centre tap, unpaired
        else
          -- high-index partner: tap (NT-1-k)
          if isq = '0' then tb := sr_re(G_NCH*(NT-1-k) - 1);
          else                tb := sr_im(G_NCH*(NT-1-k) - 1); end if;
          pa(k) <= resize(ta, PA_W) + resize(tb, PA_W);
        end if;
      end loop;
    end if;
  end process p_preadd;

  -- -------------------------------------------------------------------
  -- 32 multiplies.  Coefficients are compile-time constants, so Vivado
  -- is free to implement the cheap ones in fabric and keep DSP48 slices
  -- for the large ones - which is strictly better than forcing 32 DSPs.
  -- -------------------------------------------------------------------
  p_mult : process (clk)
  begin
    if rising_edge(clk) then
      if rst = '1' then
        pr_val <= '0';
      else
        pr_val <= pa_val;
      end if;
      pr_isq <= pa_isq;
      pr_ch  <= pa_ch;
      for k in 0 to NMULT-1 loop
        pr(k) <= resize(pa(k) * to_signed(FIR_COEF(k), W_COEF), ACC_W);
      end loop;
    end if;
  end process p_mult;

  -- -------------------------------------------------------------------
  -- Pipelined adder tree, 32 -> 16 -> 8 -> 4 -> 2 -> 1.
  -- Full precision throughout: the model sums every product before
  -- rounding, so an intermediate rounding here would break bit-exactness
  -- by up to 1 LSB and would also re-introduce the accumulator bias that
  -- stage 4 exists to avoid.
  -- -------------------------------------------------------------------
  p_tree : process (clk)
  begin
    if rising_edge(clk) then
      if rst = '1' then
        t1_val <= '0'; t2_val <= '0'; t3_val <= '0'; t4_val <= '0';
        acc_val <= '0';
      else
        t1_val <= pr_val; t2_val <= t1_val; t3_val <= t2_val;
        t4_val <= t3_val; acc_val <= t4_val;
      end if;
      t1_isq <= pr_isq; t2_isq <= t1_isq; t3_isq <= t2_isq;
      t4_isq <= t3_isq; acc_isq <= t4_isq;
      t1_ch  <= pr_ch;  t2_ch  <= t1_ch;  t3_ch  <= t2_ch;
      t4_ch  <= t3_ch;  acc_ch <= t4_ch;

      for k in 0 to 15 loop
        t1(k) <= pr(2*k) + pr(2*k+1);
      end loop;
      for k in 0 to 7 loop
        t2(k) <= t1(2*k) + t1(2*k+1);
      end loop;
      for k in 0 to 3 loop
        t3(k) <= t2(2*k) + t2(2*k+1);
      end loop;
      for k in 0 to 1 loop
        t4(k) <= t3(2*k) + t3(2*k+1);
      end loop;
      acc <= t4(0) + t4(1);
    end if;
  end process p_tree;

  -- -------------------------------------------------------------------
  -- Round, saturate, pair I with Q, and suppress the 31-sample warm-up
  -- so the emitted stream aligns with the model's centred convolution.
  -- -------------------------------------------------------------------
  p_out : process (clk)
    variable q : signed(W_DAT-1 downto 0);
  begin
    if rising_edge(clk) then
      if rst = '1' then
        out_val <= '0';
        warm    <= (others => (others => '0'));
        hold_re <= (others => '0');
      else
        out_val <= '0';
        if acc_val = '1' then
          q := shift_round_sat(acc, SH_FIR, W_DAT);
          if acc_isq = '0' then
            hold_re <= q;                         -- I arrives first
          else
            if warm(to_integer(acc_ch)) < WARMUP then
              warm(to_integer(acc_ch)) <= warm(to_integer(acc_ch)) + 1;
            else
              out_val <= '1';
              out_re  <= hold_re;
              out_im  <= q;
              out_ch  <= acc_ch;
            end if;
          end if;
        end if;
      end if;
    end if;
  end process p_out;

  m_valid <= out_val;
  m_chan  <= std_logic_vector(out_ch);
  m_re    <= std_logic_vector(out_re);
  m_im    <= std_logic_vector(out_im);

end architecture rtl;
