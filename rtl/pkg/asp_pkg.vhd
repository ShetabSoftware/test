-- =====================================================================
--  asp_pkg  -  fixed-point primitives shared by every block.
--
--  These functions are the VHDL twins of the MATLAB primitives at the
--  bottom of asp_golden_model.m.  They are the contract: if any one of
--  them differs from its model counterpart by a single LSB, every
--  downstream comparison fails and the failure looks like a datapath bug
--  in whichever block happens to be under test.  They are therefore
--  verified on their own (tb_asp_pkg) before any datapath is simulated.
--
--  All functions are pure and combinational.  None of them infers a
--  latch, a divider, or anything that is not a shift, an add, a compare
--  or a mux.  Where a function is used on a timing-critical path the
--  instantiating entity registers its output; the functions themselves
--  add no registers so the caller keeps control of the pipeline.
--
--  VHDL-93 compatible on purpose: it synthesises identically under
--  Vivado XSynth and simulates under GHDL, so the same source is used
--  for co-simulation against the MATLAB vectors and for implementation.
-- =====================================================================
library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

package asp_pkg is

  -- -------------------------------------------------------------------
  -- Saturating resize.  Hard clamp to [-2^(w-1), 2^(w-1)-1], NEVER
  -- wraparound.  A wrapped sample is a full-scale transient that spreads
  -- across the whole band; a saturated one is merely clipped.
  -- -------------------------------------------------------------------
  function clamp_s (x : signed; w : positive) return signed;

  -- -------------------------------------------------------------------
  -- Convergent (round-half-to-even) arithmetic right shift.
  --   q = round_even(x / 2^sh)
  -- Zero mean error, unlike truncation (-0.5 LSB) or round-half-away
  -- (non-zero for one-sided data).  Costs one OR gate over truncation.
  -- Requires sh < x'length.
  -- -------------------------------------------------------------------
  function conv_round_shr (x : signed; sh : natural) return signed;

  -- -------------------------------------------------------------------
  -- Round half AWAY from zero, arithmetic right shift.  This is MATLAB
  -- round(), and it is used ONLY where the model uses round() rather
  -- than its convRound() - specifically inside bfp_scale and the
  -- reciprocal-sqrt range reduction.  Keeping the two distinct matters:
  -- substituting one for the other passes most vectors and fails on
  -- exact ties, which is the hardest class of mismatch to find later.
  -- -------------------------------------------------------------------
  function round_away_shr (x : signed; sh : natural) return signed;

  -- -------------------------------------------------------------------
  -- Truncate TOWARD ZERO, arithmetic right shift.  This is MATLAB fix(),
  -- and it is what the Newton reciprocal-sqrt iteration uses.  It is not
  -- the same as an arithmetic shift for negative operands - a shift
  -- floors, fix() truncates - and the difference is one LSB, which is
  -- one LSB of a quantity that then multiplies the whole covariance.
  -- -------------------------------------------------------------------
  function trunc_shr (x : signed; sh : natural) return signed;

  -- -------------------------------------------------------------------
  -- Shift (either direction), convergent round, saturate.  The single
  -- most used primitive in the design: every multiply in the datapath is
  -- followed by exactly one of these.
  -- -------------------------------------------------------------------
  function shift_round_sat (x : signed; sh : integer; w : positive) return signed;

  -- -------------------------------------------------------------------
  -- Variable-shift versions of the two above.  The fixed-shift forms
  -- slice the operand, which is minimal logic when the shift is a
  -- constant but is not synthesisable when it is not.  Stage 5 shifts by
  -- (SH_EVD_BASE + k_i + k_j), which is data dependent, so it needs a
  -- barrel shifter and a computed guard/sticky pair instead:
  --     inc = half AND (sticky OR lsb_of_result)
  -- which is exactly round-half-to-even written without slices.
  -- -------------------------------------------------------------------
  function conv_round_shr_v (x : signed; sh : natural) return signed;
  function shift_round_sat_v (x : signed; sh : integer; w : positive) return signed;

  -- -------------------------------------------------------------------
  -- ceil(log2(x)) for x >= 1, i.e. the number of bits needed to hold x.
  -- Used for the block normalisation of the covariance before whitening.
  -- -------------------------------------------------------------------
  function ceil_log2_u (x : unsigned) return natural;

  -- -------------------------------------------------------------------
  -- Number of leading redundant sign bits, i.e. how far x can be shifted
  -- left before it overflows its own width.  Combinational priority
  -- encoder; the basis of every block-floating-point scale in the design
  -- and the reason no normalisation anywhere needs a divider.
  -- -------------------------------------------------------------------
  function lead_sign_count (x : signed) return natural;

  -- -------------------------------------------------------------------
  -- Absolute value with one extra bit, so abs(-2^(n-1)) is representable
  -- and does not silently wrap to itself.
  -- -------------------------------------------------------------------
  function abs_ext (x : signed) return unsigned;

  -- -------------------------------------------------------------------
  -- PORT CONVENTION
  -- ---------------
  -- Multi-element buses on entity ports are FLATTENED std_logic_vector,
  -- never arrays of arrays.  Element k of an M-element bus of width W
  -- occupies bits ((k+1)*W-1 downto k*W).  Two reasons, both practical:
  -- unconstrained array elements need VHDL-2008, which the Vivado IP
  -- packager still handles poorly; and a flat vector is what an AXI4-
  -- Stream TDATA field has to be anyway, so the boundary needs no
  -- conversion layer.  Inside an architecture the data is unpacked into
  -- signed arrays immediately, where it is readable.
  -- -------------------------------------------------------------------

end package asp_pkg;


package body asp_pkg is

  -- -----------------------------------------------------------------
  function clamp_s (x : signed; w : positive) return signed is
    variable xv   : signed(x'length-1 downto 0);
    variable res  : signed(w-1 downto 0);
    variable fits : boolean;
  begin
    xv := x;
    if x'length <= w then
      return resize(xv, w);
    end if;
    -- x is representable in w bits exactly when every bit above w-1
    -- equals the sign bit.
    fits := true;
    for i in w-1 to x'length-1 loop
      if xv(i) /= xv(x'length-1) then
        fits := false;
      end if;
    end loop;
    if fits then
      res := xv(w-1 downto 0);
    elsif xv(x'length-1) = '1' then
      res := (others => '0');            -- most negative: 100..0
      res(w-1) := '1';
    else
      res := (others => '1');            -- most positive: 011..1
      res(w-1) := '0';
    end if;
    return res;
  end function;

  -- -----------------------------------------------------------------
  function conv_round_shr (x : signed; sh : natural) return signed is
    variable xv    : signed(x'length-1 downto 0);
    variable f     : signed(x'length-1 downto 0);
    variable rem_u : unsigned(sh downto 0);
    variable half  : unsigned(sh downto 0);
    variable inc   : boolean;
  begin
    xv := x;
    if sh = 0 then
      return xv;
    end if;
    assert sh < x'length
      report "conv_round_shr: shift >= operand width" severity failure;

    f := shift_right(xv, sh);            -- arithmetic: this is floor()
    rem_u := (others => '0');
    rem_u(sh-1 downto 0) := unsigned(xv(sh-1 downto 0));
    half  := (others => '0');
    half(sh-1) := '1';

    if rem_u > half then
      inc := true;                       -- fraction > 1/2
    elsif rem_u = half then
      inc := (f(0) = '1');               -- exact tie: round to even
    else
      inc := false;
    end if;

    if inc then
      return f + 1;
    else
      return f;
    end if;
  end function;

  -- -----------------------------------------------------------------
  function round_away_shr (x : signed; sh : natural) return signed is
    variable xv   : signed(x'length downto 0);
    variable mag  : unsigned(x'length downto 0);
    variable q    : unsigned(x'length downto 0);
    variable half : unsigned(x'length downto 0);
    variable neg  : boolean;
  begin
    if sh = 0 then
      return x;
    end if;
    xv  := resize(x, x'length+1);
    neg := (xv(xv'high) = '1');
    if neg then
      mag := unsigned(-xv);
    else
      mag := unsigned(xv);
    end if;
    half := (others => '0');
    half(sh-1) := '1';
    q := shift_right(mag + half, sh);    -- floor(|x|/2^sh + 1/2)
    if neg then
      return resize(-signed(q), x'length);
    else
      return resize(signed(q), x'length);
    end if;
  end function;

  -- -----------------------------------------------------------------
  function trunc_shr (x : signed; sh : natural) return signed is
    variable xv : signed(x'length downto 0);
    variable m  : unsigned(x'length downto 0);
  begin
    if sh = 0 then
      return x;
    end if;
    xv := resize(x, x'length+1);
    if xv(xv'high) = '1' then
      m := shift_right(unsigned(-xv), sh);
      return resize(-signed(m), x'length);
    else
      m := shift_right(unsigned(xv), sh);
      return resize(signed(m), x'length);
    end if;
  end function;

  -- -----------------------------------------------------------------
  function shift_round_sat (x : signed; sh : integer; w : positive) return signed is
    variable ext : signed(x'length + (-sh) - 1 downto 0);
  begin
    if sh > 0 then
      return clamp_s(conv_round_shr(x, sh), w);
    elsif sh < 0 then
      ext := shift_left(resize(x, x'length + (-sh)), -sh);
      return clamp_s(ext, w);
    else
      return clamp_s(x, w);
    end if;
  end function;

  -- -----------------------------------------------------------------
  function conv_round_shr_v (x : signed; sh : natural) return signed is
    variable xv     : signed(x'length-1 downto 0);
    variable f      : signed(x'length-1 downto 0);
    variable mask   : unsigned(x'length-1 downto 0);
    variable half   : std_logic;
    variable sticky : std_logic;
  begin
    xv := x;
    if sh = 0 then
      return xv;
    end if;
    f    := shift_right(xv, sh);          -- barrel shifter, arithmetic
    half := xv(sh-1);                     -- indexed read: a multiplexer
    -- sticky = OR of bits (sh-2 .. 0), built with a shift rather than a
    -- variable-width slice
    if sh = 1 then
      sticky := '0';
    else
      mask := shift_left(to_unsigned(1, x'length), sh-1) - 1;
      if (unsigned(xv) and mask) = 0 then
        sticky := '0';
      else
        sticky := '1';
      end if;
    end if;
    if half = '1' and (sticky = '1' or f(0) = '1') then
      return f + 1;
    else
      return f;
    end if;
  end function;

  -- -----------------------------------------------------------------
  function shift_round_sat_v (x : signed; sh : integer; w : positive) return signed is
    variable ext : signed(x'length + 63 downto 0);
  begin
    if sh > 0 then
      return clamp_s(conv_round_shr_v(x, sh), w);
    elsif sh < 0 then
      ext := shift_left(resize(x, ext'length), -sh);
      return clamp_s(ext, w);
    else
      return clamp_s(x, w);
    end if;
  end function;

  -- -----------------------------------------------------------------
  function ceil_log2_u (x : unsigned) return natural is
    variable xv  : unsigned(x'length-1 downto 0);
    variable msb : integer := -1;
    variable pow : boolean;
    variable n   : natural := 0;
  begin
    xv := x;
    for i in 0 to xv'high loop
      if xv(i) = '1' then
        msb := i;
        n   := n + 1;
      end if;
    end loop;
    if msb < 0 then
      return 0;                            -- log2(0) treated as 0
    end if;
    pow := (n = 1);                        -- exact power of two
    if pow then
      return msb;
    else
      return msb + 1;
    end if;
  end function;

  -- -----------------------------------------------------------------
  function lead_sign_count (x : signed) return natural is
    variable xv : signed(x'length-1 downto 0);
    variable n  : natural := 0;
  begin
    xv := x;
    for i in xv'high-1 downto 0 loop
      if xv(i) = xv(xv'high) then
        n := n + 1;
      else
        exit;
      end if;
    end loop;
    return n;
  end function;

  -- -----------------------------------------------------------------
  function abs_ext (x : signed) return unsigned is
    variable xv : signed(x'length downto 0);
  begin
    xv := resize(x, x'length+1);
    if xv(xv'high) = '1' then
      return unsigned(-xv);
    else
      return unsigned(xv);
    end if;
  end function;

end package body asp_pkg;
