-- =====================================================================
--  asp_tb_pkg  -  shared helpers for the co-simulation testbenches.
--
--  Every stage testbench in this project follows the same pattern: read
--  the stage input vector the MATLAB golden model dumped, drive it into
--  the DUT, and compare the DUT output against the stage output vector
--  the same model dumped - BIT EXACTLY, with no tolerance.  Any
--  difference is a bug, not noise.  These helpers exist so that pattern
--  is written once rather than ten times.
--
--  FILE FORMAT (identical for every file the model emits)
--    signed decimal integers, one per line;
--    "# ..." lines are comments and are skipped;
--    complex data is interleaved I,Q;
--    multi-channel streaming data is column-major over the N_ANT x Nsamp
--    matrix, i.e. ch0..ch3 for sample 0, then ch0..ch3 for sample 1.
--
--  SIMULATION ONLY.  Nothing here is synthesisable and nothing here is
--  instantiated by the design.
-- =====================================================================
library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

library std;
use std.textio.all;

package asp_tb_pkg is

  -- Read the next integer, skipping blank and '#' comment lines.
  -- ok is false at end of file.
  procedure read_int (file f : text; variable v : out integer; variable ok : out boolean);

  -- Compare and report.  Increments err on mismatch and prints the first
  -- G_MAX_REPORT differences with enough context to locate them.
  procedure check_int (constant name  : in    string;
                       constant idx   : in    integer;
                       constant got   : in    integer;
                       constant want  : in    integer;
                       variable err   : inout integer;
                       variable shown : inout integer);

  procedure summarise (constant name : in string;
                       constant n    : in integer;
                       constant err  : in integer);

  -- -----------------------------------------------------------------
  -- 64-bit variants.  The covariance accumulators are 48 bits and their
  -- values reach ~1.7e11, well past INTEGER'HIGH.  Comparing them as
  -- VHDL integers overflows the simulator, so the wide datapath needs
  -- its own reader, comparator and decimal formatter.
  -- -----------------------------------------------------------------
  subtype big_t is signed(63 downto 0);

  procedure read_big (file f : text; variable v : out big_t; variable ok : out boolean);

  procedure check_big (constant name  : in    string;
                       constant idx   : in    integer;
                       constant got   : in    big_t;
                       constant want  : in    big_t;
                       variable err   : inout integer;
                       variable shown : inout integer);

  function to_dec (x : signed) return string;

  constant TB_MAX_REPORT : integer := 12;

end package asp_tb_pkg;


package body asp_tb_pkg is

  -- Local helper: digits are generated least-significant first.
  function reverse_str (s : string) return string is
    variable r : string(s'range);
  begin
    for i in s'range loop
      r(s'high - i + s'low) := s(i);
    end loop;
    return r;
  end function;

  procedure read_int (file f : text; variable v : out integer; variable ok : out boolean) is
    variable l     : line;
    variable good  : boolean;
    variable first : character;
  begin
    v  := 0;
    ok := false;
    while not endfile(f) loop
      readline(f, l);
      if l /= null then
        if l.all'length > 0 then
          -- Inspect the first non-blank character WITHOUT consuming it,
          -- so that a data line can then be parsed by the standard
          -- TEXTIO integer read rather than by hand.
          first := ' ';
          for i in l.all'range loop
            if l.all(i) /= ' ' and l.all(i) /= HT then
              first := l.all(i);
              exit;
            end if;
          end loop;
          if first /= '#' and first /= ' ' then
            read(l, v, good);
            if good then
              ok := true;
              return;
            end if;
          end if;
        end if;
      end if;
    end loop;
  end procedure read_int;

  procedure check_int (constant name  : in    string;
                       constant idx   : in    integer;
                       constant got   : in    integer;
                       constant want  : in    integer;
                       variable err   : inout integer;
                       variable shown : inout integer) is
  begin
    if got /= want then
      err := err + 1;
      if shown < TB_MAX_REPORT then
        shown := shown + 1;
        report name & ": MISMATCH at index " & integer'image(idx) &
               "  rtl=" & integer'image(got) &
               "  model=" & integer'image(want) &
               "  delta=" & integer'image(got - want)
          severity error;
      end if;
    end if;
  end procedure check_int;

  -- -----------------------------------------------------------------
  function to_dec (x : signed) return string is
    variable v    : signed(x'length downto 0);
    variable neg  : boolean;
    variable buf  : string(1 to 24);
    variable n    : integer := 0;
    variable d    : integer;
    variable ten  : signed(x'length downto 0);
  begin
    v := resize(x, x'length+1);
    if v = 0 then
      return "0";
    end if;
    neg := (v < 0);
    if neg then
      v := -v;
    end if;
    ten := to_signed(10, v'length);
    while v /= 0 loop
      d := to_integer(v mod ten);
      n := n + 1;
      buf(n) := character'val(character'pos('0') + d);
      v := v / ten;
    end loop;
    -- reverse
    if neg then
      return "-" & reverse_str(buf(1 to n));
    else
      return reverse_str(buf(1 to n));
    end if;
  end function;

  -- -----------------------------------------------------------------
  procedure read_big (file f : text; variable v : out big_t; variable ok : out boolean) is
    variable l     : line;
    variable acc   : big_t;
    variable neg   : boolean;
    variable seen  : boolean;
    variable c     : character;
    variable first : character;
  begin
    v  := (others => '0');
    ok := false;
    while not endfile(f) loop
      readline(f, l);
      if l /= null then
        if l.all'length > 0 then
          first := ' ';
          for i in l.all'range loop
            if l.all(i) /= ' ' and l.all(i) /= HT then
              first := l.all(i);
              exit;
            end if;
          end loop;
          if first /= '#' and first /= ' ' then
            acc  := (others => '0');
            neg  := false;
            seen := false;
            for i in l.all'range loop
              c := l.all(i);
              if c = '-' and not seen then
                neg := true;
              elsif c >= '0' and c <= '9' then
                seen := true;
                acc := resize(acc * to_signed(10, 8), 64) +
                       to_signed(character'pos(c) - character'pos('0'), 64);
              elsif seen then
                exit;
              end if;
            end loop;
            if seen then
              if neg then
                v := -acc;
              else
                v := acc;
              end if;
              ok := true;
              return;
            end if;
          end if;
        end if;
      end if;
    end loop;
  end procedure read_big;

  -- -----------------------------------------------------------------
  procedure check_big (constant name  : in    string;
                       constant idx   : in    integer;
                       constant got   : in    big_t;
                       constant want  : in    big_t;
                       variable err   : inout integer;
                       variable shown : inout integer) is
  begin
    if got /= want then
      err := err + 1;
      if shown < TB_MAX_REPORT then
        shown := shown + 1;
        report name & ": MISMATCH at index " & integer'image(idx) &
               "  rtl=" & to_dec(got) &
               "  model=" & to_dec(want) &
               "  delta=" & to_dec(got - want)
          severity error;
      end if;
    end if;
  end procedure check_big;

  procedure summarise (constant name : in string;
                       constant n    : in integer;
                       constant err  : in integer) is
  begin
    if err = 0 then
      report name & ": PASS - " & integer'image(n) &
             " values bit-identical to the MATLAB golden model"
        severity note;
    else
      report name & ": FAIL - " & integer'image(err) & " of " &
             integer'image(n) & " values differ"
        severity failure;
    end if;
  end procedure summarise;

end package body asp_tb_pkg;
