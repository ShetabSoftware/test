function oct = asp_ca_first10_octal(prn)
%ASP_CA_FIRST10_OCTAL IS-GPS-200 reference: first 10 chips of each C/A code.
%
%   oct = ASP_CA_FIRST10_OCTAL(PRN) returns the tabulated octal value from
%   IS-GPS-200 Table 3-Ia, column "First 10 Chips (Octal)".
%
%   Encoding convention of that column: the leading digit contributes ONE
%   bit and the remaining three digits contribute three bits each, so the
%   value 1440 octal means the bit pattern 1 100 100 000 = 1100100000, which
%   is the first ten C/A chips of PRN 1 in logic (not +/-1) form.
%
%   Used only by verify/test_ca_code.m.

table = [1440 1620 1710 1744 1133 1455 1131 1454 1626 1504 ...
         1642 1750 1764 1772 1775 1776 1156 1467 1633 1715 ...
         1746 1763 1063 1706 1743 1761 1770 1774 1127 1453 ...
         1625 1712];

if ~isscalar(prn) || prn < 1 || prn > 32
    error('asp_ca_first10_octal:prn', 'PRN must be an integer in 1..32.');
end

oct = table(prn);

end
