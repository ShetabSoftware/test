function code = asp_ca_code(prn)
%ASP_CA_CODE Generate one period of a GPS L1 C/A Gold code, +/-1 valued.
%
%   code = ASP_CA_CODE(PRN) returns a 1x1023 vector of +1/-1 for PRN 1..32.
%
%   Mapping: logic 1 -> -1, logic 0 -> +1 (so that XOR becomes multiply).
%
%   Implements IS-GPS-200 section 3.3.2.3:
%       G1: 1 + x^3 + x^10
%       G2: 1 + x^2 + x^3 + x^6 + x^8 + x^9 + x^10
%       C/A(PRN) = G1 xor (G2 delayed by the PRN-specific phase selection)
%   The phase selection is realised, as in the ICD, by tapping two G2 stages
%   and XOR-ing them, which is equivalent to the tabulated code phase delay.
%
%   Correctness of this generator is verified against the IS-GPS-200
%   "first 10 chips (octal)" column by verify/test_ca_code.m.  A GNSS
%   product that ships an unverified code generator will produce a receiver
%   that silently fails on a subset of PRNs; the check costs nothing and
%   must be in the regression suite.
%
%   NOTE ON THE ORIGINAL SCRIPTS
%   ----------------------------
%   GPS_SPOOFER_V9 and FPGA_FixedPoint_NullSteering used
%       prn_codes = 2*randi([0,1],N_auth,K) - 1;
%   i.e. independent random +/-1 sequences at the SAMPLE rate, not Gold
%   codes at the CHIP rate.  Three consequences:
%     1. The spectrum is flat to fs/2 instead of a 1.023 MHz sinc^2, so the
%        modelled front-end bandwidth and decimation behaviour are wrong.
%     2. Cross-correlation between "PRNs" is ~-40 dB (random, length 1e4)
%        instead of the -24 dB worst case of real Gold codes, so the paper's
%        central claim - that a strong spoofer raises the receiver noise
%        floor through cross-correlation - cannot be reproduced at all.
%     3. Code Doppler and chip-boundary effects vanish.

persistent g2Taps
if isempty(g2Taps)
    g2Taps = [ 2  6;  3  7;  4  8;  5  9;  1  9;  2 10;  1  8;  2  9; ...
               3 10;  2  3;  3  4;  5  6;  6  7;  7  8;  8  9;  9 10; ...
               1  4;  2  5;  3  6;  4  7;  5  8;  6  9;  1  3;  4  6; ...
               5  7;  6  8;  7  9;  8 10;  1  6;  2  7;  3  8;  4  9];
end

if ~isscalar(prn) || prn < 1 || prn > 32 || prn ~= floor(prn)
    error('asp_ca_code:prn', 'PRN must be an integer in 1..32.');
end

n = 1023;
g1 = -ones(1,10);   % all registers initialised to logic 1
g2 = -ones(1,10);
code = zeros(1,n);
taps = g2Taps(prn,:);

for k = 1:n
    g1Out = g1(10);
    g2Out = g2(taps(1)) * g2(taps(2));
    code(k) = g1Out * g2Out;

    g1Fb = g1(3) * g1(10);
    g2Fb = g2(2) * g2(3) * g2(6) * g2(8) * g2(9) * g2(10);

    g1 = [g1Fb, g1(1:9)];
    g2 = [g2Fb, g2(1:9)];
end

end
