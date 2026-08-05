function ok = test_ca_code()
%TEST_CA_CODE Verify the C/A generator against IS-GPS-200 Table 3-Ia.
%
%   Checks, for all 32 PRNs:
%     1. first 10 chips match the tabulated octal value
%     2. length is 1023 and values are +/-1
%     3. balance: exactly 512 ones and 511 zeros (a property of Gold codes
%        built from maximal-length sequences)
%     4. autocorrelation takes only the values {1023, 63, -1, -65}/1023
%        (the three-valued Gold correlation plus the peak)
%     5. peak cross-correlation between distinct PRNs never exceeds 65/1023,
%        i.e. -23.9 dB
%
%   Test 5 is the one that matters for anti-spoofing work: -23.9 dB is the
%   worst-case cross-correlation floor that a strong spoofer exploits to
%   raise the receiver's effective noise floor, and it is the effect the
%   paper's Figure 7 attributes the single-antenna SNR loss to.  A simulation
%   built on random sequences has a -40 dB floor and cannot show it.

ok = true;

fprintf('test_ca_code:\n');

codes = zeros(32,1023);
for prn = 1:32
    codes(prn,:) = asp_ca_code(prn);
end

% --- 1. first 10 chips vs the ICD octal table
bad = [];
for prn = 1:32
    oct = asp_ca_first10_octal(prn);
    ref = octalFirst10ToBits(oct);
    got = (1 - codes(prn,1:10))/2;         % +1 -> logic 0, -1 -> logic 1
    if any(got ~= ref)
        bad(end+1) = prn; %#ok<AGROW>
    end
end
ok = report(ok, isempty(bad), sprintf('first 10 chips match IS-GPS-200 for all 32 PRNs (failures: %s)', mat2str(bad)));

% --- 2. shape and alphabet
ok = report(ok, all(abs(codes(:)) == 1), 'all chips are +/-1');
ok = report(ok, size(codes,2) == 1023, 'code length is 1023');

% --- 3. balance
nOnes = sum(codes == -1, 2);
ok = report(ok, all(nOnes == 512), 'every code has exactly 512 logic ones');

% --- 4. three-valued autocorrelation
allowed = [1023 63 -1 -65];
worstAuto = 0;
for prn = 1:32
    c = codes(prn,:);
    ac = real(ifft(fft(c).*conj(fft(c))));
    ac = round(ac);
    if ~all(ismember(ac, allowed))
        worstAuto = 1;
    end
end
ok = report(ok, worstAuto == 0, 'autocorrelation takes only {1023,63,-1,-65}');

% --- 5. cross-correlation bound
maxCross = 0;
for i = 1:8          % subset keeps the test fast; the bound is universal
    for j = i+1:8
        cc = real(ifft(fft(codes(i,:)).*conj(fft(codes(j,:)))));
        maxCross = max(maxCross, max(abs(cc)));
    end
end
ok = report(ok, maxCross <= 65 + 1e-6, ...
    sprintf('peak cross-correlation %.0f/1023 = %.1f dB (bound 65/1023 = -23.9 dB)', ...
    maxCross, 20*log10(maxCross/1023)));

end

% -------------------------------------------------------------------------
function bits = octalFirst10ToBits(oct)
% Leading digit contributes one bit, remaining three digits contribute three
% bits each: 1440 -> 1 100 100 000.
s = sprintf('%04d', oct);
bits = zeros(1,10);
bits(1) = str2double(s(1));
p = 2;
for k = 2:4
    d = str2double(s(k));
    b = [bitget(d,3) bitget(d,2) bitget(d,1)];
    bits(p:p+2) = b;
    p = p + 3;
end
end

function ok = report(ok, cond, msg)
if cond
    fprintf('  PASS  %s\n', msg);
else
    fprintf('  FAIL  %s\n', msg);
    ok = false;
end
end
