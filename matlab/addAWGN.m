function y = addAWGN(iq, snrDb)
%ADDAWGN  Add complex AWGN at a target SNR (over the full sampled bandwidth).
%   Y = ADDAWGN(IQ, SNRDB) adds complex white Gaussian noise to IQ so that the
%   ratio of measured signal power to noise power equals SNRDB decibels.

    x = iq(:).';
    signalPower = mean(abs(x) .^ 2);
    if signalPower <= 0
        y = x;
        return;
    end
    snrLinear = 10 ^ (snrDb / 10);
    noisePower = signalPower / snrLinear;
    noise = sqrt(noisePower / 2) * ...
        (randn(size(x)) + 1i * randn(size(x)));
    y = x + noise;
end
