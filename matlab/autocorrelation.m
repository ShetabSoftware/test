function acf = autocorrelation(x, maxLag)
%AUTOCORRELATION  Biased autocorrelation of a real sequence via the FFT.
%   ACF = AUTOCORRELATION(X) returns the non-negative-lag part of the biased
%   autocorrelation of X, starting at lag 0 (so ACF(1) is the zero-lag value
%   and ACF(TAU+1) is the autocorrelation at lag TAU).
%
%   ACF = AUTOCORRELATION(X, MAXLAG) truncates the result to lags 0..MAXLAG.

    x = x(:).';
    x = x - mean(x);
    n = numel(x);
    if n == 0
        acf = [];
        return;
    end
    nfft = 2 ^ nextpow2(2 * n);
    X = fft(x, nfft);
    ac = real(ifft(X .* conj(X)));
    acf = ac(1:n);
    if nargin > 1 && ~isempty(maxLag)
        acf = acf(1:min(maxLag + 1, n));
    end
end
