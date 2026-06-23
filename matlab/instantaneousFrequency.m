function f = instantaneousFrequency(iq, fs)
%INSTANTANEOUSFREQUENCY  Instantaneous frequency (Hz) of a complex signal.
%   F = INSTANTANEOUSFREQUENCY(IQ, FS) returns the instantaneous frequency of
%   the complex baseband signal IQ sampled at FS hertz. The result is the time
%   derivative of the phase, approximated by a first difference, so for an input
%   of length N the output has length N-1.
%
%   The frequency is obtained from the angle of x[n]*conj(x[n-1]), which is
%   numerically equivalent to diff(unwrap(angle(x))) but robust to large hops.

    if fs <= 0
        error('instantaneousFrequency:fs', 'fs must be positive');
    end
    x = iq(:).';
    if numel(x) < 2
        error('instantaneousFrequency:len', 'at least two samples are required');
    end
    prod = x(2:end) .* conj(x(1:end-1));
    f = angle(prod) * (fs / (2 * pi));
end
