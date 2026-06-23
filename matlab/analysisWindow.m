function w = analysisWindow(name, n)
%ANALYSISWINDOW  Construct an analysis window of length N (no toolbox needed).
%   W = ANALYSISWINDOW(NAME, N) returns a length-N window. NAME is one of
%   'none'/'rect'/'boxcar', 'hann', 'hamming' or 'blackman'.

    name = lower(name);
    if n <= 1
        w = ones(1, max(n, 0));
        return;
    end
    m = (0:n-1);
    switch name
        case {'none', 'rect', 'boxcar'}
            w = ones(1, n);
        case 'hann'
            w = 0.5 - 0.5 * cos(2 * pi * m / (n - 1));
        case 'hamming'
            w = 0.54 - 0.46 * cos(2 * pi * m / (n - 1));
        case 'blackman'
            w = 0.42 - 0.5 * cos(2 * pi * m / (n - 1)) ...
                + 0.08 * cos(4 * pi * m / (n - 1));
        otherwise
            error('analysisWindow:name', 'unknown window: %s', name);
    end
end
