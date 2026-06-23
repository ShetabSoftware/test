function idx = localMaxima(y)
%LOCALMAXIMA  Indices of (non-strict) local maxima of a 1-D array.
%   IDX = LOCALMAXIMA(Y) returns the 1-based indices of the local maxima of Y.

    y = y(:).';
    n = numel(y);
    if n < 3
        idx = 1:n;
        return;
    end
    greaterLeft  = y(2:end-1) >= y(1:end-2);
    greaterRight = y(2:end-1) >= y(3:end);
    idx = find(greaterLeft & greaterRight) + 1;
end
