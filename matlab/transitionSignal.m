function trans = transitionSignal(iq, fs, transform)
%TRANSITIONSIGNAL  Symbol-transition pulse train from instantaneous frequency.
%   TRANS = TRANSITIONSIGNAL(IQ, FS, TRANSFORM) differences the mean-removed
%   instantaneous frequency once (emphasising transitions and removing any
%   constant carrier offset) and passes it through a memoryless nonlinearity:
%   'square' (default) or 'abs'. The DC component of the result is removed.

    if nargin < 3 || isempty(transform)
        transform = 'square';
    end
    f = instantaneousFrequency(iq, fs);
    f = f - mean(f);
    d = diff(f);
    switch lower(transform)
        case 'square'
            trans = d .^ 2;
        case 'abs'
            trans = abs(d);
        otherwise
            error('transitionSignal:transform', ...
                'transform must be ''square'' or ''abs''');
    end
    trans = trans - mean(trans);
end
