function opts = parseOptions(defaults, args)
%PARSEOPTIONS  Simple name-value option parser.
%   OPTS = PARSEOPTIONS(DEFAULTS, ARGS) starts from the struct DEFAULTS and
%   overrides any fields named in the cell array ARGS of name-value pairs.
%   Unknown names raise an error.

    opts = defaults;
    if isempty(args)
        return;
    end
    if mod(numel(args), 2) ~= 0
        error('parseOptions:pairs', 'options must be name-value pairs');
    end
    for i = 1:2:numel(args)
        name = args{i};
        if ~(ischar(name) || (exist('isstring', 'builtin') && isstring(name)))
            error('parseOptions:name', 'option name must be a string');
        end
        name = char(name);
        if ~isfield(opts, name)
            error('parseOptions:unknown', 'unknown option: %s', name);
        end
        opts.(name) = args{i + 1};
    end
end
