function demoFskBaud()
%DEMOFSKBAUD  End-to-end demonstration of blind FSK symbol-rate estimation.
%   Run with:  demoFskBaud
%
%   Synthesises several FSK signals (different rates, orders and SNRs, with a
%   residual carrier offset) and prints the blindly estimated symbol rate
%   alongside the true value.

    fs = 192000;
    configs = { ...
        struct('symbolRate', 2400, 'order', 2, 'snrDb', []),  ...
        struct('symbolRate', 4800, 'order', 2, 'snrDb', 10),  ...
        struct('symbolRate', 9600, 'order', 4, 'snrDb', 12),  ...
        struct('symbolRate', 1200, 'order', 2, 'snrDb', 8)};

    fprintf('%9s  %5s  %6s  %10s  %9s  %s\n', ...
        'true Rs', 'order', 'SNR', 'estimate', 'error %', 'method');
    fprintf('%s\n', repmat('-', 1, 62));

    for i = 1:numel(configs)
        c = configs{i};
        sig = generateFSK(c.symbolRate, fs, 6000, ...
            'order', c.order, 'snrDb', c.snrDb, ...
            'carrierOffset', 1500, 'seed', 0);
        est = estimateSymbolRate(sig.iq, fs, 'method', 'auto');
        err = (est.symbolRate - sig.symbolRate) / sig.symbolRate * 100;
        if isempty(c.snrDb)
            snrStr = 'clean';
        else
            snrStr = sprintf('%ddB', c.snrDb);
        end
        fprintf('%9.0f  %5d  %6s  %10.2f  %+9.3f  %s\n', ...
            sig.symbolRate, c.order, snrStr, est.symbolRate, err, est.method);
    end
end
