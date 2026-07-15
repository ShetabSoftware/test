function testBlindRsHF()
%TESTBLINDRSHF  Validate blindSymbolRateFSK on an HF 2-FSK Watterson channel.
%   Reproduces the user's scenario (Fs=9600, tones 600/900 Hz, HF multipath +
%   Doppler) and reports blind symbol-rate accuracy vs Eb/N0 for several true
%   symbol rates -- including Rs=300 which equals Delta f (=|f2-f1|), the case
%   that breaks a "notch Delta f" estimator.

    Fs = 9600;
    f1 = 600; f2 = 900;
    Rs_min = 150; Rs_max = 4800;
    channel_type = 'poor';
    Nsym = 4000;            % smaller than 1e5 for a fast self-test
    EbN0dB_vec = 0:2:14;
    trials = 8;

    RsList = [300 400 600 1200 2400];   % 300 == Delta f (the hard case)

    fprintf('HF 2-FSK blind symbol-rate test  (Fs=%d, f=[%d %d], %s channel)\n', ...
        Fs, f1, f2, channel_type);
    fprintf('%-8s %-8s | mean Rs_hat   RMSE(Hz)   meanRelErr\n', 'true Rs', 'Eb/N0');
    fprintf('%s\n', repmat('-', 1, 60));

    worst = 0;
    for Rs = RsList
        sps = Fs / Rs;
        for EbN0dB = EbN0dB_vec
            est = zeros(1, trials);
            for tr = 1:trials
                bits = randi([0 1], Nsym, 1);
                tx = modulate2FSK(bits, sps, Fs, f1, f2);
                rxCh = watterson(tx, Fs, channel_type);
                rx = awgnMeasured(rxCh, EbN0dB);
                est(tr) = blindSymbolRateFSK(rx, Fs, Rs_min, Rs_max);
            end
            relErr = mean(abs(est - Rs) / Rs);
            rmse = sqrt(mean((est - Rs) .^ 2));
            if EbN0dB == 0 || EbN0dB == 14
                fprintf('%-8d %-8d | %10.2f   %8.2f   %8.4f\n', ...
                    Rs, EbN0dB, mean(est), rmse, relErr);
            end
            worst = max(worst, relErr);
        end
        fprintf('%s\n', repmat('-', 1, 60));
    end
    fprintf('worst-case mean relative error across all cases: %.4f\n', worst);
end

function tx = modulate2FSK(bits, sps, Fs, f1, f2)
    Nsym = numel(bits);
    t = (0:Nsym * sps - 1) / Fs;
    fsym = f1 * ones(1, Nsym); fsym(bits == 1) = f2;
    fsamp = reshape(repmat(fsym, sps, 1), 1, []);
    tx = exp(1j * 2 * pi * fsamp .* t);   % (matches the reference generator)
end

function rx = watterson(tx, Fs, channel_type)
    switch channel_type
        case 'excellent'
            fd = 1;  pathDelays = [0 2.0e-3];      gains_db = [0 -15];
        case 'good'
            fd = 1;  pathDelays = [0 1.5e-3];      gains_db = [0 -10];
        otherwise    % 'poor'
            fd = 10; pathDelays = [0 1.5e-3 3.0e-3]; gains_db = [0 -10 -15];
    end
    d = round(pathDelays * Fs);
    hlen = max(d) + 1;
    h0 = zeros(hlen, 1); h0(1) = 1;
    for p = 2:numel(d)
        h0(d(p) + 1) = 10 ^ (gains_db(p) / 20) * exp(1j * 2 * pi * rand());
    end
    n = numel(tx);
    rx = zeros(1, n);
    idx = 0:n - 1;
    for p = 1:hlen
        if h0(p) == 0, continue; end
        dd = p - 1;
        modu = h0(p) * exp(1j * 2 * pi * fd * idx * dd / Fs);
        shifted = [zeros(1, dd), tx(1:n - dd)];
        rx = rx + shifted .* modu;
    end
end

function y = awgnMeasured(x, snrDb)
    P = mean(abs(x) .^ 2);
    N = P / 10 ^ (snrDb / 10);
    y = x + sqrt(N / 2) * (randn(size(x)) + 1j * randn(size(x)));
end
