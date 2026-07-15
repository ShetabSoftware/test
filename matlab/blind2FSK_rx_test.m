clc;
clear;
close all;

%% ====================================================================
%  Blind 2FSK Receiver + Blind Symbol-Rate Estimation
%  (HF Doppler / Watterson channel)
%
%  The blind symbol-rate estimator has been factored out into the reusable
%  function  blindSymbolRateFSK.m  (toolbox-free; runs in MATLAB and Octave).
%  This script keeps the original modulation / channel / BER / Monte-Carlo
%  structure and simply CALLS that function, so you can test the estimator in
%  isolation and tune it via name-value options.
%
%  Quick standalone test of just the estimator:
%     Rs_hat = blindSymbolRateFSK(rx, Fs, Rs_min, Rs_max)
% ====================================================================

% --- Parameters ---
channel_type = 'poor';      % 'excellent', 'good', 'poor'
EbN0dB_vec = 0:2:14;
Rs = 400;                   % 300/600/1200/2400/4800  (300 == Delta f: see note)
Nsym = 1e5;
main_loop = 100;
M = 2;
k = 1;
Fs = 9600;
sps = Fs / Rs;

f1 = 600;                   % Frequency for bit 0 (Hz)
f2 = 900;                   % Frequency for bit 1 (Hz)

% Blind Rs search band (no knowledge of true Rs)
Rs_min = 150;
Rs_max = 4800;

% NOTE on Rs == Delta f: the tone spacing here is Delta f = |f2-f1| = 300 Hz.
% The estimator notches the Delta f "tone-beat" line by default. If you test
% Rs = 300 (i.e. Rs == Delta f) set notch_tone_spacing = false, otherwise the
% true clock line is notched (an inherent blind ambiguity between a 300-baud
% clock and a 300-Hz tone beat).
notch_tone_spacing = (abs(Rs - abs(f2 - f1)) > 1e-6);

nEb = length(EbN0dB_vec);
total_ber = zeros(1, nEb);
Rs_hat_all = zeros(main_loop, nEb);
rel_err_all = zeros(main_loop, nEb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nloop = 1:main_loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Generate random bits ---
dataBits = randi([0 1], Nsym, 1);

% --- 2FSK Modulation (Oversampled) ---
t = (0:Nsym*sps-1) / Fs;
txSig = zeros(1, Nsym * sps);

for i = 1:Nsym
    if dataBits(i) == 0
        txSig((i-1)*sps+1:i*sps) = exp(1j*2*pi*f1*t((i-1)*sps+1:i*sps));
    else
        txSig((i-1)*sps+1:i*sps) = exp(1j*2*pi*f2*t((i-1)*sps+1:i*sps));
    end
end

%% ================================================
% --- HF Channel with Doppler (Watterson model) ---
% =================================================
switch channel_type
    case 'excellent'
        fd = 1;
        pathDelays = [0 2e-3];
        gains_db = [0 -15];
    case 'good'
        fd = 1;
        pathDelays = [0 1.5e-3];
        gains_db = [0 -10];
    case 'poor'
        fd = 10;
        pathDelays = [0 1.5e-3 3.0e-3];
        gains_db = [0 -10 -15];
end
pathDelays_samples = round(pathDelays * Fs);
h_len = max(pathDelays_samples) + 1;
h_initial = zeros(h_len, 1);
h_initial(1) = 1;
for p = 2:numel(pathDelays_samples)
    h_initial(pathDelays_samples(p)+1) = 10^(gains_db(p)/20) * exp(1j*2*pi*rand());
end

rxChan = zeros(size(txSig));
for i = 1:length(txSig)
    phase_factor = exp(1j * 2 * pi * fd * i / Fs);
    h_time = h_initial .* (phase_factor .^ (0:h_len-1)');
    for p = 1:h_len
        if i - p + 1 > 0
            rxChan(i) = rxChan(i) + txSig(i - p + 1) * h_time(p);
        end
    end
end
rxChan = rxChan(1:length(txSig));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- SNR sweep ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BER_mlse = zeros(1, nEb);

fprintf('\n=== Loop %d/%d | Blind 2FSK + Blind Rs Estimation ===\n', nloop, main_loop);
fprintf('True Rs = %d Hz, Fs = %d Hz | search [%d, %d] Hz\n', Rs, Fs, Rs_min, Rs_max);
fprintf('Eb/N0 (dB)    BER         Rs_hat (Hz)    RelErr\n');
fprintf('--------------------------------------------------------\n');

for idx = 1:nEb
    EbN0dB = EbN0dB_vec(idx);

    % --- Add AWGN ---
    rxChanNoisy = awgn(rxChan, EbN0dB, 'measured');

    %% ===== Blind Symbol-Rate Estimation (factored-out function) =====
    Rs_estimated = blindSymbolRateFSK(rxChanNoisy, Fs, Rs_min, Rs_max, ...
        'notchToneSpacing', notch_tone_spacing);

    Rs_hat_all(nloop, idx) = Rs_estimated;
    rel_err_all(nloop, idx) = abs(Rs_estimated - Rs) / Rs;

    %% ===== Instantaneous Frequency Detection (known sps for BER) =====
    rxBits = zeros(Nsym, 1);
    threshold = (f1 + f2) / 2;

    for i = 1:Nsym
        startIdx = (i-1)*sps + 1;
        endIdx = i*sps;
        if endIdx > length(rxChanNoisy)
            break;
        end
        sym = rxChanNoisy(startIdx:endIdx);
        phase_diff = angle(sym(2:end) .* conj(sym(1:end-1)));
        avg_freq = mean(phase_diff * Fs / (2*pi));
        if avg_freq < threshold
            rxBits(i) = 0;
        else
            rxBits(i) = 1;
        end
    end

    %% ===== Phase Ambiguity Resolution =====
    [~, ber1] = biterr(dataBits(1:length(rxBits)), rxBits);
    [~, ber2] = biterr(dataBits(1:length(rxBits)), 1 - rxBits);
    BER_mlse(idx) = min(ber1, ber2);

    fprintf('   %2d        %.5f     %8.2f       %.4f\n', ...
        EbN0dB, BER_mlse(idx), Rs_estimated, rel_err_all(nloop, idx));
end

total_ber = total_ber + BER_mlse;

end %% main loop

%% =======================
% Aggregate Monte Carlo results
%% =======================
avg_ber = total_ber / main_loop;
mean_Rs_hat = mean(Rs_hat_all, 1);
rmse_Rs = sqrt(mean((Rs_hat_all - Rs).^2, 1));
mean_rel_err = mean(rel_err_all, 1);

fprintf('\n=== Monte Carlo Summary (N=%d) | True Rs = %d Hz ===\n', main_loop, Rs);
fprintf('Eb/N0 (dB)    Avg BER     Mean Rs_hat    RMSE (Hz)    Mean RelErr\n');
fprintf('--------------------------------------------------------------------\n');
for idx = 1:nEb
    fprintf('   %2d        %.5f     %8.2f      %8.2f      %.4f\n', ...
        EbN0dB_vec(idx), avg_ber(idx), mean_Rs_hat(idx), rmse_Rs(idx), mean_rel_err(idx));
end

%% =======================
% Plot Results
%% =======================
figure('Position', [100 100 1400 500]);

% 1. BER Plot
subplot(1, 2, 1);
semilogy(EbN0dB_vec, avg_ber, 'bo-', 'LineWidth', 2.5, 'MarkerSize', 10, 'MarkerFaceColor', 'b');
hold on;
EbN0_theory = linspace(0, 16, 100);
EbN0_linear = 10.^(EbN0_theory/10);
BER_nc_fsk = 0.5 * exp(-EbN0_linear / 2);
semilogy(EbN0_theory, BER_nc_fsk, 'r--', 'LineWidth', 2);
grid on;
xlabel('E_b/N_0 (dB)', 'FontSize', 12);
ylabel('Bit Error Rate', 'FontSize', 12);
title('2FSK Receiver BER (known sps)', 'FontSize', 13, 'FontWeight', 'bold');
legend('Simulated 2FSK', 'Theoretical (Non-Coherent)', 'Location', 'southwest');
xlim([0 16]); ylim([1e-5 1]);
set(gca, 'FontSize', 11);

% 2. Blind Rs estimation accuracy
subplot(1, 2, 2);
yyaxis left
plot(EbN0dB_vec, mean_rel_err, 's-', 'LineWidth', 2.5, 'MarkerSize', 8);
ylabel('Mean |R_s^{hat}-R_s|/R_s');
yyaxis right
plot(EbN0dB_vec, rmse_Rs, 'd--', 'LineWidth', 2.5, 'MarkerSize', 8);
ylabel('RMSE (Hz)');
grid on;
xlabel('E_b/N_0 (dB)', 'FontSize', 12);
title(sprintf('Blind R_s Estimation (true R_s = %d Hz)', Rs), 'FontSize', 13, 'FontWeight', 'bold');
legend('Mean Rel. Error', 'RMSE', 'Location', 'northeast');
xlim([0 16]);
set(gca, 'FontSize', 11);

figure;
plot(EbN0dB_vec, mean_Rs_hat, 'ko-', 'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', 'k');
hold on;
plot(EbN0dB_vec([1 end]), [Rs Rs], 'r--', 'LineWidth', 2);
grid on;
xlabel('E_b/N_0 (dB)', 'FontSize', 12);
ylabel('Estimated R_s (Hz)', 'FontSize', 12);
title('Mean Blind Symbol-Rate Estimate vs SNR', 'FontSize', 13, 'FontWeight', 'bold');
legend('Mean R_s^{hat}', sprintf('True R_s = %d Hz', Rs), 'Location', 'best');
xlim([0 16]);
set(gca, 'FontSize', 11);
