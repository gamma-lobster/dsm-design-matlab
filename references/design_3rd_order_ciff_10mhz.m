%% design_3rd_order_ciff_10mhz.m
%% 3rd-Order CIFF DSM Design - Native MATLAB Workflow
%% Spec: order=3, fs=10MHz, OSR=32, H_inf=4.0, 4-bit, CIFF
%%
%% This example implements a lowpass 3rd-order DSM with optimized zeros
%% and a 4-bit quantizer. No DAC mismatch is modeled.

clear all; close all; clc;

%% Add Delta Sigma Toolbox to path (relative to this script)
addpath(fullfile(fileparts(mfilename('fullpath')), 'dstoolbox'));

fprintf('============================================================\n');
fprintf('  3rd-Order CIFF DSM Design - Native MATLAB Workflow\n');
fprintf('  fs=10MHz, OSR=32, H_inf=4.0, 4-bit quantizer\n');
fprintf('============================================================\n\n');

%% Step 1: Synthesize NTF
fprintf('[Step 1] Synthesize NTF\n');
fprintf('----------------------------------------\n');

order = 3;
OSR = 32;
H_inf = 4.0;
opt = 1;
f0 = 0;

ntf = synthesizeNTF(order, OSR, opt, H_inf, f0);

fprintf('NTF: ');
zpk(ntf)
fprintf('\n');

%% Step 2: Realize Coefficients
fprintf('[Step 2] Realize NTF to Coefficients\n');
fprintf('----------------------------------------\n');

form = 'CIFF';
[a, g, b, c] = realizeNTF(ntf, form);

fprintf('Coefficients:\n');
fprintf('  a = ['); fprintf('%.6f ', a); fprintf('] (feedback)\n');
fprintf('  g = ['); fprintf('%.6f ', g); fprintf('] (resonator)\n');
fprintf('  b = ['); fprintf('%.6f ', b); fprintf('] (input)\n');
fprintf('  c = ['); fprintf('%.6f ', c); fprintf('] (inter-stage)\n\n');

%% Step 3: Build ABCD Matrix
fprintf('[Step 3] Build ABCD Matrix\n');
fprintf('----------------------------------------\n');

ABCD = stuffABCD(a, g, b, c, form);

fprintf('ABCD Matrix (%dx%d):\n', size(ABCD, 1), size(ABCD, 2));
for i = 1:size(ABCD, 1)
    fprintf('  ['); fprintf('%10.6f ', ABCD(i, :)); fprintf(']\n');
end
fprintf('\n');

%% Step 4: Simulation Parameters
fprintf('[Step 4] Simulation Parameters\n');
fprintf('----------------------------------------\n');

fs = 10e6;
fB = fs / (2*OSR);
N = 8192;

n_bits = 4;
n_levels = 2^n_bits;
V_fs = 1.0;
LSB = 2*V_fs / (n_levels - 1);
q_levels = linspace(-V_fs, V_fs, n_levels);

fprintf('  fs = %.1f MHz\n', fs/1e6);
fprintf('  OSR = %d (BW = %.3f kHz)\n', OSR, fB/1000);
fprintf('  Quantizer: %d-bit (%d levels)\n', n_bits, n_levels);
fprintf('  Full scale: +/- %.1f V\n', V_fs);
fprintf('  LSB = %.4f V\n\n', LSB);

%% Step 5: Run Simulation
fprintf('[Step 5] Run Simulation\n');
fprintf('----------------------------------------\n');

f_bin = round(sqrt(1/7) * N / (2*OSR));
A_in = 0.5;

fprintf('  Input: %.2f V sine wave\n', A_in);
fprintf('  FFT bins: %d\n', N);
fprintf('  Signal bin: %d\n', f_bin);
fprintf('  Signal frequency: %.3f kHz\n\n', f_bin/N * fs / 1000);

u = A_in * sin(2*pi*f_bin*(0:N-1)/N);
n = 0:N-1;

n_states = size(ABCD, 1) - 1;
A_mat = ABCD(1:n_states, 1:n_states);
B_mat = ABCD(1:n_states, n_states+1:end);
C_mat = ABCD(n_states+1, 1:n_states);
D_mat = ABCD(n_states+1, n_states+1:end);

x = zeros(n_states, 1);
x_max = zeros(n_states, 1);
y = zeros(1, N);
v = zeros(1, N);

for i = 1:N
    if i == 1
        v_prev = 0;
    else
        v_prev = v(i-1);
    end

    y(i) = C_mat * x + D_mat(1)*u(i) + D_mat(2)*v_prev;

    y_clip = max(-V_fs, min(V_fs, y(i)));
    idx = round((y_clip - (-V_fs)) / LSB) + 1;
    idx = max(1, min(n_levels, idx));
    v(i) = q_levels(idx);

    x = A_mat*x + B_mat(:,1)*u(i) + B_mat(:,2)*v(i);
    x_max = max(x_max, abs(x));
end

fprintf('  Simulation complete.\n');
fprintf('  Output range: [%.4f, %.4f]\n', min(v), max(v));
fprintf('  State maxima: '); fprintf('%.4f ', x_max); fprintf('\n\n');

%% Step 6: Calculate SNR
fprintf('[Step 6] Calculate SNR\n');
fprintf('----------------------------------------\n');

[SNR, ENOB, V_out_mag, sig_bin, signal_power, noise_power] = ...
    calculate_snr_metrics(v, N, OSR);
freqs = (0:N/2)/N*fs/1000;

fprintf('  Signal power: %.4e\n', signal_power);
fprintf('  Noise power (excl. harmonics): %.4e\n', noise_power);
fprintf('  SNR: %.2f dB\n', SNR);
if isfinite(ENOB)
    fprintf('  ENOB: %.2f bits\n', ENOB);
end
fprintf('\n');

%% Step 7: Plots
fprintf('[Step 7] Generate Plots\n');
fprintf('----------------------------------------\n');

freqs_hz = (0:N/2)/N*fs;
f_min = max(fs/N, 1e3);
log_idx = freqs_hz >= f_min;
freqs_log_khz = freqs_hz(log_idx) / 1000;
V_out_dB = 20*log10(V_out_mag(1:N/2+1) + eps);

N_ntf = 8192;
freqs_ntf_hz = logspace(log10(f_min), log10(fs/2), 2000);
z_eval = exp(2*pi*1j*freqs_ntf_hz/fs);
NTF_resp = evalTF(ntf, z_eval);
NTF_dB = 20*log10(abs(NTF_resp) + eps);

fprintf('[Step 8] Run Amplitude Sweep\n');
fprintf('----------------------------------------\n');

A_sweep = 0.05:0.05:0.95;
SNR_sweep = zeros(size(A_sweep));
ENOB_sweep = zeros(size(A_sweep));
state_peak_sweep = zeros(size(A_sweep));

for k = 1:length(A_sweep)
    A_test = A_sweep(k);
    u_test = A_test * sin(2*pi*f_bin*(0:N-1)/N);
    x_test = zeros(n_states, 1);
    x_test_max = zeros(n_states, 1);
    v_test = zeros(1, N);

    for i = 1:N
        if i == 1
            v_prev = 0;
        else
            v_prev = v_test(i-1);
        end

        y_test = C_mat * x_test + D_mat(1)*u_test(i) + D_mat(2)*v_prev;
        y_clip = max(-V_fs, min(V_fs, y_test));
        idx = round((y_clip - (-V_fs)) / LSB) + 1;
        idx = max(1, min(n_levels, idx));
        v_test(i) = q_levels(idx);

        x_test = A_mat*x_test + B_mat(:,1)*u_test(i) + B_mat(:,2)*v_test(i);
        x_test_max = max(x_test_max, abs(x_test));
    end

    [SNR_sweep(k), ENOB_sweep(k)] = calculate_snr_metrics(v_test, N, OSR);
    state_peak_sweep(k) = max(x_test_max);
end

[peak_snr, peak_idx] = max(SNR_sweep);
peak_amp = A_sweep(peak_idx);

fprintf('  Peak sweep SNR: %.2f dB at amplitude %.2f V\n', peak_snr, peak_amp);
fprintf('  Sweep completed across %d amplitudes.\n\n', length(A_sweep));

fig1 = figure('Name', '3rd-Order DSM Summary');
set(fig1, 'Position', [100 100 1350 900]);

subplot(2, 2, 1);
semilogx(freqs_log_khz, V_out_dB(log_idx), 'b-', 'LineWidth', 0.75);
hold on;
semilogx([fB/1000 fB/1000], [-150 10], 'r--', 'LineWidth', 1);
hold off;
xlabel('Frequency (kHz)');
ylabel('dBFS');
title(sprintf('Output Spectrum (log f, SNR = %.1f dB)', SNR));
legend('Spectrum', 'BW Edge', 'Location', 'best');
grid on;
axis([f_min/1000 fs/2000 -140 10]);

subplot(2, 2, 2);
semilogx(freqs_ntf_hz/1000, NTF_dB, 'm-', 'LineWidth', 1.25);
hold on;
semilogx([fB/1000 fB/1000], [-120 40], 'r--', 'LineWidth', 1);
hold off;
xlabel('Frequency (kHz)');
ylabel('Magnitude (dB)');
title('NTF Magnitude (log f)');
legend('NTF', 'BW Edge', 'Location', 'best');
grid on;
axis([f_min/1000 fs/2000 -120 40]);

subplot(2, 2, 3);
stairs((0:399)/fs*1e6, v(1:400), 'b-', 'LineWidth', 1);
hold on;
plot((0:399)/fs*1e6, u(1:400), 'g-', 'LineWidth', 1.25);
hold off;
xlabel('Time (\mus)');
ylabel('Amplitude');
title('Input And Quantized Output');
legend('DSM Output', 'Input', 'Location', 'best');
grid on;

subplot(2, 2, 4);
plot(A_sweep, SNR_sweep, 'bo-', 'LineWidth', 1.25, 'MarkerSize', 4);
hold on;
plot(peak_amp, peak_snr, 'rs', 'MarkerSize', 8, 'LineWidth', 1.5);
hold off;
xlabel('Input Amplitude (V)');
ylabel('SNR (dB)');
title('Amplitude Sweep');
legend('Sweep', sprintf('Peak: %.2f V, %.1f dB', peak_amp, peak_snr), ...
    'Location', 'best');
grid on;

saveas(fig1, 'dsm_3rd_order_4bit_osr32_plots.png');
fprintf('Saved: dsm_3rd_order_4bit_osr32_plots.png\n');

%% Step 9: Save Data
save('dsm_3rd_order_4bit_osr32_results.mat', 'a', 'g', 'b', 'c', 'ABCD', ...
     'SNR', 'ENOB', 'V_out_mag', 'freqs', 'u', 'v', 'x_max', ...
     'order', 'OSR', 'H_inf', 'opt', 'form', 'fs', 'n_bits', 'A_in', ...
     'A_sweep', 'SNR_sweep', 'ENOB_sweep', 'state_peak_sweep', ...
     'peak_amp', 'peak_snr', 'freqs_ntf_hz', 'NTF_dB');
fprintf('Saved: dsm_3rd_order_4bit_osr32_results.mat\n');

if max(x_max) < 100
    status_str = 'STABLE';
else
    status_str = 'UNSTABLE';
end

fprintf('\nSummary:\n');
fprintf('  Status: %s\n', status_str);
fprintf('  No DAC mismatch modeled.\n');
fprintf('  Peak sweep SNR: %.2f dB at amplitude %.2f V\n', peak_snr, peak_amp);
fprintf('Done!\n');

function [SNR, ENOB, V_out_mag, sig_bin, signal_power, noise_power] = calculate_snr_metrics(v, N, OSR)
w = 0.5 * (1 - cos(2*pi*(0:N-1)/N));
V_out = fft(v .* w) / (N/4);
V_out_mag = abs(V_out);

fB_bins = ceil(N / (2*OSR));
[~, sig_idx] = max(V_out_mag(2:fB_bins));
sig_bin = sig_idx + 1;

sig_bins = sig_bin-1:sig_bin+1;
sig_bins = sig_bins(sig_bins >= 2 & sig_bins <= fB_bins);
signal_power = sum(V_out_mag(sig_bins).^2);

harmonic_bins = [];
for h = 2:7
    harmonic_bin = h * sig_bin;
    if harmonic_bin <= fB_bins
        harmonic_bins = [harmonic_bins, harmonic_bin-1:harmonic_bin+1];
    end
end
harmonic_bins = unique(harmonic_bins);
harmonic_bins = harmonic_bins(harmonic_bins >= 2 & harmonic_bins <= fB_bins);

exclude_bins = unique([sig_bins, harmonic_bins]);
noise_bins = setdiff(2:fB_bins, exclude_bins);
noise_power = sum(V_out_mag(noise_bins).^2);

if noise_power > 0
    SNR = 10*log10(signal_power / noise_power);
    ENOB = (SNR - 1.76) / 6.02;
else
    SNR = inf;
    ENOB = inf;
end
end
