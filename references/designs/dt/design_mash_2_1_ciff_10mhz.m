%% design_mash_2_1_ciff_10mhz.m
% Discrete-time 2-1 MASH DSM using a 2nd-order CIFF first stage and a
% 1st-order residue-shaping second stage. Both stages use 4-bit quantizers.

clearvars; close all; clc;

script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(fileparts(script_dir)), 'common'));
refs_root = setup_references_paths();
results_dir = fullfile(refs_root, 'results', 'dt');

fprintf('============================================================\n');
fprintf('  2-1 MASH DT DSM Design\n');
fprintf('  Stage 1: 2nd-order CIFF, Stage 2: 1st-order residue shaper\n');
fprintf('  fs=10MHz, OSR=32, 4-bit quantizer in both stages\n');
fprintf('============================================================\n\n');

%% Step 1: Stage NTFs
fprintf('[Step 1] Synthesize stage NTFs\n');
fprintf('----------------------------------------\n');

fs = 10e6;
OSR = 32;
f0 = 0;

stage1 = struct('order', 2, 'H_inf', 2.0, 'opt', 0, 'form', 'CIFF');
stage2 = struct('order', 1, 'H_inf', 1.5, 'opt', 0, 'form', 'CIFF');

stage1.ntf = synthesizeNTF(stage1.order, OSR, stage1.opt, stage1.H_inf, f0);
stage2.ntf = synthesizeNTF(stage2.order, OSR, stage2.opt, stage2.H_inf, f0);

[stage1.a, stage1.g, stage1.b, stage1.c] = realizeNTF(stage1.ntf, stage1.form);
[stage2.a, stage2.g, stage2.b, stage2.c] = realizeNTF(stage2.ntf, stage2.form);

stage1.ABCD = stuffABCD(stage1.a, stage1.g, stage1.b, stage1.c, stage1.form);
stage2.ABCD = stuffABCD(stage2.a, stage2.g, stage2.b, stage2.c, stage2.form);

fprintf('  Stage 1 NTF zeros at z=1 for textbook-style MASH cancellation.\n');
fprintf('  Stage 1 a = ['); fprintf('%.6f ', stage1.a); fprintf(']\n');
fprintf('  Stage 2 a = ['); fprintf('%.6f ', stage2.a); fprintf(']\n\n');

%% Step 2: Simulation setup
fprintf('[Step 2] Simulation parameters\n');
fprintf('----------------------------------------\n');

N = 8192;
n = 0:N-1;
f_bin = round(sqrt(1/7) * N / (2*OSR));
f_in = f_bin / N * fs;
A_in = 0.5;

n_bits = 4;
n_levels = 2^n_bits;
V_fs = 1.0;
q_levels = linspace(-V_fs, V_fs, n_levels);
lsb = 2 * V_fs / (n_levels - 1);

fprintf('  Input amplitude: %.2f V\n', A_in);
fprintf('  Input frequency: %.3f kHz (bin %d)\n', f_in/1e3, f_bin);
fprintf('  Quantizer levels: %d (%d-bit)\n\n', n_levels, n_bits);

u = A_in * sin(2*pi*f_bin*n/N);

%% Step 3: Run MASH simulation
fprintf('[Step 3] Run MASH simulation\n');
fprintf('----------------------------------------\n');

sim_data = simulate_mash_21(u, stage1.ABCD, stage2.ABCD, V_fs, n_levels, q_levels, lsb);

[cancel_num, cancel_den] = tfdata(stage1.ntf, 'v');
v_mash = sim_data.v1 + filter(cancel_num, cancel_den, sim_data.v2);

[SNR, ENOB, V_out_mag, sig_bin, signal_power, noise_power] = ...
    calculate_snr_metrics(v_mash, N, OSR);

fprintf('  Combined output range: [%.4f, %.4f]\n', min(v_mash), max(v_mash));
fprintf('  Stage 1 peak states: '); fprintf('%.4f ', sim_data.x1_max); fprintf('\n');
fprintf('  Stage 2 peak state: %.4f\n', sim_data.x2_max);
fprintf('  SNR: %.2f dB\n', SNR);
fprintf('  ENOB: %.2f bits\n\n', ENOB);

%% Step 4: Amplitude sweep
fprintf('[Step 4] Run amplitude sweep\n');
fprintf('----------------------------------------\n');

A_sweep = 0.05:0.05:0.95;
SNR_sweep = zeros(size(A_sweep));
ENOB_sweep = zeros(size(A_sweep));
stage1_peak_sweep = zeros(size(A_sweep));
stage2_peak_sweep = zeros(size(A_sweep));

for k = 1:numel(A_sweep)
    u_test = A_sweep(k) * sin(2*pi*f_bin*n/N);
    sim_test = simulate_mash_21(u_test, stage1.ABCD, stage2.ABCD, V_fs, n_levels, q_levels, lsb);
    v_test = sim_test.v1 + filter(cancel_num, cancel_den, sim_test.v2);
    [SNR_sweep(k), ENOB_sweep(k)] = calculate_snr_metrics(v_test, N, OSR);
    stage1_peak_sweep(k) = max(sim_test.x1_max);
    stage2_peak_sweep(k) = sim_test.x2_max;
end

[peak_snr, peak_idx] = max(SNR_sweep);
peak_amp = A_sweep(peak_idx);

fprintf('  Peak sweep SNR: %.2f dB at %.2f V\n\n', peak_snr, peak_amp);

%% Step 5: Plots
fprintf('[Step 5] Generate plots\n');
fprintf('----------------------------------------\n');

fB = fs / (2 * OSR);
freqs_hz = (0:N/2) / N * fs;
f_min = max(fs / N, 1e3);
log_idx = freqs_hz >= f_min;
V_out_dB = 20 * log10(V_out_mag(1:N/2+1) + eps);

freqs_ntf_hz = logspace(log10(f_min), log10(fs/2), 2000);
z_eval = exp(2*pi*1j*freqs_ntf_hz/fs);
NTF1_dB = 20 * log10(abs(evalTF(stage1.ntf, z_eval)) + eps);
NTF2_dB = 20 * log10(abs(evalTF(stage2.ntf, z_eval)) + eps);
NTF_mash_dB = 20 * log10(abs(evalTF(stage1.ntf, z_eval) .* evalTF(stage2.ntf, z_eval)) + eps);

fig1 = figure('Name', '2-1 MASH DSM Summary');
set(fig1, 'Position', [100 100 1400 900]);

subplot(2, 2, 1);
semilogx(freqs_hz(log_idx)/1e3, V_out_dB(log_idx), 'b-', 'LineWidth', 0.75);
hold on;
semilogx([fB/1e3 fB/1e3], [-150 10], 'r--', 'LineWidth', 1);
hold off;
xlabel('Frequency (kHz)');
ylabel('dBFS');
title(sprintf('MASH Output Spectrum (SNR = %.1f dB)', SNR));
legend('Spectrum', 'BW Edge', 'Location', 'best');
grid on;
axis([f_min/1e3 fs/2e3 -140 10]);

subplot(2, 2, 2);
semilogx(freqs_ntf_hz/1e3, NTF1_dB, 'm-', 'LineWidth', 1.0);
hold on;
semilogx(freqs_ntf_hz/1e3, NTF2_dB, 'c-', 'LineWidth', 1.0);
semilogx(freqs_ntf_hz/1e3, NTF_mash_dB, 'k-', 'LineWidth', 1.5);
semilogx([fB/1e3 fB/1e3], [-120 50], 'r--', 'LineWidth', 1);
hold off;
xlabel('Frequency (kHz)');
ylabel('Magnitude (dB)');
title('Stage And Combined NTFs');
legend('Stage 1', 'Stage 2', 'Combined', 'BW Edge', 'Location', 'best');
grid on;
axis([f_min/1e3 fs/2e3 -120 50]);

subplot(2, 2, 3);
n_plot = min(350, N);
stairs((0:n_plot-1)/fs*1e6, sim_data.v1(1:n_plot), 'b-', 'LineWidth', 1);
hold on;
plot((0:n_plot-1)/fs*1e6, sim_data.residue1(1:n_plot), 'r-', 'LineWidth', 1);
stairs((0:n_plot-1)/fs*1e6, sim_data.v2(1:n_plot), 'k-', 'LineWidth', 1);
plot((0:n_plot-1)/fs*1e6, u(1:n_plot), 'g-', 'LineWidth', 1.25);
hold off;
xlabel('Time (\mus)');
ylabel('Amplitude');
title('Stage Signals');
legend('v_1', 'Residue y_1-v_1', 'v_2', 'Input', 'Location', 'best');
grid on;

subplot(2, 2, 4);
plot(A_sweep, SNR_sweep, 'bo-', 'LineWidth', 1.25, 'MarkerSize', 4);
hold on;
plot(peak_amp, peak_snr, 'rs', 'MarkerSize', 8, 'LineWidth', 1.5);
hold off;
xlabel('Input Amplitude (V)');
ylabel('SNR (dB)');
title('Amplitude Sweep');
legend('Sweep', sprintf('Peak: %.2f V, %.1f dB', peak_amp, peak_snr), 'Location', 'best');
grid on;

plot_path = fullfile(results_dir, 'dsm_mash_2_1_ciff_10mhz_plots.png');
saveas(fig1, plot_path);

results_path = fullfile(results_dir, 'dsm_mash_2_1_ciff_10mhz_results.mat');
save(results_path, 'stage1', 'stage2', 'cancel_num', 'cancel_den', 'fs', 'OSR', ...
    'N', 'A_in', 'f_in', 'u', 'sim_data', 'v_mash', 'SNR', 'ENOB', ...
    'V_out_mag', 'sig_bin', 'signal_power', 'noise_power', 'A_sweep', ...
    'SNR_sweep', 'ENOB_sweep', 'stage1_peak_sweep', 'stage2_peak_sweep', ...
    'peak_snr', 'peak_amp', 'freqs_ntf_hz', 'NTF1_dB', 'NTF2_dB', 'NTF_mash_dB');

fprintf('Saved: %s\n', plot_path);
fprintf('Saved: %s\n', results_path);

function sim_data = simulate_mash_21(u, ABCD1, ABCD2, V_fs, n_levels, q_levels, lsb)
n_samples = numel(u);

n1 = size(ABCD1, 1) - 1;
A1 = ABCD1(1:n1, 1:n1);
B1 = ABCD1(1:n1, n1+1:end);
C1 = ABCD1(n1+1, 1:n1);
D1 = ABCD1(n1+1, n1+1:end);

n2 = size(ABCD2, 1) - 1;
A2 = ABCD2(1:n2, 1:n2);
B2 = ABCD2(1:n2, n2+1:end);
C2 = ABCD2(n2+1, 1:n2);
D2 = ABCD2(n2+1, n2+1:end);

x1 = zeros(n1, 1);
x2 = zeros(n2, 1);

y1 = zeros(1, n_samples);
v1 = zeros(1, n_samples);
residue1 = zeros(1, n_samples);
y2 = zeros(1, n_samples);
v2 = zeros(1, n_samples);
x1_max = zeros(n1, 1);
x2_max = zeros(n2, 1);

for i = 1:n_samples
    y1(i) = C1 * x1 + D1(1) * u(i);
    v1(i) = quantize_level(y1(i), V_fs, n_levels, q_levels, lsb);
    residue1(i) = y1(i) - v1(i);

    y2(i) = C2 * x2 + D2(1) * residue1(i);
    v2(i) = quantize_level(y2(i), V_fs, n_levels, q_levels, lsb);

    x1 = A1 * x1 + B1(:, 1) * u(i) + B1(:, 2) * v1(i);
    x2 = A2 * x2 + B2(:, 1) * residue1(i) + B2(:, 2) * v2(i);

    x1_max = max(x1_max, abs(x1));
    x2_max = max(x2_max, abs(x2));
end

sim_data = struct( ...
    'y1', y1, ...
    'v1', v1, ...
    'residue1', residue1, ...
    'y2', y2, ...
    'v2', v2, ...
    'x1_max', x1_max, ...
    'x2_max', max(x2_max));
end

function v = quantize_level(y, V_fs, n_levels, q_levels, lsb)
y_clip = min(max(y, -V_fs), V_fs);
idx = round((y_clip + V_fs) / lsb) + 1;
idx = min(max(idx, 1), n_levels);
v = q_levels(idx);
end

function [SNR, ENOB, V_out_mag, sig_bin, signal_power, noise_power] = calculate_snr_metrics(v, N, OSR)
w = 0.5 * (1 - cos(2*pi*(0:N-1)/N));
V_out = fft(v(:).' .* w) / (N/4);
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
        harmonic_bins = [harmonic_bins, harmonic_bin-1:harmonic_bin+1]; %#ok<AGROW>
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
