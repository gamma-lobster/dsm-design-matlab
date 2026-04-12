%% design_3rd_order_ct_dsm_10mhz.m
%% 3rd-Order Continuous-Time DSM Design - Native MATLAB Workflow
%% Spec: order=3, fs=10MHz, OSR=32, H_inf=4.0, 4-bit, CT feedback form
%%
%% This example starts from the same 3rd-order lowpass NTF used by the
%% discrete-time reference design, then maps that NTF into a continuous-time
%% loop filter with realizeNTF_ct(). The sampled CT loop is converted back
%% to a discrete-time equivalent with mapCtoD() for sample-time validation.

clearvars; close all; clc;

%% Add Delta Sigma Toolbox to path (relative to this script)
script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(fileparts(script_dir)), 'common'));
refs_root = setup_references_paths();
results_dir = fullfile(refs_root, 'results', 'ct', 'base');

fprintf('============================================================\n');
fprintf('  3rd-Order Continuous-Time DSM Design - Native MATLAB\n');
fprintf('  fs=10MHz, OSR=32, H_inf=4.0, 4-bit quantizer, tdac=[0 1]\n');
fprintf('============================================================\n\n');

%% Step 1: Synthesize the target NTF
fprintf('[Step 1] Synthesize Target NTF\n');
fprintf('----------------------------------------\n');

order = 3;
OSR = 32;
H_inf = 4.0;
opt = 1;
f0 = 0;

ntf = synthesizeNTF(order, OSR, opt, H_inf, f0);

fprintf('Target NTF: ');
zpk(ntf)
fprintf('\n');

%% Step 2: Map the NTF to a continuous-time loop filter
fprintf('[Step 2] Map NTF To Continuous-Time Loop Filter\n');
fprintf('----------------------------------------\n');

ct_form = 'FF';
tdac = [0 1];

[ABCDc, tdac2] = realizeNTF_ct(ntf, ct_form, tdac);
[Ac, Bc, Cc, Dc] = partitionABCD(ABCDc);
sys_c = ss(Ac, Bc, Cc, Dc);

fprintf('CT topology: %s\n', ct_form);
fprintf('DAC timing: [%.2f %.2f]\n', tdac(1), tdac(2));
fprintf('ABCDc Matrix (%dx%d):\n', size(ABCDc, 1), size(ABCDc, 2));
for i = 1:size(ABCDc, 1)
    fprintf('  ['); fprintf('%11.6f ', ABCDc(i, :)); fprintf(']\n');
end
fprintf('\n');

%% Step 3: Convert the sampled CT system back to a DT equivalent
fprintf('[Step 3] Build Sampled DT Equivalent\n');
fprintf('----------------------------------------\n');

sys_d = mapCtoD(sys_c, tdac2);
ABCD = [sys_d.a sys_d.b; sys_d.c sys_d.d];
[ntf_eq, stf_eq] = calculateTF(ABCD);
ntf_eq = cancelPZ(ntf_eq);

fprintf('Mapped DT-equivalent ABCD (%dx%d):\n', size(ABCD, 1), size(ABCD, 2));
for i = 1:size(ABCD, 1)
    fprintf('  ['); fprintf('%11.6f ', ABCD(i, :)); fprintf(']\n');
end
fprintf('\n');

%% Step 4: Simulation parameters
fprintf('[Step 4] Simulation Parameters\n');
fprintf('----------------------------------------\n');

fs = 10e6;
fB = fs / (2*OSR);
N = 8192;

n_bits = 4;
n_levels = 2^n_bits;
A_in = 0.5;
f_bin = round(sqrt(1/7) * N / (2*OSR));
u = A_in * sin(2*pi*f_bin*(0:N-1)/N);

fprintf('  fs = %.1f MHz\n', fs/1e6);
fprintf('  OSR = %d (BW = %.3f kHz)\n', OSR, fB/1000);
fprintf('  Quantizer: %d-bit (%d levels)\n', n_bits, n_levels);
fprintf('  Input amplitude: %.2f V\n', A_in);
fprintf('  Signal frequency: %.3f kHz\n\n', f_bin/N * fs / 1000);

%% Step 5: Check the CT pulse-response match
fprintf('[Step 5] Check Pulse-Response Match\n');
fprintf('----------------------------------------\n');

n_imp = 12;
dt_plot = 1/32;
imp_dt = -impL1(ntf, n_imp);
pulse_ct = -pulse(sys_c, [0 0; tdac2(2:end, :)], dt_plot, n_imp);
t_pulse = 0:dt_plot:n_imp;
pulse_ct_samples = -pulse(sys_c, [0 0; tdac2(2:end, :)], 1, n_imp);
pulse_err = max(abs(pulse_ct_samples(:) - imp_dt(:)));

fprintf('  Max sampled pulse/impulse mismatch: %.3e\n\n', pulse_err);

%% Step 6: Sample-time simulation using the DT equivalent
fprintf('[Step 6] Run Sample-Time Validation Simulation\n');
fprintf('----------------------------------------\n');

[v, xn, xmax] = simulateDSM(u, ABCD, n_levels);
[SNR, ENOB, V_out_mag, sig_bin, signal_power, noise_power] = ...
    calculate_snr_metrics(v, N, OSR);

fprintf('  Simulation complete.\n');
fprintf('  Output range: [%.4f, %.4f]\n', min(v), max(v));
fprintf('  State maxima: '); fprintf('%.4f ', xmax(:)); fprintf('\n');
fprintf('  Signal power: %.4e\n', signal_power);
fprintf('  Noise power (excl. harmonics): %.4e\n', noise_power);
fprintf('  SNR: %.2f dB\n', SNR);
fprintf('  ENOB: %.2f bits\n\n', ENOB);

%% Step 7: Amplitude sweep on the mapped equivalent
fprintf('[Step 7] Run Amplitude Sweep\n');
fprintf('----------------------------------------\n');

A_sweep = 0.05:0.05:0.95;
SNR_sweep = zeros(size(A_sweep));
ENOB_sweep = zeros(size(A_sweep));
state_peak_sweep = zeros(size(A_sweep));

for k = 1:length(A_sweep)
    u_test = A_sweep(k) * sin(2*pi*f_bin*(0:N-1)/N);
    [v_test, ~, xmax_test] = simulateDSM(u_test, ABCD, n_levels);
    [SNR_sweep(k), ENOB_sweep(k)] = calculate_snr_metrics(v_test, N, OSR);
    state_peak_sweep(k) = max(abs(xmax_test(:)));
end

[peak_snr, peak_idx] = max(SNR_sweep);
peak_amp = A_sweep(peak_idx);

fprintf('  Peak sweep SNR: %.2f dB at amplitude %.2f V\n\n', peak_snr, peak_amp);

%% Step 8: Generate plots
fprintf('[Step 8] Generate Plots\n');
fprintf('----------------------------------------\n');

freqs_hz = (0:N/2)/N * fs;
f_min = max(fs/N, 1e3);
log_idx = freqs_hz >= f_min;
freqs_log_khz = freqs_hz(log_idx) / 1000;
V_out_dB = 20*log10(V_out_mag(1:N/2+1) + eps);

freqs_ntf_hz = logspace(log10(f_min), log10(fs/2), 2000);
z_eval = exp(2*pi*1j*freqs_ntf_hz/fs);
NTF_target_dB = 20*log10(abs(evalTF(ntf, z_eval)) + eps);
NTF_eq_dB = 20*log10(abs(evalTF(ntf_eq, z_eval)) + eps);

fig1 = figure('Name', '3rd-Order CT DSM Summary');
set(fig1, 'Position', [100 100 1400 900]);

subplot(2, 2, 1);
semilogx(freqs_log_khz, V_out_dB(log_idx), 'b-', 'LineWidth', 0.75);
hold on;
semilogx([fB/1000 fB/1000], [-150 10], 'r--', 'LineWidth', 1);
hold off;
xlabel('Frequency (kHz)');
ylabel('dBFS');
title(sprintf('Mapped DT Output Spectrum (SNR = %.1f dB)', SNR));
legend('Spectrum', 'BW Edge', 'Location', 'best');
grid on;
axis([f_min/1000 fs/2000 -140 10]);

subplot(2, 2, 2);
semilogx(freqs_ntf_hz/1000, NTF_target_dB, 'm-', 'LineWidth', 1.25);
hold on;
semilogx(freqs_ntf_hz/1000, NTF_eq_dB, 'c--', 'LineWidth', 1.1);
semilogx([fB/1000 fB/1000], [-120 40], 'r--', 'LineWidth', 1);
hold off;
xlabel('Frequency (kHz)');
ylabel('Magnitude (dB)');
title('Target NTF vs Mapped DT Equivalent');
legend('Target NTF', 'CT->DT Equivalent', 'BW Edge', 'Location', 'best');
grid on;
axis([f_min/1000 fs/2000 -120 40]);

subplot(2, 2, 3);
stairs(0:n_imp, imp_dt, 'bo-', 'LineWidth', 1);
hold on;
plot(t_pulse, pulse_ct, 'g-', 'LineWidth', 1.1);
hold off;
xlabel('Normalized Time (samples)');
ylabel('Amplitude');
title('Loop Filter Impulse/Pulse Match');
legend('Target DT impulse', 'CT pulse response', 'Location', 'best');
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

plot_path = fullfile(results_dir, 'dsm_3rd_order_ct_10mhz_plots.png');
saveas(fig1, plot_path);
fprintf('Saved: %s\n', plot_path);

%% Step 9: Save data
results_path = fullfile(results_dir, 'dsm_3rd_order_ct_10mhz_results.mat');
save(results_path, 'order', 'OSR', 'H_inf', 'opt', 'f0', 'fs', 'fB', ...
    'n_bits', 'n_levels', 'A_in', 'A_sweep', 'SNR', 'ENOB', ...
    'SNR_sweep', 'ENOB_sweep', 'state_peak_sweep', 'peak_amp', ...
    'peak_snr', 'ntf', 'ntf_eq', 'stf_eq', 'ABCDc', 'ABCD', 'tdac', ...
    'tdac2', 'Ac', 'Bc', 'Cc', 'Dc', 'u', 'v', 'xn', 'xmax', ...
    'V_out_mag', 'sig_bin', 'signal_power', 'noise_power', ...
    'freqs_ntf_hz', 'NTF_target_dB', 'NTF_eq_dB', 'imp_dt', ...
    'pulse_ct', 'pulse_ct_samples', 't_pulse', 'pulse_err');
fprintf('Saved: %s\n', results_path);

if max(abs(xmax(:))) < 100
    status_str = 'STABLE';
else
    status_str = 'UNSTABLE';
end

fprintf('\nSummary:\n');
fprintf('  Status: %s\n', status_str);
fprintf('  CT loop filter realized with realizeNTF_ct() and tdac = [0 1].\n');
fprintf('  Validation here is at the sampling instants via mapCtoD().\n');
fprintf('  Peak sweep SNR: %.2f dB at amplitude %.2f V\n', peak_snr, peak_amp);
fprintf('Done!\n');

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
