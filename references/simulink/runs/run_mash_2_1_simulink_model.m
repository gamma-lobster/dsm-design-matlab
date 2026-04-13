%% run_mash_2_1_simulink_model.m
% Build and run the DT 2-1 MASH DSM Simulink model.

clearvars; close all; clc;

script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(fileparts(script_dir)), 'common'));
refs_root = setup_references_paths();
results_dir = fullfile(refs_root, 'results', 'simulink', 'dt');

fs = 10e6;
OSR = 32;
f0 = 0;
V_fs = 1.0;
A_in = 0.5;
N = 8192;
stop_time = (N - 1) / fs;

stage1_ntf = synthesizeNTF(2, OSR, 0, 2.0, f0);
stage2_ntf = synthesizeNTF(1, OSR, 0, 1.5, f0);

[model_name, model_path, info] = build_mash_dsm_simulink_model( ...
    'model_name', 'dsm_mash_2_1_ciff_10mhz_topology_model', ...
    'output_dir', results_dir, ...
    'stages', [ ...
        struct('name', 'Stage1', 'ntf', stage1_ntf, 'form', 'CIFF', 'n_bits', 4, 'source', 'external'), ...
        struct('name', 'Stage2', 'ntf', stage2_ntf, 'form', 'CIFF', 'n_bits', 4, 'source', 'residue1')], ...
    'fs', fs, ...
    'OSR', OSR, ...
    'V_fs', V_fs, ...
    'input_amplitude', A_in, ...
    'stop_time', stop_time, ...
    'open_model', false);

fprintf('Running Simulink model: %s\n', model_name);
load_system(model_path);
simOut = sim(model_name, 'StopTime', num2str(stop_time, 16));

u = extract_signal(simOut, 'u_mash');
y1 = extract_signal(simOut, 'y1_mash');
v1 = extract_signal(simOut, 'v1_mash');
r1 = extract_signal(simOut, 'r1_mash');
y2 = extract_signal(simOut, 'y2_mash');
v2 = extract_signal(simOut, 'v2_mash');
v_mash = extract_signal(simOut, 'v_mash');

n_avail = min([numel(u), numel(y1), numel(v1), numel(r1), numel(y2), numel(v2), numel(v_mash)]);
u = u(1:n_avail);
y1 = y1(1:n_avail);
v1 = v1(1:n_avail);
r1 = r1(1:n_avail);
y2 = y2(1:n_avail);
v2 = v2(1:n_avail);
v_mash = v_mash(1:n_avail);
t = (0:n_avail-1) / fs;

[SNR, ENOB, V_out_mag, sig_bin, signal_power, noise_power] = ...
    calculate_snr_metrics(v_mash, n_avail, OSR);

fB = fs / (2 * OSR);
freqs_hz = (0:n_avail/2) / n_avail * fs;
f_min = max(fs / n_avail, 1e3);
log_idx = freqs_hz >= f_min;
V_out_dB = 20 * log10(V_out_mag(1:n_avail/2+1) + eps);

freqs_ntf_hz = logspace(log10(f_min), log10(fs/2), 2000);
z_eval = exp(2*pi*1j*freqs_ntf_hz/fs);
NTF1_dB = 20 * log10(abs(evalTF(info.stages(1).ntf, z_eval)) + eps);
NTF2_dB = 20 * log10(abs(evalTF(info.stages(2).ntf, z_eval)) + eps);
NTF_mash_dB = 20 * log10(abs(evalTF(info.stages(1).ntf, z_eval) .* evalTF(info.stages(2).ntf, z_eval)) + eps);

fig1 = figure('Name', '2-1 MASH Simulink Summary');
set(fig1, 'Position', [100 100 1400 900]);

subplot(2, 2, 1);
semilogx(freqs_hz(log_idx)/1e3, V_out_dB(log_idx), 'b-', 'LineWidth', 0.75);
hold on;
semilogx([fB/1e3 fB/1e3], [-150 10], 'r--', 'LineWidth', 1);
hold off;
xlabel('Frequency (kHz)');
ylabel('dBFS');
title(sprintf('Simulink MASH Spectrum (SNR = %.1f dB)', SNR));
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
n_plot = min(350, n_avail);
stairs(t(1:n_plot)*1e6, v1(1:n_plot), 'b-', 'LineWidth', 1);
hold on;
plot(t(1:n_plot)*1e6, r1(1:n_plot), 'r-', 'LineWidth', 1);
stairs(t(1:n_plot)*1e6, v2(1:n_plot), 'k-', 'LineWidth', 1);
plot(t(1:n_plot)*1e6, u(1:n_plot), 'g-', 'LineWidth', 1.25);
hold off;
xlabel('Time (\mus)');
ylabel('Amplitude');
title('Stage Signals');
legend('v_1', 'Residue y_1-v_1', 'v_2', 'Input', 'Location', 'best');
grid on;

subplot(2, 2, 4);
summary_text = sprintf(['Samples: %d\n' ...
    'Signal bin: %d\n' ...
    'Signal power: %.3e\n' ...
    'Noise power: %.3e\n' ...
    'ENOB: %.2f bits\n' ...
    'Model: %s'], ...
    n_avail, sig_bin, signal_power, noise_power, ENOB, model_name);
text(0.05, 0.55, summary_text, 'FontSize', 12, 'Interpreter', 'none');
axis off;
title('Simulation Summary');

plot_path = fullfile(results_dir, 'dsm_mash_2_1_ciff_10mhz_simulink_plots.png');
saveas(fig1, plot_path);

results_path = fullfile(results_dir, 'dsm_mash_2_1_ciff_10mhz_simulink_results.mat');
save(results_path, 'u', 'y1', 'v1', 'r1', 'y2', 'v2', 'v_mash', 't', 'SNR', 'ENOB', ...
    'V_out_mag', 'sig_bin', 'signal_power', 'noise_power', 'info', 'fs', ...
    'OSR', 'N', 'A_in', 'freqs_ntf_hz', 'NTF1_dB', 'NTF2_dB', 'NTF_mash_dB', 'model_path');

fprintf('Saved: %s\n', plot_path);
fprintf('Saved: %s\n', results_path);
close_system(model_name, 0);

function sig = extract_signal(simOut, name)
if isprop(simOut, name)
    raw = simOut.(name);
else
    raw = simOut.get(name);
end

if isa(raw, 'timeseries')
    sig = raw.Data(:);
elseif isstruct(raw) && isfield(raw, 'signals')
    sig = raw.signals.values(:);
elseif isnumeric(raw)
    if isvector(raw)
        sig = raw(:);
    elseif size(raw, 2) >= 2
        sig = raw(:, end);
    else
        sig = raw(:);
    end
else
    error('Unsupported signal format for "%s".', name);
end
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
