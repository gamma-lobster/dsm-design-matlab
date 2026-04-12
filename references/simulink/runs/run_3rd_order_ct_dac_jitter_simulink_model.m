%% run_3rd_order_ct_dac_jitter_simulink_model.m
% Sweep DAC clock jitter for the 3rd-order CTDSM with compensated half-cycle delay.

clearvars; close all; clc;

script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(fileparts(script_dir)), 'common'));
refs_root = setup_references_paths();
results_dir = fullfile(refs_root, 'results', 'ct', 'jitter');

order = 3;
OSR = 32;
H_inf = 4.0;
opt = 1;
f0 = 0;
ct_form = 'FF';
tdac = [0.5 1.5];
fs = 10e6;
n_bits = 4;
V_fs = 1.0;
A_in = 0.5;
N = 8192;
stop_time = (N - 1) / fs;

% RMS DAC clock jitter values in seconds.
% For the 10 MHz clock used here (Ts = 100 ns), sub-ns jitter is a more
% meaningful sweep than deep-ps values with the current fixed-step model.
jitter_rms_values = [0, 0.5e-9, 1e-9, 2e-9, 5e-9];

ntf = synthesizeNTF(order, OSR, opt, H_inf, f0);
[ABCDc, tdac2] = realizeNTF_ct(ntf, ct_form, tdac); %#ok<ASGLU>

results = struct([]);

for k = 1:numel(jitter_rms_values)
    jitter_rms = jitter_rms_values(k);
    model_name = sprintf('dsm_3rd_order_ct_jitter_%dps_topology_model', round(jitter_rms*1e12));

    [built_name, model_path, info] = build_ct_dsm_simulink_model( ...
        'model_name', model_name, ...
        'output_dir', results_dir, ...
        'ABCDc', ABCDc, ...
        'ct_form', ct_form, ...
        'tdac', tdac, ...
        'fs', fs, ...
        'OSR', OSR, ...
        'n_bits', n_bits, ...
        'V_fs', V_fs, ...
        'input_amplitude', A_in, ...
        'stop_time', stop_time, ...
        'ct_steps_per_sample', 32, ...
        'jitter_rms', jitter_rms, ...
        'jitter_seed', 12345, ...
        'open_model', false);

    fprintf('Running jitter case %.0f ps rms: %s\n', jitter_rms*1e12, built_name);
    simOut = sim(built_name, 'StopTime', num2str(stop_time, 16));

    [u_ct, t_u] = extract_signal_with_time(simOut, 'u_dsm');
    [y_ct, ~] = extract_signal_with_time(simOut, 'y_ct_dsm');
    [y, t_y] = extract_signal_with_time(simOut, 'y_dsm');
    [v, ~] = extract_signal_with_time(simOut, 'v_dsm');

    n_avail = min([numel(y), numel(v), N]);
    y = y(1:n_avail);
    v = v(1:n_avail);
    if isempty(t_y)
        t = (0:n_avail-1).' / fs;
    else
        t = t_y(1:n_avail);
    end

    if isempty(t_u)
        u = A_in * sin(2*pi*info.input_frequency*t);
    else
        u = interp1(t_u, u_ct, t, 'previous', 'extrap');
    end

    [SNR, ENOB, V_out_mag, sig_bin, signal_power, noise_power] = ...
        calculate_snr_metrics(v, n_avail, OSR);

    results(k).jitter_rms = jitter_rms; %#ok<SAGROW>
    results(k).label = sprintf('%.0f ps', jitter_rms*1e12);
    results(k).model_name = built_name;
    results(k).model_path = model_path;
    results(k).info = info;
    results(k).u = u;
    results(k).y = y;
    results(k).y_ct = y_ct;
    results(k).v = v;
    results(k).t = t;
    results(k).SNR = SNR;
    results(k).ENOB = ENOB;
    results(k).V_out_mag = V_out_mag;
    results(k).sig_bin = sig_bin;
    results(k).signal_power = signal_power;
    results(k).noise_power = noise_power;
end

fB = fs / (2 * OSR);
freqs_hz = (0:N/2) / N * fs;
f_min = max(fs / N, 1e3);
log_idx = freqs_hz >= f_min;
freqs_log_khz = freqs_hz(log_idx) / 1000;

fig1 = figure('Name', '3rd-Order CTDSM DAC Jitter Sweep');
set(fig1, 'Position', [100 100 1450 950]);

subplot(2, 2, 1);
legend_labels = cell(1, numel(results) + 1);
for k = 1:numel(results)
    V_out_dB = 20 * log10(results(k).V_out_mag(1:N/2+1) + eps);
    semilogx(freqs_log_khz, V_out_dB(log_idx), 'LineWidth', 0.85);
    hold on;
    legend_labels{k} = sprintf('%s (%.1f dB)', results(k).label, results(k).SNR);
end
semilogx([fB/1000 fB/1000], [-150 10], 'k--', 'LineWidth', 1);
hold off;
xlabel('Frequency (kHz)');
ylabel('dBFS');
title('Output Spectrum With DAC Clock Jitter');
legend_labels{end} = 'BW Edge';
legend(legend_labels{:}, 'Location', 'best');
grid on;
axis([f_min/1000 fs/2000 -140 10]);

subplot(2, 2, 2);
snr_vals = [results.SNR];
plot(jitter_rms_values*1e12, snr_vals, 'bo-', 'LineWidth', 1.2, 'MarkerSize', 5);
xlabel('DAC Clock Jitter RMS (ps)');
ylabel('SNR (dB)');
title('SNR Versus DAC Clock Jitter');
grid on;

subplot(2, 2, 3);
n_plot = min(400, numel(results(1).t));
plot(results(1).t(1:n_plot)*1e6, results(1).u(1:n_plot), 'g-', 'LineWidth', 1.0);
hold on;
for k = 1:min(numel(results), 3)
    stairs(results(k).t(1:n_plot)*1e6, results(k).v(1:n_plot), 'LineWidth', 1.0);
end
hold off;
xlabel('Time (\mus)');
ylabel('Amplitude');
title('Sampled Output For First Few Jitter Cases');
legend_entries = [{'Input'}, arrayfun(@(r) r.label, results(1:min(numel(results),3)), 'UniformOutput', false)];
legend(legend_entries{:}, 'Location', 'best');
grid on;

subplot(2, 2, 4);
summary_lines = cell(numel(results)+1, 1);
summary_lines{1} = sprintf('Compensated CTDSM, tdac = [%.1f %.1f] Ts', tdac(1), tdac(2));
for k = 1:numel(results)
    summary_lines{k+1} = sprintf('%s rms: SNR = %.2f dB, ENOB = %.2f', ...
        results(k).label, results(k).SNR, results(k).ENOB);
end
text(0.05, 0.95, strjoin(summary_lines, '\n'), 'FontSize', 11, 'VerticalAlignment', 'top');
axis off;
title('Simulation Summary');

plot_path = fullfile(results_dir, 'dsm_3rd_order_ct_dac_jitter_simulink_plots.png');
saveas(fig1, plot_path);

results_path = fullfile(results_dir, 'dsm_3rd_order_ct_dac_jitter_simulink_results.mat');
save(results_path, 'results', 'plot_path', 'jitter_rms_values', 'tdac', ...
    'ABCDc', 'fs', 'OSR', 'N', 'A_in');

fprintf('Saved: %s\n', plot_path);
fprintf('Saved: %s\n', results_path);

function [sig, t] = extract_signal_with_time(simOut, name)
if isprop(simOut, name)
    raw = simOut.(name);
else
    raw = simOut.get(name);
end

t = [];
if isa(raw, 'timeseries')
    sig = raw.Data(:);
    t = raw.Time(:);
elseif isstruct(raw) && isfield(raw, 'signals')
    sig = raw.signals.values(:);
    if isfield(raw, 'time')
        t = raw.time(:);
    end
elseif isnumeric(raw)
    if isvector(raw)
        sig = raw(:);
    elseif size(raw, 2) >= 2
        t = raw(:, 1);
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
