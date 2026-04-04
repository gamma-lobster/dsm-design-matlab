%% run_3rd_order_simulink_model.m
% Run the generated 3rd-order CIFF Simulink model and plot its response.

clearvars; close all; clc;

script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'dstoolbox'));

order = 3;
OSR = 32;
H_inf = 4.0;
opt = 1;
f0 = 0;
form = 'CIFF';
fs = 10e6;
n_bits = 4;
V_fs = 1.0;
A_in = 0.5;
N = 8192;
stop_time = (N - 1) / fs;

ntf = synthesizeNTF(order, OSR, opt, H_inf, f0);
[a, g, b, c] = realizeNTF(ntf, form);

[model_name, model_path, info] = build_dsm_simulink_model( ...
    'model_name', 'dsm_3rd_order_ciff_10mhz_topology_model', ...
    'output_dir', script_dir, ...
    'a', a, ...
    'g', g, ...
    'b', b, ...
    'c', c, ...
    'form', form, ...
    'fs', fs, ...
    'OSR', OSR, ...
    'n_bits', n_bits, ...
    'V_fs', V_fs, ...
    'input_amplitude', A_in, ...
    'stop_time', stop_time, ...
    'implementation', 'expanded', ...
    'open_model', false);

fprintf('Running Simulink model: %s\n', model_name);
simOut = sim(model_name, 'StopTime', num2str(stop_time, 16));

u = extract_signal(simOut, 'u_dsm');
y = extract_signal(simOut, 'y_dsm');
v = extract_signal(simOut, 'v_dsm');

n_avail = min([numel(u), numel(y), numel(v)]);
u = u(1:n_avail);
y = y(1:n_avail);
v = v(1:n_avail);
t = (0:n_avail-1) / fs;

x_max = estimate_state_max(model_name, simOut);
[SNR, ENOB, V_out_mag, sig_bin, signal_power, noise_power] = ...
    calculate_snr_metrics(v, n_avail, OSR);

fB = fs / (2 * OSR);
freqs_hz = (0:n_avail/2) / n_avail * fs;
f_min = max(fs / n_avail, 1e3);
log_idx = freqs_hz >= f_min;
freqs_log_khz = freqs_hz(log_idx) / 1000;
V_out_dB = 20 * log10(V_out_mag(1:n_avail/2+1) + eps);

freqs_ntf_hz = logspace(log10(f_min), log10(fs/2), 2000);
z_eval = exp(2*pi*1j*freqs_ntf_hz/fs);
NTF_resp = evalTF(ntf, z_eval);
NTF_dB = 20 * log10(abs(NTF_resp) + eps);

fig1 = figure('Name', '3rd-Order Simulink DSM Summary');
set(fig1, 'Position', [100 100 1350 900]);

subplot(2, 2, 1);
semilogx(freqs_log_khz, V_out_dB(log_idx), 'b-', 'LineWidth', 0.75);
hold on;
semilogx([fB/1000 fB/1000], [-150 10], 'r--', 'LineWidth', 1);
hold off;
xlabel('Frequency (kHz)');
ylabel('dBFS');
title(sprintf('Simulink Output Spectrum (SNR = %.1f dB)', SNR));
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
n_plot = min(400, n_avail);
stairs(t(1:n_plot)*1e6, v(1:n_plot), 'b-', 'LineWidth', 1);
hold on;
plot(t(1:n_plot)*1e6, u(1:n_plot), 'g-', 'LineWidth', 1.25);
hold off;
xlabel('Time (\mus)');
ylabel('Amplitude');
title('Input And Quantized Output');
legend('DSM Output', 'Input', 'Location', 'best');
grid on;

subplot(2, 2, 4);
summary_text = sprintf(['Samples: %d\n' ...
    'Signal bin: %d\n' ...
    'Signal power: %.3e\n' ...
    'Noise power: %.3e\n' ...
    'ENOB: %.2f bits\n' ...
    'Peak state est.: %.3f'], ...
    n_avail, sig_bin, signal_power, noise_power, ENOB, x_max);
text(0.05, 0.65, summary_text, 'FontSize', 12);
text(0.05, 0.25, sprintf('Model: %s', model_name), 'Interpreter', 'none', 'FontSize', 10);
axis off;
title('Simulation Summary');

plot_path = fullfile(script_dir, 'dsm_3rd_order_ciff_10mhz_simulink_plots.png');
saveas(fig1, plot_path);

results_path = fullfile(script_dir, 'dsm_3rd_order_ciff_10mhz_simulink_results.mat');
save(results_path, 'u', 'y', 'v', 't', 'SNR', 'ENOB', 'V_out_mag', 'sig_bin', ...
    'signal_power', 'noise_power', 'info', 'fs', 'OSR', 'N', 'A_in', ...
    'freqs_ntf_hz', 'NTF_dB', 'x_max', 'model_path');

fprintf('Saved: %s\n', plot_path);
fprintf('Saved: %s\n', results_path);

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

function x_max = estimate_state_max(model_name, simOut)
try
    mw = get_param(model_name, 'ModelWorkspace');
    A = mw.evalin('A_dsm');
    B = mw.evalin('B_dsm');
catch
    x_max = NaN;
    return;
end

u = extract_signal(simOut, 'u_dsm');
v = extract_signal(simOut, 'v_dsm');
n = min(numel(u), numel(v));
x = zeros(size(A, 1), 1);
x_max_vec = zeros(size(A, 1), 1);
for i = 1:n
    x = A*x + B(:,1)*u(i) + B(:,2)*v(i);
    x_max_vec = max(x_max_vec, abs(x));
end
x_max = max(x_max_vec);
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
