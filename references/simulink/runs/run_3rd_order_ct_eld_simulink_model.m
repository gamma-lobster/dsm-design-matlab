%% run_3rd_order_ct_eld_simulink_model.m
% Compare uncompensated and compensated Ts/2 excess-loop-delay cases
% for the 3rd-order CTDSM Simulink model.

clearvars; close all; clc;

script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(fileparts(script_dir)), 'common'));
refs_root = setup_references_paths();
results_dir = fullfile(refs_root, 'results', 'ct', 'eld');

order = 3;
OSR = 32;
H_inf = 4.0;
opt = 1;
f0 = 0;
ct_form = 'FF';
tdac_nom = [0 1];
tdac_eld = [0.5 1.5];
fs = 10e6;
n_bits = 4;
V_fs = 1.0;
A_in = 0.5;
N = 8192;
stop_time = (N - 1) / fs;

ntf = synthesizeNTF(order, OSR, opt, H_inf, f0);
[ABCDc_nom, tdac2_nom] = realizeNTF_ct(ntf, ct_form, tdac_nom); %#ok<ASGLU>
[ABCDc_comp, tdac2_comp] = realizeNTF_ct(ntf, ct_form, tdac_eld); %#ok<ASGLU>

cases = { ...
    struct( ...
        'label', 'Uncompensated ELD', ...
        'model_name', 'dsm_3rd_order_ct_eld_uncomp_topology_model', ...
        'ABCDc', ABCDc_nom, ...
        'tdac', tdac_eld), ...
    struct( ...
        'label', 'Compensated ELD', ...
        'model_name', 'dsm_3rd_order_ct_eld_comp_topology_model', ...
        'ABCDc', ABCDc_comp, ...
        'tdac', tdac_eld) ...
    };

results = struct([]);

for k = 1:numel(cases)
    cfg = cases{k};
    [model_name, model_path, info] = build_ct_dsm_simulink_model( ...
        'model_name', cfg.model_name, ...
        'output_dir', results_dir, ...
        'ABCDc', cfg.ABCDc, ...
        'ct_form', ct_form, ...
        'tdac', cfg.tdac, ...
        'fs', fs, ...
        'OSR', OSR, ...
        'n_bits', n_bits, ...
        'V_fs', V_fs, ...
        'input_amplitude', A_in, ...
        'stop_time', stop_time, ...
        'ct_steps_per_sample', 32, ...
        'open_model', false);

    fprintf('Running %s: %s\n', cfg.label, model_name);
    simOut = sim(model_name, 'StopTime', num2str(stop_time, 16));

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

    results(k).label = cfg.label; %#ok<SAGROW>
    results(k).model_name = model_name;
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

sys_c_nom = ss(ABCDc_nom(1:order, 1:order), ABCDc_nom(1:order, order+1:end), ...
    ABCDc_nom(order+1, 1:order), ABCDc_nom(order+1, order+1:end));
sys_d_nom = mapCtoD(sys_c_nom, [-1 -1; tdac_eld]);
[ntf_uncomp, ~] = calculateTF([sys_d_nom.a sys_d_nom.b; sys_d_nom.c sys_d_nom.d]);
ntf_uncomp = cancelPZ(ntf_uncomp);

sys_c_comp = ss(ABCDc_comp(1:order, 1:order), ABCDc_comp(1:order, order+1:end), ...
    ABCDc_comp(order+1, 1:order), ABCDc_comp(order+1, order+1:end));
sys_d_comp = mapCtoD(sys_c_comp, tdac2_comp);
[ntf_comp, ~] = calculateTF([sys_d_comp.a sys_d_comp.b; sys_d_comp.c sys_d_comp.d]);
ntf_comp = cancelPZ(ntf_comp);

freqs_ntf_hz = logspace(log10(f_min), log10(fs/2), 2000);
z_eval = exp(2*pi*1j*freqs_ntf_hz/fs);
NTF_target_dB = 20 * log10(abs(evalTF(ntf, z_eval)) + eps);
NTF_uncomp_dB = 20 * log10(abs(evalTF(ntf_uncomp, z_eval)) + eps);
NTF_comp_dB = 20 * log10(abs(evalTF(ntf_comp, z_eval)) + eps);

fig1 = figure('Name', '3rd-Order CTDSM ELD Comparison');
set(fig1, 'Position', [100 100 1450 950]);

subplot(2, 2, 1);
for k = 1:numel(results)
    V_out_dB = 20 * log10(results(k).V_out_mag(1:N/2+1) + eps);
    semilogx(freqs_log_khz, V_out_dB(log_idx), 'LineWidth', 0.9);
    hold on;
end
semilogx([fB/1000 fB/1000], [-150 10], 'k--', 'LineWidth', 1);
hold off;
xlabel('Frequency (kHz)');
ylabel('dBFS');
title('Output Spectrum With Ts/2 ELD');
legend(sprintf('%s (%.1f dB)', results(1).label, results(1).SNR), ...
    sprintf('%s (%.1f dB)', results(2).label, results(2).SNR), ...
    'BW Edge', 'Location', 'best');
grid on;
axis([f_min/1000 fs/2000 -140 10]);

subplot(2, 2, 2);
semilogx(freqs_ntf_hz/1000, NTF_target_dB, 'k-', 'LineWidth', 1.2);
hold on;
semilogx(freqs_ntf_hz/1000, NTF_uncomp_dB, 'r--', 'LineWidth', 1.0);
semilogx(freqs_ntf_hz/1000, NTF_comp_dB, 'b-.', 'LineWidth', 1.0);
semilogx([fB/1000 fB/1000], [-120 40], 'k:', 'LineWidth', 1);
hold off;
xlabel('Frequency (kHz)');
ylabel('Magnitude (dB)');
title('NTF Comparison');
legend('Target NTF', 'Uncompensated ELD', 'Compensated ELD', 'BW Edge', ...
    'Location', 'best');
grid on;
axis([f_min/1000 fs/2000 -120 40]);

subplot(2, 2, 3);
n_plot = min(400, numel(results(1).t));
plot(results(1).t(1:n_plot)*1e6, results(1).u(1:n_plot), 'g-', 'LineWidth', 1.0);
hold on;
stairs(results(1).t(1:n_plot)*1e6, results(1).v(1:n_plot), 'r-', 'LineWidth', 1.0);
stairs(results(2).t(1:n_plot)*1e6, results(2).v(1:n_plot), 'b-', 'LineWidth', 1.0);
hold off;
xlabel('Time (\mus)');
ylabel('Amplitude');
title('Sampled Output Comparison');
legend('Input', 'Uncompensated ELD', 'Compensated ELD', 'Location', 'best');
grid on;

subplot(2, 2, 4);
summary_text = sprintf(['ELD timing: [%.2f %.2f] Ts\n' ...
    '%s: SNR = %.2f dB, ENOB = %.2f\n' ...
    '%s: SNR = %.2f dB, ENOB = %.2f'], ...
    tdac_eld(1), tdac_eld(2), ...
    results(1).label, results(1).SNR, results(1).ENOB, ...
    results(2).label, results(2).SNR, results(2).ENOB);
text(0.05, 0.6, summary_text, 'FontSize', 12);
axis off;
title('Simulation Summary');

plot_path = fullfile(results_dir, 'dsm_3rd_order_ct_eld_simulink_plots.png');
saveas(fig1, plot_path);

results_path = fullfile(results_dir, 'dsm_3rd_order_ct_eld_simulink_results.mat');
save(results_path, 'results', 'plot_path', 'freqs_ntf_hz', 'NTF_target_dB', ...
    'NTF_uncomp_dB', 'NTF_comp_dB', 'tdac_nom', 'tdac_eld', 'ABCDc_nom', ...
    'ABCDc_comp', 'tdac2_comp', 'fs', 'OSR', 'N', 'A_in');

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
