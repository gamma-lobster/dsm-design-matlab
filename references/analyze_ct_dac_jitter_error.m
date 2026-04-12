%% analyze_ct_dac_jitter_error.m
% Extract the equivalent DAC jitter error sequence from the CTDSM Simulink
% model, then plot its waveform and spectrum.

clearvars; close all; clc;

script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'dstoolbox'));

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
jitter_rms = 0.5e-9;
jitter_seed = 12345;

% Use a longer record and analyze the second half after the loop settles.
N_analyze = 65536;
N_total = 2 * N_analyze;
stop_time = (N_total - 1) / fs;

% Keep the same default coherent-ish tone choice used elsewhere in the repo.
f_bin = round(sqrt(1/7) * 8192 / (2*OSR));
f_in = f_bin / 8192 * fs;

ntf = synthesizeNTF(order, OSR, opt, H_inf, f0);
[ABCDc, ~] = realizeNTF_ct(ntf, ct_form, tdac);

model_name = 'ct_dac_jitter_error_analysis_model';
[built_name, model_path, info] = build_ct_dsm_simulink_model( ...
    'model_name', model_name, ...
    'output_dir', script_dir, ...
    'ABCDc', ABCDc, ...
    'ct_form', ct_form, ...
    'tdac', tdac, ...
    'fs', fs, ...
    'OSR', OSR, ...
    'n_bits', n_bits, ...
    'V_fs', V_fs, ...
    'input_amplitude', A_in, ...
    'input_frequency', f_in, ...
    'stop_time', stop_time, ...
    'ct_steps_per_sample', 32, ...
    'jitter_rms', jitter_rms, ...
    'jitter_mode', 'equivalent', ...
    'jitter_seed', jitter_seed, ...
    'open_model', false);

load_system(built_name);
add_internal_signal_logging(built_name);
save_system(built_name, model_path);

simOut = sim(built_name, 'StopTime', num2str(stop_time, 16));

[e_jit, t_e] = extract_signal_with_time(simOut, 'e_jit_dsm');
[dv, t_dv] = extract_signal_with_time(simOut, 'dv_dsm');
[dt_over_T, t_dt] = extract_signal_with_time(simOut, 'dt_over_T_dsm');

if isempty(t_e)
    t_e = (0:numel(e_jit)-1).' / fs;
end
if isempty(t_dv)
    t_dv = (0:numel(dv)-1).' / fs;
end
if isempty(t_dt)
    t_dt = (0:numel(dt_over_T)-1).' / fs;
end

e_jit = e_jit(:);
dv = dv(:);
dt_over_T = dt_over_T(:);

e_ss = e_jit(end-N_analyze+1:end);
dv_ss = dv(end-N_analyze+1:end);
dt_ss = dt_over_T(end-N_analyze+1:end);
t_ss = ((0:N_analyze-1).' / fs);

w = 0.5 * (1 - cos(2*pi*(0:N_analyze-1)/N_analyze));
E = fft(e_ss(:).' .* w) / (N_analyze/4);
E_mag = abs(E(1:N_analyze/2+1));
freqs_hz = (0:N_analyze/2) / N_analyze * fs;

% Use a simple log-log fit on the midband PSD to quantify whiteness roughly.
fit_idx = freqs_hz >= 2e3 & freqs_hz <= fs/4;
p = polyfit(log10(freqs_hz(fit_idx)), 20*log10(E_mag(fit_idx) + eps), 1);
slope_db_per_dec = p(1);

fig = figure('Name', 'CTDSM Equivalent DAC Jitter Error');
set(fig, 'Position', [120 120 1450 920]);

subplot(2, 2, 1);
n_plot = 500;
plot(t_ss(1:n_plot)*1e6, e_ss(1:n_plot), 'b-', 'LineWidth', 1.0);
xlabel('Time (\mus)');
ylabel('e_j[n]');
title('Equivalent DAC Jitter Error Waveform');
grid on;

subplot(2, 2, 2);
plot(t_ss(1:n_plot)*1e6, dv_ss(1:n_plot), 'r-', 'LineWidth', 1.0);
hold on;
plot(t_ss(1:n_plot)*1e6, dt_ss(1:n_plot), 'k-', 'LineWidth', 0.9);
hold off;
xlabel('Time (\mus)');
ylabel('Amplitude');
title('\Delta v[n] And \Delta t[n]/T');
legend('\Delta v[n]', '\Delta t[n]/T', 'Location', 'best');
grid on;

subplot(2, 2, 3);
semilogx(freqs_hz(2:end)/1e3, 20*log10(E_mag(2:end) + eps), 'b-', 'LineWidth', 1.0);
xlabel('Frequency (kHz)');
ylabel('Magnitude (dB)');
title('Equivalent DAC Jitter Error Spectrum');
grid on;
axis([1 fs/2000 -180 20]);

subplot(2, 2, 4);
summary_lines = {
    sprintf('jitter mode      = equivalent')
    sprintf('jitter rms       = %.0f ps', jitter_rms*1e12)
    sprintf('input frequency  = %.6f kHz', f_in/1e3)
    sprintf('analyzed samples = %d', N_analyze)
    sprintf('rms(e_j)         = %.6g', rms(e_ss))
    sprintf('rms(\\Delta v)     = %.6g', rms(dv_ss))
    sprintf('rms(\\Delta t/T)   = %.6g', rms(dt_ss))
    sprintf('PSD slope        = %.3f dB/dec', slope_db_per_dec)
    };
text(0.05, 0.95, strjoin(summary_lines, '\n'), 'FontSize', 11, 'VerticalAlignment', 'top');
axis off;
title('Summary');

plot_path = fullfile(script_dir, 'ct_dac_jitter_error_spectrum.png');
saveas(fig, plot_path);

results_path = fullfile(script_dir, 'ct_dac_jitter_error_spectrum.mat');
save(results_path, 'e_jit', 'dv', 'dt_over_T', 'e_ss', 'dv_ss', 'dt_ss', ...
    't_e', 't_dv', 't_dt', 't_ss', 'E_mag', 'freqs_hz', 'plot_path', ...
    'jitter_rms', 'jitter_seed', 'f_in', 'slope_db_per_dec', 'info');

fprintf('Saved: %s\n', plot_path);
fprintf('Saved: %s\n', results_path);
fprintf('Approximate PSD slope: %.3f dB/dec\n', slope_db_per_dec);

close_system(built_name, 0);

function add_internal_signal_logging(model)
signals = {
    'dac_jitter_err_1', 'e_jit_out', 'e_jit_dsm', [1540 178 1615 198]
    'dac_step_diff_1', 'dv_out', 'dv_dsm', [1320 150 1390 170]
    'dac_jitter_scale_1', 'dt_out', 'dt_over_T_dsm', [1395 235 1475 255]
    };

for i = 1:size(signals, 1)
    src = signals{i, 1};
    blk = signals{i, 2};
    var = signals{i, 3};
    pos = signals{i, 4};
    if ~isempty(find_system(model, 'SearchDepth', 1, 'Name', blk))
        continue;
    end
    add_block('simulink/Sinks/To Workspace', [model '/' blk], ...
        'Position', pos, ...
        'VariableName', var, ...
        'SaveFormat', 'Array');
    add_line(model, [src '/1'], [blk '/1'], 'autorouting', 'on');
end
end

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
