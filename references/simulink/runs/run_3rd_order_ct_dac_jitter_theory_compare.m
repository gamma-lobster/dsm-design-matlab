%% run_3rd_order_ct_dac_jitter_theory_compare.m
% Compare waveform-level and theory-faithful DAC jitter models against an
% Eq. 9.17-style SJNR estimate. The NTF integral is evaluated numerically
% in normalized frequency f and then converted to the w-domain form used by
% the textbook equation via dw = 2*pi*df.

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
Ts = 1/fs;
n_bits = 4;
n_levels = 2^n_bits;
V_fs = 1.0;
A_in = 0.5;
N = 8192;
stop_time = (N - 1) / fs;
jitter_rms = 0.5e-9;

ntf = synthesizeNTF(order, OSR, opt, H_inf, f0);
[ABCDc, ~] = realizeNTF_ct(ntf, ct_form, tdac);

cases = struct([]);
cases(1).label = 'Ideal';
cases(1).jitter_rms = 0;
cases(1).jitter_mode = 'edge';
cases(2).label = 'Edge jitter';
cases(2).jitter_rms = jitter_rms;
cases(2).jitter_mode = 'edge';
cases(3).label = 'Equivalent jitter';
cases(3).jitter_rms = jitter_rms;
cases(3).jitter_mode = 'equivalent';

for k = 1:numel(cases)
    model_name = sprintf('dsm_3rd_order_ct_jitter_%s_500ps_model', lower(strrep(cases(k).jitter_mode, ' ', '_')));
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
        'jitter_rms', cases(k).jitter_rms, ...
        'jitter_mode', cases(k).jitter_mode, ...
        'jitter_seed', 12345, ...
        'open_model', false);

    fprintf('Running %s case: %s\n', cases(k).label, built_name);
    simOut = sim(built_name, 'StopTime', num2str(stop_time, 16));

    [v, ~] = extract_signal_with_time(simOut, 'v_dsm');
    [u_ct, t_u] = extract_signal_with_time(simOut, 'u_dsm');
    [y, t_y] = extract_signal_with_time(simOut, 'y_dsm');

    n_avail = min([numel(v), numel(y), N]);
    v = v(1:n_avail);
    y = y(1:n_avail);
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

    cases(k).model_name = built_name;
    cases(k).model_path = model_path;
    cases(k).info = info;
    cases(k).t = t;
    cases(k).u = u;
    cases(k).y = y;
    cases(k).v = v;
    cases(k).SNR = SNR;
    cases(k).ENOB = ENOB;
    cases(k).V_out_mag = V_out_mag;
    cases(k).sig_bin = sig_bin;
    cases(k).signal_power = signal_power;
    cases(k).noise_power = noise_power;
end

freqs_hz = (0:N/2) / N * fs;
f_min = max(fs / N, 1e3);
log_idx = freqs_hz >= f_min;
freqs_log_khz = freqs_hz(log_idx) / 1000;
fB = fs / (2 * OSR);

% Full-band f-domain integral used by the standalone jitter calculator:
%   int_0^(1/2) |(1 - z^-1)NTF(z)|^2 df
f_int = linspace(0, 0.5, 200000);
z_int = exp(1j * 2*pi * f_int);
H_int = evalTF(ntf, z_int);
int_term_f = trapz(f_int, abs((1 - z_int.^(-1)) .* H_int).^2);
int_term_w = 2*pi * int_term_f;
delta = 2 * V_fs / (n_levels - 1);
J_theory = (jitter_rms^2 / Ts^2) * (1 / OSR) * (1 / (3*pi)) * int_term_w * (delta/2)^2;
P_sig = A_in^2 / 2;
SJNR_theory = 10 * log10(P_sig / J_theory);

fig = figure('Name', 'CTDSM DAC Jitter Theory Comparison');
set(fig, 'Position', [120 120 1450 900]);

subplot(2, 2, 1);
for k = 1:numel(cases)
    V_out_dB = 20 * log10(cases(k).V_out_mag(1:N/2+1) + eps);
    semilogx(freqs_log_khz, V_out_dB(log_idx), 'LineWidth', 0.9);
    hold on;
end
semilogx([fB/1000 fB/1000], [-150 10], 'k--', 'LineWidth', 1);
hold off;
xlabel('Frequency (kHz)');
ylabel('dBFS');
title('Output Spectrum Comparison');
legend(arrayfun(@(c) sprintf('%s (%.1f dB)', c.label, c.SNR), cases, 'UniformOutput', false), 'Location', 'best');
grid on;
axis([f_min/1000 fs/2000 -140 10]);

subplot(2, 2, 2);
bar(categorical({'Ideal', 'Edge', 'Equivalent', 'Theory'}), ...
    [cases(1).SNR, cases(2).SNR, cases(3).SNR, SJNR_theory]);
ylabel('SNR (dB)');
title('500 ps DAC Jitter SNR');
grid on;

subplot(2, 2, 3);
n_plot = min(400, numel(cases(1).t));
plot(cases(1).t(1:n_plot)*1e6, cases(1).u(1:n_plot), 'k-', 'LineWidth', 1.0);
hold on;
stairs(cases(1).t(1:n_plot)*1e6, cases(1).v(1:n_plot), 'LineWidth', 1.0);
stairs(cases(2).t(1:n_plot)*1e6, cases(2).v(1:n_plot), 'LineWidth', 1.0);
stairs(cases(3).t(1:n_plot)*1e6, cases(3).v(1:n_plot), 'LineWidth', 1.0);
hold off;
xlabel('Time (\mus)');
ylabel('Amplitude');
title('Sampled Output Comparison');
legend('Input', 'Ideal', 'Edge jitter', 'Equivalent jitter', 'Location', 'best');
grid on;

subplot(2, 2, 4);
summary_lines = {
    sprintf('jitter_{rms} = %.0f ps, tdac = [%.1f %.1f] Ts', jitter_rms*1e12, tdac(1), tdac(2))
    sprintf('Ideal SNR       = %.2f dB', cases(1).SNR)
    sprintf('Edge jitter SNR = %.2f dB', cases(2).SNR)
    sprintf('Eqv. jitter SNR = %.2f dB', cases(3).SNR)
    sprintf('Theory SJNR     = %.2f dB', SJNR_theory)
    sprintf('Theory J        = %.6g', J_theory)
    sprintf('Signal power    = %.6g', P_sig)
    sprintf('Integral f-term = %.6f', int_term_f)
    sprintf('Integral w-term = %.6f', int_term_w)
    };
text(0.05, 0.95, strjoin(summary_lines, '\n'), 'FontSize', 11, 'VerticalAlignment', 'top');
axis off;
title('Summary');

plot_path = fullfile(results_dir, 'dsm_3rd_order_ct_dac_jitter_theory_compare.png');
saveas(fig, plot_path);

results_path = fullfile(results_dir, 'dsm_3rd_order_ct_dac_jitter_theory_compare.mat');
save(results_path, 'cases', 'plot_path', 'jitter_rms', 'tdac', 'ABCDc', ...
    'fs', 'OSR', 'N', 'A_in', 'J_theory', 'SJNR_theory', 'P_sig', 'int_term_f', 'int_term_w', 'delta');

fprintf('Theory J = %.9g\n', J_theory);
fprintf('Theory SJNR = %.4f dB\n', SJNR_theory);
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

noise_bins = setdiff(2:fB_bins, [sig_bins harmonic_bins]);
noise_power = sum(V_out_mag(noise_bins).^2);
SNR = 10 * log10(signal_power / max(noise_power, eps));
ENOB = (SNR - 1.76) / 6.02;
end
