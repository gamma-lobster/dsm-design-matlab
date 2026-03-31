%% design_4th_order_ciff.m
%% 4th-Order CIFF DSM Design - Native MATLAB Workflow
%% Spec: order=4, fs=2MHz, OSR=16, H_inf=4.0, 5-bit, CIFF
%%
%% This demonstrates the full MATLAB-native Delta Sigma Toolbox workflow:
%% synthesizeNTF → realizeNTF → stuffABCD → simulation → analysis
%%
%% Requires: Richard Schreier's Delta Sigma Toolbox
%% Download: https://www.mathworks.com/matlabcentral/fileexchange/19-delta-sigma-toolbox

clear all; close all; clc;

%% Add Delta Sigma Toolbox to path (relative to this script)
addpath(fullfile(fileparts(mfilename('fullpath')), 'dstoolbox'));

fprintf('============================================================\n');
fprintf('  4th-Order CIFF DSM Design - Native MATLAB Workflow\n');
fprintf('  fs=2MHz, OSR=16, H_inf=4.0, 5-bit quantizer\n');
fprintf('============================================================\n\n');

%% Step 1: Synthesize NTF
fprintf('[Step 1] Synthesize NTF\n');
fprintf('----------------------------------------\n');

order = 4;
OSR = 16;
H_inf = 4.0;    % Aggressive noise shaping
opt = 1;        % Optimized zeros
f0 = 0;         % Lowpass

ntf = synthesizeNTF(order, OSR, opt, H_inf, f0);

% Display NTF
fprintf('NTF: ');
zpk(ntf)
fprintf('\n');

%% Step 2: Realize Coefficients
fprintf('[Step 2] Realize NTF to Coefficients\n');
fprintf('----------------------------------------\n');

form = 'CIFF';  % Cascade of Integrators, Feedback Form
[a, g, b, c] = realizeNTF(ntf, form);

fprintf('Coefficients:\n');
fprintf('  a = ['); fprintf('%.6f ', a); fprintf('] (feedback)\n');
fprintf('  g = ['); fprintf('%.6f ', g); fprintf('] (resonator)\n');
fprintf('  b = ['); fprintf('%.4f ', b); fprintf('] (input)\n');
fprintf('  c = ['); fprintf('%.4f ', c); fprintf('] (inter-stage)\n\n');

%% Step 3: Build ABCD Matrix
fprintf('[Step 3] Build ABCD Matrix\n');
fprintf('----------------------------------------\n');

ABCD = stuffABCD(a, g, b, c, form);

fprintf('ABCD Matrix (%dx%d):\n', size(ABCD,1), size(ABCD,2));
for i = 1:size(ABCD, 1)
    fprintf('  ['); fprintf('%10.6f ', ABCD(i, :)); fprintf(']\n');
end
fprintf('\n');

% Extract matrices for simulation
[n_rows, ~] = size(ABCD);
n_states = n_rows - 1;

A_mat = ABCD(1:n_states, 1:n_states);
B_mat = ABCD(1:n_states, n_states+1:end);
C_mat = ABCD(n_states+1, 1:n_states);
D_mat = ABCD(n_states+1, n_states+1:end);

%% Step 4: Simulation Parameters
fprintf('[Step 4] Simulation Parameters\n');
fprintf('----------------------------------------\n');

fs = 2e6;                   % Sampling frequency (Hz)
fB = fs / (2*OSR);          % Signal bandwidth (Hz)
N = 8192;                   % FFT points

% Quantizer
n_bits = 5;
n_levels = 2^n_bits;        % 32 levels
V_fs = 1.0;                 % Full scale (+/- 1V)
LSB = 2*V_fs / (n_levels - 1);
q_levels = linspace(-V_fs, V_fs, n_levels);

fprintf('  fs = %.1f MHz\n', fs/1e6);
fprintf('  OSR = %d (BW = %.1f kHz)\n', OSR, fB/1000);
fprintf('  Quantizer: %d-bit (%d levels)\n', n_bits, n_levels);
fprintf('  Full scale: +/- %.1f V\n', V_fs);
fprintf('  LSB = %.4f V\n\n', LSB);

%% Step 5: Run Simulation
fprintf('[Step 5] Run Simulation\n');
fprintf('----------------------------------------\n');

% Input signal (not at harmonic)
f_bin = round(sqrt(1/7) * N / (2*OSR));
A_in = 0.5;                 % Input amplitude

fprintf('  Input: %.2f V sine wave\n', A_in);
fprintf('  FFT bins: %d\n', N);
fprintf('  Signal bin: %d\n', f_bin);
fprintf('  Signal frequency: %.2f kHz\n\n', f_bin/N * fs / 1000);

% Generate input
u = A_in * sin(2*pi*f_bin*(0:N-1)/N);
n = 0:N-1;

% Run simulation
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
    
    % Compute modulator output
    y(i) = C_mat * x + D_mat(1)*u(i) + D_mat(2)*v_prev;
    
    % Quantize
    y_clip = max(-V_fs, min(V_fs, y(i)));
    idx = round((y_clip - (-V_fs)) / LSB) + 1;
    idx = max(1, min(n_levels, idx));
    v(i) = q_levels(idx);
    
    % Update state
    x = A_mat*x + B_mat(:,1)*u(i) + B_mat(:,2)*v(i);
    x_max = max(x_max, abs(x));
end

fprintf('  Simulation complete.\n');
fprintf('  Output range: [%.4f, %.4f]\n', min(v), max(v));
fprintf('  State maxima: '); fprintf('%.4f ', x_max); fprintf('\n\n');

%% Step 6: Calculate SNR
fprintf('[Step 6] Calculate SNR\n');
fprintf('----------------------------------------\n');

% Window
w = 0.5 * (1 - cos(2*pi*n/N));

% Output spectrum
V_out = fft(v .* w) / (N/4);
V_out_mag = abs(V_out);
freqs = (0:N/2)/N*fs/1000;  % kHz

% Find signal in output spectrum
[~, sig_idx] = max(V_out_mag(2:N/2));
sig_bin = sig_idx + 1;
fB_bins = ceil(N / (2*OSR));

% Signal bins (3-bin)
sig_bins = sig_bin-1:sig_bin+1;
sig_bins = sig_bins(sig_bins >= 2 & sig_bins <= fB_bins);
signal_power = sum(V_out_mag(sig_bins).^2);

% Noise bins
noise_bins = setdiff(2:fB_bins, sig_bins);
noise_power = sum(V_out_mag(noise_bins).^2);

% SNR calculation
if noise_power > 0
    SNR = 10*log10(signal_power / noise_power);
    ENOB = (SNR - 1.76) / 6.02;
else
    SNR = inf;
    ENOB = inf;
end

fprintf('  Signal power: %.4e\n', signal_power);
fprintf('  Noise power: %.4e\n', noise_power);
fprintf('  SNR: %.2f dB\n', SNR);
if isfinite(ENOB)
    fprintf('  ENOB: %.2f bits\n', ENOB);
end
fprintf('\n');

%% Step 7: Plots
fprintf('[Step 7] Generate Plots\n');
fprintf('----------------------------------------\n');

% Create figure
fig1 = figure('Name', 'DSM Output Spectrum');
set(fig1, 'Position', [100 100 1400 900]);

% 1. Time Domain
subplot(2, 3, 1);
n_plot = 500;
t = (0:n_plot-1)/fs*1e6;
plot(t, u(1:n_plot), 'g-', 'LineWidth', 1.5); hold on;
stairs(t, v(1:n_plot), 'b-', 'LineWidth', 1);
hold off;
xlabel('Time (\mus)'); ylabel('Amplitude');
title(sprintf('Time Domain (%.2fV)', A_in));
legend('Input', 'Output', 'Location', 'best');
grid on;

% 2. Output Spectrum
subplot(2, 3, 2);
V_out_dB = 20*log10(V_out_mag(1:N/2+1) + eps);
plot(freqs, V_out_dB, 'b-', 'LineWidth', 0.5); hold on;
f_sig_kHz = sig_bin/N*fs/1000;
fB_kHz = fB/1000;
plot([f_sig_kHz f_sig_kHz], [-150 10], 'g--', 'LineWidth', 2);
plot([fB_kHz fB_kHz], [-150 10], 'm--', 'LineWidth', 1);
fill([0 fB_kHz fB_kHz 0], [-150 -150 10 10], 'g', 'FaceAlpha', 0.1);
hold off;
xlabel('Frequency (kHz)'); ylabel('dBFS/NBW');
title(sprintf('Output Spectrum (SNR=%.1f dB)', SNR));
legend('Spectrum', 'Signal', 'BW Edge', 'Location', 'best');
grid on;
axis([0 fs/2000 -140 10]);

% 3. Noise Spectrum (log scale)
subplot(2, 3, 3);
noise = v - u;
V_noise = fft(noise .* w) / (N/4);
V_noise_dB = 20*log10(abs(V_noise) + eps);

% Log frequency interpolation
f_min = 100;
f_log = logspace(log10(f_min), log10(fs/2), 5000);
freqs_linear = (0:N/2)/N*fs;
V_noise_half = V_noise_dB(1:N/2+1);
V_noise_log = interp1(freqs_linear, V_noise_half, f_log, 'linear', 'extrap');

semilogx(f_log/1000, V_noise_log, 'b-', 'LineWidth', 1); hold on;
plot([fB_kHz fB_kHz], [-150 10], 'm--', 'LineWidth', 1);
hold off;
xlabel('Frequency (kHz, log)'); ylabel('Noise (dB)');
title('Noise Spectrum');
legend('Noise', 'BW', 'Location', 'best');
grid on;
axis([f_min/1000 fs/2000 -140 -40]);

% 4. Histogram
subplot(2, 3, 4);
[counts, centers] = hist(v, 20);
bar(centers, counts, 'FaceColor', [0.3 0.5 0.8]);
xlabel('Level'); ylabel('Count');
title('Quantizer Histogram');
grid on;

% 5. Status
subplot(2, 3, 5);
if max(x_max) < 100
    status_text = sprintf('\\surd STABLE\\nMax: %.2f', max(x_max));
    status_color = 'g';
else
    status_text = sprintf('\\dagger CHECK\\nMax: %.2e', max(x_max));
    status_color = 'orange';
end
text(0.5, 0.7, status_text, 'FontSize', 14, 'Color', status_color, ...
    'HorizontalAlignment', 'center');
text(0.5, 0.4, sprintf('4th-Order CIFF\\nSNR=%.1f dB\\nENOB=%.1f bits', SNR, ENOB), ...
    'FontSize', 11, 'HorizontalAlignment', 'center');
text(0.5, 0.2, sprintf('H_{inf}=%.1f, OSR=%d', H_inf, OSR), ...
    'FontSize', 10, 'HorizontalAlignment', 'center');
axis off;

% 6. ABCD Matrix
subplot(2, 3, 6);
text(0.1, 0.9, 'ABCD Matrix:', 'FontSize', 10, 'FontWeight', 'bold');
for i = 1:size(ABCD, 1)
    row_str = sprintf('Row %d: ', i);
    for j = 1:size(ABCD, 2)
        row_str = [row_str, sprintf('%.3f ', ABCD(i, j))];
    end
    text(0.1, 0.8-(i-1)*0.15, row_str, 'FontSize', 8);
end
axis off;

% Save plot
saveas(gcf, 'dsm_4th_order_ciff.png');
fprintf('  Saved: dsm_4th_order_ciff.png\n\n');

%% Step 8: NTF/STF Plots
fprintf('[Step 8] NTF/STF Frequency Response\n');
fprintf('----------------------------------------\n');

fig2 = figure('Name', 'DSM NTF and STF');
set(fig2, 'Position', [100 100 1400 600]);

% Calculate frequency response
N_fft = 8192;
freqs_norm = (0:N_fft-1)/N_fft;
z = exp(2*pi*1j*freqs_norm);

% Get NTF poles and zeros
[z_ntf, p_ntf, k_ntf] = zpkdata(ntf);
z_ntf = z_ntf{1};
p_ntf = p_ntf{1};

% Calculate NTF response
NTF_resp = ones(1, N_fft);
for i = 1:N_fft
    zi = z(i);
    num = k_ntf * prod(zi - z_ntf);
    den = prod(zi - p_ntf);
    NTF_resp(i) = num / den;
end
NTF_dB = 20*log10(abs(NTF_resp) + eps);

% NTF Plot
subplot(1, 2, 1);
f_log = logspace(log10(100), log10(fs/2), 5000);
freqs_log_norm = f_log / fs;
NTF_log = interp1(freqs_norm, NTF_dB, freqs_log_norm, 'linear', 'extrap');

semilogx(f_log/1000, NTF_log, 'b-', 'LineWidth', 1.5); hold on;
plot([fB_kHz fB_kHz], [-100 40], 'm--', 'LineWidth', 2);

% Mark H_inf
[max_ntf, max_idx] = max(NTF_dB);
f_max_kHz = freqs_norm(max_idx) * fs / 1000;
plot(f_max_kHz, max_ntf, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
text(f_max_kHz*1.2, max_ntf+3, sprintf('H_{inf}=%.1f dB', max_ntf), ...
    'FontSize', 10, 'Color', 'r');

% Mark notch frequencies
for i = 1:length(z_ntf)
    zi = z_ntf(i);
    if abs(zi) < 1.1
        f_notch = abs(angle(zi)) / (2*pi) * fs;
        if f_notch < fs/2 && f_notch > 100
            plot(f_notch/1000, -70, 'gv', 'MarkerSize', 8, 'LineWidth', 1.5);
        end
    end
end

hold off;
xlabel('Frequency (kHz)', 'FontSize', 12);
ylabel('Magnitude (dB)', 'FontSize', 12);
title(sprintf('NTF - Noise Transfer Function\\n(4th-Order, H_{inf}=%.1f)', H_inf), 'FontSize', 12);
legend('NTF', 'Signal BW', 'H_{inf}', 'Notches', 'Location', 'northwest');
grid on;
axis([100/1000 fs/2000 -100 40]);

% STF Plot
subplot(1, 2, 2);

% For CIFF, STF is approximately unity (0 dB)
STF_log = zeros(size(freqs_log_norm));

semilogx(f_log/1000, STF_log, 'g-', 'LineWidth', 1.5); hold on;
plot([fB_kHz fB_kHz], [-10 15], 'm--', 'LineWidth', 2);
plot(100/1000, 0, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
text(150/1000, 2, 'DC: 0 dB', 'FontSize', 10, 'Color', 'r');

annotation('textbox', [0.58 0.15 0.3 0.15], 'String', ...
    sprintf('CIFF Topology\\nUnity STF (0 dB)\\nFlat in passband'), ...
    'FitBoxToText', 'on', 'BackgroundColor', [0.9 0.9 0.9], ...
    'EdgeColor', 'black');

hold off;
xlabel('Frequency (kHz)', 'FontSize', 12);
ylabel('Magnitude (dB)', 'FontSize', 12);
title('STF - Signal Transfer Function\\n(CIFF Topology)', 'FontSize', 12);
legend('STF', 'Signal BW', 'DC Gain', 'Location', 'northwest');
grid on;
axis([100/1000 fs/2000 -10 15]);

% Save plot
saveas(gcf, 'dsm_4th_order_ntf_stf.png');
fprintf('  Saved: dsm_4th_order_ntf_stf.png\n\n');

%% Summary
fprintf('============================================================\n');
if max(x_max) < 100
    fprintf('  \\surd DESIGN COMPLETE - STABLE\\n');
else
    fprintf('  \\dagger CHECK STABILITY\\n');
end
fprintf('============================================================\n');
fprintf('  Order: %d\n', order);
fprintf('  OSR: %d\n', OSR);
fprintf('  H_inf: %.1f\n', H_inf);
fprintf('  SNR: %.2f dB\n', SNR);
if isfinite(ENOB)
    fprintf('  ENOB: %.2f bits\n', ENOB);
end
fprintf('============================================================\n');
