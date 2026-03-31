%% dsm_4th_order_simple.m
%% Simplified 4th-Order CIFF DSM - Headless Version
%% Runs simulation and saves results without display

clear all; clc;

%% Add Delta Sigma Toolbox to path (relative to this script)
addpath(fullfile(fileparts(mfilename('fullpath')), 'dstoolbox'));

fprintf('4th-Order CIFF DSM Simulation\n');
fprintf('=============================\n\n');

%% Specifications
order = 4;
OSR = 16;
H_inf = 4.0;
opt = 1;
fs = 2e6;
form = 'CIFF';

%% Step 1-3: Design
ntf = synthesizeNTF(order, OSR, opt, H_inf, 0);
[a, g, b, c] = realizeNTF(ntf, form);
ABCD = stuffABCD(a, g, b, c, form);

%% Step 4-5: Simulate
N = 8192;
fB = fs / (2*OSR);
f_bin = round(sqrt(1/7) * N / (2*OSR));
A_in = 0.5;
u = A_in * sin(2*pi*f_bin*(0:N-1)/N);

n_bits = 5;
n_levels = 2^n_bits;
V_fs = 1.0;
LSB = 2*V_fs / (n_levels - 1);
q_levels = linspace(-V_fs, V_fs, n_levels);

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
    if i == 1, v_prev = 0; else, v_prev = v(i-1); end
    y(i) = C_mat * x + D_mat(1)*u(i) + D_mat(2)*v_prev;
    y_clip = max(-V_fs, min(V_fs, y(i)));
    idx = round((y_clip - (-V_fs)) / LSB) + 1;
    idx = max(1, min(n_levels, idx));
    v(i) = q_levels(idx);
    x = A_mat*x + B_mat(:,1)*u(i) + B_mat(:,2)*v(i);
    x_max = max(x_max, abs(x));
end

%% Step 6: SNR
w = 0.5 * (1 - cos(2*pi*(0:N-1)/N));
V_out = fft(v .* w) / (N/4);
V_out_mag = abs(V_out);

[~, sig_idx] = max(V_out_mag(2:N/2));
sig_bin = sig_idx + 1;
fB_bins = ceil(N / (2*OSR));

sig_bins = sig_bin-1:sig_bin+1;
sig_bins = sig_bins(sig_bins >= 2 & sig_bins <= fB_bins);
signal_power = sum(V_out_mag(sig_bins).^2);
noise_bins = setdiff(2:fB_bins, sig_bins);
noise_power = sum(V_out_mag(noise_bins).^2);

SNR = 10*log10(signal_power / noise_power);
ENOB = (SNR - 1.76) / 6.02;

%% Results
fprintf('Coefficients:\n');
fprintf('  a = ['); fprintf('%.6f ', a); fprintf(']\n');
fprintf('  g = ['); fprintf('%.6f ', g); fprintf(']\n');
fprintf('  b = ['); fprintf('%.6f ', b); fprintf(']\n');
fprintf('  c = ['); fprintf('%.6f ', c); fprintf(']\n\n');

fprintf('ABCD Matrix:\n');
for i = 1:size(ABCD, 1)
    fprintf('  ['); fprintf('%10.6f ', ABCD(i, :)); fprintf(']\n');
end
fprintf('\n');

fprintf('Simulation Results:\n');
fprintf('  Input: %.2f V @ %.2f kHz\n', A_in, f_bin/N*fs/1000);
fprintf('  Output range: [%.4f, %.4f]\n', min(v), max(v));
fprintf('  State maxima: '); fprintf('%.4f ', x_max); fprintf('\n');
fprintf('  SNR: %.2f dB\n', SNR);
fprintf('  ENOB: %.2f bits\n', ENOB);
if max(x_max) < 100, status_str = 'STABLE'; else, status_str = 'UNSTABLE'; end
fprintf('  Status: %s\n', status_str);

%% Save data
freqs = (0:N/2)/N*fs;
save('dsm_4th_order_results.mat', 'a', 'g', 'b', 'c', 'ABCD', 'SNR', 'ENOB', ...
     'V_out_mag', 'freqs', 'u', 'v', 'x_max');
fprintf('\nSaved: dsm_4th_order_results.mat\n');

%% Simple ASCII spectrum plot
fprintf('\nOutput Spectrum (first 20 bins):\n');
fprintf('Bin    Freq(kHz)   dBFS\n');
fprintf('------------------------\n');
for i = 2:21
    db = 20*log10(V_out_mag(i) + eps);
    marker = '';
    if i == sig_bin, marker = ' <-- SIGNAL'; end
    fprintf('%3d    %8.2f   %6.1f%s\n', i, freqs(i)/1000, db, marker);
end

fprintf('\nDone!\n');
