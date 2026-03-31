%% dsm_quick_design.m
%% Quick-start template for Delta-Sigma Modulator Design
%% Native MATLAB workflow using Delta Sigma Toolbox
%%
%% Edit the parameters below and run!

clear all; close all; clc;

%% Add Delta Sigma Toolbox to path (relative to this script)
addpath(fullfile(fileparts(mfilename('fullpath')), 'dstoolbox'));

%% ============== EDIT THESE PARAMETERS ==============

% Modulator order
order = 3;

% Oversampling ratio
OSR = 64;

% Maximum out-of-band NTF gain (Lee's rule: keep < 2.0 for stability)
H_inf = 1.5;

% Zero optimization (0 = none, 1 = optimized)
opt = 1;

% Topology: 'CIFF', 'CIFB', 'CRFB', 'CRFF'
form = 'CIFF';

% Sampling frequency (Hz)
fs = 1e6;

% Quantizer bits
n_bits = 5;

% Input amplitude (0.0 to 1.0)
A_in = 0.5;

% FFT points for simulation
N = 8192;

%% ============== DESIGN FLOW ==============

fprintf('Delta-Sigma Modulator Design\n');
fprintf('===========================\n\n');

% Step 1: Synthesize NTF
fprintf('Synthesizing NTF...\n');
ntf = synthesizeNTF(order, OSR, opt, H_inf, 0);

% Step 2: Realize coefficients
fprintf('Realizing coefficients...\n');
[a, g, b, c] = realizeNTF(ntf, form);

% Step 3: Build ABCD
fprintf('Building ABCD matrix...\n');
ABCD = stuffABCD(a, g, b, c, form);

% Display results
fprintf('\n=== Results ===\n');
fprintf('Coefficients:\n');
fprintf('  a = ['); fprintf('%.6f ', a); fprintf(']\n');
fprintf('  g = ['); fprintf('%.6f ', g); fprintf(']\n');
fprintf('  b = ['); fprintf('%.6f ', b); fprintf(']\n');
fprintf('  c = ['); fprintf('%.6f ', c); fprintf(']\n');

% Step 4: Simulate
fprintf('\nRunning simulation...\n');

% Extract state matrices
[n_rows, ~] = size(ABCD);
n_states = n_rows - 1;
A_mat = ABCD(1:n_states, 1:n_states);
B_mat = ABCD(1:n_states, n_states+1:end);
C_mat = ABCD(n_states+1, 1:n_states);
D_mat = ABCD(n_states+1, n_states+1:end);

% Simulation parameters
fB = fs / (2*OSR);
f_bin = round(sqrt(1/7) * N / (2*OSR));
u = A_in * sin(2*pi*f_bin*(0:N-1)/N);

n_levels = 2^n_bits;
V_fs = 1.0;
LSB = 2*V_fs / (n_levels - 1);
q_levels = linspace(-V_fs, V_fs, n_levels);

% Run
x = zeros(n_states, 1);
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
end

% Step 5: Calculate SNR
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

fprintf('\n=== Performance ===\n');
fprintf('SNR:  %.2f dB\n', SNR);
fprintf('ENOB: %.2f bits\n', ENOB);

% Quick plot
figure;

% Spectrum
subplot(1, 2, 1);
freqs = (0:N/2)/N*fs/1000;
V_out_dB = 20*log10(V_out_mag(1:N/2+1) + eps);
plot(freqs, V_out_dB);
hold on;
plot([fB/1000 fB/1000], [-140 10], 'r--');
hold off;
xlabel('Frequency (kHz)');
ylabel('dBFS');
title(sprintf('Output Spectrum (SNR=%.1f dB)', SNR));
legend('Spectrum', 'Signal BW');
grid on;
axis([0 fs/2000 -140 10]);

% NTF
subplot(1, 2, 2);
f = logspace(log10(100), log10(fs/2), 1000);
[~, ~, k] = zpkdata(ntf);
z_ntf = zpkdata(ntf, 'z');
p_ntf = zpkdata(ntf, 'p');

N_fft = 8192;
freqs_norm = (0:N_fft-1)/N_fft;
z_exp = exp(2*pi*1j*freqs_norm);
NTF_resp = ones(1, N_fft);
for i = 1:N_fft
    zi = z_exp(i);
    num = k * prod(zi - cell2mat(z_ntf));
    den = prod(zi - cell2mat(p_ntf));
    NTF_resp(i) = num / den;
end
NTF_dB = 20*log10(abs(NTF_resp) + eps);
freqs_plot = freqs_norm * fs / 1000;

semilogx(freqs_plot, NTF_dB);
hold on;
plot([fB/1000 fB/1000], [-100 40], 'r--');
hold off;
xlabel('Frequency (kHz)');
ylabel('Magnitude (dB)');
title('NTF Response');
legend('NTF', 'Signal BW');
grid on;

fprintf('\nDone!\n');
