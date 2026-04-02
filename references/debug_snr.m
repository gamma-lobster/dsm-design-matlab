%% DEBUG SNR CALCULATION
% Quick script to verify SNR calculation against visual noise floor

clear; close all; clc;

%% Parameters matching the DSM simulation
N = 8192;
OSR = 16;
fs = 2e6;
fB = fs / (2*OSR);  % Signal bandwidth

%% Create a simple test signal
f_signal = round(sqrt(1/7) * N / (2*OSR));  % Same as DSM
A_signal = 0.5;  % Same amplitude as DSM

n = 0:N-1;
u = A_signal * sin(2*pi*f_signal*n/N);

% Add some white noise to simulate the expected noise floor
% Target: noise floor at -110 dB with signal at -6 dB
noise_floor_dB = -110;
signal_dB = 20*log10(A_signal);  % Should be -6 dB for A=0.5

% Noise power per bin
noise_power_per_bin = 10^(noise_floor_dB/10);
noise_total_power = noise_power_per_bin * (N/2);  % Total noise power
noise_std = sqrt(noise_total_power);

% Generate noise
noise = noise_std * randn(1, N);

% Create output = signal + noise
v = u + noise;

%% Calculate SNR using our function
w = 0.5 * (1 - cos(2*pi*(0:N-1)/N));
V_out = fft(v .* w) / (N/4);
V_out_mag = abs(V_out);

% Find signal
[~, sig_idx] = max(V_out_mag(2:N/2));
sig_bin = sig_idx + 1;
fB_bins = ceil(N / (2*OSR));

fprintf('========================================\n');
fprintf('DEBUG SNR CALCULATION\n');
fprintf('========================================\n\n');
fprintf('Parameters:\n');
fprintf('  N = %d (FFT points)\n', N);
fprintf('  OSR = %d\n', OSR);
fprintf('  fB_bins = %d (bins in signal band)\n', fB_bins);
fprintf('  Signal bin = %d\n', sig_bin);
fprintf('  Signal frequency = %.2f kHz\n\n', sig_bin/N * fs / 1000);

% Signal power (3 bins)
sig_bins = sig_bin-1:sig_bin+1;
sig_bins = sig_bins(sig_bins >= 2 & sig_bins <= fB_bins);
signal_power = sum(V_out_mag(sig_bins).^2);
fprintf('Signal power (3 bins): %.4e\n', signal_power);
fprintf('Signal magnitude: %.4f (%.2f dBFS)\n', ...
    mean(V_out_mag(sig_bins)), 20*log10(mean(V_out_mag(sig_bins)) + eps));

% Noise power
harmonic_bins = [];
for h = 2:7
    harmonic_bin = h * sig_bin;
    if harmonic_bin <= fB_bins
        harmonic_bins = [harmonic_bins, harmonic_bin-1:harmonic_bin+1];
    end
end
harmonic_bins = unique(harmonic_bins);
harmonic_bins = harmonic_bins(harmonic_bins >= 2 & harmonic_bins <= fB_bins);

exclude_bins = unique([sig_bins, harmonic_bins]);
noise_bins = setdiff(2:fB_bins, exclude_bins);
noise_power = sum(V_out_mag(noise_bins).^2);

fprintf('\nNoise power (%d bins): %.4e\n', length(noise_bins), noise_power);
fprintf('Average noise per bin: %.4e (%.2f dBFS)\n', ...
    noise_power/length(noise_bins), 10*log10(noise_power/length(noise_bins) + eps));

SNR = 10*log10(signal_power / noise_power);
fprintf('\nCalculated SNR: %.2f dB\n', SNR);

%% Visual check
figure('Name', 'Debug SNR');
freqs = (0:N/2)/N*fs/1000;
semilogx(freqs, 20*log10(abs(V_out(1:N/2+1)) + eps), 'b-', 'LineWidth', 0.5);
hold on;
xline(fB/1000, 'r--', 'Signal BW');
xline(sig_bin/N*fs/1000, 'g--', 'Signal');
hold off;
xlabel('Frequency (kHz)');
ylabel('Magnitude (dBFS)');
title(sprintf('Debug: Expected SNR ~70dB, Calculated %.2f dB', SNR));
grid on;
axis([100 fs/2000 -140 10]);