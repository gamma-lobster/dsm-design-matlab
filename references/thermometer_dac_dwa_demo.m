%% THERMOMETER_DAC_DWA_DEMO
% Demonstration of Dynamic Weighted Averaging (DWA) for DAC mismatch shaping
%
% This demo shows how DWA rotates unit cell selection to average out
% mismatch errors over time, compared to static selection.

%% Clear workspace
clc;
clear;
close all;

%% Parameters
v_fs = 1.0;        % +/- 1V full scale
n_bits = 4;        % 4-bit resolution (15 cells) - more cells for better DWA effect
n_levels = 2^n_bits;
n_comparators = n_levels - 1;
fs = 1e6;          % 1 MHz sampling
duration = 0.001;  % 1ms duration
f_signal = 5e3;    % 5 kHz sine wave
amplitude = 0.8;   % Input amplitude

% Mismatch parameters
mismatch_pct = 0.02;  % 2% standard deviation
rng(42);              % Fixed seed for reproducibility

%% Calculate derived values
num_samples = round(fs * duration);
t = (0:num_samples-1) / fs;

%% Generate unit cell weights with mismatch
nominal_weight = (2 * v_fs) / n_comparators;
bit_weights = nominal_weight * (1 + mismatch_pct * randn(1, n_comparators));

fprintf('========================================\n');
fprintf('DWA DAC DEMO\n');
fprintf('========================================\n\n');

fprintf('Configuration:\n');
fprintf('  Resolution: %d bits (%d levels)\n', n_bits, n_levels);
fprintf('  Unit cells: %d\n', n_comparators);
fprintf('  Nominal weight: %.4fV\n', nominal_weight);
fprintf('  Mismatch: %.1f%% std dev\n\n', mismatch_pct * 100);

%% Generate input signal (ramp for clear visualization)
v_in = amplitude * sin(2*pi*f_signal*t);

%% ADC: Quantize with flash ADC
fprintf('Running Flash ADC...\n');
thermometer = zeros(num_samples, n_comparators);
for i = 1:num_samples
    [therm, ~, ~] = flash_adc_quantizer(v_in(i), v_fs, n_bits);
    thermometer(i, :) = therm;
end

%% DAC: Compare three modes
fprintf('Running DAC simulations...\n\n');

% 1. Ideal DAC (no mismatch, no DWA)
v_ideal = thermometer_dac(thermometer, v_fs, n_bits);

% 2. Static DAC with mismatch (conventional thermometer)
% In static mode, cells 1:k are always used for code k
v_static_mismatch = zeros(1, num_samples);
for i = 1:num_samples
    k = sum(thermometer(i, :));
    contribution = 0;
    for j = 1:n_comparators
        if j <= k
            contribution = contribution + bit_weights(j);  % Selected: +weight
        else
            contribution = contribution - bit_weights(j);  % Not selected: -weight
        end
    end
    v_static_mismatch(i) = contribution * 0.5;
end

% 3. DWA DAC with mismatch
v_dwa_mismatch = thermometer_dac_dwa(thermometer, v_fs, n_bits, bit_weights, 1);

%% Analyze cell usage patterns
fprintf('Analyzing cell usage patterns...\n\n');

% For static DAC
cell_usage_static = zeros(1, n_comparators);
for i = 1:num_samples
    k = sum(thermometer(i, :));
    cell_usage_static(1:k) = cell_usage_static(1:k) + 1;
end

% For DWA DAC
cell_usage_dwa = zeros(1, n_comparators);
start_idx = 1;
for i = 1:num_samples
    k = sum(thermometer(i, :));
    if k > 0 && k < n_comparators
        for j = 0:k-1
            idx = mod(start_idx - 1 + j, n_comparators) + 1;
            cell_usage_dwa(idx) = cell_usage_dwa(idx) + 1;
        end
        start_idx = mod(start_idx - 1 + k, n_comparators) + 1;
    elseif k == n_comparators
        cell_usage_dwa = cell_usage_dwa + 1;
    end
end

%% Display results
fprintf('========================================\n');
fprintf('CELL USAGE COMPARISON\n');
fprintf('========================================\n');
fprintf('Cell#  Weight(%%)    Static Count    DWA Count\n');
fprintf('-----  ----------   ------------    ----------\n');
for i = 1:n_comparators
    weight_dev = (bit_weights(i) - nominal_weight) / nominal_weight * 100;
    fprintf('%3d    %+6.2f       %8d        %8d\n', ...
        i, weight_dev, cell_usage_static(i), cell_usage_dwa(i));
end
fprintf('\n');
fprintf('Static: Cells 1-%d used heavily, cells %d-%d barely used\n', ...
    round(n_comparators/2), round(n_comparators/2)+1, n_comparators);
fprintf('DWA:    All cells used more uniformly\n\n');

%% Calculate spectrum and SNDR for mismatch cases
N = 2^nextpow2(num_samples);
w = hann(num_samples)';

% Static mismatch spectrum
V_static = fft(v_static_mismatch .* w, N) / (num_samples/2);
freqs = (0:N/2-1) * fs / N;

% DWA spectrum
V_dwa = fft(v_dwa_mismatch .* w, N) / (num_samples/2);

% Calculate SNDR for both (signal at f_signal)
sig_bin = round(f_signal * N / fs) + 1;
noise_bins = setdiff(2:N/2, sig_bin-2:sig_bin+2);  % Exclude signal and nearby

% Static SNDR
sig_power_static = sum(abs(V_static(sig_bin-2:sig_bin+2)).^2);
noise_power_static = sum(abs(V_static(noise_bins)).^2);
SNDR_static = 10*log10(sig_power_static / noise_power_static);

% DWA SNDR
sig_power_dwa = sum(abs(V_dwa(sig_bin-2:sig_bin+2)).^2);
noise_power_dwa = sum(abs(V_dwa(noise_bins)).^2);
SNDR_dwa = 10*log10(sig_power_dwa / noise_power_dwa);

fprintf('========================================\n');
fprintf('SNDR COMPARISON (with %.1f%% mismatch)\n', mismatch_pct*100);
fprintf('========================================\n');
fprintf('Static (no DWA): %.2f dB\n', SNDR_static);
fprintf('DWA:             %.2f dB\n', SNDR_dwa);
fprintf('Improvement:     %.2f dB\n\n', SNDR_dwa - SNDR_static);

%% Visualizations
figure('Name', 'DWA Demo - Cell Usage and Spectrum', 'Position', [50 50 1600 1200]);

% Plot 1: Unit cell weights
subplot(3, 3, 1);
bar(1:n_comparators, bit_weights, 'FaceColor', [0.4 0.6 0.8]);
hold on;
yline(nominal_weight, 'r--', 'LineWidth', 2, 'Label', 'Nominal');
hold off;
xlabel('Cell Index');
ylabel('Weight (V)');
title(sprintf('Unit Cell Weights (%.1f%% mismatch)', mismatch_pct*100));
grid on;

% Plot 2: Cell usage - Static
subplot(3, 3, 2);
bar(1:n_comparators, cell_usage_static, 'FaceColor', [0.8 0.4 0.4]);
xlabel('Cell Index');
ylabel('Usage Count');
title('Cell Usage: Static (Conventional)');
grid on;

% Plot 3: Cell usage - DWA
subplot(3, 3, 3);
bar(1:n_comparators, cell_usage_dwa, 'FaceColor', [0.4 0.8 0.4]);
xlabel('Cell Index');
ylabel('Usage Count');
title('Cell Usage: DWA (Rotated)');
grid on;

% Plot 4: Time domain - Input
n_plot = min(500, num_samples);
subplot(3, 3, 4);
plot(t(1:n_plot)*1e6, v_in(1:n_plot), 'b-', 'LineWidth', 1);
xlabel('Time (µs)');
ylabel('Voltage (V)');
title('Input Signal');
grid on;

% Plot 5: Time domain - Static mismatch
subplot(3, 3, 5);
plot(t(1:n_plot)*1e6, v_static_mismatch(1:n_plot), 'r-', 'LineWidth', 0.8);
xlabel('Time (µs)');
ylabel('Voltage (V)');
title(sprintf('Static DAC Output (SNDR=%.1f dB)', SNDR_static));
grid on;

% Plot 6: Time domain - DWA
subplot(3, 3, 6);
plot(t(1:n_plot)*1e6, v_dwa_mismatch(1:n_plot), 'g-', 'LineWidth', 0.8);
xlabel('Time (µs)');
ylabel('Voltage (V)');
title(sprintf('DWA DAC Output (SNDR=%.1f dB)', SNDR_dwa));
grid on;

% Plot 7: Spectrum comparison
subplot(3, 3, 7);
semilogx(freqs(2:end)/1e3, 20*log10(abs(V_static(2:N/2)) + eps), 'r-', 'LineWidth', 0.8);
hold on;
semilogx(freqs(2:end)/1e3, 20*log10(abs(V_dwa(2:N/2)) + eps), 'g-', 'LineWidth', 0.8);
xline(f_signal/1e3, 'm--', 'Signal', 'LineWidth', 1);
hold off;
xlabel('Frequency (kHz)');
ylabel('Magnitude (dB)');
title('Spectrum Comparison');
legend('Static', 'DWA', 'Location', 'southwest');
grid on;
axis([100 fs/2e3 -140 -20]);

% Plot 8: DWA rotation visualization (first 20 samples)
subplot(3, 3, 8);
n_viz = min(20, num_samples);
start_idx = 1;
dwa_matrix = zeros(n_viz, n_comparators);
for i = 1:n_viz
    k = sum(thermometer(i, :));
    if k > 0 && k < n_comparators
        for j = 0:k-1
            idx = mod(start_idx - 1 + j, n_comparators) + 1;
            dwa_matrix(i, idx) = 1;
        end
        start_idx = mod(start_idx - 1 + k, n_comparators) + 1;
    elseif k == n_comparators
        dwa_matrix(i, :) = 1;
    end
end
imagesc(1:n_viz, 1:n_comparators, dwa_matrix');
colormap(gca, [1 0.8 0.8; 0.2 0.7 0.2]);
set(gca, 'YDir', 'normal');
xlabel('Sample #');
ylabel('Cell #');
title('DWA Cell Selection Pattern (First 20 samples)');
colorbar('Ticks', [0.25, 0.75], 'TickLabels', {'OFF', 'ON'});

% Plot 9: Summary
subplot(3, 3, 9);
axis off;
text(0.1, 0.9, 'DWA SUMMARY', 'FontSize', 14, 'FontWeight', 'bold');
text(0.1, 0.75, sprintf('Resolution: %d bits (%d cells)', n_bits, n_comparators), 'FontSize', 11);
text(0.1, 0.65, sprintf('Mismatch: %.1f%%', mismatch_pct*100), 'FontSize', 11);
text(0.1, 0.50, sprintf('SNDR Improvement:\n  Static: %.1f dB\n  DWA: %.1f dB\n  Gain: +%.1f dB', ...
    SNDR_static, SNDR_dwa, SNDR_dwa - SNDR_static), 'FontSize', 11, 'Color', 'b');
text(0.1, 0.25, 'DWA averages mismatch by\nrotating cell selection', 'FontSize', 10);
text(0.1, 0.10, 'Harmonics pushed to fs/m\nwhere m = # of cells', 'FontSize', 10);

%% Save
saveas(gcf, 'dwa_dac_demo.png');
fprintf('Saved: dwa_dac_demo.png\n\n');

fprintf('========================================\n');
fprintf('DEMO COMPLETE\n');
fprintf('========================================\n');
fprintf('DWA effectively:\n');
fprintf('  1. Equalizes cell usage across all cells\n');
fprintf('  2. Shapes mismatch noise out of signal band\n');
fprintf('  3. Pushes distortion to higher frequencies\n');
fprintf('  4. Improves SNDR by %.1f dB in this example\n', SNDR_dwa - SNDR_static);
fprintf('========================================\n');