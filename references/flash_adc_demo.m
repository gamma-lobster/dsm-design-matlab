%% FLASH_ADC_DEMO
% Demonstration of Flash ADC Quantizer with Thermometer Code Output
%
% Run this script to see the flash ADC in action with plots

%% Clear workspace
clc;
clear;
close all;

%% Parameters
v_fs = 1.0;        % +/- 1V full scale
n_bits = 3;        % 3-bit resolution
fs = 1000;         % Sampling frequency (Hz)
duration = 0.02;   % 20ms duration
f_signal = 50;     % 50 Hz sine wave
amplitude = 0.9;   % Input amplitude

%% Calculate derived values
n_levels = 2^n_bits;
n_comparators = n_levels - 1;
num_samples = round(fs * duration);

%% Generate time vector and input signal
t = (0:num_samples-1) / fs;
v_in = amplitude * sin(2*pi*f_signal*t);
% Add a third harmonic for interest
v_in = v_in + 0.2 * sin(2*pi*3*f_signal*t);

%% Quantize with flash ADC
fprintf('Quantizing %d samples...\n', num_samples);
[thermometer, binary, thresholds] = flash_adc_quantizer(v_in, v_fs, n_bits);

%% Display Configuration
fprintf('\n========================================\n');
fprintf('FLASH ADC CONFIGURATION\n');
fprintf('========================================\n');
fprintf('Full Scale: +/- %.2fV\n', v_fs);
fprintf('Resolution: %d bits\n', n_bits);
fprintf('Quantization Levels: %d\n', n_levels);
fprintf('Number of Comparators: %d\n', n_comparators);
fprintf('LSB: %.4fV\n', (2*v_fs)/(n_levels-1));

fprintf('\nComparator Thresholds:\n');
for i = 1:n_comparators
    fprintf('  Comp %d: Vth = %+.4fV\n', i-1, thresholds(i));
end

output_levels = linspace(-v_fs, v_fs, n_levels);
fprintf('\nOutput Levels:\n');
disp(output_levels);

%% Display Sample Values
fprintf('\n========================================\n');
fprintf('SAMPLE OUTPUT VALUES\n');
fprintf('========================================\n');
fprintf('%-12s %-12s %-20s %-12s %-6s\n', 'Time (ms)', 'Vin (V)', 'Thermometer', 'Vout (V)', 'Level');
fprintf('%-12s %-12s %-20s %-12s %-6s\n', '----------', '----------', '--------------------', '----------', '------');

% Show samples at evenly spaced points
num_display = min(20, num_samples);
step = max(1, floor(num_samples / num_display));
for i = 1:step:num_samples
    therm_str = sprintf('%d', thermometer(i, :));
    level = sum(thermometer(i, :));
    fprintf('%-12.3f %-12.4f %-20s %-12.4f %-6d\n', t(i)*1000, v_in(i), therm_str, binary(i), level);
end

%% Create Visualizations
figure('Name', 'Flash ADC Thermometer Code Output', 'Position', [100 100 1200 900]);

% Plot 1: Input signal with comparator thresholds
subplot(4, 1, 1);
plot(t*1000, v_in, 'b-', 'LineWidth', 1.5);
hold on;
colors = jet(n_comparators);
for i = 1:n_comparators
    yline(thresholds(i), '--', 'Color', colors(i,:), 'Alpha', 0.6, 'LineWidth', 0.8);
end
hold off;
ylabel('Voltage (V)');
title(sprintf('Flash ADC: Input Signal and Comparator Thresholds (%d-bit, +/-%.1fV)', n_bits, v_fs));
legend('Input Signal', 'Location', 'eastoutside');
grid on;
ylim([-v_fs-0.1, v_fs+0.1]);

% Plot 2: Thermometer code as 2D matrix
subplot(4, 1, 2);
imagesc(t*1000, 0:n_comparators-1, thermometer');
colormap(gca, [1 0.3 0.3; 0.2 1 0.2]);  % Red for 0, Green for 1
set(gca, 'YDir', 'normal');
ylabel('Comparator #');
title('Thermometer Code Output (2D Matrix) - Green = HIGH (vin > vth), Red = LOW');
xlabel('Time (ms)');
colorbar('Ticks', [0.25, 0.75], 'TickLabels', {'LOW', 'HIGH'});

% Plot 3: Binary output vs input
subplot(4, 1, 3);
plot(t*1000, v_in, 'b-', 'LineWidth', 1, 'DisplayName', 'Input');
hold on;
stairs(t*1000, binary, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Quantized Output');
hold off;
ylabel('Voltage (V)');
title('Binary Quantized Output vs Input');
legend('Location', 'best');
grid on;

% Plot 4: Quantization error
subplot(4, 1, 4);
error = v_in - binary;
plot(t*1000, error, 'g-', 'LineWidth', 0.8);
ylabel('Error (V)');
xlabel('Time (ms)');
title('Quantization Error');
grid on;
yline(0, 'k--', 'Alpha', 0.3);

%% Thermometer Code Explanation
fprintf('\n========================================\n');
fprintf('THERMOMETER CODE EXPLANATION\n');
fprintf('========================================\n');
fprintf('In a flash ADC, all comparators operate in parallel:\n');
fprintf('  * Each comparator has a different threshold voltage\n');
fprintf('  * If Vin > Vth, that comparator outputs 1 (HIGH)\n');
fprintf('  * If Vin <= Vth, that comparator outputs 0 (LOW)\n');
fprintf('  * The pattern of 1s and 0s is ''thermometer code''\n');
fprintf('  * The binary output = count the number of 1s\n\n');

fprintf('Thermometer code examples:\n');
for level = 0:n_comparators
    therm = [ones(1, level), zeros(1, n_comparators-level)];
    therm_str = sprintf('%d', therm);
    fprintf('  %s → %d comparators HIGH → Output = %+.2fV (Level %d)\n', ...
        therm_str, level, output_levels(level+1), level);
end

fprintf('\nKey advantage: Bubble immunity\n');
fprintf('  A single comparator error (a "bubble") is easily detected\n');
fprintf('  because valid thermometer code never has 010 or 101 patterns\n');

fprintf('\nDemo complete!\n');
