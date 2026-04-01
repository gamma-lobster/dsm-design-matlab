%% THERMOMETER_DAC_DEMO
% Demonstration of Flash ADC followed by Thermometer DAC
%
% This shows that the DAC properly reconstructs the quantized signal
% from the thermometer code output of the ADC.

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

%% ADC: Quantize with flash ADC
fprintf('Running Flash ADC...\n');
[thermometer, binary_adc, thresholds] = flash_adc_quantizer(v_in, v_fs, n_bits);

%% DAC: Convert thermometer code back to analog
fprintf('Running Thermometer DAC...\n');
binary_dac = thermometer_dac(thermometer, v_fs, n_bits);

%% Verify DAC matches ADC output
fprintf('\n========================================\n');
fprintf('ADC-DAC VERIFICATION\n');
fprintf('========================================\n');
fprintf('Checking if DAC output matches ADC binary output...\n');
max_error = max(abs(binary_adc(:) - binary_dac(:)));
fprintf('Maximum difference: %e\n', max_error);
if max_error < 1e-10
    fprintf('✓ DAC output matches ADC output (within numerical precision)\n');
else
    fprintf('✗ Warning: Difference detected!\n');
end

%% Display configuration
fprintf('\n========================================\n');
fprintf('CONFIGURATION\n');
fprintf('========================================\n');
fprintf('Full Scale: +/- %.2fV\n', v_fs);
fprintf('Resolution: %d bits\n', n_bits);
fprintf('Quantization Levels: %d\n', n_levels);
fprintf('Number of Comparators: %d\n', n_comparators);
fprintf('Bit Weight: %.4fV\n', (2*v_fs)/n_comparators);
fprintf('\nDAC Operation:\n');
fprintf('  Each thermometer bit has weight = ±%.4fV\n', (2*v_fs)/n_comparators);
fprintf('  Bit = 1: contribute +weight\n');
fprintf('  Bit = 0: contribute -weight\n');
fprintf('  Output scaled to +/-%.2fV range\n', v_fs);

%% Sample values showing ADC-DAC chain
fprintf('\n========================================\n');
fprintf('ADC-DAC CHAIN SAMPLE VALUES\n');
fprintf('========================================\n');
fprintf('%-12s %-12s %-20s %-12s %-12s\n', ...
    'Time (ms)', 'Vin (V)', 'Thermometer', 'ADC Out', 'DAC Out');
fprintf('%-12s %-12s %-20s %-12s %-12s\n', ...
    '----------', '----------', '--------------------', '----------', '----------');

num_display = min(16, num_samples);
step = max(1, floor(num_samples / num_display));
for i = 1:step:num_samples
    therm_str = sprintf('%d', thermometer(i, :));
    fprintf('%-12.3f %-12.4f %-20s %-12.4f %-12.4f\n', ...
        t(i)*1000, v_in(i), therm_str, binary_adc(i), binary_dac(i));
end

%% Create Visualizations
figure('Name', 'Flash ADC + Thermometer DAC Chain', 'Position', [100 100 1400 1000]);

% Plot 1: Input signal with comparator thresholds
subplot(5, 1, 1);
plot(t*1000, v_in, 'b-', 'LineWidth', 1.5);
hold on;
colors = jet(n_comparators);
for i = 1:n_comparators
    yline(thresholds(i), '--', 'Color', colors(i,:), 'Alpha', 0.5, 'LineWidth', 0.7);
end
hold off;
ylabel('Voltage (V)');
title(sprintf('Input Signal with Comparator Thresholds (%d-bit Flash ADC)', n_bits));
grid on;
ylim([-v_fs-0.15, v_fs+0.15]);

% Plot 2: Thermometer code as 2D matrix
subplot(5, 1, 2);
imagesc(t*1000, 0:n_comparators-1, thermometer');
colormap(gca, [1 0.3 0.3; 0.2 1 0.2]);
set(gca, 'YDir', 'normal');
ylabel('Comparator #');
title('Thermometer Code Output from Flash ADC');
colorbar('Ticks', [0.25, 0.75], 'TickLabels', {'LOW', 'HIGH'});

% Plot 3: ADC output vs Input
subplot(5, 1, 3);
plot(t*1000, v_in, 'b-', 'LineWidth', 1, 'DisplayName', 'Input');
hold on;
stairs(t*1000, binary_adc, 'r-', 'LineWidth', 1.5, 'DisplayName', 'ADC Output');
hold off;
ylabel('Voltage (V)');
title('Flash ADC: Quantized Output vs Input');
legend('Location', 'best');
grid on;

% Plot 4: DAC output
subplot(5, 1, 4);
plot(t*1000, v_in, 'b-', 'LineWidth', 1, 'DisplayName', 'Original Input');
hold on;
stairs(t*1000, binary_dac, 'm-', 'LineWidth', 1.5, 'DisplayName', 'DAC Output');
stairs(t*1000, binary_adc, 'r--', 'LineWidth', 1, 'DisplayName', 'ADC Output');
hold off;
ylabel('Voltage (V)');
title('Thermometer DAC: Reconstructed Output');
legend('Location', 'best');
grid on;

% Plot 5: ADC-DAC difference
subplot(5, 1, 5);
adc_error = v_in - binary_adc;
dac_error = v_in - binary_dac;
plot(t*1000, adc_error, 'g-', 'LineWidth', 0.8, 'DisplayName', 'ADC Quantization Error');
hold on;
plot(t*1000, dac_error, 'c--', 'LineWidth', 0.8, 'DisplayName', 'DAC Error');
hold off;
ylabel('Error (V)');
xlabel('Time (ms)');
title('Quantization Error: Input - Quantized Output');
legend('Location', 'best');
grid on;
yline(0, 'k--', 'Alpha', 0.3);

%% Explanation
fprintf('\n========================================\n');
fprintf('DAC OPERATION EXPLANATION\n');
fprintf('========================================\n');
fprintf('The thermometer DAC works as follows:\n');
fprintf('  1. Count the number of 1s and 0s in the thermometer code\n');
fprintf('  2. Each bit has weight = %.4fV\n', (2*v_fs)/n_comparators);
fprintf('  3. Contribution = (num_ones - num_zeros) * weight\n');
fprintf('  4. Output = Contribution / 2 (to scale to +/-%.1fV)\n\n', v_fs);

fprintf('Example calculations:\n');
for level = [0, 2, 4, 7]
    ones_count = level;
    zeros_count = n_comparators - level;
    weight = (2*v_fs)/n_comparators;
    raw_sum = (ones_count - zeros_count) * weight;
    output = raw_sum / 2;
    therm = [ones(1, ones_count), zeros(1, zeros_count)];
    fprintf('  Thermometer %s:\n', sprintf('%d', therm));
    fprintf('    %d ones, %d zeros → (%d - %d) * %.4f = %.4f → /2 = %+.4fV\n', ...
        ones_count, zeros_count, ones_count, zeros_count, weight, raw_sum, output);
end

fprintf('\n========================================\n');
fprintf('DEMO COMPLETE\n');
fprintf('========================================\n');
