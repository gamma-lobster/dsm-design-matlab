%% THERMOMETER_DAC_MISMATCH_DEMO
% Demonstration of DAC Unit Cell Mismatch Effects
%
% Shows how random variations in unit cell weights affect DAC linearity
% and introduce distortion/noise in the reconstructed signal.

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

%% Mismatch parameters
mismatch_pct = 0.02;  % 2% unit cell mismatch
seed = 42;            % Random seed for reproducibility

%% Calculate derived values
n_levels = 2^n_bits;
n_comparators = n_levels - 1;
num_samples = round(fs * duration);

%% Generate time vector and input signal
t = (0:num_samples-1) / fs;
v_in = amplitude * sin(2*pi*f_signal*t);

%% ADC: Quantize with flash ADC
fprintf('Running Flash ADC...\n');
[thermometer, binary_ideal, thresholds] = flash_adc_quantizer(v_in, v_fs, n_bits);

%% DAC: Ideal (no mismatch)
fprintf('Running Ideal DAC...\n');
binary_ideal_dac = thermometer_dac(thermometer, v_fs, n_bits);

%% DAC: With mismatch
fprintf('Running DAC with %.1f%% unit cell mismatch...\n', mismatch_pct*100);
[binary_mismatch, bit_weights] = thermometer_dac_mismatch(thermometer, v_fs, n_bits, mismatch_pct, seed);

%% Display bit weights
fprintf('\n========================================\n');
fprintf('UNIT CELL WEIGHTS\n');
fprintf('========================================\n');
nominal_weight = (2 * v_fs) / n_comparators;
fprintf('Nominal weight per cell: %.6fV\n', nominal_weight);
fprintf('Actual weights (with %.1f%% mismatch):\n', mismatch_pct*100);
for i = 1:n_comparators
    deviation_pct = (bit_weights(i) - nominal_weight) / nominal_weight * 100;
    fprintf('  Cell %d: %.6fV (%.3f%% deviation)\n', i, bit_weights(i), deviation_pct);
end

%% Calculate errors
error_ideal = binary_ideal - binary_ideal_dac;
error_mismatch = binary_ideal - binary_mismatch;
inldnl = calculate_inl_dnl(bit_weights, v_fs, n_bits);

%% Display INL/DNL
fprintf('\n========================================\n');
fprintf('INTEGRAL/DIFFERENTIAL NONLINEARITY\n');
fprintf('========================================\n');
fprintf('INL (max): %.4f LSB\n', max(abs(inldnl.inl)));
fprintf('DNL (max): %.4f LSB\n', max(abs(inldnl.dnl)));

%% Create Visualizations
figure('Name', 'DAC Unit Cell Mismatch Effects', 'Position', [100 100 1400 1200]);

% Plot 1: Bit weights
subplot(5, 1, 1);
b = bar(0:n_comparators-1, bit_weights);
hold on;
yline(nominal_weight, 'r--', 'LineWidth', 2, 'Label', sprintf('Nominal: %.4fV', nominal_weight));
hold off;
xlabel('Unit Cell Index');
ylabel('Weight (V)');
title(sprintf('Unit Cell Weights (%.1f%% Mismatch)', mismatch_pct*100));
grid on;

% Plot 2: Input vs Ideal DAC vs Mismatched DAC
subplot(5, 1, 2);
plot(t*1000, v_in, 'b-', 'LineWidth', 1, 'DisplayName', 'Input');
hold on;
stairs(t*1000, binary_ideal_dac, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Ideal DAC');
stairs(t*1000, binary_mismatch, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Mismatched DAC');
hold off;
ylabel('Voltage (V)');
title('DAC Output Comparison');
legend('Location', 'best');
grid on;

% Plot 3: Quantization error (Ideal vs Mismatch)
subplot(5, 1, 3);
plot(t*1000, v_in - binary_ideal_dac, 'g-', 'LineWidth', 0.8, 'DisplayName', 'Ideal Error');
hold on;
plot(t*1000, v_in - binary_mismatch, 'r-', 'LineWidth', 0.8, 'DisplayName', 'Mismatch Error');
hold off;
ylabel('Error (V)');
title('Reconstruction Error');
legend('Location', 'best');
grid on;

% Plot 4: Mismatch-induced error only
subplot(5, 1, 4);
mismatch_error_only = binary_ideal_dac - binary_mismatch;
plot(t*1000, mismatch_error_only, 'm-', 'LineWidth', 1);
ylabel('Error (V)');
title(sprintf('Mismatch-Induced Error Only (%.1f%% mismatch)', mismatch_pct*100));
grid on;
yline(0, 'k--', 'Alpha', 0.3);

% Plot 5: INL
subplot(5, 1, 5);
stem(inldnl.code, inldnl.inl, 'filled');
xlabel('DAC Code');
ylabel('INL (LSB)');
title('Integral Nonlinearity');
grid on;
yline(0, 'k--', 'Alpha', 0.3);

%% Summary
fprintf('\n========================================\n');
fprintf('SUMMARY\n');
fprintf('========================================\n');
fprintf('Max ideal quantization error: %.4fV\n', max(abs(error_ideal)));
fprintf('Max mismatch-induced error: %.4fV\n', max(abs(mismatch_error_only)));
fprintf('SNR degradation due to mismatch: %.2f dB\n', ...
    20*log10(std(v_in - binary_ideal_dac) / std(mismatch_error_only)));

fprintf('\nMismatch Analysis:\n');
fprintf('  Unit cell variation causes:\n');
fprintf('    - Non-uniform step sizes (DNL errors)\n');
fprintf('    - Deviation from ideal transfer curve (INL errors)\n');
fprintf('    - Harmonic distortion in output spectrum\n');
fprintf('    - Reduced SFDR (spurious-free dynamic range)\n');

fprintf('\nDemo complete!\n');

%% Helper function: Calculate INL and DNL
function inldnl = calculate_inl_dnl(bit_weights, v_fs, n_bits)
    n_levels = 2^n_bits;
    n_comparators = n_levels - 1;
    nominal_weight = (2 * v_fs) / n_comparators;
    
    % Calculate actual step sizes and output levels
    % For thermometer DAC, output at code k depends on sum of first k weights
    actual_levels = zeros(1, n_levels);
    for k = 0:n_levels-1
        % For code k: k ones followed by (n_comparators-k) zeros
        contribution = 0;
        for j = 1:n_comparators
            if j <= k
                contribution = contribution + bit_weights(j);
            else
                contribution = contribution - bit_weights(j);
            end
        end
        actual_levels(k+1) = contribution * 0.5;
    end
    
    % Ideal levels
    ideal_levels = linspace(-v_fs, v_fs, n_levels);
    
    % LSB size
    lsb = (2 * v_fs) / (n_levels - 1);
    
    % INL: deviation from ideal in LSBs
    inl = (actual_levels - ideal_levels) / lsb;
    
    % DNL: deviation of step size from ideal in LSBs
    dnl = diff(actual_levels) / lsb - 1;
    dnl = [0, dnl];  % DNL at code 0 is 0
    
    inldnl.inl = inl;
    inldnl.dnl = dnl;
    inldnl.code = 0:n_levels-1;
end
