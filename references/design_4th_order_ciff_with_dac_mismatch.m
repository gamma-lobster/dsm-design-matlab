%% design_4th_order_ciff_with_dac_mismatch.m
%% 4th-Order CIFF DSM with Flash ADC and Thermometer DAC
%% Compares: No mismatch vs 0.5% DAC unit cell mismatch
%%
%% This integrates the flash ADC quantizer and thermometer DAC models
%% into the DSM simulation loop.

clear all; close all; clc;

%% Add paths
addpath(fullfile(fileparts(mfilename('fullpath')), 'dstoolbox'));

fprintf('============================================================\n');
fprintf('  4th-Order CIFF DSM with Flash ADC + Thermometer DAC\n');
fprintf('  Comparing: No Mismatch vs 0.5%% DAC Unit Cell Mismatch\n');
fprintf('============================================================\n\n');

%% DSM Design Parameters
order = 4;
OSR = 16;
H_inf = 4.0;
opt = 1;
f0 = 0;

%% Step 1: Synthesize NTF
fprintf('[Step 1] Synthesize NTF\n');
ntf = synthesizeNTF(order, OSR, opt, H_inf, f0);

%% Step 2: Realize Coefficients
form = 'CIFF';
[a, g, b, c] = realizeNTF(ntf, form);

%% Step 3: Build ABCD Matrix
ABCD = stuffABCD(a, g, b, c, form);

[n_rows, ~] = size(ABCD);
n_states = n_rows - 1;

A_mat = ABCD(1:n_states, 1:n_states);
B_mat = ABCD(1:n_states, n_states+1:end);
C_mat = ABCD(n_states+1, 1:n_states);
D_mat = ABCD(n_states+1, n_states+1:end);

%% Step 4: Simulation Parameters
fs = 2e6;                   % Sampling frequency (Hz)
fB = fs / (2*OSR);          % Signal bandwidth (Hz)
N = 8192;                   % FFT points

% Flash ADC Parameters (5-bit, 32 levels)
n_bits = 5;
n_levels = 2^n_bits;        % 32 levels
V_fs = 1.0;                 % Full scale (+/- 1V)

fprintf('  fs = %.1f MHz\n', fs/1e6);
fprintf('  OSR = %d (BW = %.1f kHz)\n', OSR, fB/1000);
fprintf('  Flash ADC: %d-bit (%d levels)\n', n_bits, n_levels);
fprintf('  DAC Mismatch: 0%% (ideal) vs 0.5%% (static) vs 0.5%% (with DWA)\n\n');

%% Step 5: Input Signal
f_bin = round(sqrt(1/7) * N / (2*OSR));
A_in = 0.5;

fprintf('  Input: %.2f V sine wave @ %.2f kHz\n\n', A_in, f_bin/N * fs / 1000);

u = A_in * sin(2*pi*f_bin*(0:N-1)/N);
n = 0:N-1;

%% Step 6: Run Simulations
fprintf('[Step 6] Run Simulations\n');
fprintf('----------------------------------------\n');

% Simulation 1: Ideal DAC (no mismatch)
fprintf('  Running: Ideal DAC (no mismatch)...\n');
[v_ideal_dac, v_ideal_adc, therm_ideal, y_ideal] = run_dsm_with_flash_adc(...
    u, A_mat, B_mat, C_mat, D_mat, n_bits, V_fs);

% Simulation 2: DAC with 0.5% mismatch (static)
mismatch_pct = 0.005;  % 0.5%
seed = 42;
fprintf('  Running: DAC with %.1f%% mismatch (static)...\n', mismatch_pct*100);
[v_mismatch_dac, v_mismatch_adc, therm_mismatch, y_mismatch, bit_weights] = run_dsm_with_flash_adc_mismatch(...
    u, A_mat, B_mat, C_mat, D_mat, n_bits, V_fs, mismatch_pct, seed);

% Simulation 3: DAC with 0.5% mismatch + DWA
fprintf('  Running: DAC with %.1f%% mismatch (with DWA)...\n', mismatch_pct*100);
[v_dwa_dac, v_dwa_adc, therm_dwa, y_dwa] = run_dsm_with_flash_adc_dwa(...
    u, A_mat, B_mat, C_mat, D_mat, n_bits, V_fs, mismatch_pct, seed);

fprintf('  Simulations complete.\n\n');

%% Step 7: Calculate SNR for both
fprintf('[Step 7] Calculate SNR\n');
fprintf('----------------------------------------\n');

% Use ADC binary output for SNR calculation (as in real DSM)
[SNR_ideal, ENOB_ideal] = calculate_snr(v_ideal_adc, u, N, OSR);
[SNR_mismatch, ENOB_mismatch] = calculate_snr(v_mismatch_adc, u, N, OSR);
[SNR_dwa, ENOB_dwa] = calculate_snr(v_dwa_adc, u, N, OSR);

% Calculate SNDR for mismatch cases (includes distortion)
SNDR_mismatch = calculate_sndr(v_mismatch_adc, u, N, OSR, f_bin);
SNDR_dwa = calculate_sndr(v_dwa_adc, u, N, OSR, f_bin);

fprintf('  Ideal DAC:\n');
fprintf('    SNR:  %.2f dB\n', SNR_ideal);
fprintf('    ENOB: %.2f bits\n\n', ENOB_ideal);

fprintf('  0.5%% Mismatch DAC (static):\n');
fprintf('    SNR:   %.2f dB\n', SNR_mismatch);
fprintf('    SNDR:  %.2f dB\n', SNDR_mismatch);
fprintf('    ENOB:  %.2f bits\n', ENOB_mismatch);
fprintf('    SNR degradation: %.2f dB\n\n', SNR_ideal - SNR_mismatch);

fprintf('  0.5%% Mismatch DAC (with DWA):\n');
fprintf('    SNR:   %.2f dB\n', SNR_dwa);
fprintf('    SNDR:  %.2f dB\n', SNDR_dwa);
fprintf('    ENOB:  %.2f bits\n', ENOB_dwa);
fprintf('    SNR improvement vs static: %.2f dB\n\n', SNR_dwa - SNR_mismatch);

%% Step 8: Analyze Mismatch Effects
fprintf('[Step 8] DAC Mismatch Analysis\n');
fprintf('----------------------------------------\n');

nominal_weight = (2 * V_fs) / (n_levels - 1);
fprintf('  Nominal unit cell weight: %.6fV\n', nominal_weight);
fprintf('  Standard deviation: %.4f%%\n\n', mismatch_pct * 100);

% Calculate actual INL/DNL from bit_weights
inldnl = calculate_inl_dnl_thermometer(bit_weights, V_fs, n_bits);
fprintf('  INL (max): %.3f LSB\n', max(abs(inldnl.inl)));
fprintf('  DNL (max): %.3f LSB\n\n', max(abs(inldnl.dnl)));

%% Step 9: Generate Plots
fprintf('[Step 9] Generate Comparison Plots\n');
fprintf('----------------------------------------\n');

generate_comparison_plots(u, v_ideal_adc, v_mismatch_adc, v_dwa_adc, v_ideal_dac, v_mismatch_dac, ...
    therm_ideal, therm_mismatch, n_bits, fs, N, OSR, SNR_ideal, SNR_mismatch, SNR_dwa, bit_weights, V_fs);

%% Summary
fprintf('\n============================================================\n');
fprintf('  SUMMARY\n');
fprintf('============================================================\n');
fprintf('  Ideal DAC:\n');
fprintf('    SNR:  %.2f dB, ENOB: %.2f bits\n', SNR_ideal, ENOB_ideal);
fprintf('  0.5%% Mismatch DAC (static):\n');
fprintf('    SNR:  %.2f dB, ENOB: %.2f bits\n', SNR_mismatch, ENOB_mismatch);
fprintf('    SNDR: %.2f dB (includes distortion)\n', SNDR_mismatch);
fprintf('    Degradation: %.2f dB\n', SNR_ideal - SNR_mismatch);
fprintf('    INL: %.3f LSB, DNL: %.3f LSB\n', max(abs(inldnl.inl)), max(abs(inldnl.dnl)));
fprintf('  0.5%% Mismatch DAC (with DWA):\n');
fprintf('    SNR:  %.2f dB, ENOB: %.2f bits\n', SNR_dwa, ENOB_dwa);
fprintf('    SNDR: %.2f dB (includes distortion)\n', SNDR_dwa);
fprintf('    Improvement vs static: %.2f dB\n', SNR_dwa - SNR_mismatch);
fprintf('============================================================\n');

%% ========================================================================
%% Helper Functions
%% ========================================================================

function [v_dac, v_adc, thermometer, y_all] = run_dsm_with_flash_adc(u, A_mat, B_mat, C_mat, D_mat, n_bits, V_fs)
    %% Run DSM simulation with Flash ADC + Ideal Thermometer DAC
    %% Returns:
    %%   v_dac - DAC output (same as ADC in ideal case)
    %%   v_adc - ADC binary output (the actual DSM output)
    N = length(u);
    n_states = size(A_mat, 1);
    n_levels = 2^n_bits;
    
    x = zeros(n_states, 1);
    v_dac = zeros(1, N);
    v_adc = zeros(1, N);
    y_all = zeros(1, N);
    thermometer = zeros(N, n_levels - 1);
    
    for i = 1:N
        if i == 1
            v_prev = 0;
        else
            v_prev = v_dac(i-1);
        end
        
        % Compute modulator output (before quantizer)
        y = C_mat * x + D_mat(1)*u(i) + D_mat(2)*v_prev;
        y_all(i) = y;
        
        % Flash ADC quantization
        [therm, binary_adc, ~] = flash_adc_quantizer(y, V_fs, n_bits);
        thermometer(i, :) = therm;
        v_adc(i) = binary_adc;  % This is the DSM output
        
        % Ideal DAC: output equals ADC output
        v_dac(i) = v_adc(i);  % For ideal DAC
        
        % Update state
        x = A_mat*x + B_mat(:,1)*u(i) + B_mat(:,2)*v_dac(i);
    end
end

function [v_dac, v_adc, thermometer, y_all, bit_weights] = run_dsm_with_flash_adc_mismatch(u, A_mat, B_mat, C_mat, D_mat, n_bits, V_fs, mismatch_pct, seed)
    %% Run DSM simulation with Flash ADC + Mismatched Thermometer DAC
    %% Returns:
    %%   v_dac - DAC output (goes to feedback loop)
    %%   v_adc - ADC binary output (the actual DSM output)
    N = length(u);
    n_states = size(A_mat, 1);
    n_levels = 2^n_bits;
    n_comparators = n_levels - 1;
    
    % Generate mismatched bit weights (once)
    if nargin >= 9 && ~isempty(seed)
        rng(seed);
    end
    
    nominal_weight = (2 * V_fs) / n_comparators;
    bit_weights = nominal_weight * (1 + mismatch_pct * randn(1, n_comparators));
    
    x = zeros(n_states, 1);
    v_dac = zeros(1, N);  % DAC output for feedback
    v_adc = zeros(1, N);  % ADC binary output (measured output)
    y_all = zeros(1, N);
    thermometer = zeros(N, n_comparators);
    
    for i = 1:N
        if i == 1
            v_prev = 0;
        else
            v_prev = v_dac(i-1);  % Feedback from DAC
        end
        
        % Compute modulator output (before quantizer)
        y = C_mat * x + D_mat(1)*u(i) + D_mat(2)*v_prev;
        y_all(i) = y;
        
        % Flash ADC quantization
        [therm, binary_adc, ~] = flash_adc_quantizer(y, V_fs, n_bits);
        thermometer(i, :) = therm;
        v_adc(i) = binary_adc;  % This is the actual DSM output
        
        % Mismatched DAC reconstruction (for feedback)
        contribution = 0;
        for j = 1:n_comparators
            if therm(j) == 1
                contribution = contribution + bit_weights(j);
            else
                contribution = contribution - bit_weights(j);
            end
        end
        v_dac(i) = contribution * 0.5;  % This goes to feedback
        
        % Update state using DAC output
        x = A_mat*x + B_mat(:,1)*u(i) + B_mat(:,2)*v_dac(i);
    end
end

function [v_dac, v_adc, thermometer, y_all, bit_weights] = run_dsm_with_flash_adc_dwa(u, A_mat, B_mat, C_mat, D_mat, n_bits, V_fs, mismatch_pct, seed)
    %% Run DSM simulation with Flash ADC + DWA Mismatched Thermometer DAC
    %% Returns:
    %%   v_dac - DAC output (goes to feedback loop)
    %%   v_adc - ADC binary output (the actual DSM output)
    N = length(u);
    n_states = size(A_mat, 1);
    n_levels = 2^n_bits;
    n_comparators = n_levels - 1;
    
    % Generate mismatched bit weights (once)
    if nargin >= 9 && ~isempty(seed)
        rng(seed);
    end
    
    nominal_weight = (2 * V_fs) / n_comparators;
    bit_weights = nominal_weight * (1 + mismatch_pct * randn(1, n_comparators));
    
    x = zeros(n_states, 1);
    v_dac = zeros(1, N);  % DAC output for feedback
    v_adc = zeros(1, N);  % ADC binary output (measured output)
    y_all = zeros(1, N);
    thermometer = zeros(N, n_comparators);
    
    % DWA state: starting index for cell selection (1-indexed)
    dwa_start_idx = 1;
    
    for i = 1:N
        if i == 1
            v_prev = 0;
        else
            v_prev = v_dac(i-1);  % Feedback from DAC
        end
        
        % Compute modulator output (before quantizer)
        y = C_mat * x + D_mat(1)*u(i) + D_mat(2)*v_prev;
        y_all(i) = y;
        
        % Flash ADC quantization
        [therm, binary_adc, ~] = flash_adc_quantizer(y, V_fs, n_bits);
        thermometer(i, :) = therm;
        v_adc(i) = binary_adc;  % This is the actual DSM output
        
        % DWA: Count number of ones (cells to activate)
        k = sum(therm);
        
        % DWA cell selection with rotation
        if k == 0
            % All zeros - no cells activated
            contribution = -sum(bit_weights);
            v_dac(i) = contribution * 0.5;
            % Start position stays same (or could advance by 1)
        elseif k == n_comparators
            % All ones - all cells activated
            contribution = sum(bit_weights);
            v_dac(i) = contribution * 0.5;
            % Advance by k (mod n_comparators)
            dwa_start_idx = mod(dwa_start_idx - 1 + k, n_comparators) + 1;
        else
            % Normal case: select k consecutive cells starting from dwa_start_idx
            contribution = 0;
            for j = 0:k-1
                idx = mod(dwa_start_idx - 1 + j, n_comparators) + 1;  % Wrap around
                if therm(idx)
                    contribution = contribution + bit_weights(idx);
                else
                    contribution = contribution - bit_weights(idx);
                end
            end
            v_dac(i) = contribution * 0.5;
            
            % Update start position for next sample
            dwa_start_idx = mod(dwa_start_idx - 1 + k, n_comparators) + 1;
        end
        
        % Update state using DAC output
        x = A_mat*x + B_mat(:,1)*u(i) + B_mat(:,2)*v_dac(i);
    end
end

function [thermometer, binary, thresholds] = flash_adc_quantizer(v_in, v_fs, n_bits)
    %% Flash ADC Quantizer with Thermometer Code Output
    n_levels = 2^n_bits;
    n_comparators = n_levels - 1;
    
    v_min = -v_fs;
    v_max = v_fs;
    lsb = (v_max - v_min) / (n_levels - 1);
    thresholds = linspace(v_min + lsb/2, v_max - lsb/2, n_comparators);
    
    % Clip input
    v_in_clipped = max(v_min, min(v_max, v_in));
    
    % Generate thermometer code
    thermometer = zeros(1, n_comparators);
    for i = 1:n_comparators
        thermometer(i) = double(v_in_clipped > thresholds(i));
    end
    
    % Calculate binary output
    level_idx = sum(thermometer);
    output_levels = linspace(v_min, v_max, n_levels);
    binary = output_levels(level_idx + 1);
end

function [SNR, ENOB] = calculate_snr(v, u, N, OSR)
    %% Calculate SNR (Signal-to-Noise Ratio)
    %% SNR = Signal Power / Noise Power
    %% Signal = fundamental only
    %% Noise = everything in-band except DC, fundamental, and harmonics
    w = 0.5 * (1 - cos(2*pi*(0:N-1)/N));
    V_out = fft(v .* w) / (N/4);
    V_out_mag = abs(V_out);
    
    % Find fundamental signal
    [~, sig_idx] = max(V_out_mag(2:N/2));
    sig_bin = sig_idx + 1;
    fB_bins = ceil(N / (2*OSR));
    
    % Signal bins (3-bin around fundamental)
    sig_bins = sig_bin-1:sig_bin+1;
    sig_bins = sig_bins(sig_bins >= 2 & sig_bins <= fB_bins);
    signal_power = sum(V_out_mag(sig_bins).^2);
    
    % Identify harmonic bins to exclude from noise (up to 7th harmonic)
    harmonic_bins = [];
    for h = 2:7  % 2nd through 7th harmonics
        harmonic_bin = h * sig_bin;
        if harmonic_bin <= fB_bins
            % Exclude 3 bins around each harmonic
            harmonic_bins = [harmonic_bins, harmonic_bin-1:harmonic_bin+1];
        end
    end
    harmonic_bins = unique(harmonic_bins);
    harmonic_bins = harmonic_bins(harmonic_bins >= 2 & harmonic_bins <= fB_bins);
    
    % Noise bins = in-band excluding DC (bin 1), signal, and harmonics
    exclude_bins = unique([sig_bins, harmonic_bins]);
    noise_bins = setdiff(2:fB_bins, exclude_bins);
    noise_power = sum(V_out_mag(noise_bins).^2);
    
    if noise_power > 0
        SNR = 10*log10(signal_power / noise_power);
        ENOB = (SNR - 1.76) / 6.02;
    else
        SNR = inf;
        ENOB = inf;
    end
end

function SNDR = calculate_sndr(v, u, N, OSR, f_bin)
    %% Calculate SNDR (Signal-to-Noise-and-Distortion Ratio)
    %% SNDR = Signal Power / (Noise + Distortion) Power
    %% Signal = fundamental only (NOT harmonics!)
    %% Noise+Distortion = everything in-band except DC and fundamental
    %%                    (this includes harmonics in the denominator)
    w = 0.5 * (1 - cos(2*pi*(0:N-1)/N));
    V_out = fft(v .* w) / (N/4);
    V_out_mag = abs(V_out);
    
    fB_bins = ceil(N / (2*OSR));
    
    % Signal = fundamental only (3 bins)
    % Note: Use actual signal bin from spectrum, not input f_bin
    [~, sig_idx] = max(V_out_mag(2:N/2));
    sig_bin = sig_idx + 1;
    sig_bins = sig_bin-1:sig_bin+1;
    sig_bins = sig_bins(sig_bins >= 2 & sig_bins <= fB_bins);
    signal_power = sum(V_out_mag(sig_bins).^2);
    
    % Noise + Distortion = everything in-band except DC and fundamental
    noise_distortion_bins = setdiff(2:fB_bins, sig_bins);
    noise_distortion_power = sum(V_out_mag(noise_distortion_bins).^2);
    
    if noise_distortion_power > 0
        SNDR = 10*log10(signal_power / noise_distortion_power);
    else
        SNDR = inf;
    end
end

function inldnl = calculate_inl_dnl_thermometer(bit_weights, v_fs, n_bits)
    %% Calculate INL and DNL for thermometer DAC
    n_levels = 2^n_bits;
    n_comparators = n_levels - 1;
    nominal_weight = (2 * v_fs) / n_comparators;
    
    % Calculate actual output levels
    actual_levels = zeros(1, n_levels);
    for k = 0:n_levels-1
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
    
    % LSB
    lsb = (2 * v_fs) / (n_levels - 1);
    
    % INL and DNL
    inl = (actual_levels - ideal_levels) / lsb;
    dnl = diff(actual_levels) / lsb - 1;
    dnl = [0, dnl];
    
    inldnl.inl = inl;
    inldnl.dnl = dnl;
    inldnl.code = 0:n_levels-1;
end

function generate_comparison_plots(u, v_ideal_adc, v_mismatch_adc, v_dwa_adc, v_ideal_dac, v_mismatch_dac, ...
    therm_ideal, therm_mismatch, n_bits, fs, N, OSR, SNR_ideal, SNR_mismatch, SNR_dwa, bit_weights, V_fs)
    %% Generate comparison plots
    n_comparators = size(therm_ideal, 2);
    n_plot = 500;
    t = (0:n_plot-1)/fs*1e6;  % microseconds
    
    figure('Name', 'DSM with Flash ADC - Ideal vs Mismatched DAC', 'Position', [50 50 1600 1200]);
    
    % 1. Time domain comparison (ADC outputs)
    subplot(3, 3, 1);
    plot(t, u(1:n_plot), 'b-', 'LineWidth', 1); hold on;
    stairs(t, v_ideal_adc(1:n_plot), 'g-', 'LineWidth', 1.5);
    stairs(t, v_mismatch_adc(1:n_plot), 'r-', 'LineWidth', 1.5);
    hold off;
    xlabel('Time (\mus)'); ylabel('Amplitude (V)');
    title('ADC Output: Time Domain');
    legend('Input', 'Ideal DAC', 'Mismatched DAC', 'Location', 'best');
    grid on;
    
    % 2. Spectra comparison (using ADC output)
    subplot(3, 3, 2);
    w = 0.5 * (1 - cos(2*pi*(0:N-1)/N));
    V_ideal = fft(v_ideal_adc .* w) / (N/4);
    V_mismatch = fft(v_mismatch_adc .* w) / (N/4);
    V_dwa = fft(v_dwa_adc .* w) / (N/4);
    freqs = (0:N/2)/N*fs/1000;
    
    semilogx(freqs, 20*log10(abs(V_ideal(1:N/2+1)) + eps), 'g-', 'LineWidth', 1.5); hold on;
    semilogx(freqs, 20*log10(abs(V_mismatch(1:N/2+1)) + eps), 'r-', 'LineWidth', 1);
    semilogx(freqs, 20*log10(abs(V_dwa(1:N/2+1)) + eps), 'b-', 'LineWidth', 1);
    fB_kHz = fs / (2*OSR) / 1000;
    plot([fB_kHz fB_kHz], [-140 10], 'm--', 'LineWidth', 2);
    hold off;
    xlabel('Frequency (kHz)'); ylabel('Magnitude (dBFS/NBW)');
    title(sprintf('ADC Output Spectrum\nG: SNR=%.1f dB, R: SNR=%.1f dB, B: SNR=%.1f dB', SNR_ideal, SNR_mismatch, SNR_dwa));
    legend('Ideal', 'Static', 'DWA', 'Signal BW', 'Location', 'northwest');
    grid on;
    axis([100 fs/2000 -140 10]);
    
    % 3. Mismatch-induced error (ADC outputs) - Static vs DWA
    subplot(3, 3, 3);
    mismatch_error = v_ideal_adc - v_mismatch_adc;
    dwa_error = v_ideal_adc - v_dwa_adc;
    plot(t, mismatch_error(1:n_plot), 'r-', 'LineWidth', 1); hold on;
    plot(t, dwa_error(1:n_plot), 'b-', 'LineWidth', 1);
    hold off;
    xlabel('Time (\mus)'); ylabel('Error (V)');
    title('Error: Ideal - Output');
    legend('Static', 'DWA', 'Location', 'best');
    grid on; yline(0, 'k--', 'Alpha', 0.3);
    
    % 4. Ideal DAC thermometer code
    subplot(3, 3, 4);
    imagesc(t, 0:n_comparators-1, therm_ideal(1:n_plot, :)');
    colormap(gca, [1 0.3 0.3; 0.2 1 0.2]);
    set(gca, 'YDir', 'normal');
    xlabel('Time (\mus)'); ylabel('Comparator #');
    title('Ideal DAC: Thermometer Code');
    
    % 5. Mismatched DAC thermometer code
    subplot(3, 3, 5);
    imagesc(t, 0:n_comparators-1, therm_mismatch(1:n_plot, :)');
    colormap(gca, [1 0.3 0.3; 0.2 1 0.2]);
    set(gca, 'YDir', 'normal');
    xlabel('Time (\mus)'); ylabel('Comparator #');
    title('Static DAC: Thermometer Code');
    
    % 6. Histogram comparison (ADC outputs)
    subplot(3, 3, 6);
    histogram(v_ideal_adc, 32, 'FaceColor', 'g', 'FaceAlpha', 0.5); hold on;
    histogram(v_mismatch_adc, 32, 'FaceColor', 'r', 'FaceAlpha', 0.5);
    histogram(v_dwa_adc, 32, 'FaceColor', 'b', 'FaceAlpha', 0.3);
    hold off;
    xlabel('Output Level (V)'); ylabel('Count');
    title('ADC Output Histogram');
    legend('Ideal', 'Static', 'DWA', 'Location', 'best');
    
    % 7. SNR comparison bar chart
    subplot(3, 3, 7);
    bar_data = [SNR_ideal, SNR_mismatch, SNR_dwa];
    bar_colors = {'g', 'r', 'b'};
    h = bar(bar_data);
    h.FaceColor = 'flat';
    h.CData(1,:) = [0.2 0.8 0.2];
    h.CData(2,:) = [0.8 0.2 0.2];
    h.CData(3,:) = [0.2 0.4 0.8];
    set(gca, 'XTickLabel', {'Ideal', 'Static', 'DWA'});
    ylabel('SNR (dB)');
    title(sprintf('SNR Comparison\nDWA Improvement: +%.1f dB', SNR_dwa - SNR_mismatch));
    grid on;
    
    % 8. INL plot
    subplot(3, 3, 8);
    inldnl = calculate_inl_dnl_thermometer(bit_weights, V_fs, n_bits);
    stem(inldnl.code, inldnl.inl, 'filled', 'MarkerFaceColor', 'b');
    xlabel('DAC Code'); ylabel('INL (LSB)');
    title(sprintf('Integral Nonlinearity\nMax: %.3f LSB', max(abs(inldnl.inl))));
    grid on; yline(0, 'k--', 'Alpha', 0.3);
    
    % 9. DNL plot
    subplot(3, 3, 9);
    stem(inldnl.code, inldnl.dnl, 'filled', 'MarkerFaceColor', 'r');
    xlabel('DAC Code'); ylabel('DNL (LSB)');
    title(sprintf('Differential Nonlinearity\nMax: %.3f LSB', max(abs(inldnl.dnl))));
    grid on; yline(0, 'k--', 'Alpha', 0.3);
    yline(1, 'g--', 'Alpha', 0.3);
    yline(-1, 'g--', 'Alpha', 0.3);
    
    % Save
    saveas(gcf, 'dsm_4th_order_dac_mismatch_comparison.png');
    fprintf('  Saved: dsm_4th_order_dac_mismatch_comparison.png\n');
end
