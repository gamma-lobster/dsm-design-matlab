%% FLASH_ADC_QUANTIZER
% Flash ADC Quantizer with Thermometer Code Output
%
% A flash ADC uses a bank of comparators, each with a slightly different
% reference voltage. The output is thermometer code: all 0s below the
% input level, all 1s above.
%
% For N-bit resolution (2^N levels), we have (2^N - 1) comparators.
% Example: 3-bit flash ADC has 7 comparators, producing 7-bit thermometer code.
%
% Usage:
%   [thermometer, binary, thresholds] = flash_adc_quantizer(v_in, v_fs, n_bits)
%
% Inputs:
%   v_in    - Input voltage(s) [scalar or vector]
%   v_fs    - Full-scale voltage (symmetric +/-v_fs), default = 1.0
%   n_bits  - Number of bits, default = 3
%
% Outputs:
%   thermometer - Matrix of thermometer code (samples x (2^n_bits - 1))
%   binary      - Quantized output values
%   thresholds  - Comparator threshold voltages
%
% Example:
%   v_in = linspace(-1, 1, 100);
%   [therm, bin, thresh] = flash_adc_quantizer(v_in, 1.0, 3);

function [thermometer, binary, thresholds] = flash_adc_quantizer(v_in, v_fs, n_bits)
    % Default parameters
    if nargin < 2
        v_fs = 1.0;
    end
    if nargin < 3
        n_bits = 3;
    end
    
    % Calculate number of levels and comparators
    n_levels = 2^n_bits;
    n_comparators = n_levels - 1;
    
    % Voltage range
    v_min = -v_fs;
    v_max = v_fs;
    
    % LSB size
    lsb = (v_max - v_min) / (n_levels - 1);
    
    % Comparator thresholds (evenly spaced)
    % For +/-1V with 8 levels: thresholds at -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75
    thresholds = linspace(v_min + lsb/2, v_max - lsb/2, n_comparators);
    
    % Ensure input is column vector
    v_in = v_in(:);
    n_samples = length(v_in);
    
    % Clip input to valid range
    v_in_clipped = max(v_min, min(v_max, v_in));
    
    % Generate thermometer code
    % Each comparator outputs 1 if v_in > threshold
    thermometer = zeros(n_samples, n_comparators);
    for i = 1:n_comparators
        thermometer(:, i) = double(v_in_clipped > thresholds(i));
    end
    
    % Calculate binary output from thermometer code
    % Count number of 1s in thermometer code for each sample
    level_idx = sum(thermometer, 2);
    
    % Map to output voltage levels
    output_levels = linspace(v_min, v_max, n_levels);
    binary = output_levels(level_idx + 1);  % +1 because MATLAB is 1-indexed
    
    % Ensure binary is a row vector (1 x n_samples) to match v_in format
    binary = binary(:)';
    
    % Transpose thresholds to row vector for cleaner display
    thresholds = thresholds(:)';
end
