%% THERMOMETER_DAC
% DAC that converts thermometer code to analog voltage
%
% The DAC takes thermometer code as input and produces an analog output.
% Each thermometer bit has weight = full_scale / (2^n_bits - 1)
% Bit = 1 contributes positive weight
% Bit = 0 contributes negative weight
%
% Usage:
%   v_out = thermometer_dac(thermometer, v_fs, n_bits)
%
% Inputs:
%   thermometer - Matrix of thermometer code (samples x (2^n_bits - 1))
%   v_fs        - Full-scale voltage, default = 1.0
%   n_bits      - Number of bits, default = 3
%
% Output:
%   v_out - Analog output voltage
%
% Example:
%   v_in = linspace(-1, 1, 100);
%   [therm, ~, ~] = flash_adc_quantizer(v_in, 1.0, 3);
%   v_dac = thermometer_dac(therm, 1.0, 3);
%   % v_dac should equal the binary quantized output

function v_out = thermometer_dac(thermometer, v_fs, n_bits)
    % Default parameters
    if nargin < 2
        v_fs = 1.0;
    end
    if nargin < 3
        n_bits = 3;
    end
    
    % Calculate derived values
    n_levels = 2^n_bits;
    n_comparators = n_levels - 1;
    
    % Weight per bit
    bit_weight = (2 * v_fs) / n_comparators;  % full_scale / (2^n_bits - 1)
    
    % For each sample: DAC output = sum over all bits
    % Bit = 1: contribute +bit_weight
    % Bit = 0: contribute -bit_weight
    %
    % Mathematically: v_out = sum( (2*thermometer - 1) * bit_weight )
    % Which equals: v_out = (2*sum(thermometer) - n_comparators) * bit_weight
    %
    % Simplifying:
    %   v_out = (2*level - n_comparators) * bit_weight
    %   v_out = (2*level - (n_levels-1)) * (2*v_fs)/(n_levels-1)
    %
    % This maps:
    %   level = 0 (all zeros)  → v_out = -(n_levels-1)*bit_weight = -2*v_fs
    %   level = n_comparators (all ones) → v_out = +n_comparators*bit_weight = +2*v_fs
    %
    % Wait, that gives range -2*v_fs to +2*v_fs...
    % 
    % Alternative interpretation: use thermometer to create weighted sum
    % that matches the ADC output levels
    
    num_samples = size(thermometer, 1);
    v_out = zeros(num_samples, 1);
    
    for i = 1:num_samples
        % Count ones and zeros
        num_ones = sum(thermometer(i, :));
        num_zeros = n_comparators - num_ones;
        
        % DAC output: positive contribution from 1s, negative from 0s
        v_out(i) = (num_ones * bit_weight) - (num_zeros * bit_weight);
        
        % This simplifies to:
        % v_out = (num_ones - num_zeros) * bit_weight
        % v_out = (num_ones - (n_comparators - num_ones)) * bit_weight
        % v_out = (2*num_ones - n_comparators) * bit_weight
    end
    
    % Scale to match expected ADC output range
    % The above gives -2*v_fs to +2*v_fs, but we want -v_fs to +v_fs
    % So scale by 0.5
    v_out = v_out * 0.5;
    
    % Ensure output is a row vector (1 x n_samples) to match v_in format
    v_out = v_out(:)';
end
