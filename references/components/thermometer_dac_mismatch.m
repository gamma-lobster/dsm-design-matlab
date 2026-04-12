%% THERMOMETER_DAC_MISMATCH
% DAC with unit cell mismatch modeling
%
% Models real-world DAC non-ideality where each unit cell has slightly
% different weight due to manufacturing variations.
%
% Mismatch model:
%   Actual weight = nominal_weight + error
%   error ~ N(0, (percentage * nominal_weight)^2)
%
% Usage:
%   [v_out, bit_weights] = thermometer_dac_mismatch(thermometer, v_fs, n_bits, mismatch_pct, seed)
%
% Inputs:
%   thermometer   - Matrix of thermometer code (samples x (2^n_bits - 1))
%   v_fs          - Full-scale voltage, default = 1.0
%   n_bits        - Number of bits, default = 3
%   mismatch_pct  - Mismatch percentage (e.g., 0.01 for 1%), default = 0
%   seed          - Random seed for reproducibility (optional)
%
% Outputs:
%   v_out       - Analog output voltage with mismatch
%   bit_weights - Actual bit weights including mismatch (for analysis)
%
% Example:
%   [therm, ~, ~] = flash_adc_quantizer(v_in, 1.0, 3);
%   [v_dac, weights] = thermometer_dac_mismatch(therm, 1.0, 3, 0.02);
%   % 2% mismatch applied to each unit cell

function [v_out, bit_weights] = thermometer_dac_mismatch(thermometer, v_fs, n_bits, mismatch_pct, seed)
    % Default parameters
    if nargin < 2
        v_fs = 1.0;
    end
    if nargin < 3
        n_bits = 3;
    end
    if nargin < 4
        mismatch_pct = 0;  % No mismatch by default
    end
    
    % Calculate derived values
    n_levels = 2^n_bits;
    n_comparators = n_levels - 1;
    
    % Nominal weight per bit (signed weight: +weight for 1, -weight for 0)
    % Each bit contributes +W if 1, -W if 0
    % Net contribution = (2*bit - 1) * W
    nominal_weight = (2 * v_fs) / n_comparators;
    
    % Set random seed if provided (for reproducibility)
    if nargin >= 5 && ~isempty(seed)
        rng(seed);
    else
        % Use default random state
    end
    
    % Generate mismatched bit weights
    % Each comparator cell has independent mismatch
    bit_weights = zeros(1, n_comparators);
    for i = 1:n_comparators
        % Nominal signed weight (same for all bits in thermometer DAC)
        nominal = nominal_weight;
        
        % Gaussian error: mean=0, std=mismatch_pct * nominal
        if mismatch_pct > 0
            error = randn() * mismatch_pct * nominal;
        else
            error = 0;
        end
        
        % Actual weight includes mismatch
        bit_weights(i) = nominal + error;
    end
    
    % Compute DAC output for each sample
    num_samples = size(thermometer, 1);
    v_out = zeros(1, num_samples);
    
    for i = 1:num_samples
        % For each bit: if thermometer=1, use +weight; if 0, use -weight
        % But the DAC output is the sum of actual contributions
        % Let's define: contribution = bit_weights(j) * (2*thermometer(i,j) - 1)
        contribution = 0;
        for j = 1:n_comparators
            if thermometer(i, j) == 1
                contribution = contribution + bit_weights(j);
            else
                contribution = contribution - bit_weights(j);
            end
        end
        
        % Scale to output voltage range
        % The raw contribution ranges from -sum(weights) to +sum(weights)
        % We want output to be approximately -v_fs to +v_fs
        % So we need to scale by (v_fs / sum_of_weights)
        v_out(i) = contribution * 0.5;
    end
    
    % Ensure row vector output
    v_out = v_out(:)';
end
