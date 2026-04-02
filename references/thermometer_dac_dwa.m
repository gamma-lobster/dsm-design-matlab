%% THERMOMETER_DAC_DWA
% DAC with Dynamic Weighted Averaging (DWA) for mismatch shaping
%
% DWA rotates the starting position of thermometer code selection
% to average out unit cell mismatch errors over time.
%
% Usage:
%   [v_out, start_idx] = thermometer_dac_dwa(thermometer, v_fs, n_bits, bit_weights, start_idx)
%
% Inputs:
%   thermometer - Matrix of thermometer code (samples x (2^n_bits - 1))
%   v_fs        - Full-scale voltage, default = 1.0
%   n_bits      - Number of bits, default = 3
%   bit_weights - Optional: custom weights for each cell (for mismatch simulation)
%   start_idx   - Optional: starting index for first sample (1-indexed)
%
% Output:
%   v_out       - Analog output voltage (with DWA rotation)
%   start_idx   - Final starting index (for continuing DWA)
%
% DWA Algorithm:
%   For each sample with thermometer code (k ones):
%     1. Select k consecutive cells starting from current position
%     2. Wrap around if needed (circular buffer)
%     3. Next start position = position after last selected cell
%
% Example with 7 cells (3-bit), code = 3 (three 1s):
%   Sample 1: cells 1,2,3 activated → next start = 4
%   Sample 2: cells 4,5,6 activated → next start = 7
%   Sample 3: cells 7,1,2 activated (wrap) → next start = 3
%   Sample 4: cells 3,4,5 activated → next start = 6
%   ... and so on

function [v_out, start_idx_out] = thermometer_dac_dwa(thermometer, v_fs, n_bits, bit_weights, start_idx_in)
    % Default parameters
    if nargin < 2
        v_fs = 1.0;
    end
    if nargin < 3
        n_bits = 3;
    end
    
    n_levels = 2^n_bits;
    n_comparators = n_levels - 1;
    
    % Default: equal weights (no mismatch)
    if nargin < 4 || isempty(bit_weights)
        nominal_weight = (2 * v_fs) / n_comparators;
        bit_weights = nominal_weight * ones(1, n_comparators);
    end
    
    % Ensure bit_weights is correct size
    if length(bit_weights) ~= n_comparators
        error('bit_weights must have length = 2^n_bits - 1 = %d', n_comparators);
    end
    
    % Starting index (1-indexed, MATLAB style)
    if nargin < 5 || isempty(start_idx_in)
        start_idx = 1;
    else
        start_idx = start_idx_in;
    end
    
    num_samples = size(thermometer, 1);
    v_out = zeros(1, num_samples);
    
    for i = 1:num_samples
        % Count number of ones in thermometer code
        k = sum(thermometer(i, :));  % k = number of cells to activate
        
        if k == 0
            % Special case: all zeros - don't activate any cells
            % Next start stays the same (or could increment by 1, depending on implementation)
            % We'll keep it same for now
            contribution = -sum(bit_weights);  % All cells contribute -weight
            v_out(i) = contribution * 0.5;
            continue;
        end
        
        if k == n_comparators
            % Special case: all ones - activate all cells
            contribution = sum(bit_weights);  % All cells contribute +weight
            v_out(i) = contribution * 0.5;
            % Next start advances by k (mod n_comparators)
            start_idx = mod(start_idx - 1 + k, n_comparators) + 1;
            continue;
        end
        
        % DWA: Select k consecutive cells starting from start_idx
        selected_cells = zeros(1, n_comparators);
        
        % Create indices with wrap-around
        for j = 0:k-1
            idx = mod(start_idx - 1 + j, n_comparators) + 1;  % 1-indexed
            selected_cells(idx) = 1;
        end
        
        % Calculate DAC output: sum of selected cells positive, unselected negative
        contribution = 0;
        for j = 1:n_comparators
            if selected_cells(j)
                contribution = contribution + bit_weights(j);  % Selected: +weight
            else
                contribution = contribution - bit_weights(j);  % Not selected: -weight
            end
        end
        
        v_out(i) = contribution * 0.5;
        
        % Update start index for next sample
        start_idx = mod(start_idx - 1 + k, n_comparators) + 1;
    end
    
    % Return final start index for continuation
    start_idx_out = start_idx;
end