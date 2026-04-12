%% test_build_3rd_order_ct_simulink_model.m
% Build and smoke-test a Simulink CTDSM model for the 3rd-order example.

clearvars; close all; clc;

script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(fileparts(script_dir)), 'common'));
refs_root = setup_references_paths();

order = 3;
OSR = 32;
H_inf = 4.0;
opt = 1;
f0 = 0;
ct_form = 'FF';
tdac = [0 1];
fs = 10e6;
n_bits = 4;
V_fs = 1.0;
stop_time = 2e-4;

ntf = synthesizeNTF(order, OSR, opt, H_inf, f0);

[model_name, model_path, info] = build_ct_dsm_simulink_model( ...
    'model_name', 'dsm_3rd_order_ct_10mhz_topology_model', ...
    'output_dir', fullfile(refs_root, 'results', 'simulink', 'ct'), ...
    'ntf', ntf, ...
    'ct_form', ct_form, ...
    'tdac', tdac, ...
    'fs', fs, ...
    'OSR', OSR, ...
    'n_bits', n_bits, ...
    'V_fs', V_fs, ...
    'input_amplitude', 0.5, ...
    'stop_time', stop_time, ...
    'ct_steps_per_sample', 32, ...
    'open_model', false);

fprintf('Built model: %s\n', model_name);
fprintf('Saved to: %s\n', model_path);
fprintf('CT form: %s\n', info.ct_form);
fprintf('tdac2:\n');
disp(info.tdac2);

load_system(model_path);
simOut = sim(model_name, 'StopTime', num2str(stop_time, 16));

u = extract_signal(simOut, 'u_dsm');
y_ct = extract_signal(simOut, 'y_ct_dsm');
y = extract_signal(simOut, 'y_dsm');
v = extract_signal(simOut, 'v_dsm');

fprintf('Captured samples:\n');
fprintf('  u:    %d\n', numel(u));
fprintf('  y_ct: %d\n', numel(y_ct));
fprintf('  y:    %d\n', numel(y));
fprintf('  v:    %d\n', numel(v));

close_system(model_name, 0);

function sig = extract_signal(simOut, name)
if isprop(simOut, name)
    raw = simOut.(name);
else
    raw = simOut.get(name);
end

if isa(raw, 'timeseries')
    sig = raw.Data(:);
elseif isstruct(raw) && isfield(raw, 'signals')
    sig = raw.signals.values(:);
elseif isnumeric(raw)
    if isvector(raw)
        sig = raw(:);
    elseif size(raw, 2) >= 2
        sig = raw(:, end);
    else
        sig = raw(:);
    end
else
    error('Unsupported signal format for "%s".', name);
end
end
