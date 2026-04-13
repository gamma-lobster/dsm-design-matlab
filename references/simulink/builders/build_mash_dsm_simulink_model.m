%% BUILD_MASH_DSM_SIMULINK_MODEL
% Build a generic discrete-time cascaded-residue MASH DSM Simulink model.
%
% Supported entry points:
%   build_mash_dsm_simulink_model(spec_struct)
%   build_mash_dsm_simulink_model('name', value, ...)
%
% The builder currently supports the common MASH pattern where:
%   - stage 1 processes the external input u
%   - each later stage processes the residue from the previous stage
%   - each stage is a single-quantizer DT loop with inputs [signal; v_k]
%   - the final digital output is a sum of filtered stage quantizer outputs
%
% By default the digital cancellation filters are:
%   H_1(z) = 1
%   H_k(z) = prod_{i=1}^{k-1} NTF_i(z),  k >= 2
%
% This default matches textbook MASH error-cancellation for cascaded stages
% when each later stage shapes the previous stage residue.

function [model_name, model_path, info] = build_mash_dsm_simulink_model(varargin)
    spec = parse_inputs(varargin{:});

    addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'common'));
    setup_references_paths();

    if isempty(spec.fs)
        error('build_mash_dsm_simulink_model:MissingFs', ...
            'Sampling frequency fs must be provided.');
    end

    if isempty(spec.stages)
        error('build_mash_dsm_simulink_model:MissingStages', ...
            'Provide a non-empty stages struct array.');
    end

    if isempty(spec.input_frequency)
        if ~isempty(spec.OSR)
            N_tone = 8192;
            f_bin = round(sqrt(1/7) * N_tone / (2*spec.OSR));
            spec.input_frequency = f_bin / N_tone * spec.fs;
        else
            spec.input_frequency = spec.fs / 128;
        end
    end

    if isempty(spec.stop_time)
        spec.stop_time = max(200/spec.fs, 1e-4);
    end

    stages = resolve_stages(spec);
    [combiner_filters, combiner_desc] = resolve_combiner_filters(spec, stages);

    if numel(combiner_filters) ~= numel(stages)
        error('build_mash_dsm_simulink_model:CombinerLengthMismatch', ...
            'Need exactly one digital combiner filter per stage.');
    end

    if bdIsLoaded(spec.model_name)
        close_system(spec.model_name, 0);
    end

    model_path = fullfile(spec.output_dir, [spec.model_name '.slx']);
    if exist(model_path, 'file')
        delete(model_path);
    end

    new_system(spec.model_name);
    set_param(spec.model_name, 'Solver', 'FixedStepDiscrete');
    set_param(spec.model_name, 'FixedStep', num2str(1/spec.fs, 16));
    set_param(spec.model_name, 'StopTime', num2str(spec.stop_time, 16));

    prepare_workspace(spec, stages);
    build_model(spec, stages, combiner_filters);

    Simulink.BlockDiagram.arrangeSystem(spec.model_name);
    save_system(spec.model_name, model_path);

    if spec.open_model
        open_system(spec.model_name);
    else
        close_system(spec.model_name);
    end

    model_name = spec.model_name;
    info = struct();
    info.stages = stages;
    info.combiner_filters = combiner_filters;
    info.combiner_description = combiner_desc;
    info.fs = spec.fs;
    info.OSR = spec.OSR;
    info.input_amplitude = spec.input_amplitude;
    info.input_frequency = spec.input_frequency;
    info.model_path = model_path;
end

function spec = parse_inputs(varargin)
    defaults = struct( ...
        'model_name', 'dsm_mash_topology_model', ...
        'output_dir', pwd, ...
        'stages', [], ...
        'digital_combiner_filters', [], ...
        'digital_combiner_mode', 'cumulative-ntf', ...
        'fs', [], ...
        'OSR', [], ...
        'V_fs', 1.0, ...
        'input_amplitude', 0.5, ...
        'input_frequency', [], ...
        'stop_time', [], ...
        'open_model', false);

    if nargin == 1 && isstruct(varargin{1})
        spec = defaults;
        user = varargin{1};
        names = fieldnames(user);
        for i = 1:numel(names)
            spec.(names{i}) = user.(names{i});
        end
        return;
    end

    p = inputParser;
    fields = fieldnames(defaults);
    for i = 1:numel(fields)
        addParameter(p, fields{i}, defaults.(fields{i}));
    end
    parse(p, varargin{:});
    spec = p.Results;
end

function stages = resolve_stages(spec)
    raw_stages = spec.stages;
    resolved = cell(1, numel(raw_stages));

    for k = 1:numel(raw_stages)
        stage = raw_stages(k);
        if ~isfield(stage, 'name') || isempty(stage.name)
            stage.name = sprintf('Stage%d', k);
        end
        if ~isfield(stage, 'form') || isempty(stage.form)
            stage.form = 'CIFF';
        end
        if ~isfield(stage, 'n_bits') || isempty(stage.n_bits)
            if isfield(stage, 'n_levels') && ~isempty(stage.n_levels)
                stage.n_bits = log2(stage.n_levels);
            else
                stage.n_bits = 1;
            end
        end
        if ~isfield(stage, 'n_levels') || isempty(stage.n_levels)
            stage.n_levels = 2^stage.n_bits;
        end
        if ~isfield(stage, 'source') || isempty(stage.source)
            if k == 1
                stage.source = 'external';
            else
                stage.source = sprintf('residue%d', k-1);
            end
        end

        resolved{k} = resolve_stage_data(stage, spec);
    end

    stages = [resolved{:}];
end

function stage = resolve_stage_data(stage, spec)
    if ~isfield(stage, 'ABCD')
        stage.ABCD = [];
    end
    if ~isfield(stage, 'ntf')
        stage.ntf = [];
    end

    if isempty(stage.ABCD)
        if isempty(stage.ntf)
            required = {'order', 'opt', 'H_inf'};
            for i = 1:numel(required)
                if ~isfield(stage, required{i}) || isempty(stage.(required{i}))
                    error('build_mash_dsm_simulink_model:MissingStageField', ...
                        'Stage "%s" is missing required field "%s".', stage.name, required{i});
                end
            end
            stage.ntf = synthesizeNTF(stage.order, spec.OSR, stage.opt, stage.H_inf, 0);
        end
        [stage.a, stage.g, stage.b, stage.c] = realizeNTF(stage.ntf, stage.form);
        stage.ABCD = stuffABCD(stage.a, stage.g, stage.b, stage.c, stage.form);
    else
        if isempty(stage.ntf)
            error('build_mash_dsm_simulink_model:MissingStageNtf', ...
                'Stage "%s" needs its matching NTF when ABCD is provided directly.', stage.name);
        end
        try
            [stage.a, stage.g, stage.b, stage.c] = mapABCD(stage.ABCD, stage.form);
        catch
            stage.a = [];
            stage.g = [];
            stage.b = [];
            stage.c = [];
        end
    end

    [stage.A, stage.B, stage.C, stage.D] = partition_abcd(stage.ABCD);
    if size(stage.B, 2) ~= 2 || size(stage.D, 2) ~= 2
        error('build_mash_dsm_simulink_model:UnsupportedStageShape', ...
            'Stage "%s" must have two inputs [signal; v_k].', stage.name);
    end
end

function [filters, description] = resolve_combiner_filters(spec, stages)
    n_stages = numel(stages);

    if ~isempty(spec.digital_combiner_filters)
        filters = spec.digital_combiner_filters;
        description = 'user-specified';
        return;
    end

    mode = lower(spec.digital_combiner_mode);
    switch mode
        case 'cumulative-ntf'
            filters = repmat({struct('num', 1, 'den', 1, 'label', '')}, 1, n_stages);
            running_num = 1;
            running_den = 1;
            filters{1}.label = '1';
            for k = 2:n_stages
                [num_km1, den_km1] = tfdata(stages(k-1).ntf, 'v');
                running_num = conv(running_num, num_km1);
                running_den = conv(running_den, den_km1);
                filters{k} = struct( ...
                    'num', running_num, ...
                    'den', running_den, ...
                    'label', sprintf('prod(NTF_1:NTF_%d)', k-1));
            end
            description = 'cumulative previous-stage NTF products';
        otherwise
            error('build_mash_dsm_simulink_model:UnsupportedCombinerMode', ...
                'Unsupported digital_combiner_mode "%s".', spec.digital_combiner_mode);
    end
end

function prepare_workspace(spec, stages)
    mw = get_param(spec.model_name, 'ModelWorkspace');
    assignin(mw, 'Ts_mash', 1/spec.fs);
    assignin(mw, 'n_stages_mash', numel(stages));

    for k = 1:numel(stages)
        assignin(mw, sprintf('A%d_mash', k), stages(k).A);
        assignin(mw, sprintf('B%d_mash', k), stages(k).B);
        assignin(mw, sprintf('C%d_mash', k), stages(k).C);
        assignin(mw, sprintf('D%d_mash', k), stages(k).D);
    end
end

function build_model(spec, stages, filters)
    n_stages = numel(stages);
    x_step = 430;
    y_main = 110;
    y_filter = 360;

    add_block('simulink/Sources/Sine Wave', [spec.model_name '/Input'], ...
        'Position', [30 y_main-15 85 y_main+15], ...
        'Amplitude', num2str(spec.input_amplitude, 16), ...
        'Frequency', num2str(2*pi*spec.input_frequency, 16), ...
        'SampleTime', num2str(1/spec.fs, 16), ...
        'Bias', '0');

    add_block('simulink/Sinks/To Workspace', [spec.model_name '/u_out'], ...
        'Position', [35 y_main+40 105 y_main+60], ...
        'VariableName', 'u_mash', ...
        'SaveFormat', 'Array');
    add_line(spec.model_name, 'Input/1', 'u_out/1');
    final_sum_inputs = '';

    for k = 1:n_stages
        x0 = 120 + (k-1) * x_step;

        build_stage_loop(spec, stages(k), k, x0, y_main);

        if k < n_stages
            add_stage_residue(spec, k, x0, y_main);
        end

        filter_name = sprintf('CombinerFilter%d', k);
        add_block('simulink/Discrete/Discrete Filter', [spec.model_name '/' filter_name], ...
            'Position', [x0+255 y_filter x0+390 y_filter+30], ...
            'Numerator', mat2str(filters{k}.num, 16), ...
            'Denominator', mat2str(filters{k}.den, 16), ...
            'InitialStates', '0');
        add_line(spec.model_name, sprintf('Quantizer%d/1', k), sprintf('%s/1', filter_name), 'autorouting', 'on');

        final_sum_inputs = [final_sum_inputs '+']; %#ok<AGROW>
    end

    sum_x = 120 + n_stages * x_step + 20;
    add_block('simulink/Math Operations/Sum', [spec.model_name '/MashOutput'], ...
        'Position', [sum_x y_filter+5 sum_x+30 y_filter+35], ...
        'Inputs', final_sum_inputs);

    for k = 1:n_stages
        add_line(spec.model_name, sprintf('CombinerFilter%d/1', k), sprintf('MashOutput/%d', k), 'autorouting', 'on');
    end

    add_block('simulink/Sinks/To Workspace', [spec.model_name '/vmash_out'], ...
        'Position', [sum_x+60 y_filter-45 sum_x+140 y_filter-25], ...
        'VariableName', 'v_mash', ...
        'SaveFormat', 'Array');
    add_line(spec.model_name, 'MashOutput/1', 'vmash_out/1', 'autorouting', 'on');

    add_block('simulink/Sinks/Scope', [spec.model_name '/Signals'], ...
        'Position', [sum_x+90 y_filter+10 sum_x+145 y_filter+10+25*(n_stages+2)], ...
        'NumInputPorts', num2str(n_stages + 2));
    add_line(spec.model_name, 'Input/1', 'Signals/1', 'autorouting', 'on');
    for k = 1:n_stages
        add_line(spec.model_name, sprintf('Quantizer%d/1', k), sprintf('Signals/%d', k+1), 'autorouting', 'on');
    end
    add_line(spec.model_name, 'MashOutput/1', sprintf('Signals/%d', n_stages+2), 'autorouting', 'on');
end

function build_stage_loop(spec, stage, k, x0, y_mid)
n_states = size(stage.A, 1);
signal_tag = sprintf('stage%d_signal', k);
v_tag = sprintf('stage%d_v', k);
vprev_tag = sprintf('stage%d_vprev', k);
state_x = x0 + 40;
state_dx = 95;
y_state = y_mid + 150;
y_state_sum = y_state + 15;
y_fb = y_state + 110;

add_stage_signal_router(spec, stage, k, x0, y_mid, signal_tag);

add_block('simulink/Discrete/Unit Delay', [spec.model_name sprintf('/v%d_prev', k)], ...
    'Position', [x0+285 y_fb x0+325 y_fb+30], ...
    'SampleTime', 'Ts_mash', ...
    'InitialCondition', '0');
add_block('simulink/Signal Routing/Goto', [spec.model_name sprintf('/v%d_prev_tag', k)], ...
    'Position', [x0+340 y_fb+5 x0+390 y_fb+25], ...
    'GotoTag', vprev_tag, ...
    'TagVisibility', 'global');
add_line(spec.model_name, sprintf('v%d_prev/1', k), sprintf('v%d_prev_tag/1', k), 'autorouting', 'on');

for i = 1:n_states
    sum_terms = nnz(stage.A(i, :) ~= 0) + (stage.B(i, 1) ~= 0) + (stage.B(i, 2) ~= 0);
    sum_terms = max(sum_terms, 1);

    add_block('simulink/Math Operations/Sum', [spec.model_name sprintf('/Sum_x%d_%d_next', k, i)], ...
        'Position', [state_x-55+(i-1)*state_dx y_state_sum state_x-25+(i-1)*state_dx y_state_sum+20], ...
        'Inputs', repmat('+', 1, sum_terms));
    add_block('simulink/Discrete/Unit Delay', [spec.model_name sprintf('/x%d_%d', k, i)], ...
        'Position', [state_x+(i-1)*state_dx y_state state_x+40+(i-1)*state_dx y_state+30], ...
        'SampleTime', 'Ts_mash', ...
        'InitialCondition', '0');
    add_line(spec.model_name, sprintf('Sum_x%d_%d_next/1', k, i), sprintf('x%d_%d/1', k, i), 'autorouting', 'on');
end

for i = 1:n_states
    input_port = 1;
    for j = 1:n_states
        if stage.A(i, j) == 0
            continue;
        end
        src = sprintf('x%d_%d/1', k, j);
        if stage.A(i, j) == 1
            add_line(spec.model_name, src, sprintf('Sum_x%d_%d_next/%d', k, i, input_port), 'autorouting', 'on');
        else
            gain_name = sprintf('Gain_A%d_%d_%d', k, i, j);
            add_block('simulink/Math Operations/Gain', [spec.model_name '/' gain_name], ...
                'Position', [state_x+(j-1)*state_dx+50 y_mid+65 state_x+(j-1)*state_dx+100 y_mid+85], ...
                'Gain', num2str(stage.A(i, j), 16));
            add_line(spec.model_name, src, sprintf('%s/1', gain_name), 'autorouting', 'on');
            add_line(spec.model_name, sprintf('%s/1', gain_name), sprintf('Sum_x%d_%d_next/%d', k, i, input_port), 'autorouting', 'on');
        end
        input_port = input_port + 1;
    end

    if stage.B(i, 1) ~= 0
        add_block('simulink/Signal Routing/From', [spec.model_name sprintf('/from_u%d_%d', k, i)], ...
            'Position', [x0-35 y_state-20+30*(i-1) x0+15 y_state+30*(i-1)], ...
            'GotoTag', signal_tag);
        if stage.B(i, 1) == 1
            add_line(spec.model_name, sprintf('from_u%d_%d/1', k, i), sprintf('Sum_x%d_%d_next/%d', k, i, input_port), 'autorouting', 'on');
        else
            gain_name = sprintf('Gain_Bu%d_%d', k, i);
            add_block('simulink/Math Operations/Gain', [spec.model_name '/' gain_name], ...
                'Position', [x0+20 y_state-20+30*(i-1) x0+70 y_state+30*(i-1)], ...
                'Gain', num2str(stage.B(i, 1), 16));
            add_line(spec.model_name, sprintf('from_u%d_%d/1', k, i), sprintf('%s/1', gain_name), 'autorouting', 'on');
            add_line(spec.model_name, sprintf('%s/1', gain_name), sprintf('Sum_x%d_%d_next/%d', k, i, input_port), 'autorouting', 'on');
        end
        input_port = input_port + 1;
    end

    if stage.B(i, 2) ~= 0
        add_block('simulink/Signal Routing/From', [spec.model_name sprintf('/from_v%d_%d', k, i)], ...
            'Position', [x0-35 y_fb+20+30*(i-1) x0+15 y_fb+40+30*(i-1)], ...
            'GotoTag', v_tag);
        if stage.B(i, 2) == 1
            add_line(spec.model_name, sprintf('from_v%d_%d/1', k, i), sprintf('Sum_x%d_%d_next/%d', k, i, input_port), 'autorouting', 'on');
        else
            gain_name = sprintf('Gain_Bv%d_%d', k, i);
            add_block('simulink/Math Operations/Gain', [spec.model_name '/' gain_name], ...
                'Position', [x0+20 y_fb+20+30*(i-1) x0+70 y_fb+40+30*(i-1)], ...
                'Gain', num2str(stage.B(i, 2), 16));
            add_line(spec.model_name, sprintf('from_v%d_%d/1', k, i), sprintf('%s/1', gain_name), 'autorouting', 'on');
            add_line(spec.model_name, sprintf('%s/1', gain_name), sprintf('Sum_x%d_%d_next/%d', k, i, input_port), 'autorouting', 'on');
        end
    end
end

y_terms = nnz(stage.C ~= 0) + (stage.D(1) ~= 0) + (stage.D(2) ~= 0);
add_block('simulink/Math Operations/Sum', [spec.model_name sprintf('/Sum_y%d', k)], ...
    'Position', [x0+220 y_mid-10 x0+250 y_mid+20], ...
    'Inputs', repmat('+', 1, max(y_terms, 1)));

input_port = 1;
for j = 1:n_states
    if stage.C(j) == 0
        continue;
    end
    if stage.C(j) == 1
        add_line(spec.model_name, sprintf('x%d_%d/1', k, j), sprintf('Sum_y%d/%d', k, input_port), 'autorouting', 'on');
    else
        gain_name = sprintf('Gain_C%d_%d', k, j);
        add_block('simulink/Math Operations/Gain', [spec.model_name '/' gain_name], ...
            'Position', [x0+105+(j-1)*55 y_mid+55 x0+155+(j-1)*55 y_mid+75], ...
            'Gain', num2str(stage.C(j), 16));
        add_line(spec.model_name, sprintf('x%d_%d/1', k, j), sprintf('%s/1', gain_name), 'autorouting', 'on');
        add_line(spec.model_name, sprintf('%s/1', gain_name), sprintf('Sum_y%d/%d', k, input_port), 'autorouting', 'on');
    end
    input_port = input_port + 1;
end

if stage.D(1) ~= 0
    add_block('simulink/Signal Routing/From', [spec.model_name sprintf('/from_uy%d', k)], ...
        'Position', [x0+95 y_mid-55 x0+145 y_mid-35], ...
        'GotoTag', signal_tag);
    if stage.D(1) == 1
        add_line(spec.model_name, sprintf('from_uy%d/1', k), sprintf('Sum_y%d/%d', k, input_port), 'autorouting', 'on');
    else
        gain_name = sprintf('Gain_Du%d', k);
        add_block('simulink/Math Operations/Gain', [spec.model_name '/' gain_name], ...
            'Position', [x0+155 y_mid-55 x0+205 y_mid-35], ...
            'Gain', num2str(stage.D(1), 16));
        add_line(spec.model_name, sprintf('from_uy%d/1', k), sprintf('%s/1', gain_name), 'autorouting', 'on');
        add_line(spec.model_name, sprintf('%s/1', gain_name), sprintf('Sum_y%d/%d', k, input_port), 'autorouting', 'on');
    end
    input_port = input_port + 1;
end

if stage.D(2) ~= 0
    add_block('simulink/Signal Routing/From', [spec.model_name sprintf('/from_vprevy%d', k)], ...
        'Position', [x0+95 y_mid+105 x0+145 y_mid+125], ...
        'GotoTag', vprev_tag);
    if stage.D(2) == 1
        add_line(spec.model_name, sprintf('from_vprevy%d/1', k), sprintf('Sum_y%d/%d', k, input_port), 'autorouting', 'on');
    else
        gain_name = sprintf('Gain_Dv%d', k);
        add_block('simulink/Math Operations/Gain', [spec.model_name '/' gain_name], ...
            'Position', [x0+155 y_mid+105 x0+205 y_mid+125], ...
            'Gain', num2str(stage.D(2), 16));
        add_line(spec.model_name, sprintf('from_vprevy%d/1', k), sprintf('%s/1', gain_name), 'autorouting', 'on');
        add_line(spec.model_name, sprintf('%s/1', gain_name), sprintf('Sum_y%d/%d', k, input_port), 'autorouting', 'on');
    end
end

add_block('simulink/User-Defined Functions/MATLAB Fcn', [spec.model_name sprintf('/Quantizer%d', k)], ...
    'Position', [x0+300 y_mid-15 x0+400 y_mid+20], ...
    'MATLABFcn', sprintf('dsm_quantizer(u, %.16g, %d)', spec.V_fs, stage.n_levels));
add_line(spec.model_name, sprintf('Sum_y%d/1', k), sprintf('Quantizer%d/1', k), 'autorouting', 'on');

add_block('simulink/Signal Routing/Goto', [spec.model_name sprintf('/v%d_tag', k)], ...
    'Position', [x0+415 y_mid-10 x0+470 y_mid+10], ...
    'GotoTag', v_tag, ...
    'TagVisibility', 'global');
add_line(spec.model_name, sprintf('Quantizer%d/1', k), sprintf('v%d_prev/1', k), 'autorouting', 'on');
add_line(spec.model_name, sprintf('Quantizer%d/1', k), sprintf('v%d_tag/1', k), 'autorouting', 'on');

add_block('simulink/Sinks/To Workspace', [spec.model_name sprintf('/y%d_out', k)], ...
    'Position', [x0+215 y_mid-90 x0+285 y_mid-70], ...
    'VariableName', sprintf('y%d_mash', k), ...
    'SaveFormat', 'Array');
add_block('simulink/Sinks/To Workspace', [spec.model_name sprintf('/v%d_out', k)], ...
    'Position', [x0+405 y_mid-90 x0+475 y_mid-70], ...
    'VariableName', sprintf('v%d_mash', k), ...
    'SaveFormat', 'Array');
add_line(spec.model_name, sprintf('Sum_y%d/1', k), sprintf('y%d_out/1', k), 'autorouting', 'on');
add_line(spec.model_name, sprintf('Quantizer%d/1', k), sprintf('v%d_out/1', k), 'autorouting', 'on');
end

function add_stage_signal_router(spec, stage, k, x0, y_mid, signal_tag)
add_block('simulink/Signal Routing/Goto', [spec.model_name sprintf('/stage%d_signal_tag', k)], ...
    'Position', [x0+5 y_mid-10 x0+60 y_mid+10], ...
    'GotoTag', signal_tag, ...
    'TagVisibility', 'global');

switch lower(stage.source)
    case 'external'
        add_line(spec.model_name, 'Input/1', sprintf('stage%d_signal_tag/1', k), 'autorouting', 'on');
    otherwise
        expected_source = sprintf('residue%d', k-1);
        if ~strcmpi(stage.source, expected_source)
            error('build_mash_dsm_simulink_model:UnsupportedSource', ...
                'Stage "%s" uses source "%s". Only cascaded residue chaining is supported now.', ...
                stage.name, stage.source);
        end
        add_line(spec.model_name, sprintf('Residue%d/1', k-1), sprintf('stage%d_signal_tag/1', k), 'autorouting', 'on');
end
end

function add_stage_residue(spec, k, x0, y_mid)
add_block('simulink/Math Operations/Sum', [spec.model_name sprintf('/Residue%d', k)], ...
    'Position', [x0+360 y_mid+80 x0+390 y_mid+110], ...
    'Inputs', '+-');
add_block('simulink/Sinks/To Workspace', [spec.model_name sprintf('/r%d_out', k)], ...
    'Position', [x0+420 y_mid+120 x0+490 y_mid+140], ...
    'VariableName', sprintf('r%d_mash', k), ...
    'SaveFormat', 'Array');
add_line(spec.model_name, sprintf('Sum_y%d/1', k), sprintf('Residue%d/1', k), 'autorouting', 'on');
add_line(spec.model_name, sprintf('Quantizer%d/1', k), sprintf('Residue%d/2', k), 'autorouting', 'on');
add_line(spec.model_name, sprintf('Residue%d/1', k), sprintf('r%d_out/1', k), 'autorouting', 'on');
end

function [A, B, C, D] = partition_abcd(ABCD)
    n_states = size(ABCD, 1) - 1;
    A = ABCD(1:n_states, 1:n_states);
    B = ABCD(1:n_states, n_states+1:end);
    C = ABCD(n_states+1, 1:n_states);
    D = ABCD(n_states+1, n_states+1:end);
end
