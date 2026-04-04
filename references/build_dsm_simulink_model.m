%% BUILD_DSM_SIMULINK_MODEL
% Build a generic single-quantizer DSM Simulink model from a toolbox design.
%
% Supported entry points:
%   build_dsm_simulink_model(spec_struct)
%   build_dsm_simulink_model('name', value, ...)
%
% The builder accepts either:
%   1. ABCD + form
%   2. a, g, b, c + form
%   3. ntf + form
%   4. order, OSR, H_inf, opt, f0, form
%
% The generated model supports two implementations:
%   state-space : compact Discrete State-Space loop filter
%   expanded    : visible state-update blocks, gains, and summing junctions
%
% The expanded mode is useful when you want something closer to the
% appendix topology figures while still being generated from the actual ABCD
% realization produced by the toolbox.

function [model_name, model_path, info] = build_dsm_simulink_model(varargin)
    spec = parse_inputs(varargin{:});

    if isempty(spec.fs)
        error('build_dsm_simulink_model:MissingFs', ...
            'Sampling frequency fs must be provided.');
    end

    if isempty(spec.form)
        error('build_dsm_simulink_model:MissingForm', ...
            'A realization form such as CIFF or CIFB must be provided.');
    end

    addpath(fullfile(fileparts(mfilename('fullpath')), 'dstoolbox'));

    if isempty(spec.n_levels)
        if ~isempty(spec.n_bits)
            spec.n_levels = 2^spec.n_bits;
        else
            spec.n_levels = 2;
        end
    end

    if isempty(spec.n_bits)
        spec.n_bits = log2(spec.n_levels);
    end

    [ABCD, coeffs, ntf] = resolve_model_data(spec);
    n_states = size(ABCD, 1) - 1;
    A = ABCD(1:n_states, 1:n_states);
    B = ABCD(1:n_states, n_states+1:end);
    C = ABCD(n_states+1, 1:n_states);
    D = ABCD(n_states+1, n_states+1:end);

    if size(B, 2) ~= 2 || size(D, 2) ~= 2
        error('build_dsm_simulink_model:UnsupportedABCD', ...
            'This builder currently supports only single-quantizer ABCD matrices with two inputs [u; v].');
    end

    if isempty(spec.input_frequency)
        if ~isempty(spec.OSR)
            N_tone = 8192;
            f_bin = round(sqrt(1/7) * N_tone / (2*spec.OSR));
            spec.input_frequency = f_bin/N_tone * spec.fs;
        else
            spec.input_frequency = spec.fs / 128;
        end
    end

    if isempty(spec.stop_time)
        spec.stop_time = max(200/spec.fs, 1e-4);
    end

    if bdIsLoaded(spec.model_name)
        close_system(spec.model_name, 0);
    end

    if exist(fullfile(spec.output_dir, [spec.model_name '.slx']), 'file')
        delete(fullfile(spec.output_dir, [spec.model_name '.slx']));
    end

    new_system(spec.model_name);

    set_param(spec.model_name, 'Solver', 'FixedStepDiscrete');
    set_param(spec.model_name, 'FixedStep', num2str(1/spec.fs, 16));
    set_param(spec.model_name, 'StopTime', num2str(spec.stop_time, 16));

    mw = get_param(spec.model_name, 'ModelWorkspace');
    assignin(mw, 'A_dsm', A);
    assignin(mw, 'B_dsm', B);
    assignin(mw, 'C_dsm', C);
    assignin(mw, 'D_dsm', D);
    assignin(mw, 'Ts_dsm', 1/spec.fs);
    assignin(mw, 'x0_dsm', zeros(n_states, 1));

    if strcmpi(spec.implementation, 'expanded')
        build_expanded_loop(spec, A, B, C, D, n_states, coeffs);
    else
        build_state_space_loop(spec);
    end

    if ~strcmpi(spec.implementation, 'expanded')
        Simulink.BlockDiagram.arrangeSystem(spec.model_name);
    end

    model_name = spec.model_name;
    model_path = fullfile(spec.output_dir, [spec.model_name '.slx']);
    save_system(spec.model_name, model_path);

    if spec.open_model
        open_system(spec.model_name);
    else
        close_system(spec.model_name);
    end

    info = struct();
    info.ABCD = ABCD;
    info.coeffs = coeffs;
    info.ntf = ntf;
    info.form = spec.form;
    info.fs = spec.fs;
    info.n_levels = spec.n_levels;
    info.input_amplitude = spec.input_amplitude;
    info.input_frequency = spec.input_frequency;
    info.implementation = spec.implementation;
    info.model_path = model_path;
end

function spec = parse_inputs(varargin)
    defaults = struct( ...
        'model_name', 'dsm_simulink_model', ...
        'output_dir', pwd, ...
        'form', '', ...
        'ABCD', [], ...
        'a', [], ...
        'g', [], ...
        'b', [], ...
        'c', [], ...
        'ntf', [], ...
        'order', [], ...
        'OSR', [], ...
        'H_inf', [], ...
        'opt', 1, ...
        'f0', 0, ...
        'fs', [], ...
        'n_bits', [], ...
        'n_levels', [], ...
        'V_fs', 1.0, ...
        'input_amplitude', 0.5, ...
        'input_frequency', [], ...
        'stop_time', [], ...
        'implementation', 'state-space', ...
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

function [ABCD, coeffs, ntf] = resolve_model_data(spec)
    coeffs = struct('a', [], 'g', [], 'b', [], 'c', []);
    ntf = [];

    if ~isempty(spec.ABCD)
        ABCD = spec.ABCD;
        if any(~isempty([spec.a, spec.g, spec.b, spec.c]))
            coeffs = struct('a', spec.a, 'g', spec.g, 'b', spec.b, 'c', spec.c);
        else
            try
                [coeffs.a, coeffs.g, coeffs.b, coeffs.c] = mapABCD(ABCD, spec.form);
            catch
                coeffs = struct('a', [], 'g', [], 'b', [], 'c', []);
            end
        end
        return;
    end

    if ~isempty(spec.a) && ~isempty(spec.b) && ~isempty(spec.c)
        coeffs = struct('a', spec.a, 'g', spec.g, 'b', spec.b, 'c', spec.c);
        ABCD = stuffABCD(coeffs.a, coeffs.g, coeffs.b, coeffs.c, spec.form);
        return;
    end

    if ~isempty(spec.ntf)
        ntf = spec.ntf;
        [coeffs.a, coeffs.g, coeffs.b, coeffs.c] = realizeNTF(ntf, spec.form);
        ABCD = stuffABCD(coeffs.a, coeffs.g, coeffs.b, coeffs.c, spec.form);
        return;
    end

    required = {'order', 'OSR', 'H_inf'};
    for i = 1:numel(required)
        if isempty(spec.(required{i}))
            error('build_dsm_simulink_model:MissingDesignData', ...
                'Need ABCD, coefficients, ntf, or full design inputs (order, OSR, H_inf, ...).');
        end
    end

    ntf = synthesizeNTF(spec.order, spec.OSR, spec.opt, spec.H_inf, spec.f0);
    [coeffs.a, coeffs.g, coeffs.b, coeffs.c] = realizeNTF(ntf, spec.form);
    ABCD = stuffABCD(coeffs.a, coeffs.g, coeffs.b, coeffs.c, spec.form);
end

function build_state_space_loop(spec)
    add_block('simulink/Sources/Sine Wave', [spec.model_name '/Input'], ...
        'Position', [40, 140, 90, 170], ...
        'Amplitude', num2str(spec.input_amplitude, 16), ...
        'Frequency', num2str(2*pi*spec.input_frequency, 16), ...
        'SampleTime', num2str(1/spec.fs, 16), ...
        'Bias', '0');

    add_block('simulink/Signal Routing/Mux', [spec.model_name '/Mux_u_v'], ...
        'Position', [145, 125, 150, 205], ...
        'Inputs', '2');

    add_block('simulink/Discrete/Discrete State-Space', ...
        [spec.model_name '/LoopFilter'], ...
        'Position', [220, 118, 350, 212], ...
        'A', 'A_dsm', ...
        'B', 'B_dsm', ...
        'C', 'C_dsm', ...
        'D', 'D_dsm', ...
        'X0', 'x0_dsm', ...
        'SampleTime', 'Ts_dsm');

    add_quantizer_and_outputs(spec, [405, 140, 520, 190], 'LoopFilter/1', 'Input/1', 'LoopFilter/1');

    add_block('simulink/Signal Routing/Goto', [spec.model_name '/v_fb'], ...
        'Position', [565, 105, 620, 125], ...
        'GotoTag', 'v_fb', ...
        'TagVisibility', 'global');

    add_block('simulink/Signal Routing/From', [spec.model_name '/v_fb_from'], ...
        'Position', [85, 175, 135, 195], ...
        'GotoTag', 'v_fb');

    add_line(spec.model_name, 'Input/1', 'Mux_u_v/1');
    add_line(spec.model_name, 'v_fb_from/1', 'Mux_u_v/2');
    add_line(spec.model_name, 'Mux_u_v/1', 'LoopFilter/1');
end

function build_expanded_loop(spec, A, B, C, D, n_states, coeffs)
    add_block('simulink/Sources/Sine Wave', [spec.model_name '/Input'], ...
        'Position', [40, 220, 90, 250], ...
        'Amplitude', num2str(spec.input_amplitude, 16), ...
        'Frequency', num2str(2*pi*spec.input_frequency, 16), ...
        'SampleTime', num2str(1/spec.fs, 16), ...
        'Bias', '0');

    add_block('simulink/Signal Routing/Goto', [spec.model_name '/u_sig'], ...
        'Position', [115, 220, 165, 240], ...
        'GotoTag', 'u_sig', ...
        'TagVisibility', 'global');
    add_line(spec.model_name, 'Input/1', 'u_sig/1');

    add_block('simulink/Signal Routing/Goto', [spec.model_name '/v_fb'], ...
        'Position', [1090, 210, 1140, 230], ...
        'GotoTag', 'v_fb', ...
        'TagVisibility', 'global');

    state_x0 = 300;
    state_dx = 170;
    state_y = 240;
    sum_y = 295;
    gain_y0 = 150;
    gain_dy = 45;

    for i = 1:n_states
        sum_terms = 0;
        for j = 1:n_states
            if A(i, j) ~= 0
                sum_terms = sum_terms + 1;
            end
        end
        if B(i, 1) ~= 0
            sum_terms = sum_terms + 1;
        end
        if B(i, 2) ~= 0
            sum_terms = sum_terms + 1;
        end
        if sum_terms == 0
            sum_terms = 1;
        end

        sum_name = sprintf('%s/Sum_x%d_next', spec.model_name, i);
        delay_name = sprintf('%s/x%d', spec.model_name, i);
        add_block('simulink/Math Operations/Sum', sum_name, ...
            'Position', [state_x0 - 55 + (i-1)*state_dx, state_y + 5, state_x0 - 25 + (i-1)*state_dx, state_y + 25], ...
            'Inputs', repmat('+', 1, sum_terms));
        add_block('simulink/Discrete/Unit Delay', delay_name, ...
            'Position', [state_x0 + (i-1)*state_dx, state_y, state_x0 + 40 + (i-1)*state_dx, state_y + 30], ...
            'SampleTime', 'Ts_dsm', ...
            'InitialCondition', sprintf('x0_dsm(%d)', i));
        add_line(spec.model_name, sprintf('Sum_x%d_next/1', i), sprintf('x%d/1', i));
    end

    for i = 1:n_states
        input_port = 1;
        for j = 1:n_states
            if A(i, j) == 0
                continue;
            end
            src = sprintf('x%d/1', j);
            if A(i, j) == 1
                add_line(spec.model_name, src, sprintf('Sum_x%d_next/%d', i, input_port), 'autorouting', 'on');
            else
                gain_name = sprintf('%s/Gain_A_%d_%d', spec.model_name, i, j);
                if i == 2 && j == 3
                    gy = state_y - 35;
                    gx = state_x0 + 2*state_dx + 55;
                else
                    gy = gain_y0 + (j-1)*gain_dy;
                    gx = state_x0 + (j-1)*state_dx + 60;
                end
                add_block('simulink/Math Operations/Gain', gain_name, ...
                    'Position', [gx, gy, gx + 50, gy + 20], ...
                    'Gain', num2str(A(i, j), 16));
                add_line(spec.model_name, src, sprintf('Gain_A_%d_%d/1', i, j), 'autorouting', 'on');
                add_line(spec.model_name, sprintf('Gain_A_%d_%d/1', i, j), sprintf('Sum_x%d_next/%d', i, input_port), 'autorouting', 'on');
            end
            input_port = input_port + 1;
        end

        if B(i, 1) ~= 0
            from_name = sprintf('%s/from_u_%d', spec.model_name, i);
            add_block('simulink/Signal Routing/From', from_name, ...
                'Position', [state_x0 - 170, sum_y - 40 + 45*(i-1), state_x0 - 120, sum_y - 20 + 45*(i-1)], ...
                'GotoTag', 'u_sig');
            if B(i, 1) == 1
                add_line(spec.model_name, sprintf('from_u_%d/1', i), sprintf('Sum_x%d_next/%d', i, input_port), 'autorouting', 'on');
            else
                label_name = input_gain_label(spec, coeffs, i);
                gain_name = sprintf('%s/%s', spec.model_name, label_name);
                add_block('simulink/Math Operations/Gain', gain_name, ...
                    'Position', [state_x0 - 110, sum_y - 40 + 45*(i-1), state_x0 - 60, sum_y - 20 + 45*(i-1)], ...
                    'Gain', num2str(B(i, 1), 16));
                add_line(spec.model_name, sprintf('from_u_%d/1', i), sprintf('%s/1', label_name), 'autorouting', 'on');
                add_line(spec.model_name, sprintf('%s/1', label_name), sprintf('Sum_x%d_next/%d', i, input_port), 'autorouting', 'on');
            end
            input_port = input_port + 1;
        end

        if B(i, 2) ~= 0
            from_name = sprintf('%s/from_v_%d', spec.model_name, i);
            add_block('simulink/Signal Routing/From', from_name, ...
                'Position', [state_x0 - 170, sum_y + 120 + 30*(i-1), state_x0 - 120, sum_y + 140 + 30*(i-1)], ...
                'GotoTag', 'v_fb');
            if B(i, 2) == 1
                add_line(spec.model_name, sprintf('from_v_%d/1', i), sprintf('Sum_x%d_next/%d', i, input_port), 'autorouting', 'on');
            else
                label_name = feedback_gain_label(spec, coeffs, i);
                gain_name = sprintf('%s/%s', spec.model_name, label_name);
                add_block('simulink/Math Operations/Gain', gain_name, ...
                    'Position', [state_x0 - 110, sum_y + 120 + 30*(i-1), state_x0 - 60, sum_y + 140 + 30*(i-1)], ...
                    'Gain', num2str(B(i, 2), 16));
                add_line(spec.model_name, sprintf('from_v_%d/1', i), sprintf('%s/1', label_name), 'autorouting', 'on');
                add_line(spec.model_name, sprintf('%s/1', label_name), sprintf('Sum_x%d_next/%d', i, input_port), 'autorouting', 'on');
            end
        end
    end

    y_terms = nnz(C ~= 0) + (D(1) ~= 0) + (D(2) ~= 0);
    add_block('simulink/Math Operations/Sum', [spec.model_name '/Sum_y'], ...
        'Position', [900, 95, 930, 115], ...
        'Inputs', repmat('+', 1, max(y_terms, 1)));

    input_port = 1;
    for j = 1:n_states
        if C(j) == 0
            continue;
        end
        if C(j) == 1
            add_line(spec.model_name, sprintf('x%d/1', j), sprintf('Sum_y/%d', input_port), 'autorouting', 'on');
        else
            label_name = output_gain_label(spec, coeffs, j);
            gain_name = sprintf('%s/%s', spec.model_name, label_name);
            gx = 820;
            gy = 60 + 45*(j-1);
            add_block('simulink/Math Operations/Gain', gain_name, ...
                'Position', [gx, gy, gx + 50, gy + 20], ...
                'Gain', num2str(C(j), 16));
            add_line(spec.model_name, sprintf('x%d/1', j), sprintf('%s/1', label_name), 'autorouting', 'on');
            add_line(spec.model_name, sprintf('%s/1', label_name), sprintf('Sum_y/%d', input_port), 'autorouting', 'on');
        end
        input_port = input_port + 1;
    end

    if D(1) ~= 0
        add_block('simulink/Signal Routing/From', [spec.model_name '/from_u_y'], ...
            'Position', [760, 120, 810, 140], ...
            'GotoTag', 'u_sig');
        if D(1) == 1
            add_line(spec.model_name, 'from_u_y/1', sprintf('Sum_y/%d', input_port), 'autorouting', 'on');
        else
            add_block('simulink/Math Operations/Gain', [spec.model_name '/b_direct'], ...
                'Position', [820, 120, 870, 140], ...
                'Gain', num2str(D(1), 16));
            add_line(spec.model_name, 'from_u_y/1', 'b_direct/1', 'autorouting', 'on');
            add_line(spec.model_name, 'b_direct/1', sprintf('Sum_y/%d', input_port), 'autorouting', 'on');
        end
        input_port = input_port + 1;
    end

    if D(2) ~= 0
        add_block('simulink/Signal Routing/From', [spec.model_name '/from_v_y'], ...
            'Position', [760, 270, 810, 290], ...
            'GotoTag', 'v_fb');
        if D(2) == 1
            add_line(spec.model_name, 'from_v_y/1', sprintf('Sum_y/%d', input_port), 'autorouting', 'on');
        else
            add_block('simulink/Math Operations/Gain', [spec.model_name '/Gain_D2'], ...
                'Position', [820, 270, 870, 290], ...
                'Gain', num2str(D(2), 16));
            add_line(spec.model_name, 'from_v_y/1', 'Gain_D2/1', 'autorouting', 'on');
            add_line(spec.model_name, 'Gain_D2/1', sprintf('Sum_y/%d', input_port), 'autorouting', 'on');
        end
    end

    add_quantizer_and_outputs(spec, [970, 80, 1075, 120], 'Sum_y/1', 'Input/1', 'Sum_y/1');
end

function add_quantizer_and_outputs(spec, quant_pos, quant_input, input_scope_src, y_scope_src)
    quantizer_expr = sprintf('dsm_quantizer(u, %.16g, %d)', spec.V_fs, spec.n_levels);

    add_block('simulink/User-Defined Functions/MATLAB Fcn', ...
        [spec.model_name '/IdealQuantizer'], ...
        'Position', quant_pos, ...
        'MATLABFcn', quantizer_expr);

    add_block('simulink/Sinks/Scope', [spec.model_name '/Signals'], ...
        'Position', [1135, 60, 1185, 150], ...
        'NumInputPorts', '3');

    add_block('simulink/Sinks/To Workspace', [spec.model_name '/v_out'], ...
        'Position', [1085, 125, 1135, 145], ...
        'VariableName', 'v_dsm', ...
        'SaveFormat', 'Array');

    add_block('simulink/Sinks/To Workspace', [spec.model_name '/y_out'], ...
        'Position', [980, 150, 1035, 170], ...
        'VariableName', 'y_dsm', ...
        'SaveFormat', 'Array');

    add_block('simulink/Sinks/To Workspace', [spec.model_name '/u_out'], ...
        'Position', [90, 180, 160, 200], ...
        'VariableName', 'u_dsm', ...
        'SaveFormat', 'Array');

    add_line(spec.model_name, input_scope_src, 'u_out/1', 'autorouting', 'on');
    add_line(spec.model_name, quant_input, 'IdealQuantizer/1', 'autorouting', 'on');
    add_line(spec.model_name, 'IdealQuantizer/1', 'v_fb/1', 'autorouting', 'on');
    add_line(spec.model_name, 'IdealQuantizer/1', 'v_out/1', 'autorouting', 'on');
    add_line(spec.model_name, y_scope_src, 'y_out/1', 'autorouting', 'on');
    add_line(spec.model_name, input_scope_src, 'Signals/1', 'autorouting', 'on');
    add_line(spec.model_name, y_scope_src, 'Signals/2', 'autorouting', 'on');
    add_line(spec.model_name, 'IdealQuantizer/1', 'Signals/3', 'autorouting', 'on');
end

function name = input_gain_label(spec, coeffs, idx)
if strcmpi(spec.form, 'CIFF') && ~isempty(coeffs.b) && numel(coeffs.b) >= idx
    name = sprintf('b%d', idx);
else
    name = sprintf('Gain_Bu_%d', idx);
end
end

function name = feedback_gain_label(spec, coeffs, idx)
if strcmpi(spec.form, 'CIFF')
    name = sprintf('vfb%d', idx);
else
    name = sprintf('Gain_Bv_%d', idx);
end
end

function name = output_gain_label(spec, coeffs, idx)
if strcmpi(spec.form, 'CIFF') && ~isempty(coeffs.a) && numel(coeffs.a) >= idx
    name = sprintf('a%d', idx);
else
    name = sprintf('Gain_C_%d', idx);
end
end
