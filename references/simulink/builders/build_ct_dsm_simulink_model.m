%% BUILD_CT_DSM_SIMULINK_MODEL
% Build a generic continuous-time single-quantizer DSM Simulink model.
%
% Supported entry points:
%   build_ct_dsm_simulink_model(spec_struct)
%   build_ct_dsm_simulink_model('name', value, ...)
%
% The builder accepts either:
%   1. ABCDc
%   2. Ac, Bc, Cc, Dc
%   3. ntf + ct_form + tdac
%   4. order, OSR, H_inf, opt, f0, ct_form, tdac
%
% The generated model uses:
%   - continuous-time integrators for the loop filter states
%   - an explicit sampler before the quantizer
%   - piecewise-constant DAC feedback branches
%
% Current timing support covers delayed DAC pulses whose width is no more
% than one sample period, i.e. 0 <= tdac(:,1) < tdac(:,2) and
% tdac(:,2) - tdac(:,1) <= 1.

function [model_name, model_path, info] = build_ct_dsm_simulink_model(varargin)
    spec = parse_inputs(varargin{:});

    addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'common'));
    setup_references_paths();

    if isempty(spec.fs)
        error('build_ct_dsm_simulink_model:MissingFs', ...
            'Sampling frequency fs must be provided.');
    end

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

    [ABCDc, tdac2, ntf] = resolve_ct_model_data(spec);
    [Ac, Bc, Cc, Dc] = partitionABCD(ABCDc);
    n_states = size(Ac, 1);
    n_inputs = size(Bc, 2);
    n_feedback = n_inputs - 1;

    if n_feedback < 1
        error('build_ct_dsm_simulink_model:UnsupportedABCDc', ...
            'ABCDc must contain one signal input and at least one DAC feedback input.');
    end

    if size(Dc, 2) ~= n_inputs
        error('build_ct_dsm_simulink_model:InvalidFeedthrough', ...
            'Bc and Dc must have the same number of inputs.');
    end

    if size(tdac2, 1) ~= n_feedback + 1
        error('build_ct_dsm_simulink_model:TdacMismatch', ...
            'The tdac2 matrix must have one row per CT input, including the input signal row.');
    end

    unsupported = tdac2(2:end, 1) < 0 | ...
        tdac2(2:end, 2) <= tdac2(2:end, 1) | ...
        (tdac2(2:end, 2) - tdac2(2:end, 1)) > 1;
    if any(unsupported)
        error('build_ct_dsm_simulink_model:UnsupportedTdac', ...
            ['This builder currently supports delayed DAC pulses with width ' ...
             'no more than one sample period: 0 <= t1 < t2 and t2-t1 <= 1.']);
    end

    if bdIsLoaded(spec.model_name)
        close_system(spec.model_name, 0);
    end

    model_path = fullfile(spec.output_dir, [spec.model_name '.slx']);
    if exist(model_path, 'file')
        delete(model_path);
    end

    new_system(spec.model_name);

    Ts = 1/spec.fs;
    ct_step = Ts / spec.ct_steps_per_sample;
    set_param(spec.model_name, 'SolverType', 'Fixed-step');
    set_param(spec.model_name, 'Solver', spec.solver);
    set_param(spec.model_name, 'FixedStep', num2str(ct_step, 16));
    set_param(spec.model_name, 'StopTime', num2str(spec.stop_time, 16));

    mw = get_param(spec.model_name, 'ModelWorkspace');
    assignin(mw, 'Ac_ct', Ac);
    assignin(mw, 'Bc_ct', Bc);
    assignin(mw, 'Cc_ct', Cc);
    assignin(mw, 'Dc_ct', Dc);
    assignin(mw, 'Ts_dsm', Ts);
    assignin(mw, 'x0_ct', zeros(n_states, 1));

    build_ct_loop(spec, Ac, Bc, Cc, Dc, tdac2);

    model_name = spec.model_name;
    save_system(spec.model_name, model_path);

    if spec.open_model
        open_system(spec.model_name);
    else
        close_system(spec.model_name);
    end

    info = struct();
    info.ABCDc = ABCDc;
    info.Ac = Ac;
    info.Bc = Bc;
    info.Cc = Cc;
    info.Dc = Dc;
    info.tdac = spec.tdac;
    info.tdac2 = tdac2;
    info.ntf = ntf;
    info.ct_form = spec.ct_form;
    info.fs = spec.fs;
    info.n_levels = spec.n_levels;
    info.input_amplitude = spec.input_amplitude;
    info.input_frequency = spec.input_frequency;
    info.ct_steps_per_sample = spec.ct_steps_per_sample;
    info.jitter_rms = spec.jitter_rms;
    info.jitter_mode = spec.jitter_mode;
    info.jitter_seed = spec.jitter_seed;
    info.model_path = model_path;
end

function spec = parse_inputs(varargin)
    defaults = struct( ...
        'model_name', 'ct_dsm_simulink_model', ...
        'output_dir', pwd, ...
        'ABCDc', [], ...
        'Ac', [], ...
        'Bc', [], ...
        'Cc', [], ...
        'Dc', [], ...
        'ntf', [], ...
        'order', [], ...
        'OSR', [], ...
        'H_inf', [], ...
        'opt', 1, ...
        'f0', 0, ...
        'ct_form', 'FF', ...
        'tdac', [0 1], ...
        'fs', [], ...
        'n_bits', [], ...
        'n_levels', [], ...
        'V_fs', 1.0, ...
        'input_amplitude', 0.5, ...
        'input_frequency', [], ...
        'stop_time', [], ...
        'ct_steps_per_sample', 32, ...
        'jitter_rms', 0, ...
        'jitter_mode', 'edge', ...
        'jitter_seed', 12345, ...
        'solver', 'ode4', ...
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

function [ABCDc, tdac2, ntf] = resolve_ct_model_data(spec)
    ntf = [];

    if ~isempty(spec.ABCDc)
        ABCDc = spec.ABCDc;
        tdac2 = normalize_tdac_for_abcdc(spec.tdac, ABCDc);
        return;
    end

    if ~isempty(spec.Ac) && ~isempty(spec.Bc) && ~isempty(spec.Cc) && ~isempty(spec.Dc)
        ABCDc = [spec.Ac spec.Bc; spec.Cc spec.Dc];
        tdac2 = normalize_tdac_for_abcdc(spec.tdac, ABCDc);
        return;
    end

    if ~isempty(spec.ntf)
        ntf = spec.ntf;
        [ABCDc, tdac2] = realizeNTF_ct(ntf, spec.ct_form, spec.tdac);
        return;
    end

    required = {'order', 'OSR', 'H_inf'};
    for i = 1:numel(required)
        if isempty(spec.(required{i}))
            error('build_ct_dsm_simulink_model:MissingDesignData', ...
                ['Need ABCDc, state-space matrices, ntf, or full design ' ...
                 'inputs (order, OSR, H_inf, ...).']);
        end
    end

    ntf = synthesizeNTF(spec.order, spec.OSR, spec.opt, spec.H_inf, spec.f0);
    [ABCDc, tdac2] = realizeNTF_ct(ntf, spec.ct_form, spec.tdac);
end

function tdac2 = normalize_tdac_for_abcdc(tdac, ABCDc)
    n_states = size(ABCDc, 1) - 1;
    n_inputs = size(ABCDc, 2) - n_states;

    if iscell(tdac)
        error('build_ct_dsm_simulink_model:UnsupportedTdacCell', ...
            'When passing ABCDc directly, provide tdac as a numeric matrix.');
    end

    if isempty(tdac)
        tdac = [0 1];
    end

    if isequal(size(tdac), [1 2])
        tdac2 = [-1 -1; repmat(tdac, n_inputs-1, 1)];
        return;
    end

    if isequal(size(tdac), [n_inputs 2])
        tdac2 = tdac;
        return;
    end

    error('build_ct_dsm_simulink_model:InvalidTdac', ...
        'tdac must be either a 1x2 vector or an nInputs x 2 matrix.');
end

function build_ct_loop(spec, Ac, Bc, Cc, Dc, tdac2)
    model = spec.model_name;
    n_states = size(Ac, 1);
    n_feedback = size(Bc, 2) - 1;

    add_block('simulink/Sources/Sine Wave', [model '/Input'], ...
        'Position', [40, 245, 90, 275], ...
        'Amplitude', num2str(spec.input_amplitude, 16), ...
        'Frequency', num2str(2*pi*spec.input_frequency, 16), ...
        'SampleTime', '0', ...
        'Bias', '0');

    add_block('simulink/Signal Routing/Goto', [model '/u_sig'], ...
        'Position', [115, 245, 165, 265], ...
        'GotoTag', 'u_sig', ...
        'TagVisibility', 'global');
    add_line(model, 'Input/1', 'u_sig/1');

    build_feedback_branches(spec, n_feedback, tdac2);
    build_state_equations(spec, Ac, Bc);
    build_output_equation(spec, Cc, Dc);
    add_sampler_quantizer_outputs(spec, n_feedback);
end

function build_feedback_branches(spec, n_feedback, tdac2)
    model = spec.model_name;

    add_block('simulink/Signal Routing/Goto', [model '/v_q'], ...
        'Position', [1120, 105, 1170, 125], ...
        'GotoTag', 'v_q', ...
        'TagVisibility', 'global');

    for k = 1:n_feedback
        t1 = tdac2(k+1, 1);
        t2 = tdac2(k+1, 2);
        branch_tag = sprintf('v_dac_%d', k);

        add_block('simulink/Signal Routing/From', [model sprintf('/from_vq_%d', k)], ...
            'Position', [1110, 180 + 70*(k-1), 1160, 200 + 70*(k-1)], ...
            'GotoTag', 'v_q');

        pulse_width = t2 - t1;
        branch_src = sprintf('from_vq_%d/1', k);

        if spec.jitter_rms > 0
            if strcmpi(spec.jitter_mode, 'equivalent')
                if abs(pulse_width - 1) > 1e-12
                    error('build_ct_dsm_simulink_model:UnsupportedEquivalentJitterPulse', ...
                        'Equivalent DAC jitter modeling currently supports only full-width NRZ-like DAC pulses.');
                end
                branch_src = add_equivalent_jitter_path(spec, k, branch_src);
            elseif ~strcmpi(spec.jitter_mode, 'edge')
                error('build_ct_dsm_simulink_model:UnsupportedJitterMode', ...
                    'jitter_mode must be either ''edge'' or ''equivalent''.');
            end
        end

        if abs(t1) < 1e-12 && abs(t2 - 1) < 1e-12
            if spec.jitter_rms > 0 && strcmpi(spec.jitter_mode, 'edge')
                error('build_ct_dsm_simulink_model:NeedPositiveDelayForJitter', ...
                    'Clock jitter modeling currently requires a positive nominal DAC delay.');
            end
            add_block('simulink/Signal Routing/Goto', [model sprintf('/%s', branch_tag)], ...
                'Position', [1200, 180 + 70*(k-1), 1260, 200 + 70*(k-1)], ...
                'GotoTag', branch_tag, ...
                'TagVisibility', 'global');
            add_line(model, branch_src, sprintf('%s/1', branch_tag));
        else
            delayed_name = sprintf('dac_delay_%d', k);
            delayed_src = branch_src;
            if abs(t1) > 1e-12 || (spec.jitter_rms > 0 && strcmpi(spec.jitter_mode, 'edge'))
                if spec.jitter_rms > 0 && strcmpi(spec.jitter_mode, 'edge') && abs(pulse_width - 1) > 1e-12
                    error('build_ct_dsm_simulink_model:UnsupportedJitterPulse', ...
                        'Clock jitter modeling is currently supported only for full-width NRZ-like DAC pulses.');
                end

                if spec.jitter_rms > 0 && strcmpi(spec.jitter_mode, 'edge')
                    add_block('simulink/Sources/Random Number', [model sprintf('/dac_jitter_%d', k)], ...
                        'Position', [1125, 225 + 70*(k-1), 1215, 255 + 70*(k-1)], ...
                        'Mean', '0', ...
                        'Variance', num2str(spec.jitter_rms^2, 16), ...
                        'SampleTime', 'Ts_dsm', ...
                        'Seed', num2str(spec.jitter_seed + k - 1));
                    add_block('simulink/Sources/Constant', [model sprintf('/dac_delay_nom_%d', k)], ...
                        'Position', [1125, 270 + 70*(k-1), 1215, 290 + 70*(k-1)], ...
                        'Value', sprintf('%.16g*Ts_dsm', t1));
                    add_block('simulink/Math Operations/Sum', [model sprintf('/dac_delay_sum_%d', k)], ...
                        'Position', [1245, 232 + 70*(k-1), 1275, 258 + 70*(k-1)], ...
                        'Inputs', '++');
                    add_block('simulink/Discontinuities/Saturation', [model sprintf('/dac_delay_clip_%d', k)], ...
                        'Position', [1305, 230 + 70*(k-1), 1355, 260 + 70*(k-1)], ...
                        'LowerLimit', '0', ...
                        'UpperLimit', 'inf');
                    add_block('simulink/Continuous/Variable Time Delay', [model sprintf('/%s', delayed_name)], ...
                        'Position', [1390, 180 + 70*(k-1), 1485, 210 + 70*(k-1)], ...
                        'MaximumDelay', sprintf('%.16g*Ts_dsm', max(t2 + 0.1, t1 + 8*spec.jitter_rms*spec.fs)));
                    add_line(model, delayed_src, sprintf('%s/1', delayed_name), 'autorouting', 'on');
                    add_line(model, sprintf('dac_jitter_%d/1', k), sprintf('dac_delay_sum_%d/1', k), 'autorouting', 'on');
                    add_line(model, sprintf('dac_delay_nom_%d/1', k), sprintf('dac_delay_sum_%d/2', k), 'autorouting', 'on');
                    add_line(model, sprintf('dac_delay_sum_%d/1', k), sprintf('dac_delay_clip_%d/1', k), 'autorouting', 'on');
                    add_line(model, sprintf('dac_delay_clip_%d/1', k), sprintf('%s/2', delayed_name), 'autorouting', 'on');
                else
                    add_block('simulink/Continuous/Transport Delay', [model sprintf('/%s', delayed_name)], ...
                        'Position', [1190, 180 + 70*(k-1), 1270, 200 + 70*(k-1)], ...
                        'DelayTime', sprintf('%.16g*Ts_dsm', t1));
                    add_line(model, delayed_src, sprintf('%s/1', delayed_name), 'autorouting', 'on');
                end
                delayed_src = sprintf('%s/1', delayed_name);
            end

            if abs(pulse_width - 1) < 1e-12
                add_block('simulink/Signal Routing/Goto', [model sprintf('/%s', branch_tag)], ...
                    'Position', [1510, 180 + 70*(k-1), 1570, 200 + 70*(k-1)], ...
                    'GotoTag', branch_tag, ...
                    'TagVisibility', 'global');
                add_line(model, delayed_src, sprintf('%s/1', branch_tag), 'autorouting', 'on');
            else
                pulse_width_pct = 100 * pulse_width;
                add_block('simulink/Sources/Pulse Generator', [model sprintf('/dac_gate_%d', k)], ...
                    'Position', [1190, 230 + 70*(k-1), 1260, 260 + 70*(k-1)], ...
                    'PulseType', 'Time based', ...
                    'Amplitude', '1', ...
                    'Period', 'Ts_dsm', ...
                    'PulseWidth', num2str(pulse_width_pct, 16), ...
                    'PhaseDelay', sprintf('%.16g*Ts_dsm', t1));
                add_block('simulink/Math Operations/Product', [model sprintf('/dac_shape_%d', k)], ...
                    'Position', [1295, 195 + 70*(k-1), 1325, 225 + 70*(k-1)]);
                add_block('simulink/Signal Routing/Goto', [model sprintf('/%s', branch_tag)], ...
                    'Position', [1360, 200 + 70*(k-1), 1420, 220 + 70*(k-1)], ...
                    'GotoTag', branch_tag, ...
                    'TagVisibility', 'global');
                add_line(model, delayed_src, sprintf('dac_shape_%d/1', k), 'autorouting', 'on');
                add_line(model, sprintf('dac_gate_%d/1', k), sprintf('dac_shape_%d/2', k), 'autorouting', 'on');
                add_line(model, sprintf('dac_shape_%d/1', k), sprintf('%s/1', branch_tag), 'autorouting', 'on');
            end
        end
    end
end

function branch_src = add_equivalent_jitter_path(spec, branch_idx, input_src)
    model = spec.model_name;
    y0 = 180 + 70*(branch_idx-1);
    branch_src = sprintf('dac_jitter_sum_%d/1', branch_idx);

    add_block('simulink/Discrete/Unit Delay', [model sprintf('/dac_prev_%d', branch_idx)], ...
        'Position', [1185, y0 - 10, 1225, y0 + 20], ...
        'SampleTime', 'Ts_dsm', ...
        'InitialCondition', '0');
    add_block('simulink/Math Operations/Sum', [model sprintf('/dac_step_diff_%d', branch_idx)], ...
        'Position', [1260, y0 - 8, 1290, y0 + 18], ...
        'Inputs', '+-');
    add_block('simulink/Sources/Random Number', [model sprintf('/dac_jitter_%d', branch_idx)], ...
        'Position', [1185, y0 + 38, 1275, y0 + 68], ...
        'Mean', '0', ...
        'Variance', num2str(spec.jitter_rms^2, 16), ...
        'SampleTime', 'Ts_dsm', ...
        'Seed', num2str(spec.jitter_seed + branch_idx - 1));
    add_block('simulink/Math Operations/Gain', [model sprintf('/dac_jitter_scale_%d', branch_idx)], ...
        'Position', [1315, y0 + 42, 1375, y0 + 64], ...
        'Gain', '1/Ts_dsm');
    add_block('simulink/Math Operations/Product', [model sprintf('/dac_jitter_err_%d', branch_idx)], ...
        'Position', [1410, y0 + 2, 1440, y0 + 32]);
    add_block('simulink/Math Operations/Sum', [model sprintf('/dac_jitter_sum_%d', branch_idx)], ...
        'Position', [1475, y0 - 8, 1505, y0 + 18], ...
        'Inputs', '++');

    add_line(model, input_src, sprintf('dac_prev_%d/1', branch_idx), 'autorouting', 'on');
    add_line(model, input_src, sprintf('dac_step_diff_%d/1', branch_idx), 'autorouting', 'on');
    add_line(model, sprintf('dac_prev_%d/1', branch_idx), sprintf('dac_step_diff_%d/2', branch_idx), 'autorouting', 'on');
    add_line(model, sprintf('dac_jitter_%d/1', branch_idx), sprintf('dac_jitter_scale_%d/1', branch_idx), 'autorouting', 'on');
    add_line(model, sprintf('dac_step_diff_%d/1', branch_idx), sprintf('dac_jitter_err_%d/1', branch_idx), 'autorouting', 'on');
    add_line(model, sprintf('dac_jitter_scale_%d/1', branch_idx), sprintf('dac_jitter_err_%d/2', branch_idx), 'autorouting', 'on');
    add_line(model, input_src, sprintf('dac_jitter_sum_%d/1', branch_idx), 'autorouting', 'on');
    add_line(model, sprintf('dac_jitter_err_%d/1', branch_idx), sprintf('dac_jitter_sum_%d/2', branch_idx), 'autorouting', 'on');
end

function build_state_equations(spec, Ac, Bc)
    model = spec.model_name;
    n_states = size(Ac, 1);
    n_feedback = size(Bc, 2) - 1;

    state_x0 = 330;
    state_dx = 170;
    state_y = 250;
    gain_y0 = 115;
    gain_dy = 38;

    for i = 1:n_states
        sum_terms = nnz(Ac(i, :) ~= 0) + (Bc(i, 1) ~= 0) + nnz(Bc(i, 2:end) ~= 0);
        sum_terms = max(sum_terms, 1);

        add_block('simulink/Math Operations/Sum', [model sprintf('/Sum_xdot_%d', i)], ...
            'Position', [state_x0 - 55 + (i-1)*state_dx, state_y + 5, ...
                         state_x0 - 25 + (i-1)*state_dx, state_y + 25], ...
            'Inputs', repmat('+', 1, sum_terms));
        add_block('simulink/Math Operations/Gain', [model sprintf('/fs_scale_%d', i)], ...
            'Position', [state_x0 - 5 + (i-1)*state_dx, state_y + 5, ...
                         state_x0 + 45 + (i-1)*state_dx, state_y + 25], ...
            'Gain', '1/Ts_dsm');
        add_block('simulink/Continuous/Integrator', [model sprintf('/x%d', i)], ...
            'Position', [state_x0 + 70 + (i-1)*state_dx, state_y, ...
                         state_x0 + 110 + (i-1)*state_dx, state_y + 30], ...
            'InitialCondition', sprintf('x0_ct(%d)', i));
        add_line(model, sprintf('Sum_xdot_%d/1', i), sprintf('fs_scale_%d/1', i));
        add_line(model, sprintf('fs_scale_%d/1', i), sprintf('x%d/1', i));
    end

    for i = 1:n_states
        input_port = 1;

        for j = 1:n_states
            if Ac(i, j) == 0
                continue;
            end
            src = sprintf('x%d/1', j);
            if Ac(i, j) == 1
                add_line(model, src, sprintf('Sum_xdot_%d/%d', i, input_port), 'autorouting', 'on');
            else
                gain_name = sprintf('Gain_Ac_%d_%d', i, j);
                gx = state_x0 + (j-1)*state_dx + 55;
                gy = gain_y0 + (j-1)*gain_dy;
                add_block('simulink/Math Operations/Gain', [model '/' gain_name], ...
                    'Position', [gx, gy, gx + 50, gy + 20], ...
                    'Gain', num2str(Ac(i, j), 16));
                add_line(model, src, [gain_name '/1'], 'autorouting', 'on');
                add_line(model, [gain_name '/1'], sprintf('Sum_xdot_%d/%d', i, input_port), 'autorouting', 'on');
            end
            input_port = input_port + 1;
        end

        if Bc(i, 1) ~= 0
            add_block('simulink/Signal Routing/From', [model sprintf('/from_u_%d', i)], ...
                'Position', [state_x0 - 170, state_y - 85 + 42*(i-1), state_x0 - 120, state_y - 65 + 42*(i-1)], ...
                'GotoTag', 'u_sig');
            if Bc(i, 1) == 1
                add_line(model, sprintf('from_u_%d/1', i), sprintf('Sum_xdot_%d/%d', i, input_port), 'autorouting', 'on');
            else
                gain_name = sprintf('Gain_Bu_%d', i);
                add_block('simulink/Math Operations/Gain', [model '/' gain_name], ...
                    'Position', [state_x0 - 110, state_y - 85 + 42*(i-1), state_x0 - 60, state_y - 65 + 42*(i-1)], ...
                    'Gain', num2str(Bc(i, 1), 16));
                add_line(model, sprintf('from_u_%d/1', i), [gain_name '/1'], 'autorouting', 'on');
                add_line(model, [gain_name '/1'], sprintf('Sum_xdot_%d/%d', i, input_port), 'autorouting', 'on');
            end
            input_port = input_port + 1;
        end

        for k = 1:n_feedback
            if Bc(i, k+1) == 0
                continue;
            end
            add_block('simulink/Signal Routing/From', [model sprintf('/from_vdac_%d_%d', i, k)], ...
                'Position', [state_x0 - 175, state_y + 120 + 32*(i-1) + 22*(k-1), ...
                             state_x0 - 120, state_y + 140 + 32*(i-1) + 22*(k-1)], ...
                'GotoTag', sprintf('v_dac_%d', k));
            if Bc(i, k+1) == 1
                add_line(model, sprintf('from_vdac_%d_%d/1', i, k), sprintf('Sum_xdot_%d/%d', i, input_port), 'autorouting', 'on');
            else
                gain_name = sprintf('Gain_Bv_%d_%d', i, k);
                add_block('simulink/Math Operations/Gain', [model '/' gain_name], ...
                    'Position', [state_x0 - 110, state_y + 120 + 32*(i-1) + 22*(k-1), ...
                                 state_x0 - 60, state_y + 140 + 32*(i-1) + 22*(k-1)], ...
                    'Gain', num2str(Bc(i, k+1), 16));
                add_line(model, sprintf('from_vdac_%d_%d/1', i, k), [gain_name '/1'], 'autorouting', 'on');
                add_line(model, [gain_name '/1'], sprintf('Sum_xdot_%d/%d', i, input_port), 'autorouting', 'on');
            end
            input_port = input_port + 1;
        end
    end
end

function build_output_equation(spec, Cc, Dc)
    model = spec.model_name;
    n_states = numel(Cc);
    n_feedback = numel(Dc) - 1;

    y_terms = nnz(Cc ~= 0) + (Dc(1) ~= 0) + nnz(Dc(2:end) ~= 0);
    y_terms = max(y_terms, 1);

    add_block('simulink/Math Operations/Sum', [model '/Sum_y'], ...
        'Position', [930, 95, 960, 115], ...
        'Inputs', repmat('+', 1, y_terms));

    input_port = 1;
    for j = 1:n_states
        if Cc(j) == 0
            continue;
        end
        if Cc(j) == 1
            add_line(model, sprintf('x%d/1', j), sprintf('Sum_y/%d', input_port), 'autorouting', 'on');
        else
            gain_name = sprintf('Gain_Cc_%d', j);
            gx = 850;
            gy = 55 + 40*(j-1);
            add_block('simulink/Math Operations/Gain', [model '/' gain_name], ...
                'Position', [gx, gy, gx + 50, gy + 20], ...
                'Gain', num2str(Cc(j), 16));
            add_line(model, sprintf('x%d/1', j), [gain_name '/1'], 'autorouting', 'on');
            add_line(model, [gain_name '/1'], sprintf('Sum_y/%d', input_port), 'autorouting', 'on');
        end
        input_port = input_port + 1;
    end

    if Dc(1) ~= 0
        add_block('simulink/Signal Routing/From', [model '/from_u_y'], ...
            'Position', [760, 120, 810, 140], ...
            'GotoTag', 'u_sig');
        if Dc(1) == 1
            add_line(model, 'from_u_y/1', sprintf('Sum_y/%d', input_port), 'autorouting', 'on');
        else
            add_block('simulink/Math Operations/Gain', [model '/Gain_Du'], ...
                'Position', [820, 120, 870, 140], ...
                'Gain', num2str(Dc(1), 16));
            add_line(model, 'from_u_y/1', 'Gain_Du/1', 'autorouting', 'on');
            add_line(model, 'Gain_Du/1', sprintf('Sum_y/%d', input_port), 'autorouting', 'on');
        end
        input_port = input_port + 1;
    end

    for k = 1:n_feedback
        if Dc(k+1) == 0
            continue;
        end
        add_block('simulink/Signal Routing/From', [model sprintf('/from_vdac_y_%d', k)], ...
            'Position', [760, 260 + 35*(k-1), 815, 280 + 35*(k-1)], ...
            'GotoTag', sprintf('v_dac_%d', k));
        if Dc(k+1) == 1
            add_line(model, sprintf('from_vdac_y_%d/1', k), sprintf('Sum_y/%d', input_port), 'autorouting', 'on');
        else
            gain_name = sprintf('Gain_Dv_%d', k);
            add_block('simulink/Math Operations/Gain', [model '/' gain_name], ...
                'Position', [825, 260 + 35*(k-1), 875, 280 + 35*(k-1)], ...
                'Gain', num2str(Dc(k+1), 16));
            add_line(model, sprintf('from_vdac_y_%d/1', k), [gain_name '/1'], 'autorouting', 'on');
            add_line(model, [gain_name '/1'], sprintf('Sum_y/%d', input_port), 'autorouting', 'on');
        end
        input_port = input_port + 1;
    end
end

function add_sampler_quantizer_outputs(spec, n_feedback)
    model = spec.model_name;
    quantizer_expr = sprintf('dsm_quantizer(u, %.16g, %d)', spec.V_fs, spec.n_levels);

    add_block('simulink/Discrete/Zero-Order Hold', [model '/Sampler'], ...
        'Position', [1000, 90, 1060, 120], ...
        'SampleTime', 'Ts_dsm');

    add_block('simulink/User-Defined Functions/MATLAB Fcn', [model '/IdealQuantizer'], ...
        'Position', [1085, 90, 1180, 120], ...
        'MATLABFcn', quantizer_expr);

    add_block('simulink/Sinks/Scope', [model '/Signals'], ...
        'Position', [1500, 70, 1555, 180], ...
        'NumInputPorts', '4');

    add_block('simulink/Sinks/To Workspace', [model '/u_out'], ...
        'Position', [95, 300, 165, 320], ...
        'VariableName', 'u_dsm', ...
        'SaveFormat', 'Array');
    add_block('simulink/Sinks/To Workspace', [model '/y_ct_out'], ...
        'Position', [970, 150, 1040, 170], ...
        'VariableName', 'y_ct_dsm', ...
        'SaveFormat', 'Array');
    add_block('simulink/Sinks/To Workspace', [model '/y_samp_out'], ...
        'Position', [1090, 145, 1160, 165], ...
        'VariableName', 'y_dsm', ...
        'SaveFormat', 'Array');
    add_block('simulink/Sinks/To Workspace', [model '/v_out'], ...
        'Position', [1205, 145, 1270, 165], ...
        'VariableName', 'v_dsm', ...
        'SaveFormat', 'Array');

    add_line(model, 'Input/1', 'u_out/1', 'autorouting', 'on');
    add_line(model, 'Sum_y/1', 'Sampler/1', 'autorouting', 'on');
    add_line(model, 'Sum_y/1', 'y_ct_out/1', 'autorouting', 'on');
    add_line(model, 'Sampler/1', 'IdealQuantizer/1', 'autorouting', 'on');
    add_line(model, 'Sampler/1', 'y_samp_out/1', 'autorouting', 'on');
    add_line(model, 'IdealQuantizer/1', 'v_q/1', 'autorouting', 'on');
    add_line(model, 'IdealQuantizer/1', 'v_out/1', 'autorouting', 'on');

    add_line(model, 'Input/1', 'Signals/1', 'autorouting', 'on');
    add_line(model, 'Sum_y/1', 'Signals/2', 'autorouting', 'on');
    add_line(model, 'Sampler/1', 'Signals/3', 'autorouting', 'on');
    add_line(model, 'IdealQuantizer/1', 'Signals/4', 'autorouting', 'on');

    if n_feedback < 1
        error('build_ct_dsm_simulink_model:MissingFeedback', ...
            'At least one feedback branch is required.');
    end
end
