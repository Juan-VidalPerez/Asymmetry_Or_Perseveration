function plot_figS6_S15(varargin)
% PLOT_FIGS6_S15 Generates plots for Figures S6 and S15: Confirmation Bias and Perseveration.
%
% USAGE:
%   plot_figS6_S15()        % Default MAP fits
%   plot_figS6_S15('MLE')   % MLE fits

    %% 1. Initialization & Data Loading
    [data, is_mle, filename] = load_figure_data(varargin{:});
    
    % Configuration
    exps = {'L1','L2','P1','P2','C1','C2','C3','C4','S1a','S1b'};
    offsets = [-0.5, 0, 0.5];
    colors = struct('bias', [234, 172, 139; 181, 101, 118; 151, 189, 161] / 255, ...
                    'pers', [181, 101, 118; 107, 142, 185; 151, 189, 161] / 255);
    
    % Stats Containers: [Experiment x Model/Comparison x Metric]
    p_vs_zero = NaN(length(exps), 3, 3); 
    p_compare = NaN(length(exps), 2, 2); % Only for metrics 1 and 2

    figure('Position', [100, 100, 1200, 950], 'Color', 'w');

    %% 2. Main Plotting & Calculation Loop
    for m_idx = 1:3 % 1=Abs Bias, 2=Norm Bias, 3=Phi
        subplot(3, 1, m_idx); hold on;
        
        % Select models and colors
        if m_idx < 3
            model_group = {data.parameters_CB, data.parameters_CBPERS, data.parameters_CBsim_CBPERSfit};
            current_colors = colors.bias;
        else
            model_group = {data.parameters_CBPERS, data.parameters_PSL, data.parameters_CBsim_CBPERSfit};
            current_colors = colors.pers;
        end

        for e_idx = 1:length(exps)
            vals = cell(1, 3);
            
            for f_idx = 1:3
                raw_params = model_group{f_idx}{e_idx};
                
                % Squeeze 3D simulation data
                if f_idx == 3 && ndims(raw_params) > 2
                    raw_params = squeeze(mean(raw_params, 3));
                end
                
                vals{f_idx} = get_metric(raw_params, m_idx);
                
                % Plotting
                x_pos = e_idx * 3 + offsets(f_idx);
                draw_point(x_pos, vals{f_idx}, current_colors(f_idx, :));
                
                % Stat: vs Zero
                [~, p_vs_zero(e_idx, f_idx, m_idx)] = ttest(vals{f_idx});
            end
            
            % Stat: Between Models (Specific to your requested comparisons)
            if m_idx < 3
                [~, p_compare(e_idx, 1, m_idx)] = ttest(vals{1}, vals{3}); % CB vs CBsim
                [~, p_compare(e_idx, 2, m_idx)] = ttest(vals{2}, vals{3}); % CBPERS vs CBsim
            end
        end
        
        style_axes(m_idx, exps, is_mle);
    end

    %% 3. Statistics Report
    print_report(p_vs_zero, p_compare, filename);
end

%% --- Helper Functions ---

function [data, is_mle, fname] = load_figure_data(varargin)
    is_mle = nargin > 0 && strcmpi(varargin{1}, 'MLE');
    if is_mle
        fname = 'MLE_fits.mat';
    elseif nargin > 0 && strcmpi(varargin{1}, 'MAP_wide')
        fname = 'MAP_fits_wide.mat';
    else
        fname = 'MAP_fits.mat';
    end

    full_path = fullfile('Data', fname);
    try
        data = load(full_path);
        fprintf('--- Generating Figure: Loading %s ---\n', full_path);
    catch ME
        error('Could not load %s. Ensure it is in the "Data" folder. Error: %s', fname, ME.message);
    end
end

function v = get_metric(p, type)
    switch type
        case 1, v = p(:, 2) - p(:, 3); % Abs Bias
        case 2, v = (p(:, 2) - p(:, 3)) ./ (p(:, 2) + p(:, 3)); % Norm Bias
        case 3, v = p(:, end); % Phi
    end
end

function draw_point(x, data, col)
    mu = mean(data);
    err = std(data) / sqrt(length(data));
    errorbar(x, mu, err, '.k', 'CapSize', 0, 'LineWidth', 1.3, 'HandleVisibility', 'off');
    plot(x, mu, 'ko', 'MarkerFaceColor', col, 'MarkerSize', 8, 'LineWidth', 1.3);
end

function style_axes(idx, exps, is_mle)
    plot(xlim, [0 0], 'k--', 'LineWidth', 1.3, 'HandleVisibility', 'off');
    set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'Box', 'off', 'TickDir', 'out');
    xlim([0, length(exps) * 3 + 3]);
    
    labels = {'\alpha_c - \alpha_d', '(\alpha_c - \alpha_d)/(\alpha_c + \alpha_d)', '\phi'};
    ylabel(labels{idx});

    if idx < 3
        xticklabels({});
        legend({'Data (LA fit)', 'Data (hybrid fit)', 'CB sim (hybrid fit)'}, ...
               'Orientation', 'horizontal', 'Location', 'northwest');
        if is_mle && idx == 1, ylim([-0.2, 0.4]); end
    else
        xticks(3:3:length(exps)*3);
        xticklabels(exps);
        xlabel('Experiment', 'FontSize', 14);
        legend({'Data (hybrid fit)', 'Data (PSL fit)', 'CB sim (hybrid fit)'}, ...
               'Orientation', 'horizontal', 'Location', 'northwest');
        if is_mle, ylim([-1, 8]); else, ylim([-0.1, 3]); end
    end
end

function print_report(p_zero, p_comp, fname)
    fprintf('\n======================================================\n');
    fprintf('STATISTICAL TESTS (Data: %s)\n', upper(fname));
    fprintf('======================================================\n');
    
    titles = {'ABS BIAS', 'NORM BIAS'};
    for i = 1:2
        fprintf('\n--- %s ---\n', titles{i});
        fprintf('T-tests vs Zero (1=CB, 2=CBPERS, 3=CBsim):\n');
        disp(p_zero(:,:,i));
        fprintf('Between-Model (1=CB vs CBsim, 2=CBPERS vs CBsim):\n');
        disp(p_comp(:,:,i));
    end
    
    fprintf('\n--- PERSEVERATION (PHI) ---\n');
    fprintf('T-tests vs Zero (1=CBPERS, 2=PSL, 3=CBsim):\n');
    disp(p_zero(:,:,3));
    fprintf('-------------------- END OF TESTS --------------------\n\n');
end