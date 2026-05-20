function plot_fig2(varargin)
% PLOT_FIG2 Generates plots for Figure 2: Confirmation Bias and Perseveration.
%
% USAGE:
%   plot_fig2()           % Default MAP fits (Figure 2)
%   plot_fig2('MLE')      % MLE fits (Figure S14)
%   plot_fig2('MAP_wide') % Wide MAP fits (Figure S7)

    %% 1. Initialization & Data Loading
    [data, is_mle, filename] = load_figure_data(varargin{:});
    
    % Configuration
    exps   = {'L1','L2','P1','P2','C1','C2','C3','C4','S1a','S1b'};
    offsets = [-0.5, 0, 0.5];
    colors = struct('bias', [234, 172, 139; 181, 101, 118; 172, 136, 187] / 255, ...
                    'pers', [181, 101, 118; 107, 142, 185; 172, 136, 187] / 255);
    
    % Stats Containers: [Experiment x Model/Comparison x Metric]
    p_vs_zero = NaN(length(exps), 3, 3); 
    p_compare = NaN(length(exps), 2, 3); 

    figure('Position', [100, 100, 1200, 950], 'Color', 'w');

    %% 2. Main Plotting & Calculation Loop
    for m_idx = 1:3 % Metric Index: 1=Abs Bias, 2=Norm Bias, 3=Phi
        subplot(3, 1, m_idx); hold on;
        
        % Assign specific models for this metric
        if m_idx < 3
            model_group = {data.parameters_CB, data.parameters_CBPERS, data.parameters_PSLsim_CBPERSfit};
            current_colors = colors.bias;
        else
            model_group = {data.parameters_CBPERS, data.parameters_PSL, data.parameters_PSLsim_CBPERSfit};
            current_colors = colors.pers;
        end

        for e_idx = 1:length(exps)
            vals = cell(1, 3);
            
            for f_idx = 1:3 % Model index within subplot
                raw_params = model_group{f_idx}{e_idx};
                
                % Squeeze simulation data (3D) to 2D if necessary
                if f_idx == 3 && ndims(raw_params) > 2
                    raw_params = squeeze(mean(raw_params, 3));
                end
                
                % Calculate metric and store for stats
                vals{f_idx} = calculate_metric_values(raw_params, m_idx);
                
                % Plot Data Point
                x_pos = e_idx * 3 + offsets(f_idx);
                render_point(x_pos, vals{f_idx}, current_colors(f_idx, :));
                
                % Statistical Test: vs Zero
                [~, p_vs_zero(e_idx, f_idx, m_idx)] = ttest(vals{f_idx});
            end
            
            % Statistical Test: Between Models (1 vs 2, 2 vs 3)
            [~, p_compare(e_idx, 1, m_idx)] = ttest(vals{1}, vals{2});
            [~, p_compare(e_idx, 2, m_idx)] = ttest(vals{2}, vals{3});
        end
        
        apply_formatting(m_idx, exps, is_mle);
    end

    %% 3. Output Detailed Statistics
    display_stats_report(p_vs_zero, p_compare, filename);
end

%% --- Support Functions ---

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

function v = calculate_metric_values(p, type)
    switch type
        case 1 % Absolute Bias: (Alpha_Confirm - Alpha_Disconfirm)
            v = p(:, 2) - p(:, 3);
        case 2 % Normalized Bias
            v = (p(:, 2) - p(:, 3)) ./ (p(:, 2) + p(:, 3));
        case 3 % Perseveration (Phi)
            v = p(:, end);
    end
end

function render_point(x, data_vec, col)
    mu = mean(data_vec);
    err = std(data_vec) / sqrt(length(data_vec));
    errorbar(x, mu, err, '.k', 'CapSize', 0, 'LineWidth', 1.3, 'HandleVisibility', 'off');
    plot(x, mu, 'ko', 'MarkerFaceColor', col, 'MarkerSize', 8, 'LineWidth', 1.3);
end

function apply_formatting(idx, exps, is_mle)
    plot(xlim, [0 0], 'k--', 'LineWidth', 1.3, 'HandleVisibility', 'off');
    set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'Box', 'off', 'TickDir', 'out');
    xlim([0, length(exps) * 3 + 3]);
    xticks(3:3:length(exps)*3);
    xticklabels(exps);

    if idx == 1
        ylabel('\alpha_c - \alpha_d');
        if is_mle, ylim([-0.2, 0.4]); end
        legend({'Data (LA fit)', 'Data (hybrid fit)', 'PSL sim (hybrid fit)'}, ...
               'Orientation', 'horizontal', 'Location', 'northwest');
    elseif idx == 2
        ylabel('(\alpha_c - \alpha_d)/(\alpha_c + \alpha_d)');
    else
        ylabel('\phi');
        if is_mle, ylim([-1, 8]); else, ylim([-0.1, 3]); end
        xlabel('Experiment', 'FontSize', 14);
        legend({'Data (hybrid fit)', 'Data (PSL fit)', 'PSL sim (hybrid fit)'}, ...
               'Orientation', 'horizontal', 'Location', 'northwest');
    end
end

function display_stats_report(p_zero, p_comp, fname)
    fprintf('\n======================================================\n');
    fprintf('STATISTICAL TESTS (Data: %s)\n', upper(fname));
    fprintf('======================================================\n');
    
    titles = {'ABS CONFIRMATION BIAS', 'NORM CONFIRMATION BIAS', 'PERSEVERATION (PHI)'};
    comp_labels = {'1 vs 2', '2 vs 3'};
    
    for i = 1:3
        fprintf('\n--- %s ---\n', titles{i});
        fprintf('T-tests vs Zero (Models 1, 2, 3):\n');
        disp(p_zero(:,:,i));
        fprintf('Between-Model T-tests (%s):\n', strjoin(comp_labels, ', '));
        disp(p_comp(:,:,i));
    end
    fprintf('-------------------- END OF TESTS --------------------\n\n');
end