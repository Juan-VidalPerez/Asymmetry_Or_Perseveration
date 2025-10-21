function plot_fig2(varargin)
% Generates plots for Figure 2, combining confirmation bias and perseveration.
%
% This function creates a single figure with three subplots and prints
% statistical test results to the command window. By default, it uses data
% from MAP fits. An optional argument can specify using MLE fits instead.
%
% USAGE:
%   plot_fig2()           % --- For Figure 2 (uses Data/MAP_fits.mat)
%   plot_fig2('MLE')      % --- For Figure S10 (uses Data/MLE_fits.mat)
%
% INPUT:
%   varargin{1} (optional) - String, 'MLE', to specify using Maximum
%                            Likelihood Estimation results. If omitted,
%                            defaults to 'MAP'.
%
% REQUIREMENTS:
%   The files 'MAP_fits.mat' and/or 'MLE_fits.mat' must be in the 'Data'
%   subdirectory relative to the location of this function.

%% Data Loading
data_path = 'Data'; % Define the data subdirectory
is_mle = false; % Flag for conditional formatting
filename = 'MAP_fits.mat';
if nargin > 0 && strcmpi(varargin{1}, 'MLE')
    filename = 'MLE_fits.mat';
    is_mle = true;
end

full_filepath = fullfile(data_path, filename); % Construct the full path
fprintf('--- Generating Figure 2/S10: Loading data from %s ---\n', full_filepath);

% Load pre-computed model fits.
try
    loaded_data = load(full_filepath); % Load into a struct
catch ME
    error('Could not load data. Ensure "%s" is in the "%s" subdirectory. Details: %s', filename, data_path, ME.message);
end

% Assemble the main parameter cell array that the script expects.
if ~isfield(loaded_data, 'parameters_CB') || ~isfield(loaded_data, 'parameters_CBPERS') || ~isfield(loaded_data, 'parameters_PSL') || ~isfield(loaded_data, 'parameters_PSLsim_CBPERSfit')
    error('Could not find required variables in "%s". The file may be corrupt or incomplete.', filename);
end
parameters = {loaded_data.parameters_CB, loaded_data.parameters_CBPERS, loaded_data.parameters_PSL, loaded_data.parameters_PSLsim_CBPERSfit};

%% Setup and Initialization
close all;
experiments = {'L1','L2','P1','P2','C1','C2','C3','C4','S1a','S1b'};
offsets = [-0.5, 0, 0.5];
colors_bias = [234, 172, 139; 181, 101, 118; 172, 136, 187] / 255;
colors_pers = [181, 101, 118; 107, 142, 185; 172, 136, 187] / 255;
bias_data = cell(length(experiments), 3, 2); % {exp, model, metric}
p_bias_vs_zero = NaN(length(experiments), 3, 2);
p_bias_between_models = NaN(length(experiments), 2, 2);
p_phi_vs_zero = NaN(length(experiments), 3);
figure('Position', [100, 100, 1200, 950]);

%% Main Plotting Loop (for 3 subplots)
for metric_idx = 1:3
    subplot(3, 1, metric_idx);
    hold on;
    if metric_idx == 1 || metric_idx == 2 % CONFIRMATION BIAS PLOTS
        params_to_plot = {parameters{1}, parameters{2}, parameters{4}};
        plot_colors = colors_bias;
    else % PERSEVERATION (PHI) PLOT
        params_to_plot = {parameters{2}, parameters{3}, parameters{4}};
        plot_colors = colors_pers;
    end
    for exp_idx = 1:length(experiments)
        for model_idx = 1:3
            params_current = params_to_plot{model_idx}{exp_idx};
            if (model_idx == 3) && (ndims(params_current) > 2)
                 params_current = squeeze(mean(params_current, 3));
            end
            if metric_idx == 1 % Absolute confirmation bias
                metric_value = params_current(:, 2) - params_current(:, 3);
                bias_data{exp_idx, model_idx, 1} = metric_value;
                [~, p_bias_vs_zero(exp_idx, model_idx, 1)] = ttest(metric_value);
            elseif metric_idx == 2 % Normalized confirmation bias
                metric_value = (params_current(:, 2) - params_current(:, 3)) ./ (params_current(:, 2) + params_current(:, 3));
                bias_data{exp_idx, model_idx, 2} = metric_value;
                [~, p_bias_vs_zero(exp_idx, model_idx, 2)] = ttest(metric_value);
            else % Perseveration weight phi
                metric_value = params_current(:, end);
                [~, p_phi_vs_zero(exp_idx, model_idx)] = ttest(metric_value);
            end
            x_pos = exp_idx * 3 + offsets(model_idx);
            sem = std(metric_value) / sqrt(length(metric_value));
            errorbar(x_pos, mean(metric_value), sem, '.k', 'CapSize', 0, 'LineWidth', 1.3);
            plot(x_pos, mean(metric_value), 'ko', 'MarkerFaceColor', plot_colors(model_idx, :), 'MarkerSize', 8, 'LineWidth', 1.3);
        end
        
        if metric_idx == 1 % Absolute bias:
            [~, p_bias_between_models(exp_idx, 1, 1)] = ttest(bias_data{exp_idx, 1, 1}, bias_data{exp_idx, 2, 1}); % CB vs CBPERS
            [~, p_bias_between_models(exp_idx, 2, 1)] = ttest(bias_data{exp_idx, 2, 1}, bias_data{exp_idx, 3, 1}); % CBPERS vs PSLsim
        elseif metric_idx == 2% Normalized bias:
            [~, p_bias_between_models(exp_idx, 1, 2)] = ttest(bias_data{exp_idx, 1, 2}, bias_data{exp_idx, 2, 2}); % CB vs CBPERS
            [~, p_bias_between_models(exp_idx, 2, 2)] = ttest(bias_data{exp_idx, 2, 2}, bias_data{exp_idx, 3, 2}); % CBPERS vs PSLsim
        end
    end
    plot([0, length(experiments) * 3 + 3], [0, 0], 'k--', 'LineWidth', 1.3);
    xlim([0, length(experiments) * 3 + 3]);
    set(gca, 'FontSize', 12, 'LineWidth', 1.2);
    if metric_idx == 1
        ylabel('\alpha_c - \alpha_d');
        if is_mle
            ylim([-0.2, 0.4]);
        end
        legend({'', 'Data (CB fit)', '', 'Data (CBPERS fit)', '', 'PSL sim (CBPERS fit)'}, ...
               'Orientation', 'horizontal', 'Location', 'northwest', 'Box', 'off');
        xticklabels({});
    elseif metric_idx == 2
        ylabel('(\alpha_c - \alpha_d)/(\alpha_c + \alpha_d)');
        xticklabels({});
    else % metric_idx == 3
        ylabel('\phi');
        if is_mle
            ylim([-1, 8]);
        else
            ylim([-0.1, 3]);
        end
        xticks([1:length(experiments)] * 3);
        xticklabels(experiments);
        xlabel('Experiment', 'FontSize', 14);
        legend({'', 'Data (CBPERS fit)', '', 'Data (PSL fit)', '', 'PSL sim (CBPERS fit)'}, ...
               'Orientation', 'horizontal', 'Location', 'northwest', 'Box', 'off');
    end
end
%% Display Statistics in Command Window 📊
fprintf('\n======================================================\n');
fprintf('STATISTICAL TESTS (using data from %s)\n', upper(filename));
fprintf('======================================================\n\n');
fprintf('--- CONFIRMATION BIAS: T-tests vs. Zero ---\n\n');
fprintf('Absolute Bias (alpha_c - alpha_d):\n');
fprintf('Rows are experiments, Columns are models (1=CB, 2=CBPERS, 3=PSLsim)\n');
disp(p_bias_vs_zero(:, :, 1));
fprintf('Normalized Bias:\n');
fprintf('Rows are experiments, Columns are models (1=CB, 2=CBPERS, 3=PSLsim)\n');
disp(p_bias_vs_zero(:, :, 2));
fprintf('\n--- CONFIRMATION BIAS: Between-Model T-tests ---\n\n');
fprintf('Absolute Bias (alpha_c - alpha_d):\n');
fprintf('Rows are experiments, Columns are comparisons (1=CB vs CBPERS, 2=CBPERS vs PSLsim)\n');
disp(p_bias_between_models(:, :, 1));
fprintf('Normalized Bias:\n');
fprintf('Rows are experiments, Columns are comparisons (1=CB vs CBPERS, 2=CBPERS vs PSLsim)\n');
disp(p_bias_between_models(:, :, 2));
fprintf('\n--- PERSEVERATION (phi): T-tests vs. Zero ---\n\n');
fprintf('Rows are experiments, Columns are models (1=CBPERS, 2=PSL, 3=PSLsim)\n');
disp(p_phi_vs_zero);
fprintf('-------------------- END OF TESTS --------------------\n\n');
end