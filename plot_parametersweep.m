function plot_parametersweep(varargin)
% Generates plots for Figure 5c,d,e, showing parameter recovery results from
% systematic parameter sweeps. It plots how fitted bias and learning rates
% change when one generative parameter is varied at a time.
%
% This function creates a single figure with three rows of plots.
%
% USAGE:
%   plot_parametersweep()       % --- Uses MLE results (Default)
%   plot_parametersweep('MAP')  % --- Uses MAP results
%
% REQUIREMENTS:
%   The file 'parameters_sweep.mat' must be in the current MATLAB path.

%% 1. Data Loading and Preparation
fprintf('--- Generating Parameter Sweep Figure: Loading data... ---\n');

% Set default to MLE and check for optional 'MAP' input.
fit_type = 'MLE';
if nargin > 0 && strcmpi(varargin{1}, 'MAP')
    fit_type = 'MAP';
end

try
    data = load('data/parameters_sweep.mat');
catch ME
    error('Could not load data. Ensure "parameters_sweep.mat" is in the MATLAB path. Details: %s', ME.message);
end
close all;

% Select the correct data based on the chosen fit type
if strcmp(fit_type, 'MAP')
    params_fitted_sweep = data.parameters_sweep_MAP; % Cell array {swept_param}[sweep_val, sim, fitted_param]
    sweep_values = data.swept_MAP;                   % Cell array {swept_param}[sweep_val]
    generative_ref = data.generative_MAP;            % Vector of reference generative params
else % MLE
    params_fitted_sweep = data.parameters_sweep_MLE;
    sweep_values = data.swept_MLE;
    generative_ref = data.generative_MLE;
end

% --- Parameter Definitions ---
swept_param_labels = {'\beta', '\alpha', '\tau', '\phi'};
num_swept_params = length(params_fitted_sweep);

%% 2. Generate Figure
figure('Position', [100, 100, 1200, 700]);
sgtitle(['Parameter Recovery Sweep (' fit_type ')'], 'FontSize', 16, 'FontWeight', 'bold');

% --- Define colors for plots ---
color_cb = [172, 136, 187] / 255; % Purple for confirmation bias
color_ac = [32, 138, 34] / 255;  % Green for alpha_c
color_ad = [192, 0, 0] / 255;    % Red for alpha_d

for col_idx = 1:num_swept_params % Loop over columns (swept parameters)
    
    % --- Extract common data for this column ---
    fitted_params_matrix = params_fitted_sweep{col_idx}; % [sweep_val, sim, fitted_param]
    x_sweep = sweep_values{col_idx};                     % Vector of swept values
    
    %% ROW 1: Absolute Confirmation Bias
    subplot(3, num_swept_params, col_idx);
    hold on;
    data_abs_cb = squeeze(fitted_params_matrix(:,2,:) - fitted_params_matrix(:,3,:));
    plot_sweep_metric(x_sweep, data_abs_cb, color_cb, generative_ref, col_idx, 'y=0', NaN);
    if col_idx == 1, ylabel('Fitted \alpha_c-\alpha_d'); end
    xticklabels({});

    %% ROW 2: Normalized Confirmation Bias
    subplot(3, num_swept_params, col_idx + num_swept_params);
    hold on;
    data_norm_cb = squeeze((fitted_params_matrix(:,2,:) - fitted_params_matrix(:,3,:)) ./ (fitted_params_matrix(:,2,:) + fitted_params_matrix(:,3,:)));
    plot_sweep_metric(x_sweep, data_norm_cb, color_cb, generative_ref, col_idx, 'y=0', NaN);
    if col_idx == 1, ylabel('Fitted Normalized CB'); end
    xticklabels({});

    %% ROW 3: Individual Learning Rates (alpha_c and alpha_d)
    ax_lr = subplot(3, num_swept_params, col_idx + 2*num_swept_params); % Get axes handle
    hold on;
    
    % Store y-data ranges to calculate combined ylim later
    all_y_min = [];
    all_y_max = [];
    
    % --- Plot alpha_c (green) ---
    data_ac = squeeze(fitted_params_matrix(:,2,:));
    [y_min_ac, y_max_ac] = plot_sweep_metric(x_sweep, data_ac, color_ac, generative_ref, col_idx, 'y=x_or_const', 2);
    all_y_min = [all_y_min, y_min_ac];
    all_y_max = [all_y_max, y_max_ac];

    % --- Plot alpha_d (red) on the same axes ---
    data_ad = squeeze(fitted_params_matrix(:,3,:));
    [y_min_ad, y_max_ad] = plot_sweep_metric(x_sweep, data_ad, color_ad, generative_ref, col_idx, 'y=x_or_const', 2);
    all_y_min = [all_y_min, y_min_ad];
    all_y_max = [all_y_max, y_max_ad];
    
    % --- Set YLim based on combined range ---
    final_y_min = min(all_y_min);
    final_y_max = max(all_y_max);
    % Add a small buffer
    y_range = final_y_max - final_y_min;
    ylim(ax_lr, [final_y_min - 0.05*y_range, final_y_max + 0.05*y_range]);
    
    % --- Finalize Subplot ---
    if col_idx == 1, ylabel('Fitted \alpha'); end
    xlabel(['Generative ' swept_param_labels{col_idx}]);
    xlim([min(x_sweep), max(x_sweep)]); % Ensure xlim is set after potential ylim changes reset it

end % End loop over columns

end % End of main function

%% -------------------- Helper Plotting Function --------------------
function [y_min_plot, y_max_plot] = plot_sweep_metric(x_axis, data, color, generative_ref, swept_idx, ground_truth_type, plot_param_idx)
    % Plots a single metric sweep and returns the min/max y-values plotted.
    
    mean_data = squeeze(mean(data, 1));
    sem_data = squeeze(std(data, 0, 1)) ./ sqrt(size(data, 1));
    
    if size(mean_data, 2) > 1; mean_data = mean_data'; end
    if size(sem_data, 2) > 1; sem_data = sem_data'; end
    if size(x_axis, 2) > 1; x_axis = x_axis'; end

    ground_truth_y = NaN(size(x_axis)); % Initialize ground truth vector
    switch ground_truth_type
        case 'y=0'
            ground_truth_y = zeros(size(x_axis));
        case 'y=x_or_const'
            if swept_idx == 2
                ground_truth_y = x_axis;
            else
                ground_truth_y = generative_ref(plot_param_idx) * ones(size(x_axis));
            end
    end
    plot(x_axis, ground_truth_y, 'k-', 'LineWidth', 1.5);

    y_upper = mean_data + sem_data;
    y_lower = mean_data - sem_data;
    y_patch = [y_upper; flipud(y_lower)];
    x_patch = [x_axis; flipud(x_axis)];
    patch(x_patch, y_patch, 'k', 'FaceColor', color, 'FaceAlpha', 0.3, 'LineStyle', 'none');
    plot(x_axis, mean_data, 'LineWidth', 2, 'Color', color);
    
    gener_param_ref_val = generative_ref(swept_idx);
    current_ylim = get(gca, 'YLim'); % Get YLim *before* plotting vertical line
    plot([gener_param_ref_val, gener_param_ref_val], current_ylim, 'k--', 'LineWidth', 1.5);
    
    % Calculate min/max y values including ground truth and SEM bounds
    y_min_plot = min([y_lower; ground_truth_y]);
    y_max_plot = max([y_upper; ground_truth_y]);
    
    xlim([min(x_axis), max(x_axis)]);
end