function plot_CBsweep(varargin)
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
    data = load('CBsweep_noPERS.mat');
catch ME
    error('Could not load data. Ensure "parameters_sweep.mat" is in the MATLAB path. Details: %s', ME.message);
end
close all;

% Select the correct data based on the chosen fit type
if strcmp(fit_type, 'MAP')
    params_fitted_sweep = data.parameters_CBPERsim_CBPERSfit_MAP; % Cell array {swept_param}[sweep_val, sim, fitted_param]
    sweep_values = data.cb;                   % Cell array {swept_param}[sweep_val]
    generative_ref = data.cb;            % Vector of reference generative params
else % MLE
    params_fitted_sweep = data.parameters_CBPERsim_CBPERSfit_MLE;
    sweep_values = data.cb;
    generative_ref = data.cb;
end

% --- Parameter Definitions ---
swept_param_labels = {'\alpha_c-\alpha_d'};
num_swept_params = length(params_fitted_sweep);

%% 2. Generate Figure
figure('Position', [100, 100, 1200, 700]);
sgtitle(['Confirmation Bias Recovery Sweep (' fit_type ')'], 'FontSize', 16, 'FontWeight', 'bold');

% --- Define colors for plots ---
color_cb = [172, 136, 187] / 255; % Purple for confirmation bias
color_ac = [32, 138, 34] / 255;  % Green for alpha_c
color_ad = [192, 0, 0] / 255;    % Red for alpha_d


    
% --- Extract common data for this column ---
fitted_params_matrix = params_fitted_sweep; % [sweep_val, sim, fitted_param]
x_sweep = sweep_values;                     % Vector of swept values

%% ROW 1: Absolute Confirmation Bias
subplot(1,3,1);
hold on;
data_abs_cb = fitted_params_matrix(:,:,2) - fitted_params_matrix(:,:,3);
plot_sweep_metric(x_sweep, data_abs_cb, color_cb, generative_ref);
ylabel('Fitted \alpha_c-\alpha_d')


% %% ROW 2: Normalized Confirmation Bias
% subplot(1,3);
% hold on;
% data_norm_cb = squeeze((fitted_params_matrix(:,2,:) - fitted_params_matrix(:,3,:)) ./ (fitted_params_matrix(:,2,:) + fitted_params_matrix(:,3,:)));
% plot_sweep_metric(x_sweep, data_norm_cb, color_cb, generative_ref, col_idx, 'y=0', NaN);
% if col_idx == 1, ylabel('Fitted Normalized CB'); end
% xticklabels({});

%% ROW 3: Individual Learning Rates (alpha_c and alpha_d)
subplot(1,3,2); % Get axes handle
data_abs_cb = fitted_params_matrix(:,:,2);
plot_sweep_metric(x_sweep, data_abs_cb, color_ac, 0.4+generative_ref/2);
ylabel('Fitted \alpha_c')

xlabel('Generative \alpha_c - \alpha_d')

subplot(1,3,3); % Get axes handle
data_abs_cb = fitted_params_matrix(:,:,3);
plot_sweep_metric(x_sweep, data_abs_cb, color_ad, 0.4-generative_ref/2);
ylabel('Fitted \alpha_d')





end % End of main function

%% -------------------- Helper Plotting Function --------------------
function [y_min_plot, y_max_plot] = plot_sweep_metric(x_axis, data, color, ground_truth)
    % Plots a single metric sweep and returns the min/max y-values plotted.
    
    mean_data = squeeze(mean(data, 2));
    sem_data = squeeze(std(data, 0, 2)) ./ sqrt(size(data, 2));
    
    if size(mean_data, 2) > 1; mean_data = mean_data'; end
    if size(sem_data, 2) > 1; sem_data = sem_data'; end
    if size(x_axis, 2) > 1; x_axis = x_axis'; end

    ground_truth_y = NaN(size(x_axis)); % Initialize ground truth vector


    ground_truth_y = ground_truth';
      
    plot(x_axis, ground_truth_y, 'k-', 'LineWidth', 1.5);
    hold on

    y_upper = mean_data + sem_data;
    y_lower = mean_data - sem_data;
    y_patch = [y_upper; flipud(y_lower)];
    x_patch = [x_axis; flipud(x_axis)];
    patch(x_patch, y_patch, 'k', 'FaceColor', color, 'FaceAlpha', 0.3, 'LineStyle', 'none');
    hold on
    plot(x_axis, mean_data, 'LineWidth', 2, 'Color', color);
    hold on
    % 

    
    % Calculate min/max y values including ground truth and SEM bounds
    y_min_plot = min([y_lower; ground_truth_y]);
    y_max_plot = max([y_upper; ground_truth_y]);
    ylim([y_min_plot,y_max_plot])
    xlim([min(x_axis), max(x_axis)]);


    gener_param_ref_val = 0;
    current_ylim = get(gca, 'YLim'); % Get YLim *before* plotting vertical line
    current_xlim = get(gca, 'XLim'); % Get YLim *before* plotting vertical line
    plot([0,0], current_ylim, '--','Color',[0.4 0.4 0.4], 'LineWidth', 1.5);
    % hold on
    % plot( current_xlim, [0,0], '--','Color',[0.4 0.4 0.4], 'LineWidth', 1.5);
end