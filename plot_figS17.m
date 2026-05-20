function plot_figS17()
% Generates Supplementary Figure S17.
%
% This function performs a parameter sweep analysis to characterize which 
% parameters of the Pure Learning Asymmetry (LA) model most strongly drive 
% the emergence of artifactual (spurious) perseveration.
%
% It systematically varies generative parameters of the LA model (beta, 
% alpha_c, alpha_d, and net bias) and plots the resulting fitted 
% perseveration parameter (phi) recovered by the full hybrid model.
%
% REQUIREMENTS:
%   The file 'figureS17_data.mat' must be in the 'data/' path.

%% 1. Data Loading and Preparation
fprintf('--- Generating Figure S17: Loading LA Parameter Sweep data... ---\n');


try
    % Loading the sweep data specifically for perseveration (phi) recovery
    data = load('data/figureS17_data.mat');
catch ME
    error('Could not load data. Ensure "figureS17_data.mat" is in the data folder. Details: %s', ME.message);
end
close all;

params_fitted_sweep = data.parameters_sweep_MLE;
sweep_values = data.swept_MLE;
generative_ref = data.generative_MLE; % Reference values from P2 MLE fits


% --- Parameter Definitions ---
% Represents the parameters of the Pure LA model being swept
swept_param_labels = {'\beta', '\alpha_c', '\alpha_d', '\alpha_c-\alpha_d'};
num_swept_params = length(params_fitted_sweep);

%% 2. Generate Figure
figure('Position', [100, 100, 1200, 400], 'Color', 'w', 'Name', 'Figure S17');

% Color for fitted perseveration (phi) - Light green/sage
color_pers = [151, 189, 161] / 255; 

for col_idx = 1:num_swept_params 
    
    % --- Extract data for current swept parameter ---
    fitted_params_matrix = params_fitted_sweep{col_idx}; % [sweep_val, sim, fitted_param]
    x_sweep = sweep_values{col_idx};                     
    
    % Create subplot for a 1x4 horizontal layout
    subplot(1, num_swept_params, col_idx);
    hold on;
    
    % Extract recovered phi (index 5 in the hybrid model fit)
    % We squeeze to obtain a [sweep_val x simulations] matrix
    data_PERS = squeeze(fitted_params_matrix(:,5,:));
    
    % Plot recovered phi against the varied generative LA parameter
    plot_sweep_metric(x_sweep, data_PERS, color_pers, generative_ref, col_idx);
    
    if col_idx == 1, ylabel('Fitted \phi (Perseveration)'); end
    xlabel(['Generative ' swept_param_labels{col_idx}]);
    
    % Clean up axes
    set(gca, 'Box', 'off', 'TickDir', 'out');
    xlim([min(x_sweep), max(x_sweep)]); 
end

% Add a horizontal legend at the bottom of the figure
legend({'Ground truth (\phi=0)', '', 'Fitted \phi', 'P2 Participant Mean'}, ...
    'Orientation', 'horizontal', 'Location', 'southoutside');

end

%% -------------------- Helper Plotting Function --------------------
function [y_min_plot, y_max_plot] = plot_sweep_metric(x_axis, data, color, generative_ref, swept_idx)
    % Plots mean and SEM of recovered parameters against generative values.
    
    % Calculate summary statistics across simulations (dimension 2)
    mean_data = mean(data, 1);
    sem_data = std(data, 0,1) ./ sqrt(size(data, 1));
    
    % Standardize vector orientations for plotting
    mean_data = mean_data(:);
    sem_data = sem_data(:);
    x_axis = x_axis(:);

    % Ground truth: Since data is generated with a Pure LA model, true phi is 0
    ground_truth_y = zeros(size(x_axis));
    plot(x_axis, ground_truth_y, 'k-', 'LineWidth', 1.5);
    
    % Shaded SEM region (Fitted value uncertainty)
    y_upper = mean_data + sem_data;
    y_lower = mean_data - sem_data;
    patch([x_axis; flipud(x_axis)], [y_upper; flipud(y_lower)], 'k', ...
        'FaceColor', color, 'FaceAlpha', 0.3, 'LineStyle', 'none');
    
    % Plot mean fitted value
    plot(x_axis, mean_data, 'LineWidth', 2, 'Color', color);
    
    % Vertical reference line: Mean parameter value of participants in experiment P2
    gener_param_ref_val = generative_ref(swept_idx);
    current_ylim = ylim;
    plot([gener_param_ref_val, gener_param_ref_val], current_ylim, 'k--', 'LineWidth', 1.5);
    
    % Return min/max for potential subplot synchronization
    y_min_plot = min([y_lower; ground_truth_y]);
    y_max_plot = max([y_upper; ground_truth_y]);
end