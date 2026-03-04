function plot_figS20() 
% PLOT_S20 Generates supplementary figure S20 (MLE only).
%
% This function shows the distribution of generative negative phi values 
% and the resulting spurious fitted confirmation bias (alpha_c - alpha_d) 
% when simulating data with negative perseveration and fitting with the 
% CBPERS model using MLE.
%
% REQUIREMENTS:
%   The file 'figure5_data.mat' must be in the 'Data' subdirectory
%   relative to this function's location.

%% 1. Data Loading and Preparation
data_path = 'Data'; 
filename = 'figureS20_data.mat';
full_filepath = fullfile(data_path, filename); 

fprintf('--- Generating Figure S20: Loading MLE data from %s ---\n', full_filepath);

try
    data = load(full_filepath); 
catch ME
    error('Could not load data. Ensure "%s" is in the "%s" subdirectory. Details: %s', filename, data_path, ME.message);
end
close all;

% Exclusively select MLE data
generative_phi = data.phi_MLE;
fitted_params = data.parameters_negphi_MLE;

% --- Calculate Fitted Bias Metrics ---
if iscell(fitted_params)
    fitted_params = vertcat(fitted_params{:});
end

% Extract alpha_c (col 2) and alpha_d (col 3)
fitted_cb_abs = fitted_params(:, 2) - fitted_params(:, 3);
fitted_cb_norm = (fitted_params(:, 2) - fitted_params(:, 3)) ./ (fitted_params(:, 2) + fitted_params(:, 3));

% Perform t-tests vs zero for fitted bias (ground truth is zero)
[~, p_abs, ~, stats_abs] = ttest(fitted_cb_abs);
[~, p_norm, ~, stats_norm] = ttest(fitted_cb_norm);

%% 2. Plotting Configuration
colors = {[107, 142, 185]/255, [172, 136, 187]/255}; % Blue for phi, Purple for CB
plot_data = {generative_phi, fitted_cb_abs, fitted_cb_norm};
plot_labels = {'Generative \phi', 'Fitted \alpha_c-\alpha_d', 'Fitted Normalized CB'};
plot_colors = {colors{1}, colors{2}, colors{2}};

%% 3. Generate Figure
figure('Position', [100, 100, 900, 400], 'Color', 'w', 'Name', 'Figure S20 (MLE)');
sgtitle('Figure S20: Negative Perseveration Analysis (MLE)', 'FontSize', 14, 'FontWeight', 'bold');

for i = 1:3
    subplot(1, 3, i);
    hold on;
    
    current_data = plot_data{i};
    current_color = plot_colors{i};
    
    % Horizontal histogram (distribution)
    [counts, bins] = hist(current_data, 15);
    barh(bins, 1 + counts / (2 * max(counts)), 'FaceColor', current_color, 'FaceAlpha', 0.5);
    
    % Scaling and bounds
    min_val = min(current_data) - 0.1 * range(current_data);
    max_val = max(current_data) + 0.1 * range(current_data);
    if isnan(min_val) || isnan(max_val) || isinf(min_val) || isinf(max_val) || range(current_data)==0
        min_val = -1; max_val = 1; 
    end
    
    % White-out mask for the bar chart background
    rectangle('Position', [-1, min_val, 2.01, max_val - min_val], 'FaceColor', [1 1 1], 'EdgeColor', [1 1 1]);
    
    % Zero reference line
    plot([-0.5, 2], [0 0], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
    
    % Summary Statistics (Mean + SEM)
    mean_val = mean(current_data);
    sem_val = std(current_data) / sqrt(length(current_data));
    errorbar(0.6, mean_val, sem_val, '.k', 'CapSize', 0, 'LineWidth', 1.3);
    plot(0.6, mean_val, 'ok', 'MarkerSize', 8, 'MarkerFaceColor', current_color, 'LineWidth', 1.2);
    
    % Aesthetics
    xticks(0.8);
    xticklabels('');
    ylabel(plot_labels{i});
    xlim([0, 1.9]);
    set(gca, 'Box', 'off', 'TickDir', 'out');
    
    % Dynamic Y-Axis Limits
    if i == 1 
         ylim_dynamic = [min_val, max(0.1, max_val)];
    else 
        ylim_dynamic = [min(-0.1, min_val), max(0.1, max_val)];
    end
    
    if diff(ylim_dynamic) < 1e-6
        ylim_dynamic = ylim_dynamic + [-0.1, 0.1];
    end
    ylim(ylim_dynamic);
end

%% 4. Statistical Analyses Output
fprintf('\n==================================================================\n');
fprintf('STATISTICAL TESTS FOR FIGURE S20 (MLE Estimation)\n');
fprintf('Comparing fitted bias metrics to 0 (Ground Truth)\n');
fprintf('==================================================================\n\n');
fprintf('--- Fitted Absolute CB vs. Zero ---\n');
fprintf('  Mean = %.4f, SEM = %.4f, t = %.3f, p = %.4f\n', ...
        mean(fitted_cb_abs), std(fitted_cb_abs)/sqrt(length(fitted_cb_abs)), stats_abs.tstat, p_abs);
fprintf('\n--- Fitted Normalized CB vs. Zero ---\n');
fprintf('  Mean = %.4f, SEM = %.4f, t = %.3f, p = %.4f\n', ...
        mean(fitted_cb_norm), std(fitted_cb_norm)/sqrt(length(fitted_cb_norm)), stats_norm.tstat, p_norm);
fprintf('\n-------------------------- END OF TESTS --------------------------\n\n');

end