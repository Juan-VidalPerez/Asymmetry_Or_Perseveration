function plot_fig5(varargin) 
% Generates plots for Figure 5c, 5d, and 5e, showing the distribution of
% generative negative phi values and the resulting fitted confirmation bias
% when simulating data with negative perseveration and fitting with CBPERS.
%
% USAGE:
%   plot_fig5()       % --- Uses MLE results (Default)
%   plot_fig5('MAP')  % --- Uses MAP results
%
% REQUIREMENTS:
%   The file 'figure5_data.mat' must be in the 'Data' subdirectory
%   relative to this function's location.

%% 1. Data Loading and Preparation
data_path = 'Data'; % Define the data subdirectory
filename = 'figure5_data.mat';
full_filepath = fullfile(data_path, filename); % Construct the full path
fprintf('--- Generating Figure 5c,d,e: Loading data from %s ---\n', full_filepath);

% Set default to MLE and check for optional 'MAP' input.
fit_type = 'MLE'; % Default changed to MLE
if nargin > 0 && strcmpi(varargin{1}, 'MAP')
    fit_type = 'MAP';
end

try
    data = load(full_filepath); % Load data from the Data subdirectory
catch ME
    error('Could not load data. Ensure "%s" is in the "%s" subdirectory. Details: %s', filename, data_path, ME.message);
end
close all;

% Select the correct data based on the chosen fit type
if strcmp(fit_type, 'MAP')
    generative_phi = data.phi_MAP;
    fitted_params = data.parameters_negphi_MAP;
else % MLE
    generative_phi = data.phi_MLE;
    fitted_params = data.parameters_negphi_MLE;
end

% --- Calculate Fitted Bias Metrics ---
if iscell(fitted_params)
    fitted_params = vertcat(fitted_params{:});
end
fitted_cb_abs = fitted_params(:, 2) - fitted_params(:, 3);
fitted_cb_norm = (fitted_params(:, 2) - fitted_params(:, 3)) ./ (fitted_params(:, 2) + fitted_params(:, 3));

% Perform t-tests vs zero for fitted bias
[~, p_abs, ~, stats_abs] = ttest(fitted_cb_abs);
[~, p_norm, ~, stats_norm] = ttest(fitted_cb_norm);

%% 2. Plotting Configuration
colors = {[107, 142, 185]/255, [172, 136, 187]/255}; % Blue for phi, Purple for CB
plot_data = {generative_phi, fitted_cb_abs, fitted_cb_norm};
plot_labels = {'Generative \phi', 'Fitted \alpha_c-\alpha_d', 'Fitted Normalized CB'};
plot_colors = {colors{1}, colors{2}, colors{2}};

%% 3. Generate Figure
figure('Position', [100, 100, 900, 400]);
sgtitle(['Figure 5c,d,e (' fit_type ')'], 'FontSize', 14, 'FontWeight', 'bold');
for i = 1:3
    subplot(1, 3, i);
    hold on;
    
    current_data = plot_data{i};
    current_color = plot_colors{i};
    
    [counts, bins] = hist(current_data, 15);
    barh(bins, 1 + counts / (2 * max(counts)), 'FaceColor', current_color, 'FaceAlpha', 0.5);
    
    min_val = min(current_data) - 0.1 * range(current_data);
    max_val = max(current_data) + 0.1 * range(current_data);
    if isnan(min_val) || isnan(max_val) || isinf(min_val) || isinf(max_val) || range(current_data)==0
        min_val = -1; max_val = 1; % Handle edge cases gracefully
    end
    rectangle('Position', [-1, min_val, 2.01, max_val - min_val], 'FaceColor', [1 1 1], 'EdgeColor', [1 1 1]);
    
    plot([-0.5, 2], [0 0], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
    
    mean_val = mean(current_data);
    sem_val = std(current_data) / sqrt(length(current_data));
    errorbar(0.6, mean_val, sem_val, '.k', 'CapSize', 0, 'LineWidth', 1.3);
    plot(0.6, mean_val, 'ok', 'MarkerSize', 8, 'MarkerFaceColor', current_color, 'LineWidth', 1.2);
    
    xticks(0.8);
    xticklabels('');
    ylabel(plot_labels{i});
    xlim([0, 1.9]);
    
    if i == 1 % Generative phi
         ylim_dynamic = [min_val, max(0.1, max_val)];
    else % Fitted CB metrics
        ylim_dynamic = [min(-0.1, min_val), max(0.1, max_val)];
    end
    % Ensure ylim doesn't become degenerate if min/max are too close
    if diff(ylim_dynamic) < 1e-6
        ylim_dynamic = ylim_dynamic + [-0.1, 0.1];
    end
    ylim(ylim_dynamic);
         
end
%% 4. Statistical Analyses Output
fprintf('\n==================================================================\n');
fprintf('STATISTICAL TESTS FOR FIGURE 5d,e (Fit Type: %s)\n', fit_type);
fprintf('T-tests comparing fitted bias metrics to zero\n');
fprintf('==================================================================\n\n');
fprintf('--- Fitted Absolute CB vs. Zero ---\n');
fprintf('  Mean = %.4f, SEM = %.4f, t = %.3f, p = %.4f\n', ...
        mean(fitted_cb_abs), std(fitted_cb_abs)/sqrt(length(fitted_cb_abs)), stats_abs.tstat, p_abs);
fprintf('\n--- Fitted Normalized CB vs. Zero ---\n');
fprintf('  Mean = %.4f, SEM = %.4f, t = %.3f, p = %.4f\n', ...
        mean(fitted_cb_norm), std(fitted_cb_norm)/sqrt(length(fitted_cb_norm)), stats_norm.tstat, p_norm);
fprintf('\n-------------------------- END OF TESTS --------------------------\n\n');
end