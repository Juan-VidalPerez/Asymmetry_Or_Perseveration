function plot_fig4(varargin)
% Generates plots for Figure 4 and supplementary figures S5 & S11, showing
% parameter recovery as a function of the number of training sessions.
%
% This function creates two figures:
% 1. Confirmation Bias parameters (alpha_c - alpha_d, alpha_c, alpha_d).
% 2. Perseveration parameters (tau, phi).
%
% USAGE:
%   plot_fig4()         % --- For Figure 4 & S5 (uses MLE fits)
%   plot_fig4('MAP')    % --- For Figure S11 (uses MAP fits)
%
% REQUIREMENTS:
%   The file 'figure4_fits.mat' must be in the current MATLAB path.

%% -------------------- Data Loading and Preparation --------------------
fprintf('--- Generating Figure 4/S5/S11: Loading data... ---\n');

% Set default to MAP and check for optional 'MLE' input.
fit_type = 'MLE';
if nargin > 0 && strcmpi(varargin{1}, 'MAP')
    fit_type = 'MAP';
end

try
    data = load('data/figure4_fits.mat');
catch ME
    error('Could not load data. Ensure "figure4_fits.mat" is in the MATLAB path. Details: %s', ME.message);
end
close all;

% Select the correct data based on the chosen fit type
if strcmp(fit_type, 'MAP')
    params_fitted = data.parameters_PSLsim_CBPERSfit_MAP;
    params_generating = data.parameters_PSL_MAP;
    session_lengths = data.session_length_MAP;
else
    params_fitted = data.parameters_PSLsim_CBPERSfit_MLE;
    params_generating = data.parameters_PSL_MLE;
    session_lengths = data.session_length_MLE;
end

% The generating parameters are the ground truth (asymptotes).
% We take the mean across simulations to get a single value for each parameter.
generating_asymptotes = [0 params_generating(2) params_generating(2) params_generating(3) params_generating(4)];
x_axis = log2(session_lengths);

%% -------------------- Plotting Configuration --------------------
% Define configurations for each parameter to be plotted
plot_configs = {
  % {ylabel, data_func, asymptote_func, fig_num, subplot_pos, color}
    {'Fitted \alpha_c - \alpha_d', @(p) p(:,:,2) - p(:,:,3), @(a) a(2) - a(3), 1, 1, [172, 136, 187]/255}, ...
    {'Fitted \alpha_c',          @(p) p(:,:,2),             @(a) a(2),          1, 2, [32, 138, 34]/255}, ...
    {'Fitted \alpha_d',          @(p) p(:,:,3),             @(a) a(3),          1, 2, [192, 0, 0]/255}, ...
    {'Fitted \tau',              @(p) p(:,:,end-1),         @(a) a(end-1),      2, 1, [172, 136, 187]/255}, ...
    {'Fitted \phi',              @(p) p(:,:,end),           @(a) a(end),        2, 2, [172, 136, 187]/255}
};

% Pre-allocate for statistical results
p_values = NaN(size(plot_configs, 1), length(x_axis));
t_stats = NaN(size(plot_configs, 1), length(x_axis));
mean_diffs = NaN(size(plot_configs, 1), length(x_axis));

%% -------------------- Figure Generation --------------------
% Create figures
fig1 = figure('Position', [100, 100, 600, 300]);
fig2 = figure('Position', [750, 100, 600, 300]);


% Loop through each parameter configuration and create the plots
for mm = 1:size(plot_configs, 2)
    % Select the correct figure and subplot
    figure(plot_configs{mm}{4});
    subplot(1,2, plot_configs{mm}{5});
    hold on;

    % Extract data and asymptote using the predefined functions
    data = plot_configs{mm}{2}(params_fitted);
    asymptote = plot_configs{mm}{3}(generating_asymptotes);
    
    % --- Plotting ---
    mean_data = mean(data, 2);
    sem = std(data, 0, 2) ./ sqrt(size(data, 2));
    
    % Create shaded area for SEM
    y_patch = [mean_data + sem; flipud(mean_data - sem)];
    patch([x_axis'; flipud(x_axis')], y_patch, 'k', 'FaceColor', plot_configs{mm}{6}, 'FaceAlpha', 0.3, 'LineStyle', 'none');
    
    % Plot mean line
    plot(x_axis, mean_data, 'LineWidth', 2, 'Color', plot_configs{mm}{6});
    
    % Plot ground truth asymptote line
    plot(x_axis, asymptote * ones(1, length(x_axis)), 'k--', 'LineWidth', 1.5);
    
    % --- Aesthetics and Labels ---
    ylabel(plot_configs{mm}{1});
    xlim([min(x_axis), max(x_axis)]);
    xticks(x_axis);
    
    % Configure x-tick labels
    if strcmp(fit_type, 'MAP')
        xticklabels(arrayfun(@(x) num2str(2^x), x_axis, 'UniformOutput', false));
    else
        xticks([0:6])
        xticklabels(arrayfun(@(x) num2str(2^x), [0:6], 'UniformOutput', false));
    end
    xlabel('Number of training sessions');


    % Add legends
    if ismember(mm,[2,3])
        legend({'', '\alpha_c', '', '', '\alpha_d'}, 'Location', 'best');
    else
        legend({'', 'Fitted', 'Generative'}, 'Orientation', 'horizontal', 'Location', 'best');
    end

    % --- Statistical Tests ---
    for ii = 1:size(data, 1)
        [~, p, ~, stats] = ttest(data(ii, :) - asymptote);
        p_values(mm, ii) = p;
        t_stats(mm, ii) = stats.tstat;
        mean_diffs(mm, ii) = mean(data(ii, :)) - asymptote;
    end
end

%% -------------------- Statistical Analyses Output --------------------
fprintf('\n==================================================================\n');
fprintf('STATISTICAL TESTS FOR FIGURE 4/S5/S11 (Fit Type: %s)\n', fit_type);
fprintf('T-tests comparing fitted parameter to generating parameter (asymptote)\n');
fprintf('==================================================================\n\n');

for mm = 1:size(plot_configs, 2)
    fprintf('--- Parameter: %s ---\n', plot_configs{mm}{1});
    for ii = 1:length(x_axis)
        fprintf('  Sessions: %-4d | Mean Diff: %+.4f, t = %+.3f, p = %.4f\n', ...
                2^x_axis(ii), mean_diffs(mm, ii), t_stats(mm, ii), p_values(mm, ii));
    end
    fprintf('\n');
end
fprintf('-------------------------- END OF TESTS --------------------------\n\n');

end