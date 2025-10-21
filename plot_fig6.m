function plot_fig6()
% Generates plots for Figure 6, which uses a bootstrap analysis to assess
% whether perseveration alone (in the PSL model) can spuriously generate
% a confirmation bias when fit with the full CBPERS model.
%
% This function creates two figures:
% 1. A plot showing the analysis for each experiment individually.
% 2. A plot showing the analysis aggregated across all experiments.
%
% USAGE:
%   plot_fig6()
%
% REQUIREMENTS:
%   The file 'MLE_fits.mat' must be in the current MATLAB path.

%% 1. Data Loading and Preparation
fprintf('--- Generating Figure 6: Loading data and setting up... ---\n');
try
    data = load('data/MLE_fits.mat');
catch ME
    error('Could not load data. Ensure "MLE_fits.mat" is in the MATLAB path. Details: %s', ME.message);
end
close all;

% Parameters from the CBPERS model fit on REAL data. This is our ground truth.
real_data_params = data.parameters_CBPERS;

% Parameters from the CBPERS model fit on data SIMULATED from the PSL model.
% This will be used to generate the null distribution.
sim_data_params = data.parameters_PSLsim_CBPERSfit;

% --- Configuration ---
bootstrap_samples = 10000;
experiments = {'L1','L2','P1','P2','C1','C2','C3','C4','S1a','S1b'};
colors = struct(...
    'sim_dist', [172, 136, 187]/255, ... % Purple for simulation distribution
    'real_sig', [181, 101, 118]/255, ... % Reddish for significant real data point
    'real_ns',  [1, 1, 1] ...             % White for non-significant
);
p_empirical_by_exp = NaN(length(experiments), 2);

%% 2. By-Experiment Bootstrap Analysis and Plotting
figure('Position', [100, 100, 800, 800]);
sgtitle('Figure 6: Bootstrap Analysis by Experiment', 'FontSize', 16, 'FontWeight', 'bold');
% Loop over the two metrics: 1 = absolute CB, 2 = normalized CB
for metric_idx = 1:2
    subplot(2, 1, metric_idx);
    hold on;
    % Loop backwards to ensure y-axis labels match the plot order
    for exp_idx = length(experiments):-1:1
        
        % --- Bootstrap the simulation data to create a null distribution ---
        sim_params_exp = sim_data_params{exp_idx};
        n_subjects = size(sim_params_exp, 1);
        n_sims_per_subject = size(sim_params_exp, 3);
        
        bootstrapped_means = NaN(bootstrap_samples, 1);
        for i = 1:bootstrap_samples
            % For each bootstrap sample, resample subjects with replacement
            resampled_subjects_idx = randsample(n_subjects, n_subjects, true);
            
            % For each resampled subject, also pick a random simulation run
            resampled_sims_idx = randsample(n_sims_per_subject, n_subjects, true);
            
            % Create the bootstrap parameter matrix
            bootstrap_param_sample = NaN(n_subjects, size(sim_params_exp, 2));
            for j = 1:n_subjects
                bootstrap_param_sample(j,:) = sim_params_exp(resampled_subjects_idx(j), :, resampled_sims_idx(j));
            end
            % Calculate the mean confirmation bias for this bootstrap sample
            mean_params = mean(bootstrap_param_sample, 1);
            if metric_idx == 1
                bootstrapped_means(i) = mean_params(2) - mean_params(3);
            else
                bootstrapped_means(i) = (mean_params(2) - mean_params(3)) / (mean_params(2) + mean_params(3));
            end
        end
        % --- Calculate the metric for the real data ---
        real_params_exp = real_data_params{exp_idx};
        if metric_idx == 1
            real_data_metric = real_params_exp(:, 2) - real_params_exp(:, 3);
        else
            real_data_metric = (real_params_exp(:, 2) - real_params_exp(:, 3)) ./ (real_params_exp(:, 2) + real_params_exp(:, 3));
        end
        mean_real_metric = mean(real_data_metric);
        % --- Calculate empirical p-value ---
        % This is the proportion of bootstrapped means (from PSL sims) that are
        % greater than the mean observed in the real data.
        p_empirical_by_exp(exp_idx, metric_idx) = mean(bootstrapped_means > mean_real_metric);
        % --- Plotting for the current experiment ---
        y_pos = exp_idx * 1.5; % This is now the X position
        
        % Plot the bootstrapped null distribution as a horizontal histogram
        [counts, bins] = hist(bootstrapped_means, 20);
        barh(bins, y_pos + counts / (1.2 * max(counts)), 'FaceColor', colors.sim_dist, 'FaceAlpha', 0.5);
        rectangle('Position', [-1, min(bins)-1, y_pos + 1, max(bins)-min(bins)+2], 'FaceColor', 'w', 'EdgeColor', 'w');
        % Plot the mean from the real data as a circle
        % Color it red if significantly > 0, otherwise white.
        [~, p_ttest] = ttest(real_data_metric, 0, 'Tail', 'right');
        marker_face_color = colors.real_ns;
        if p_ttest < 0.05
            marker_face_color = colors.real_sig;
        end
        plot(y_pos, mean_real_metric, 'ok', 'MarkerSize', 8, 'MarkerFaceColor', marker_face_color, 'LineWidth', 1.2);
    end
    
    % --- Finalize Subplot Aesthetics ---
    plot([0.5, 1.5 * length(experiments) + 1], [0, 0], 'k--'); % Zero line
    xticks(1.5:1.5:(1.5 * length(experiments)));
    xticklabels(fliplr(experiments)); % Flipped to match plotting order
    xlim([0.5, 1.5 * length(experiments) + 1]);
    
    if metric_idx == 1
        ylabel('\alpha_c - \alpha_d');
        ylim([-0.25, 0.5]);
    else
        ylabel('Normalized CB');
        ylim([-0.5, .8]);
    end
end
xlabel('Experiment');

%% 3. Aggregated Bootstrap Analysis and Plotting
figure('Position', [950, 100, 400, 800]);
sgtitle('Figure 6: Aggregated Analysis', 'FontSize', 16, 'FontWeight', 'bold');
p_empirical_agg = NaN(1, 2);

for metric_idx = 1:2
    subplot(2, 1, metric_idx);
    hold on;

    % --- Calculate the grand mean metric for real data ---
    real_means_per_exp = NaN(1, length(experiments));
    for exp_idx = 1:length(experiments)
        real_params_exp = real_data_params{exp_idx};
        if metric_idx == 1
            real_means_per_exp(exp_idx) = mean(real_params_exp(:, 2) - real_params_exp(:, 3));
        else
            real_means_per_exp(exp_idx) = mean((real_params_exp(:, 2) - real_params_exp(:, 3)) ./ (real_params_exp(:, 2) + real_params_exp(:, 3)));
        end
    end
    grand_mean_real = mean(real_means_per_exp);

    % --- Bootstrap the aggregated simulation data ---
    bootstrapped_grand_means = NaN(bootstrap_samples, 1);
    for i = 1:bootstrap_samples
        mean_across_exps = NaN(1, length(experiments));
        for exp_idx = 1:length(experiments)
            sim_params_exp = sim_data_params{exp_idx};
            n_subjects = size(sim_params_exp, 1);
            n_sims_per_subject = size(sim_params_exp, 3);
            
            resampled_subjects_idx = randsample(n_subjects, n_subjects, true);
            resampled_sims_idx = randsample(n_sims_per_subject, n_subjects, true);

            bootstrap_param_sample = NaN(n_subjects, size(sim_params_exp, 2));
            for j = 1:n_subjects
                 bootstrap_param_sample(j, :) = sim_params_exp(resampled_subjects_idx(j), :, resampled_sims_idx(j));
            end
            
            mean_params = mean(bootstrap_param_sample, 1);
            if metric_idx == 1
                mean_across_exps(exp_idx) = mean_params(2) - mean_params(3);
            else
                mean_across_exps(exp_idx) = (mean_params(2) - mean_params(3)) / (mean_params(2) + mean_params(3));
            end
        end
        bootstrapped_grand_means(i) = mean(mean_across_exps);
    end

    % --- Calculate empirical p-value for the aggregated data ---
    p_empirical_agg(metric_idx) = mean(bootstrapped_grand_means > grand_mean_real);

    % --- Plotting (axes swapped to match Part 2) ---
    x_pos = 1.5;
    
    [counts, bins] = hist(bootstrapped_grand_means, 20);
    barh(bins, x_pos + counts / (1.2 * max(counts)), 'FaceColor', colors.sim_dist, 'FaceAlpha', 0.5);
    rectangle('Position', [-1, min(bins)-1, x_pos + 1, max(bins)-min(bins)+2], 'FaceColor', 'w', 'EdgeColor', 'w');

    [~, p_ttest] = ttest(real_means_per_exp, 0, 'Tail', 'right');
    marker_face_color = colors.real_ns;
    if p_ttest < 0.05
        marker_face_color = colors.real_sig;
    end
    plot(x_pos, grand_mean_real, 'ok', 'MarkerSize', 8, 'MarkerFaceColor', marker_face_color, 'LineWidth', 1.2);

    % --- Finalize Subplot Aesthetics ---
    plot([0.5, 3], [0, 0], 'k--');
    xticks(1.5); xticklabels('All Experiments');
    xlim([0.5, 3]);
    
    if metric_idx == 1
        ylabel('\alpha_c - \alpha_d');
        ylim([-0.1, 0.2]);
    else
        ylabel('Normalized CB');
        ylim([-0.2, 0.35]);
    end
end

%% 4. Statistical Analyses Output
fprintf('\n======================================================\n');
fprintf('STATISTICAL TESTS FOR FIGURE 6\n');
fprintf('Empirical p-value = Proportion of bootstrapped means (from PSL sims)\n');
fprintf('that are greater than the mean observed in the real data.\n');
fprintf('======================================================\n\n');
fprintf('--- By-Experiment Empirical P-Values ---\n');
fprintf('%-12s | %-20s | %-20s\n', 'Experiment', 'Absolute CB (p)', 'Normalized CB (p)');
fprintf('-----------------------------------------------------\n');
for i = 1:length(experiments)
    fprintf('%-12s | p = %-18.4f | p = %-18.4f\n', experiments{i}, p_empirical_by_exp(i,1), p_empirical_by_exp(i,2));
end
fprintf('\n--- Aggregated Empirical P-Values ---\n');
fprintf('  Absolute CB:   p = %.4f\n', p_empirical_agg(1));
fprintf('  Normalized CB: p = %.4f\n', p_empirical_agg(2));
fprintf('\n-------------------- END OF TESTS --------------------\n\n');

end