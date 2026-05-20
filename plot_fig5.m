function plot_fig5()
% PLOT_FIG5 Generates plots for Figure 5 and Supplementary Figure S18.
% This function uses a bootstrap analysis to assess whether perseveration 
% alone (PSL model) can spuriously generate a confirmation bias.
%
%   - Figure 5: Absolute CB (Metric 1: alpha_c - alpha_d)
%   - Figure S20: Normalized CB (Metric 2: (ac-ad)/(ac+ad))
%
% REQUIREMENTS:
%   The file 'figure5_data.mat' must be in the 'data/' subdirectory.

%% 1. Data Loading and Preparation
fprintf('--- Generating Figure 5 & S18: Loading data... ---\n');
try
    % Path updated to data/ subdirectory as requested
    data = load('data/figure5_data.mat');
catch ME
    error('Could not load data. Ensure "figure5_data.mat" is in the data/ folder. Details: %s', ME.message);
end
close all;

% real_data_params: CBPERS model fit on real experimental data.
real_data_params = data.parameters_CBPERS;
% sim_data_params: CBPERS model fit on data simulated from the PSL model (Null distribution).
sim_data_params = data.parameters_CBPERSsim_CBPERSfit;

% --- Configuration ---
bootstrap_samples = 10000;
experiments = {'L1','L2','P1','P2','C1','C2','C3','C4','S1a','S1b'};
colors = struct(...
    'sim_dist', [181, 101, 118]/255, ... % Color for the PSL-simulated null distribution
    'real_sig', [181, 101, 118]/255, ... % Color for real data markers if empirical p < 0.05
    'real_ns',  [1 1 1] ...             % White for non-significant markers
);

p_empirical_by_exp = NaN(length(experiments), 2);

%% 2. By-Experiment Bootstrap Analysis and Plotting
% Generates a single figure containing both the absolute and normalized metrics.
figure('Position', [100, 100, 800, 800], 'Name', 'Figure 5 & S18: Experiment-wise');
sgtitle('Hybrid Bootstrap Analysis', 'FontSize', 16, 'FontWeight', 'bold');

% metric_idx 1 = Fig 5 (Absolute); metric_idx 2 = S18 (Normalized)
for metric_idx = 1:2
    subplot(2, 1, metric_idx);
    hold on;
    
    % Loop through experiments in reverse to maintain top-to-bottom visual order
    for exp_idx = length(experiments):-1:1
        
        % --- Bootstrap simulation data ---
        sim_params_exp = sim_data_params{exp_idx};
        n_subjects = size(sim_params_exp, 1);
        n_sims_per_subject = size(sim_params_exp, 4);
        
        bootstrapped_means = NaN(bootstrap_samples, 1);
        for i = 1:bootstrap_samples
            % Resample subjects with replacement and select a random simulation run
            resampled_sims_idx = randsample(n_sims_per_subject, n_subjects, true);
            resampled_anti_idx = randsample(2, n_subjects, true);
            
            bootstrap_param_sample = NaN(n_subjects, size(sim_params_exp, 2));
            for j = 1:n_subjects
                bootstrap_param_sample(j,:) = sim_params_exp(j,:,resampled_anti_idx(j), resampled_sims_idx(j));
            end
            
            % Compute mean for this bootstrap iteration
            mean_params = mean(bootstrap_param_sample, 1);
            if metric_idx == 1
                bootstrapped_means(i) = mean_params(2) - mean_params(3);
            else
                bootstrapped_means(i) = (mean_params(2) - mean_params(3)) / (mean_params(2) + mean_params(3));
            end
        end
        
        % --- Real data calculation ---
        real_params_exp = real_data_params{exp_idx};
        if metric_idx == 1
            real_data_metric = real_params_exp(:, 2) - real_params_exp(:, 3);
        else
            real_data_metric = (real_params_exp(:, 2) - real_params_exp(:, 3)) ./ (real_params_exp(:, 2) + real_params_exp(:, 3));
        end
        mean_real_metric = mean(real_data_metric);
        
        % --- Empirical p-value ---
        % Proportions of simulation-based means larger than the observed real mean.
        p_empirical_by_exp(exp_idx, metric_idx) = mean(bootstrapped_means > mean_real_metric);
        
        % --- Visual Plotting ---
        y_pos = exp_idx * 1.5; 
        
        % Plot simulated distribution as horizontal bar histogram
        [counts, bins] = hist(bootstrapped_means, 20);
        barh(bins, y_pos + counts / (1.2 * max(counts)), 'FaceColor', colors.sim_dist, 'FaceAlpha', 0.5);
        
        % Mask background for overlapping clarity
        rectangle('Position', [-1, min(bins)-1, y_pos + 1, max(bins)-min(bins)+2], 'FaceColor', 'w', 'EdgeColor', 'w');
        
        % Plot the observed real data mean
        marker_face_color = colors.real_ns;
        if p_empirical_by_exp(exp_idx, metric_idx) < 0.05
            marker_face_color = colors.real_sig;
        end
        plot(y_pos, mean_real_metric, 'ok', 'MarkerSize', 8, 'MarkerFaceColor', marker_face_color, 'LineWidth', 1.2);
    end
    
    % Final aesthetics for experiment plots
    plot([0.5, 1.5 * length(experiments) + 1], [0, 0], 'k--');
    xticks(1.5:1.5:(1.5 * length(experiments)));
    xticklabels(experiments);
    xlim([0.5, 1.5 * length(experiments) + 1]);
    
    if metric_idx == 1
        ylabel('Absolute Bias (\alpha_c - \alpha_d)'); % Figure 5
        ylim([-0.25, 0.5]);
    else
        ylabel('Normalized Bias'); % S18
        ylim([-0.5, .8]);
    end
end
xlabel('Experiment');

%% 3. Aggregated Bootstrap Analysis and Plotting
% Generates a second figure showing the data pooled across all experiments.
figure('Position', [950, 100, 400, 800], 'Name', 'Aggregated Bootstrap');
sgtitle('Aggregated Analysis', 'FontSize', 16, 'FontWeight', 'bold');
p_empirical_agg = NaN(1, 2);

for metric_idx = 1:2
    subplot(2, 1, metric_idx);
    hold on;
    
    % Calculate grand mean from real data
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
    
    % Bootstrapping the pooled mean
    bootstrapped_grand_means = NaN(bootstrap_samples, 1);
    for i = 1:bootstrap_samples
        mean_across_exps = NaN(1, length(experiments));
        for exp_idx = 1:length(experiments)
            sim_params_exp = sim_data_params{exp_idx};
            n_subjects = size(sim_params_exp, 1);
            n_sims_per_subject = size(sim_params_exp, 4);
            
            res_sim_idx = randsample(n_sims_per_subject, n_subjects, true);
            res_anti_idx = randsample(2, n_subjects, true);
            bootstrap_param_sample = NaN(n_subjects, size(sim_params_exp, 2));
            for j = 1:n_subjects
                 bootstrap_param_sample(j, :) = sim_params_exp(j,:,res_anti_idx(j), res_sim_idx(j));
            end
            
            mu = mean(bootstrap_param_sample, 1);
            if metric_idx == 1
                mean_across_exps(exp_idx) = mu(2) - mu(3);
            else
                mean_across_exps(exp_idx) = (mu(2) - mu(3)) / (mu(2) + mu(3));
            end
        end
        bootstrapped_grand_means(i) = mean(mean_across_exps);
    end
    
    p_empirical_agg(metric_idx) = mean(bootstrapped_grand_means > grand_mean_real);
    
    % Pooled Plotting
    x_pos = 1.5;
    [counts, bins] = hist(bootstrapped_grand_means, 20);
    barh(bins, x_pos + counts / (1.2 * max(counts)), 'FaceColor', colors.sim_dist, 'FaceAlpha', 0.5);
    rectangle('Position', [-1, min(bins)-1, x_pos + 1, max(bins)-min(bins)+2], 'FaceColor', 'w', 'EdgeColor', 'w');
    
    m_color = colors.real_ns;
    if p_empirical_agg(metric_idx) < 0.05, m_color = colors.real_sig; end
    plot(x_pos, grand_mean_real, 'ok', 'MarkerSize', 8, 'MarkerFaceColor', m_color, 'LineWidth', 1.2);
    
    plot([0.5, 3], [0, 0], 'k--');
    xticks(1.5); xticklabels('All Experiments');
    xlim([0.5, 3]);
    
    if metric_idx == 1
        ylabel('Absolute Bias (\alpha_c - \alpha_d)'); % Figure 5
        ylim([-0.1, 0.2]);
    else
        ylabel('Normalized Bias'); % S18
        ylim([-0.2, 0.35]);
    end
end

%% 4. Statistical Output
fprintf('\n======================================================\n');
fprintf('STATISTICAL RESULTS FOR FIGURE 5 & S18\n');
fprintf('======================================================\n\n');
for i = 1:length(experiments)
    fprintf('%-12s | Abs p = %.4f | Norm p = %.4f\n', experiments{i}, p_empirical_by_exp(i,1), p_empirical_by_exp(i,2));
end
fprintf('------------------------------------------------------\n');
fprintf('Pooled Results | Abs p = %.4f | Norm p = %.4f\n', p_empirical_agg(1), p_empirical_agg(2));
end