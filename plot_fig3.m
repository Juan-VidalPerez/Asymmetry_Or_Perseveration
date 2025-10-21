function plot_fig3()
% Generates all plots for Figure 3: Model and parameter recovery analysis.
%
% This function creates four figures:
% 1. Figure 3a: Parameter Recovery - Scatter plots of generating vs. fitted parameters
%    for the CBPERS model.
% 2. Figure 3b: Recovery Error - Summary of the difference (fitted - generating).
% 3. Figure 3c: Aggregated MAP vs. MLE Comparison on Real Data.
% 4. Figure 3d: Aggregated MAP vs. MLE Comparison on Data Simulated from PSL model.
%
% USAGE:
%   plot_fig3()
%
% REQUIREMENTS:
%   Files 'MAP_fits.mat', 'MLE_fits.mat', and 'CBPERS_recoveries.mat' must be
%   in the 'Data' subdirectory relative to this function's location.

%% -------------------- Data Loading and Preparation --------------------
data_path = 'Data'; % Define the data subdirectory
fprintf('--- Generating Figure 3: Loading data from %s directory... ---\n', data_path);
try
    map_data = load(fullfile(data_path, 'MAP_fits.mat'));
    mle_data = load(fullfile(data_path, 'MLE_fits.mat'));
    recovery_data = load(fullfile(data_path, 'CBPERS_recoveries.mat'));
catch ME
    error('Could not load data. Ensure "MAP_fits.mat", "MLE_fits.mat", and "CBPERS_recoveries.mat" are in the "%s" subdirectory. Details: %s', data_path, ME.message);
end
close all;

% --- Prepare data for Recovery Analysis (Fig 3a & 3b) from CBPERS_recoveries.mat ---
generated_params = {recovery_data.parameters_CBPERSgener_MAP, recovery_data.parameters_CBPERSgener_MLE};
fitted_params = {recovery_data.parameters_CBPERSfitted_MAP, recovery_data.parameters_CBPERSfitted_MLE};
fit_types = {'MAP', 'MLE'};
recovery_errors = cell(2, 2); % To store recovery errors {fit_type, metric}
phi_generating_all = cell(1,2);

% --- Prepare data for Aggregated Comparisons (Fig 3c & 3d) from main fit files ---
params_agg_real_fitted = {map_data.parameters_CBPERS, mle_data.parameters_CBPERS};
params_agg_real_generating = {map_data.parameters_PSL, mle_data.parameters_PSL}; % Used in Fig 3c, metric 3
params_agg_sim_fitted = {map_data.parameters_PSLsim_CBPERSfit, mle_data.parameters_PSLsim_CBPERSfit};
params_agg_sim_generating = {map_data.parameters_PSL, mle_data.parameters_PSL}; % Used in Fig 3d, metric 3

%% -------------------- Figure 3a: Parameter Recovery Scatter Plots --------------------
figure('Position', [100, 100, 800, 700]);
sgtitle('Figure 3a: Parameter Recovery (CBPERS model)', 'FontSize', 16, 'FontWeight', 'bold');
color_scatter = {[107, 142, 185]/255, [172, 136, 187]/255}; % Blue for phi, Purple for CB

for i = 1:2 % Loop over MAP (1) and MLE (2)
    gen_phi_vec = []; fit_phi_vec = [];
    gen_cb_vec = [];  fit_cb_vec = [];
    
    for exp_idx = 1:numel(generated_params{i})
        gen_matrix_exp = generated_params{i}{exp_idx};
        fit_matrix_exp = fitted_params{i}{exp_idx};
        
        gen_phi_vec = [gen_phi_vec; gen_matrix_exp(:, end)];
        fit_phi_vec = [fit_phi_vec; fit_matrix_exp(:, end)];
        gen_cb_vec = [gen_cb_vec; gen_matrix_exp(:, 2) - gen_matrix_exp(:, 3)];
        fit_cb_vec = [fit_cb_vec; fit_matrix_exp(:, 2) - fit_matrix_exp(:, 3)];
    end
    
    recovery_errors{i, 1} = fit_phi_vec - gen_phi_vec;
    recovery_errors{i, 2} = fit_cb_vec - gen_cb_vec;
    phi_generating_all{i} = gen_phi_vec;
    
    subplot(2, 2, i);
    scatter(gen_phi_vec, fit_phi_vec, 10, color_scatter{1}, 'filled', 'MarkerFaceAlpha', 0.15);
    hold on;
    p_min = min([gen_phi_vec; fit_phi_vec]); p_max = max([gen_phi_vec; fit_phi_vec]);
    plot([p_min, p_max], [p_min, p_max], 'k--', 'LineWidth', 1.5);
    xlabel('Generating \phi'); ylabel('Fitted \phi');
    title(fit_types{i}); xlim([p_min, p_max]); ylim([p_min, p_max]);
    
    subplot(2, 2, 2 + i);
    scatter(gen_cb_vec, fit_cb_vec, 10, color_scatter{2}, 'filled', 'MarkerFaceAlpha', 0.15);
    hold on;
    p_min = min([gen_cb_vec; fit_cb_vec]); p_max = max([gen_cb_vec; fit_cb_vec]);
    plot([p_min, p_max], [p_min, p_max], 'k--', 'LineWidth', 1.5);
    xlabel('Generating \alpha_c - \alpha_d');
    ylabel('Fitted \alpha_c - \alpha_d');
    xlim([p_min, p_max]); ylim([p_min, p_max]);
end

%% -------------------- Figure 3b: Recovery Error Summary --------------------
figure('Position', [300, 200, 500, 600]);
sgtitle('Figure 3b: Recovery Error (Fitted - Generating)', 'FontSize', 16, 'FontWeight', 'bold');
marker_shapes = {'o', 's'};
color_summary = {[107, 142, 185]/255, [172, 136, 187]/255}; % Blue for phi, Purple for CB

subplot(2, 1, 1); hold on;
for i = 1:2 % MAP, MLE
    mean_err = mean(recovery_errors{i, 1});
    sem_err = std(recovery_errors{i, 1}) / sqrt(length(recovery_errors{i, 1}));
    errorbar(i-1, mean_err, sem_err, '.k', 'CapSize', 0, 'LineWidth', 1.3);
    plot(i-1, mean_err, marker_shapes{i}, 'MarkerFaceColor', color_summary{1}, 'MarkerEdgeColor', 'k', 'MarkerSize', 10, 'LineWidth', 1.3);
end
plot([-0.5, 1.5], [0, 0], '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2);
xlim([-0.5, 1.5]); ylim([-0.4, 0.1]);
xticks([0, 1]); xticklabels({'MAP', 'MLE'});
ylabel('Error in \phi');

subplot(2, 1, 2); hold on;
for i = 1:2 % MAP, MLE
    mean_err = mean(recovery_errors{i, 2});
    sem_err = std(recovery_errors{i, 2}) / sqrt(length(recovery_errors{i, 2}));
    errorbar(i-1, mean_err, sem_err, '.k', 'CapSize', 0, 'LineWidth', 1.3);
    plot(i-1, mean_err, marker_shapes{i}, 'MarkerFaceColor', color_summary{2}, 'MarkerEdgeColor', 'k', 'MarkerSize', 10, 'LineWidth', 1.3);
end
plot([-0.5, 1.5], [0, 0], '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2);
xlim([-0.5, 1.5]); ylim([-0.04, 0.07]);
xticks([0, 1]); xticklabels({'MAP', 'MLE'});
ylabel('Error in \alpha_c - \alpha_d');

%% -------------------- Figure 3c: Aggregated MAP vs. MLE Comparison (Real Data) --------------------
figure('Position', [500, 300, 900, 400]);
sgtitle('Figure 3c: Aggregated Comparison on Real Data', 'FontSize', 16, 'FontWeight', 'bold');
marker_shapes_agg = {'o', 's'};
colors_c = {[181, 101, 118]/255, [181, 101, 118]/255};
p_all_c = NaN(1, 3);
t_all_c = NaN(1, 3);
mean_diff_c = NaN(1, 3);

for metric_idx = 1:3
    subplot(1, 3, metric_idx); hold on;
    mean_vals_across_exp = NaN(10, 2);
    
    for fit_idx = 1:2 % 1=MAP, 2=MLE
        vbi = [];
        for exp_idx = 1:10
            fitted = params_agg_real_fitted{fit_idx}{exp_idx};
            generating = params_agg_real_generating{fit_idx}{exp_idx}; % Only used for metric_idx==3
            
            if metric_idx == 1
                metric_val = fitted(:, 2) - fitted(:, 3);
            elseif metric_idx == 2
                metric_val = (fitted(:, 2) - fitted(:, 3)) ./ (fitted(:, 2) + fitted(:, 3));
            else % Plot fitted phi directly (not the error vs generating)
                metric_val = fitted(:, end);
            end
            vbi= [vbi; metric_val]; % Aggregate all subjects across experiments
            mean_vals_across_exp(exp_idx, fit_idx) = mean(metric_val); % Store mean per experiment for t-test
        end
       
        sem = std(vbi) / sqrt(length(vbi)); % SEM across all subjects
        errorbar(fit_idx, mean(vbi), sem, '.k', 'CapSize', 0, 'LineWidth', 1.3);
        plot(fit_idx, mean(vbi), marker_shapes_agg{fit_idx}, 'MarkerFaceColor', colors_c{fit_idx}, 'MarkerEdgeColor', 'k', 'MarkerSize', 10, 'LineWidth', 1.3);
    end
    
    % T-test compares the *means across experiments*
    [~, p_all_c(metric_idx), ~, stats_c] = ttest(mean_vals_across_exp(:, 1), mean_vals_across_exp(:, 2));
    t_all_c(metric_idx) = stats_c.tstat;
    mean_diff_c(metric_idx) = mean(mean_vals_across_exp(:, 1) - mean_vals_across_exp(:, 2));
    
    plot([0, 3], [0, 0], 'k--');
    xlim([0.5, 2.5]); xticks([1, 2]); xticklabels({'MAP', 'MLE'});
    if metric_idx == 1
        ylabel('\alpha_c - \alpha_d'); ylim([-0.2, 0.4]);
        legend({'','MAP (with priors)','','MLE (without priors)'}, 'Orientation', 'horizontal', 'Location', 'southoutside');
    elseif metric_idx == 2
        ylabel('Normalized CB'); ylim([-0.2, 0.8]);
    else
        ylabel('\phi'); ylim([-.5, 7]); % Adjusted label
    end
end

%% -------------------- Figure 3d: Aggregated MAP vs. MLE Comparison (Simulated Data) --------------------
figure('Position', [700, 400, 900, 400]);
sgtitle('Figure 3d: Aggregated Comparison on Simulated Data (PSLsim fits)', 'FontSize', 16, 'FontWeight', 'bold');
colors_d = {[172, 136, 187]/255, [172, 136, 187]/255};
p_all_d = NaN(1, 3);
t_all_d = NaN(1, 3);
mean_diff_d = NaN(1, 3);

for metric_idx = 1:3
    subplot(1, 3, metric_idx); hold on;
    mean_vals_across_exp = NaN(10, 2);
    
    for fit_idx = 1:2 % 1=MAP, 2=MLE
        vbi = [];
        for exp_idx = 1:10
            fitted = params_agg_sim_fitted{fit_idx}{exp_idx};
            generating = params_agg_sim_generating{fit_idx}{exp_idx};
            
            % Average across simulations first
            if ndims(fitted) > 2 && size(fitted,3) > 1
                fitted_avg = squeeze(mean(fitted, 3));
            else
                fitted_avg = fitted;
            end
            
            if metric_idx == 1
                metric_val = fitted_avg(:, 2) - fitted_avg(:, 3);
            elseif metric_idx == 2
                metric_val = (fitted_avg(:, 2) - fitted_avg(:, 3)) ./ (fitted_avg(:, 2) + fitted_avg(:, 3));
            else % Plot error vs generating phi
                metric_val = fitted_avg(:, end) - generating(:, end);
            end
            vbi = [vbi;metric_val]; % Aggregate all subjects across experiments
            mean_vals_across_exp(exp_idx, fit_idx) = mean(metric_val); % Store mean per experiment for t-test
        end
        
        sem = std(vbi) / sqrt(length(vbi)); % SEM across all subjects
        errorbar(fit_idx, mean(vbi), sem, '.k', 'CapSize', 0, 'LineWidth', 1.3);
        plot(fit_idx, mean(vbi), marker_shapes_agg{fit_idx}, 'MarkerFaceColor', colors_d{fit_idx}, 'MarkerEdgeColor', 'k', 'MarkerSize', 10, 'LineWidth', 1.3);
    end
    
    % T-test compares the *means across experiments*
    [~, p_all_d(metric_idx), ~, stats_d] = ttest(mean_vals_across_exp(:, 1), mean_vals_across_exp(:, 2));
    t_all_d(metric_idx) = stats_d.tstat;
    mean_diff_d(metric_idx) = mean(mean_vals_across_exp(:, 1) - mean_vals_across_exp(:, 2));
    
    plot([0, 3], [0, 0], 'k--');
    xlim([0.5, 2.5]); xticks([1, 2]); xticklabels({'MAP', 'MLE'});
    if metric_idx == 1
        ylabel('\alpha_c - \alpha_d'); ylim([-0.2, 0.4]);
        legend({'','MAP (with priors)','','MLE (without priors)'}, 'Orientation', 'horizontal', 'Location', 'southoutside');
    elseif metric_idx == 2
        ylabel('Normalized CB'); ylim([-0.2, 0.8]);
    else
        ylabel('\phi_{fitted} - \phi_{generating}'); ylim([-3, 1]);
    end
end

%% -------------------- Statistical Analyses Output --------------------
fprintf('\n======================================================\n');
fprintf('STATISTICAL TESTS FOR FIGURE 3\n');
fprintf('======================================================\n\n');
fprintf('--- Part A/B: Recovery Error T-Tests (vs. Zero) ---\n');
[~, p_map_phi_err, ~, stats_map_phi] = ttest(recovery_errors{1, 1});
[~, p_mle_phi_err, ~, stats_mle_phi] = ttest(recovery_errors{2, 1});
[~, p_map_cb_err, ~, stats_map_cb] = ttest(recovery_errors{1, 2});
[~, p_mle_cb_err, ~, stats_mle_cb] = ttest(recovery_errors{2, 2});
fprintf('  MAP Phi Error vs 0: mean = %.3f, t = %.3f, p = %.4f\n', mean(recovery_errors{1,1}), stats_map_phi.tstat, p_map_phi_err);
fprintf('  MLE Phi Error vs 0: mean = %.3f, t = %.3f, p = %.4f\n', mean(recovery_errors{2,1}), stats_mle_phi.tstat, p_mle_phi_err);
fprintf('  MAP CB Error vs 0:  mean = %.3f, t = %.3f, p = %.4f\n', mean(recovery_errors{1,2}), stats_map_cb.tstat, p_map_cb_err);
fprintf('  MLE CB Error vs 0:  mean = %.3f, t = %.3f, p = %.4f\n', mean(recovery_errors{2,2}), stats_mle_cb.tstat, p_mle_cb_err);
fprintf('\n--- Part B: Recovery Error T-Tests (MAP vs. MLE) ---\n');
[~, p_pers, ~, stats_pers] = ttest2(recovery_errors{1, 1}, recovery_errors{2, 1});
[~, p_cb, ~, stats_cb] = ttest2(recovery_errors{1, 2}, recovery_errors{2, 2});
fprintf('  Phi Error (MAP vs MLE): mean diff = %.3f, t = %.3f, p = %.4f\n', mean(recovery_errors{1,1})-mean(recovery_errors{2,1}), stats_pers.tstat, p_pers);
fprintf('  CB Error (MAP vs MLE):  mean diff = %.3f, t = %.3f, p = %.4f\n', mean(recovery_errors{1,2})-mean(recovery_errors{2,2}), stats_cb.tstat, p_cb);
fprintf('\n--- Part A: Correlation between generating phi and recovery error (Spearman) ---\n');
[r_map_phi, p_map_phi_corr] = corr(phi_generating_all{1}, recovery_errors{1, 1}, 'Type', 'Spearman');
[r_map_cb, p_map_cb_corr] = corr(phi_generating_all{1}, recovery_errors{1, 2}, 'Type', 'Spearman');
[r_mle_phi, p_mle_phi_corr] = corr(phi_generating_all{2}, recovery_errors{2, 1}, 'Type', 'Spearman');
[r_mle_cb, p_mle_cb_corr] = corr(phi_generating_all{2}, recovery_errors{2, 2}, 'Type', 'Spearman');
fprintf('  MAP: Generating phi vs. phi error: r = %.3f, p = %.4f\n', r_map_phi, p_map_phi_corr);
fprintf('  MAP: Generating phi vs. CB error:  r = %.3f, p = %.4f\n', r_map_cb, p_map_cb_corr);
fprintf('  MLE: Generating phi vs. phi error: r = %.3f, p = %.4f\n', r_mle_phi, p_mle_phi_corr);
fprintf('  MLE: Generating phi vs. CB error:  r = %.3f, p = %.4f\n', r_mle_cb, p_mle_cb_corr);
fprintf('\n--- Part C: Aggregated MAP vs. MLE t-tests on REAL data ---\n');
fprintf('  CB difference:        mean diff = %.3f, t = %.3f, p = %.4f\n', mean_diff_c(1), t_all_c(1), p_all_c(1));
fprintf('  Normalized CB difference: mean diff = %.3f, t = %.3f, p = %.4f\n', mean_diff_c(2), t_all_c(2), p_all_c(2));
fprintf('  Phi difference:       mean diff = %.3f, t = %.3f, p = %.4f\n', mean_diff_c(3), t_all_c(3), p_all_c(3));
fprintf('\n--- Part D: Aggregated MAP vs. MLE t-tests on SIMULATED data ---\n');
fprintf('  CB difference:        mean diff = %.3f, t = %.3f, p = %.4f\n', mean_diff_d(1), t_all_d(1), p_all_d(1));
fprintf('  Normalized CB difference: mean diff = %.3f, t = %.3f, p = %.4f\n', mean_diff_d(2), t_all_d(2), p_all_d(2));
fprintf('  Phi error difference: mean diff = %.3f, t = %.3f, p = %.4f\n', mean_diff_d(3), t_all_d(3), p_all_d(3));
fprintf('\n-------------------- END OF TESTS --------------------\n\n');
end