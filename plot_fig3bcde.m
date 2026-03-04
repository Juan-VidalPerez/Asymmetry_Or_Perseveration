function plot_fig3bcde()
% Generates all plots for Figure 3: Model and parameter recovery analysis.
%
% This function creates four figures:
% 1. Figure 3b: Parameter Recovery - Scatter plots of generating vs. fitted parameters
%    for the CBPERS model.
% 2. Figure 3c: Recovery Error - Summary of the difference (fitted - generating) 
%    on Hybrid simulation fits.
% 3. Figure 3d: Aggregated MAP vs. MLE Comparison on Real Behavioral Data.
% 4. Figure 3e: Aggregated MAP vs. MLE Comparison on Data Simulated from PSL model.
%
% USAGE:
%   plot_fig3bcde()
%
% REQUIREMENTS:
%   Files 'MAP_fits.mat', 'MLE_fits.mat', and 'figure3bc_data.mat' must be
%   in the 'Data' subdirectory relative to this function's location.

%% -------------------- Data Loading and Preparation --------------------
data_path = 'Data'; 
fprintf('--- Generating Figure 3: Loading data from %s directory... ---\n', data_path);
try
    map_data = load(fullfile(data_path, 'MAP_fits.mat'));
    mle_data = load(fullfile(data_path, 'MLE_fits.mat'));
    recovery_data = load(fullfile(data_path, 'figure3bc_data.mat'));
catch ME
    error('Could not load data. Ensure "MAP_fits.mat", "MLE_fits.mat", and "figure3bc_data.mat" are in the "%s" subdirectory. Details: %s', data_path, ME.message);
end
close all;

% --- Prepare data for Recovery Analysis (Fig 3b) ---
generated_params = {recovery_data.parameters_CBPERSgener_MAP, recovery_data.parameters_CBPERSgener_MLE};
fitted_params = {recovery_data.parameters_CBPERSfitted_MAP, recovery_data.parameters_CBPERSfitted_MLE};
fit_types = {'MAP', 'MLE'};
recovery_errors = cell(2, 2); 
phi_generating_all = cell(1,2);

% --- Prepare data for Aggregated Comparisons (Fig 3d & 3e) ---
params_agg_real_fitted = {map_data.parameters_CBPERS, mle_data.parameters_CBPERS};
params_agg_real_generating = {map_data.parameters_PSL, mle_data.parameters_PSL}; 
params_agg_sim_fitted = {map_data.parameters_PSLsim_CBPERSfit, mle_data.parameters_PSLsim_CBPERSfit};
params_agg_sim_generating = {map_data.parameters_PSL, mle_data.parameters_PSL}; 

%% -------------------- Figure 3b: Parameter Recovery Scatter Plots --------------------
figure('Position', [100, 100, 800, 700], 'Name', 'Figure 3b');
sgtitle('Figure 3b: Parameter Recovery (CBPERS model)', 'FontSize', 16, 'FontWeight', 'bold');
color_scatter = {[107, 142, 185]/255, [172, 136, 187]/255}; 

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

%% -------------------- Figure 3c: Recovery Error Summary (Hybrid simulation) --------------------
figure('Position', [500, 300, 900, 400], 'Name', 'Figure 3c');
sgtitle('Figure 3c: Aggregated Recovery Error (Hybrid simulation)', 'FontSize', 16, 'FontWeight', 'bold');
marker_shapes_agg = {'o', 's'};
colors_d = {[172, 136, 187]/255, [107, 142, 185]/255, [172, 136, 187]/255};

for metric_idx = 1:3
    subplot(1, 3, metric_idx); hold on;
    mean_vals_across_exp = NaN(10, 2);
    
    for fit_idx = 1:2 % 1=MAP, 2=MLE
        vbi = [];
        for exp_idx = 1:10
            fitted = fitted_params{fit_idx}{exp_idx};
            generating = generated_params{fit_idx}{exp_idx};
            
            if ndims(fitted) > 2 && size(fitted,3) > 1
                fitted_avg = squeeze(mean(fitted, 3));
            else
                fitted_avg = fitted;
            end
            
            if metric_idx == 1
                metric_val = fitted_avg(:, 2) - fitted_avg(:, 3) - (generating(:, 2) - generating(:, 3));
            elseif metric_idx == 3
                metric_val = (fitted_avg(:, 2) - fitted_avg(:, 3)) ./ (fitted_avg(:, 2) + fitted_avg(:, 3)) - (generating(:, 2) - generating(:, 3)) ./ (generating(:, 2) + generating(:, 3));
            else 
                metric_val = fitted_avg(:, end) - generating(:, end);
            end
            vbi = [vbi; metric_val]; 
            mean_vals_across_exp(exp_idx, fit_idx) = mean(metric_val); 
        end
        
        sem = std(vbi) / sqrt(length(vbi)); 
        errorbar(fit_idx, mean(vbi), sem, '.k', 'CapSize', 0, 'LineWidth', 1.3);
        plot(fit_idx, mean(vbi), marker_shapes_agg{fit_idx}, 'MarkerFaceColor', colors_d{metric_idx}, 'MarkerEdgeColor', 'k', 'MarkerSize', 10, 'LineWidth', 1.3);
    end
    
    plot([0, 3], [0, 0], 'k--');
    xlim([0.5, 2.5]); xticks([1, 2]); xticklabels({'MAP', 'MLE'});
    if metric_idx == 1
        ylabel('CB_{fitted} - CB_{generating}'); ylim([-0.04, 0.06]);
        legend({'','MAP','','MLE'}, 'Orientation', 'horizontal', 'Location', 'southoutside');
    elseif metric_idx == 3
        ylabel('Norm. CB Error'); ylim([-0.2, 0.8]);
    else
        ylabel('\phi_{fitted} - \phi_{generating}'); ylim([-.4, .1]);
    end
end

%% -------------------- Figure 3d: Aggregated MAP vs. MLE Comparison (Real Data) --------------------
figure('Position', [500, 300, 900, 400], 'Name', 'Figure 3d');
sgtitle('Figure 3d: Aggregated Comparison on Real Data', 'FontSize', 16, 'FontWeight', 'bold');
colors_c = {[172, 136, 187]/255, [107, 142, 185]/255, [172, 136, 187]/255};
p_all_c = NaN(1, 3); t_all_c = NaN(1, 3); mean_diff_c = NaN(1, 3);

for metric_idx = 1:3
    subplot(1, 3, metric_idx); hold on;
    mean_vals_across_exp = NaN(10, 2);
    
    for fit_idx = 1:2 % 1=MAP, 2=MLE
        vbi = [];
        for exp_idx = 1:10
            fitted = params_agg_real_fitted{fit_idx}{exp_idx};
            if metric_idx == 1
                metric_val = fitted(:, 2) - fitted(:, 3);
            elseif metric_idx == 3
                metric_val = (fitted(:, 2) - fitted(:, 3)) ./ (fitted(:, 2) + fitted(:, 3));
            else 
                metric_val = fitted(:, end);
            end
            vbi= [vbi; metric_val]; 
            mean_vals_across_exp(exp_idx, fit_idx) = mean(metric_val); 
        end
       
        sem = std(vbi) / sqrt(length(vbi)); 
        errorbar(fit_idx, mean(vbi), sem, '.k', 'CapSize', 0, 'LineWidth', 1.3);
        plot(fit_idx, mean(vbi), marker_shapes_agg{fit_idx}, 'MarkerFaceColor', colors_c{metric_idx}, 'MarkerEdgeColor', 'k', 'MarkerSize', 10, 'LineWidth', 1.3);
    end
    
    [~, p_all_c(metric_idx), ~, stats_c] = ttest(mean_vals_across_exp(:, 1), mean_vals_across_exp(:, 2));
    t_all_c(metric_idx) = stats_c.tstat;
    mean_diff_c(metric_idx) = mean(mean_vals_across_exp(:, 1) - mean_vals_across_exp(:, 2));
    
    plot([0, 3], [0, 0], 'k--');
    xlim([0.5, 2.5]); xticks([1, 2]); xticklabels({'MAP', 'MLE'});
    if metric_idx == 1
        ylabel('\alpha_c - \alpha_d'); ylim([-0.2, 0.4]);
        legend({'','MAP','','MLE'}, 'Orientation', 'horizontal', 'Location', 'southoutside');
    elseif metric_idx == 3
        ylabel('Normalized CB'); ylim([-0.2, 0.8]);
    else
        ylabel('\phi'); ylim([-.5, 7]); 
    end
end

%% -------------------- Figure 3e: Aggregated MAP vs. MLE Comparison (Simulated PSL) --------------------
figure('Position', [700, 400, 900, 400], 'Name', 'Figure 3e');
sgtitle('Figure 3e: Aggregated Comparison on Simulated Data (PSL simulation)', 'FontSize', 16, 'FontWeight', 'bold');
p_all_e = NaN(1, 3); t_all_e = NaN(1, 3); mean_diff_e = NaN(1, 3);

for metric_idx = 1:3
    subplot(1, 3, metric_idx); hold on;
    mean_vals_across_exp = NaN(10, 2);
    
    for fit_idx = 1:2 % 1=MAP, 2=MLE
        vbi = [];
        for exp_idx = 1:10
            fitted = params_agg_sim_fitted{fit_idx}{exp_idx};
            generating = params_agg_sim_generating{fit_idx}{exp_idx};
            
            if ndims(fitted) > 2 && size(fitted,3) > 1
                fitted_avg = squeeze(mean(fitted, 3));
            else
                fitted_avg = fitted;
            end
            
            if metric_idx == 1
                metric_val = fitted_avg(:, 2) - fitted_avg(:, 3);
            elseif metric_idx == 3
                metric_val = (fitted_avg(:, 2) - fitted_avg(:, 3)) ./ (fitted_avg(:, 2) + fitted_avg(:, 3));
            else 
                metric_val = fitted_avg(:, end) - generating(:, end);
            end
            vbi = [vbi; metric_val]; 
            mean_vals_across_exp(exp_idx, fit_idx) = mean(metric_val); 
        end
        
        sem = std(vbi) / sqrt(length(vbi)); 
        errorbar(fit_idx, mean(vbi), sem, '.k', 'CapSize', 0, 'LineWidth', 1.3);
        plot(fit_idx, mean(vbi), marker_shapes_agg{fit_idx}, 'MarkerFaceColor', colors_d{metric_idx}, 'MarkerEdgeColor', 'k', 'MarkerSize', 10, 'LineWidth', 1.3);
    end
    
    [~, p_all_e(metric_idx), ~, stats_e] = ttest(mean_vals_across_exp(:, 1), mean_vals_across_exp(:, 2));
    t_all_e(metric_idx) = stats_e.tstat;
    mean_diff_e(metric_idx) = mean(mean_vals_across_exp(:, 1) - mean_vals_across_exp(:, 2));
    
    plot([0, 3], [0, 0], 'k--');
    xlim([0.5, 2.5]); xticks([1, 2]); xticklabels({'MAP', 'MLE'});
    if metric_idx == 1
        ylabel('CB_{fitted} (CB_{generating}=0)'); ylim([-0.2, 0.4]);
        legend({'','MAP','','MLE'}, 'Orientation', 'horizontal', 'Location', 'southoutside');
    elseif metric_idx == 3
        ylabel('Normalized CB'); ylim([-0.2, 0.8]);
    else
        ylabel('\phi_{fit} - \phi_{gen}'); ylim([-3, 1]);
    end
end

%% -------------------- Statistical Analyses Output --------------------
fprintf('\n======================================================\n');
fprintf('STATISTICAL TESTS FOR FIGURE 3\n');
fprintf('======================================================\n\n');

fprintf('--- Part 3b-c: Recovery Error T-Tests (vs. Zero) ---\n');
[~, p_map_phi_err, ~, stats_map_phi] = ttest(recovery_errors{1, 1});
[~, p_mle_phi_err, ~, stats_mle_phi] = ttest(recovery_errors{2, 1});
[~, p_map_cb_err, ~, stats_map_cb] = ttest(recovery_errors{1, 2});
[~, p_mle_cb_err, ~, stats_mle_cb] = ttest(recovery_errors{2, 2});
fprintf('  MAP Phi Error vs 0: mean = %.3f, t = %.3f, p = %.4f\n', mean(recovery_errors{1,1}), stats_map_phi.tstat, p_map_phi_err);
fprintf('  MLE Phi Error vs 0: mean = %.3f, t = %.3f, p = %.4f\n', mean(recovery_errors{2,1}), stats_mle_phi.tstat, p_mle_phi_err);
fprintf('  MAP CB Error vs 0:  mean = %.3f, t = %.3f, p = %.4f\n', mean(recovery_errors{1,2}), stats_map_cb.tstat, p_map_cb_err);
fprintf('  MLE CB Error vs 0:  mean = %.3f, t = %.3f, p = %.4f\n', mean(recovery_errors{2,2}), stats_mle_cb.tstat, p_mle_cb_err);

fprintf('\n--- Part 3b-c: Recovery Error Comparison (MAP vs. MLE) ---\n');
[~, p_pers, ~, stats_pers] = ttest2(recovery_errors{1, 1}, recovery_errors{2, 1});
[~, p_cb, ~, stats_cb] = ttest2(recovery_errors{1, 2}, recovery_errors{2, 2});
fprintf('  Phi Error (MAP vs MLE): mean diff = %.3f, t = %.3f, p = %.4f\n', mean(recovery_errors{1,1})-mean(recovery_errors{2,1}), stats_pers.tstat, p_pers);
fprintf('  CB Error (MAP vs MLE):  mean diff = %.3f, t = %.3f, p = %.4f\n', mean(recovery_errors{1,2})-mean(recovery_errors{2,2}), stats_cb.tstat, p_cb);

fprintf('\n--- Part 3d: Aggregated MAP vs. MLE t-tests on REAL data ---\n');
fprintf('  CB difference:        mean diff = %.3f, t = %.3f, p = %.4f\n', mean_diff_c(1), t_all_c(1), p_all_c(1));
fprintf('  Normalized CB diff:   mean diff = %.3f, t = %.3f, p = %.4f\n', mean_diff_c(2), t_all_c(2), p_all_c(2));
fprintf('  Phi difference:       mean diff = %.3f, t = %.3f, p = %.4f\n', mean_diff_c(3), t_all_c(3), p_all_c(3));

fprintf('\n--- Part 3e: Aggregated MAP vs. MLE t-tests on SIMULATED data ---\n');
fprintf('  CB difference:        mean diff = %.3f, t = %.3f, p = %.4f\n', mean_diff_e(1), t_all_e(1), p_all_e(1));
fprintf('  Normalized CB diff:   mean diff = %.3f, t = %.3f, p = %.4f\n', mean_diff_e(2), t_all_e(2), p_all_e(2));
fprintf('  Phi error difference: mean diff = %.3f, t = %.3f, p = %.4f\n', mean_diff_e(3), t_all_e(3), p_all_e(3));

fprintf('\n-------------------- END OF TESTS --------------------\n\n');
end