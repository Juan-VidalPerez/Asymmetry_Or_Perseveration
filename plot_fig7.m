function plot_fig7(n_sims)
% PLOT_FIG7 Generates plots for Figure 7.
% Shows average choice rates in specific 'New' task conditions 
% (Low-mean bandits vs High-mean bandits) as confirmation bias or
% perseveration parameters are varied.
%
% USAGE:
%   plot_fig7(10000); (as in paper) 
%
% INPUTS:
%   n_sims - Number of simulated subjects per parameter step.

%% 1. Setup & Data Loading
close all;
fprintf('--- Generating Figure 7: Loading base parameters from MAP_fits.mat ---\n');
try
    map_data = load('data/MAP_fits.mat');
catch ME
    error('Could not load data. Ensure "MAP_fits.mat" is in the MATLAB path. Details: %s', ME.message);
end

% --- Extract average parameters from Experiment P2 (index 4) ---
avg_params_cb_p2 = mean(map_data.parameters_CB{4}, 1);
base_beta_cb = avg_params_cb_p2(1);
base_mean_alpha_cb = 0.4; % Fixed per original instructions

avg_params_psl_p2 = mean(map_data.parameters_PSL{4}, 1);
base_beta_psl = avg_params_psl_p2(1);
base_alpha_psl = avg_params_psl_p2(2);
base_tau_psl = avg_params_psl_p2(3);

% --- Simulation & Parameter Configurations ---
sim_config = struct(...
    'CB', struct(...
        'base_params', [base_beta_cb, base_mean_alpha_cb], ...
        'values', -0.4:0.1:0.4, ...
        'param_label', '\alpha_c - \alpha_d', ...
        'model_idx', 1, ...
        'color_map_func', @(n) colormap(slanCM('berlin', n)) ...
    ), ...
    'PERS', struct(...
        'base_params', [base_beta_psl, base_alpha_psl, base_tau_psl], ...
        'values', -3:0.5:3, ...
        'param_label', '\phi', ...
        'model_idx', 2, ...
        'color_map_func', @(n) colormap(slanCM('vanimo', n)) ...
    ) ...
);
param_types = {'CB', 'PERS'};

conditions_in_sim = [3, 4];   
plot_indices = [1, 2];      
choice_mapping_plot = [2, 3]; 
titles_plot = {'Low-mean bandits', 'High-mean bandits'}; 
ylab_rate_plot = {'Wide bandit choice rate', 'Wide bandit choice rate'}; 
marker = 'o';

%% 2. Initialize Figure
figure('Position', [100, 100, 800, 600], 'Color', 'w');

%% 3. Simulation and Plotting Loop
for type_idx = 1:length(param_types)
    current_type = param_types{type_idx};
    cfg = sim_config.(current_type);
    
    % Generate colors for the parameter sweep
    colors = cfg.color_map_func(length(cfg.values));
    
    for row_idx = 1:length(plot_indices)
        cond_idx_sim = conditions_in_sim(row_idx); 
        subplot_idx = (row_idx - 1) * 2 + type_idx;
        subplot(2, 2, subplot_idx);
        hold on;
        
        all_mean_rates = NaN(1, length(cfg.values));
        all_sem_rates = NaN(1, length(cfg.values));
        
        for val_idx = 1:length(cfg.values)
            current_value = cfg.values(val_idx);
            
            % --- Construct parameter set ---
            if strcmp(current_type, 'CB')
                beta_sim = cfg.base_params(1);
                mean_alpha = cfg.base_params(2);
                lr1 = mean_alpha + current_value / 2; lr2 = mean_alpha - current_value / 2;
                lr1 = max(0, min(1, lr1)); lr2 = max(0, min(1, lr2));
                sim_params = repmat([beta_sim, lr1, lr2], [n_sims, 1]);
            else % PERS
                sim_params = repmat([cfg.base_params(1:3), current_value], [n_sims, 1]);
            end
            
            % --- Run simulation ---
            [sim_data] = simulate_newtask(sim_params, cfg.model_idx);
            
            % --- Calculate average choice rate ---
            sim_rates = NaN(n_sims, 1);
            for ss = 1:n_sims
                choices_in_cond = sim_data{2}{ss}(sim_data{1}{ss} == cond_idx_sim);
                if ~isempty(choices_in_cond)
                    sim_rates(ss) = mean(choices_in_cond == choice_mapping_plot(row_idx));
                end
            end
            
            all_mean_rates(val_idx) = nanmean(sim_rates);
            all_sem_rates(val_idx) = nanstd(sim_rates) / sqrt(sum(~isnan(sim_rates)));
            
            % --- Plot individual point with error bar ---
            errorbar(current_value, all_mean_rates(val_idx), all_sem_rates(val_idx), ...
                     'Color', colors(val_idx, :), 'LineWidth', 1.25, 'CapSize', 0);
            plot(current_value, all_mean_rates(val_idx), marker, 'Color', colors(val_idx, :), ...
                 'LineWidth', 1.25, 'MarkerSize', 8, 'MarkerFaceColor', colors(val_idx, :), ...
                 'MarkerEdgeColor', 'black');
        end 
        
        % --- Dynamic Y-Axis Adjustment ---
        % Calculate limits based on mean +/- SEM with a 10% buffer
        data_min = min(all_mean_rates - all_sem_rates);
        data_max = max(all_mean_rates + all_sem_rates);
        data_range = data_max - data_min;
        if data_range == 0, data_range = 0.1; end % Avoid degenerate limits
        
        ylim([data_min - 0.1*data_range, data_max + 0.1*data_range]);
        
        % --- Finalize Subplot ---
        plot(cfg.values([1, end]), [0.5, 0.5], '-', 'LineWidth', 1, 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off'); % Chance line
        xlim([min(cfg.values), max(cfg.values)]);
        
        title([titles_plot{row_idx} ' (' current_type ')']);
        ylabel(ylab_rate_plot{row_idx});
        
        if row_idx == length(plot_indices)
            xlabel(cfg.param_label);
        else
             xticklabels({}); 
        end
        grid on;
    end 
end
end