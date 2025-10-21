function plot_fig8(n_sims)
% Generates plots for Figure 8, showing average choice rates in specific
% 'New' task conditions (HL vs LL and HH vs LH) as confirmation bias or
% perseveration parameters are varied.
% Creates a single figure with a 2x2 layout (Condition x Parameter Type).
%
% USAGE:
%   plot_fig8(100); % Specify the number of simulations per parameter value
%
% INPUTS:
%   n_sims - Number of simulated subjects.
%
% REQUIREMENTS:
%   The function 'simulate_newtask.m' (containing the simulation logic for
%   the 'New' task) and the file 'MAP_fits.mat' must be in the MATLAB path.

%% 1. Setup & Data Loading
close all;
fprintf('--- Generating Figure 8: Loading base parameters from MAP_fits.mat ---\n');

try
    map_data = load('data/MAP_fits.mat');
catch ME
    error('Could not load data. Ensure "MAP_fits.mat" is in the MATLAB path. Details: %s', ME.message);
end

% --- Extract average parameters from Experiment P2 (index 4) ---
avg_params_cb_p2 = mean(map_data.parameters_CB{4}, 1);
base_beta_cb = avg_params_cb_p2(1);
base_mean_alpha_cb = 0.4; % Fixed as per instructions

avg_params_psl_p2 = mean(map_data.parameters_PSL{4}, 1);
base_beta_psl = avg_params_psl_p2(1);
base_alpha_psl = avg_params_psl_p2(2);
base_tau_psl = avg_params_psl_p2(3);

% --- Simulation & Parameter Configurations ---
sim_config = struct(...
    'CB', struct(...
        'base_params', [base_beta_cb, base_mean_alpha_cb], ...
        'values', -0.6:0.1:0.6, ...
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

% --- Task/Plotting Configurations (Simplified for Conditions 3 & 4 only) ---
conditions_in_sim = [3, 4];   % The original condition numbers to extract from sim_data
plot_indices = [1, 2];      % Indices for plotting loop (rows)
choice_mapping_plot = [2, 3]; % Target choices for conditions 3 and 4
titles_plot = {'Low-mean bandits', 'High-mean bandits'}; % Titles for conditions 3 and 4
ylab_rate_plot = {'Wide bandit choice rate', 'Wide bandit choice rate'}; % Y-labels for conditions 3 and 4
marker = 'o';

%% 2. Initialize Figure
figure('Position', [100, 100, 800, 600]);

%% 3. Simulation and Plotting Loop
% Loop through parameter types (CB, PERS) to create columns
for type_idx = 1:length(param_types)
    current_type = param_types{type_idx};
    cfg = sim_config.(current_type);
    colors = cfg.color_map_func(length(cfg.values));

    % Loop through the selected conditions (using plot_indices 1 and 2)
    for row_idx = 1:length(plot_indices)
        cond_idx_sim = conditions_in_sim(row_idx); % Actual condition index from simulation (3 or 4)

        % Subplot indexing for 2x2 layout
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

            % --- Calculate average choice rate for the target choice in this condition ---
            sim_rates = NaN(n_sims, 1);
            for ss = 1:n_sims
                % Use cond_idx_sim to extract data from simulation output
                choices_in_cond = sim_data{2}{ss}(sim_data{1}{ss} == cond_idx_sim);
                if isempty(choices_in_cond)
                    warning('No trials found for condition %d in simulation %d. Check simulation logic.', cond_idx_sim, ss);
                    sim_rates(ss) = NaN;
                else
                    % Use row_idx to get correct target choice from the simplified array
                    sim_rates(ss) = mean(choices_in_cond == choice_mapping_plot(row_idx));
                end
            end
            all_mean_rates(val_idx) = nanmean(sim_rates); % Use nanmean
            all_sem_rates(val_idx) = nanstd(sim_rates) / sqrt(sum(~isnan(sim_rates)));

            % --- Plot individual point with error bar ---
            errorbar(current_value, all_mean_rates(val_idx), all_sem_rates(val_idx), ...
                     'Color', colors(val_idx, :), 'LineWidth', 1.25, 'CapSize', 0);
            plot(current_value, all_mean_rates(val_idx), marker, 'Color', colors(val_idx, :), ...
                 'LineWidth', 1.25, 'MarkerSize', 8, 'MarkerFaceColor', colors(val_idx, :), ...
                 'MarkerEdgeColor', 'black');

        end % End loop over values

        % --- Finalize Subplot ---
        plot(cfg.values([1, end]), [0.5, 0.5], '-', 'LineWidth', 1, 'Color', [0.5 0.5 0.5]); % Chance line
        xlim([min(cfg.values), max(cfg.values)]);
        ylim([0, 1]); % Probability scale
        % Use row_idx to get title and ylabel from simplified arrays
        title([titles_plot{row_idx} ' (' current_type ')']);
        ylabel(ylab_rate_plot{row_idx});
        
        % Only add x-label to bottom plots (row_idx == 2)
        if row_idx == length(plot_indices)
            xlabel(cfg.param_label);
        else
             xticklabels({}); % Remove x-ticks for top row
        end

    end % End loop over conditions (rows)
end % End loop over parameter types (columns)

end % End of main function