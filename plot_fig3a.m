function plot_fig3a()
% PLOT_FIG3A Generates Figure 3a: Assessing the effect of prior width (SD).
%
% This function assesses how the standard deviation of the Gaussian 
% perseveration parameter (phi) prior affects the detection of spurious 
% confirmation bias when fitting PSL-generated data with a hybrid model.
%
% USAGE:
%   plot_fig3a()
%
% REQUIREMENTS:
%   The file 'data/figureSD_data.mat' must exist relative to this script.

    %% 1. Data Loading and Preparation
    fprintf('--- Generating Figure 3a: Assessing Prior Width Effect ---\n');
    
    data_file = 'data/figure3a_data.mat';
    if ~exist(data_file, 'file')
        error('Could not find %s. Please ensure it is in the data/ folder.', data_file);
    end
    
    loaded = load(data_file);
    
    % Extract core variables from the simulation data
    params_fitted = loaded.parameters_PSL_sd; % [SD_steps x Subjects x Parameters]
    params_gen    = loaded.parameters_P2;     % Generating parameters (Empirical P2)
    phi_sd_axis   = loaded.sd;                % X-axis: Prior SD values (0.5 to 10)
    
    % Define ground truth asymptotes based on the PSL generating model
    % Mapping: [beta, lr_pos, lr_neg, tau, phi]
    % Note: Since it is a PSL model, lr_pos = lr_neg (Ground truth bias = 0).
    gen_asymptotes = [0, params_gen(2), params_gen(2), params_gen(3), params_gen(4)];

    %% 2. Plotting Configuration
    % Configuration format: {Label, DataFunc, AsymptoteValue, Fig#, Subplot#, Color}
    configs = { ...
        {'Spurious \alpha_c - \alpha_d', @(p) p(:,:,2) - p(:,:,3), gen_asymptotes(2)-gen_asymptotes(3), 1, 1, [172, 136, 187]/255}; ...
        {'Fitted \alpha_c',          @(p) p(:,:,2),             gen_asymptotes(2),              1, 2, [32, 138, 34]/255}; ...
        {'Fitted \alpha_d',          @(p) p(:,:,3),             gen_asymptotes(3),              1, 2, [192, 0, 0]/255}; ...
        {'Fitted \tau',              @(p) p(:,:,end-1),         gen_asymptotes(4),              2, 1, [172, 136, 187]/255}; ...
        {'Fitted \phi',              @(p) p(:,:,end),           gen_asymptotes(5),              2, 2, [172, 136, 187]/255}  ...
    };

    % Initialize Stat Storage
    n_metrics = size(configs, 1);
    n_steps   = length(phi_sd_axis);
    stats_out = struct('p', NaN(n_metrics, n_steps), 't', NaN(n_metrics, n_steps), 'diff', NaN(n_metrics, n_steps));

    %% 3. Figure Generation
    close all;
    fig_handles(1) = figure('Name', 'Confirmation Bias Metrics', 'Position', [100, 100, 800, 400], 'Color', 'w');
    fig_handles(2) = figure('Name', 'Perseveration Metrics', 'Position', [950, 100, 800, 400], 'Color', 'w');

    for m = 1:n_metrics
        cfg = configs{m};
        
        % Select target figure and subplot
        figure(fig_handles(cfg{4}));
        subplot(1, 2, cfg{5}); hold on;
        
        % Process data
        raw_data = cfg{2}(params_fitted); % Dimensions: [SD_steps x Subjects]
        mu  = mean(raw_data, 2);
        sem = std(raw_data, 0, 2) ./ sqrt(size(raw_data, 2));
        
        % Plot Shaded SEM Area
        fill_x = [phi_sd_axis(:); flipud(phi_sd_axis(:))];
        fill_y = [mu + sem; flipud(mu - sem)];
        patch(fill_x, fill_y, cfg{6}, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        
        % Plot Mean Line and Ground Truth Asymptote
        plot(phi_sd_axis, mu, 'Color', cfg{6}, 'LineWidth', 2.5, 'DisplayName', cfg{1});
        plot(xlim, [cfg{3}, cfg{3}], 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
        
        % Aesthetics
        xlabel('Standard deviation of \phi prior');
        ylabel(cfg{1});
        set(gca, 'Box', 'off', 'TickDir', 'out', 'FontSize', 11);
        grid on;
        
        % Specific Legend for Alpha Comparison Subplot
        if m == 2 || m == 3
            legend('show', 'Location', 'best');
        end
        
        % Perform Statistical Tests across the SD axis
        for s = 1:n_steps
            [~, p, ~, s_res] = ttest(raw_data(s, :) - cfg{3});
            stats_out.p(m, s)    = p;
            stats_out.t(m, s)    = s_res.tstat;
            stats_out.diff(m, s) = mu(s) - cfg{3};
        end
    end

    %% 4. Command Window Statistical Report
    display_report(configs, phi_sd_axis, stats_out);
end

%% --- Helper: Statistical Report Generator ---
function display_report(configs, x_axis, s)
    fprintf('\n======================================================\n');
    fprintf('STATISTICAL ANALYSIS: Spurious Bias vs. Prior Width\n');
    fprintf('======================================================\n');
    
    for m = 1:size(configs, 1)
        fprintf('\nMetric: %s\n', configs{m}{1});
        fprintf('%-10s | %-10s | %-8s | %-8s\n', 'Prior SD', 'Mean Diff', 't-stat', 'p-val');
        fprintf('------------------------------------------------------\n');
        for i = 1:length(x_axis)
            fprintf('%-10.1f | %-+10.4f | %-+8.3f | %.4f\n', ...
                x_axis(i), s.diff(m, i), s.t(m, i), s.p(m, i));
        end
    end
    fprintf('\n------------------ END OF REPORT ------------------\n');
end