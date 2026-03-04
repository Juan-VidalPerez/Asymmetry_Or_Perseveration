function plot_fig4a(varargin)
% PLOT_FIG4A Generates plots for Figure 4a and supplementary figures S10 and S17.
%
% This function assesses parameter recovery as a function of training sessions.
% The x-axis is scaled such that ticks only appear at powers of 2.
%
% USAGE:
%   plot_fig4a()         % Default: MLE fits (Fig 4a, S17)
%   plot_fig4a('MAP')    % MAP fits (Fig S10)

%% -------------------- Data Loading and Preparation --------------------
fit_type = 'MLE';
if nargin > 0 && strcmpi(varargin{1}, 'MAP')
    fit_type = 'MAP';
end

try
    data = load('data/figure4a_data.mat');
catch
    error('Ensure "figure4a_data.mat" is in the data/ folder.');
end
close all;

if strcmp(fit_type, 'MAP')
    params_fitted     = data.parameters_PSLsim_CBPERSfit_MAP;
    params_generating = data.parameters_PSL_MAP;
    session_lengths   = data.session_length_MAP;
else
    params_fitted     = data.parameters_PSLsim_CBPERSfit_MLE;
    params_generating = data.parameters_PSL_MLE;
    session_lengths   = data.session_length_MLE;
end

% Convert sessions to log2 space for the x-axis
x_axis = log2(session_lengths); 
gen_asymptotes = [0, params_generating(2), params_generating(2), ...
                     params_generating(3), params_generating(4)];

%% -------------------- Plotting Configuration --------------------
plot_configs = { ...
    {'Fitted \alpha_c - \alpha_d', @(p) p(:,:,2) - p(:,:,3), @(a) a(2) - a(3), 1, 1, [172, 136, 187]/255}; ...
    {'Fitted \alpha_c',          @(p) p(:,:,2),             @(a) a(2),          1, 2, [32, 138, 34]/255}; ...
    {'Fitted \alpha_d',          @(p) p(:,:,3),             @(a) a(3),          1, 2, [192, 0, 0]/255}; ...
    {'Fitted \tau',              @(p) p(:,:,end-1),         @(a) a(end-1),      2, 1, [172, 136, 187]/255}; ...
    {'Fitted \phi',              @(p) p(:,:,end),           @(a) a(end),        2, 2, [172, 136, 187]/255}  ...
};

%% -------------------- Figure Generation --------------------
fig_handles(1) = figure('Position', [100, 100, 600, 300], 'Color', 'w');
fig_handles(2) = figure('Position', [750, 100, 600, 300], 'Color', 'w');

for mm = 1:size(plot_configs, 1)
    cfg = plot_configs{mm};
    figure(fig_handles(cfg{4}));
    subplot(1, 2, cfg{5}); hold on;

    % Data processing
    y_data = cfg{2}(params_fitted);
    mu  = mean(y_data, 2);
    sem = std(y_data, 0, 2) ./ sqrt(size(y_data, 2));
    
    % Shaded Error Bars
    fill_x = [x_axis(:); flipud(x_axis(:))];
    fill_y = [mu + sem; flipud(mu - sem)];
    patch(fill_x, fill_y, cfg{6}, 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    
    % Plots
    plot(x_axis, mu, 'LineWidth', 2, 'Color', cfg{6}, 'DisplayName', cfg{1});
    plot(xlim, [cfg{3}(gen_asymptotes), cfg{3}(gen_asymptotes)], 'k--', 'HandleVisibility', 'off');
    
    % --- X-AXIS FIX ---
    % Force ticks to only occur at the actual session points (log2 values)
    plot_x_axis=0:round(max(x_axis));
    xticks(plot_x_axis); 
    % Label them as the original session numbers (2^x)
    xticklabels(cellstr(num2str(2.^plot_x_axis'))); 
    
    xlabel('Number of training sessions');
    ylabel(cfg{1});
    set(gca, 'Box', 'off', 'TickDir', 'out', 'FontSize', 10);
    xlim([min(x_axis) max(x_axis)]);
    
    if ismember(mm, [2, 3]), legend('show', 'Location', 'best'); end
end


end