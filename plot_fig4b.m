function plot_fig4b(varargin)
% PLOT_FIG4B Generates Figure 4b and supplementary figures S11 and S19.
%
% This function performs a "Parameter Sweep" analysis. It plots how 
% fitted parameters (Bias and Learning Rates) change when one generative 
% parameter (beta, alpha, tau, or phi) is varied systematically.
%
% FIGURES PRODUCED:
%   - Figure 4b: Main parameter sweep (MLE fits) showing spurious bias induction.
%   - Figure S11: MAP version of the parameter sweep.
%   - Figure S19: Detailed MLE sweep across all generative parameters.
%
% USAGE:
%   plot_fig4b()         % Default: MLE fits (Fig 4b, S21)
%   plot_fig4b('MAP')    % MAP fits (Fig S12)

%% 1. Data Loading and Preparation
fit_type = 'MLE';
if nargin > 0 && strcmpi(varargin{1}, 'MAP')
    fit_type = 'MAP';
end

fprintf('--- Generating Fig 4b/S11/S19: Loading %s sweep data... ---\n', fit_type);

try
    data = load('data/parameters_sweep.mat');
catch
    error('Ensure "parameters_sweep.mat" is in the data/ folder.');
end
close all;

% Select data based on fit type
if strcmp(fit_type, 'MAP')
    params_fitted_sweep = data.parameters_sweep_MAP; 
    sweep_values = data.swept_MAP;                   
    generative_ref = data.generative_MAP;            
else 
    params_fitted_sweep = data.parameters_sweep_MLE;
    sweep_values = data.swept_MLE;
    generative_ref = data.generative_MLE;
end

swept_param_labels = {'\beta', '\alpha', '\tau', '\phi'};
num_swept_params = length(params_fitted_sweep);

%% 2. Generate Figure
figure('Position', [50, 50, 1200, 800], 'Color', 'w', 'Name', ['Fig 4b: ' fit_type]);
sgtitle(['Parameter Recovery Sweep: ' fit_type ' Estimation'], 'FontSize', 16, 'FontWeight', 'bold');

% Define visual palette
color_cb = [172, 136, 187] / 255; % Purple: Bias
color_ac = [32, 138, 34] / 255;  % Green: alpha_c
color_ad = [192, 0, 0] / 255;    % Red: alpha_d

for col_idx = 1:num_swept_params 
    
    % Extract data for current swept parameter column
    % matrix shape: [sweep_val_idx, subject_idx, fitted_param_idx]
    fitted_params_matrix = params_fitted_sweep{col_idx}; 
    x_sweep = sweep_values{col_idx};                     
    
    %% ROW 1: Absolute Spurious Confirmation Bias (alpha_c - alpha_d)
    subplot(3, num_swept_params, col_idx);
    hold on;
    data_abs_cb = squeeze(fitted_params_matrix(:,2,:) - fitted_params_matrix(:,3,:));
    plot_sweep_metric(x_sweep, data_abs_cb, color_cb, generative_ref, col_idx, 'y=0', 0);
    if col_idx == 1, ylabel('Fitted \alpha_c - \alpha_d'); end
    title(['Sweeping ' swept_param_labels{col_idx}]);
    set(gca, 'XTickLabel', []);

    %% ROW 2: Normalized Confirmation Bias
    subplot(3, num_swept_params, col_idx + num_swept_params);
    hold on;
    % (ac - ad) / (ac + ad)
    data_norm_cb = squeeze((fitted_params_matrix(:,2,:) - fitted_params_matrix(:,3,:)) ./ ...
                           (fitted_params_matrix(:,2,:) + fitted_params_matrix(:,3,:)));
    plot_sweep_metric(x_sweep, data_norm_cb, color_cb, generative_ref, col_idx, 'y=0', 0);
    if col_idx == 1, ylabel('Norm. Bias'); end
    set(gca, 'XTickLabel', []);

    %% ROW 3: Individual Learning Rates (alpha_c vs alpha_d)
    subplot(3, num_swept_params, col_idx + 2*num_swept_params);
    hold on;
    
    % Plot alpha_c (Positive/Confirmatory)
    data_ac = squeeze(fitted_params_matrix(:,2,:));
    plot_sweep_metric(x_sweep, data_ac, color_ac, generative_ref, col_idx, 'y=x_or_const', 2);

    % Plot alpha_d (Negative/Disconfirmatory)
    data_ad = squeeze(fitted_params_matrix(:,3,:));
    plot_sweep_metric(x_sweep, data_ad, color_ad, generative_ref, col_idx, 'y=x_or_const', 2);
    
    if col_idx == 1, ylabel('Fitted \alpha_{c,d}'); end
    xlabel(['Generative ' swept_param_labels{col_idx}]);
    grid on;
end 

end

%% -------------------- Helper Plotting Function --------------------
function plot_sweep_metric(x_axis, data, color, generative_ref, swept_idx, ground_truth_type, param_idx)
    % Processes mean/SEM and plots against the ground truth.
    
    mu  = mean(data, 1);
    sem = std(data, 0, 1) ./ sqrt(size(data, 1));
    
    % Ensure orientation is correct for patch
    x_axis = x_axis(:); mu = mu(:); sem = sem(:);

    % Define Ground Truth line
    if strcmp(ground_truth_type, 'y=0')
        gt_y = zeros(size(x_axis));
    else
        % If we are sweeping alpha itself (col 2), gt is the x-axis. 
        % Otherwise, it's the fixed reference value for alpha.
        if swept_idx == 2
            gt_y = x_axis;
        else
            gt_y = generative_ref(2) * ones(size(x_axis));
        end
    end
    plot(x_axis, gt_y, 'k-', 'LineWidth', 1.2, 'HandleVisibility', 'off');

    % Shaded Error Region
    patch([x_axis; flipud(x_axis)], [mu+sem; flipud(mu-sem)], color, ...
          'EdgeColor', 'none', 'FaceAlpha', 0.25, 'HandleVisibility', 'off');
    
    % Mean Line
    plot(x_axis, mu, 'Color', color, 'LineWidth', 2);
    
    % Reference line (Vertical) showing the "Standard" parameter value
    ref_val = generative_ref(swept_idx);
    yl = ylim;
    plot([ref_val, ref_val], yl, 'k:', 'LineWidth', 1.1, 'HandleVisibility', 'off');
    
    set(gca, 'Box', 'off', 'TickDir', 'out');
    xlim([min(x_axis), max(x_axis)]);
end