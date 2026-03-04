function [parametersLPP, LPP] = fit_models(data, model, fit_procedure, varargin)
% FIT_MODELS Fits reinforcement learning models to behavioral data using
% Maximum a Posteriori (MAP) or Maximum Likelihood Estimation (MLE).
% It is adapted from Palminteri et al., (2023), Nat Rev Neurosci.
%
% INPUTS:
%   data          - Can be a string to load a specific dataset or a cell 
%                   array containing the behavioral data.
%                   * As a string, valid options are: 'S1a', 'S1b', 'P2', 
%                     'P1', 'L1', 'L2', 'C1', 'C2', 'C3', 'C4'.
%                   * As a cell array, must be in the format {sta, cho, out, cou}:
%                       - sta: {1 x nSubs} cell array of state vectors.
%                       - cho: {1 x nSubs} cell array of choice vectors (1 or 2).
%                       - out: {1 x nSubs} cell array of outcome vectors (-1 or 1).
%                       - cou: {1 x nSubs} cell array of counterfactual outcome vectors.
%                       OPTIONAL: if the simulation follows the structure
%                       of the experiments C1, C2 or C3 (with forced
%                       choices), one must include an extra cell {sta, cho,
%                       out, cou, obs}:
%                       - obs: {1 x nSubs} cell array of trial type vectors (1=free, 0=forced).
%
%   model         - An integer specifying the model to fit:
%                   1: Simple Rescorla-Wagner (RW) model (beta, lr1) .
%                   2: RW model with asymmetric learning rates (LA) (beta, lr1, lr2).
%                   3: RW model with perseveration (PSL) (beta, lr1, tau, phi).
%                   4: Full model with asymmetric learning and
%                   perseveration (HYBRID)
%                      (beta, lr1, lr2, tau, phi).
%                   NOTE: all fits based on experimental design C1, C2, C3
%                   also include an l3 parameter.
%
%   fit_procedure - A string specifying the fitting method:
%                   'MAP': Maximum a Posteriori estimation (uses priors).
%                   'MLE': Maximum Likelihood Estimation (no priors).
%
%   varargin{1} (optional) - Double, specifying the standard deviation of the 
%               prior for the perseveration parameter (phi) when fitting MAP 'MAP'.
%
% OUTPUTS:
%   parametersLPP - [nSubs x nParams] matrix of the fitted parameters for
%                   each subject. The column order depends on the model:
%                   * Model 1: [beta, lr1] (C1-3: [beta, lr1, lr3])
%                   * Model 2: [beta, lr1, lr2] (C1-3: [beta, lr1, lr2, lr3])
%                   * Model 3: [beta, lr1, tau, phi] (C1-3: [beta, lr1, lr3, tau, phi])
%                   * Model 4: [beta, lr1, lr2, tau, phi] (C1-3: [beta, lr1, lr2, lr3, tau, phi])
%
%   LPP           - [nSubs x 1] vector of the log posterior probability
%                   (or log likelihood adjusted for Laplace approximation)
%                   for each subject's best fit.
%
% USAGE:
%   [params, lpp] = fit_models('L1', 2, 'MAP')
%   [params, lpp] = fit_models(data_cell, 4, 'MLE')
%   [params, lpp] = fit_models('P1', 3, 'MAP', 3) % MAP with phi prior SD = 3
%
% REQUIREMENTS:
%   - MATLAB Optimization Toolbox (for fmincon).
%   - MATLAB Parallel Computing Toolbox (for parfor loop).

    %% 1. Input Validation
    if ~(ischar(data) || isstring(data) || iscell(data))
        error('fit_models:InvalidInput', "Input 'data' must be a string or a cell array.");
    end
    
    if ~isnumeric(model) || ~isscalar(model) || ~ismember(model, 1:4)
        error('fit_models:InvalidInput', "Input 'model' must be an integer from 1 to 4.");
    end
    
    if ~ischar(fit_procedure) || ~ismember(upper(fit_procedure), {'MAP', 'MLE'})
        error('fit_models:InvalidInput', "Input 'fit_procedure' must be either 'MAP' or 'MLE'.");
    end

    %% 2. Routing Logic
    % We separate data into two categories: 
    % 1. "C-design" (C1, C2, C3) which includes forced-choice observation trials.
    % 2. "Standard-design" (everything else).
    
    is_C_type = false;

    if iscell(data)
        if length(data) == 5
            is_C_type = true;
        elseif length(data) ~= 4
            error('fit_models:InvalidCell', "Cell array must contain 4 elements (standard) or 5 (C1-3 design).");
        end
    else
        % Data is a string experiment ID
        if any(strcmpi(data, {'C1', 'C2', 'C3'}))
            is_C_type = true;
        end
    end

    %% 3. Execute Fitting
    % Call the appropriate sub-function based on the design type.
    if is_C_type
        [parametersLPP, LPP] = fit_models_C123(data, model, fit_procedure, varargin{:});
    else
        [parametersLPP, LPP] = fit_models_noC123(data, model, fit_procedure, varargin{:});
    end

end