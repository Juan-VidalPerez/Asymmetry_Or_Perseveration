function [sim_data] = simulate_models(parameters, model, experiment, varargin)
% SIMULATE_MODELS Simulates behavioral data for reinforcement learning models.
%
% This function routes simulation to specific sub-functions based on the 
% experimental design (handling 'C1-C3' as Chambon-style designs).
%
% USAGE:
%   sim_data = simulate_models(params, 4, 'P1');      % Default sessions
%   sim_data = simulate_models(params, 4, 'P1', 10);  % 10 sessions
%
% INPUTS:
%   parameters  - [nSubs x nParams] matrix where each row contains the
%                 parameters for one simulated subject.
%   model       - An integer specifying the model to simulate:
%                 1: Simple Rescorla-Wagner (RW) model (symmetric learning).
%                 2: Learning assymetry (LA) model (asymmetric learning).
%                 3: Perseveration (PSL) model (symmetric learning).
%                 4: Full (HYBRID) model (asymmetric learning + perseveration).
%   experiment  - A string specifying the experimental design to simulate:
%                 'P2', 'S1b', 'P1', 'S1a', 'L1', 'L2', 'C1','C2','C3','C4'.
%   varargin{1} (optional) - Double, the number of sessions to simulate. 
%                            If omitted, the experiment default is used.
%
% OUTPUT:
%   sim_data    - A cell array containing the simulated data:
%                 * Standard: 1x4 {states, choices, outcomes, counterfactuals}
%                 * C1, C2, C3: 1x5 {states, choices, outcomes, counterfactuals, free}

    %% 1. Validation
    if ~isnumeric(model) || ~isscalar(model) || ~ismember(model, 1:4)
        error('simulate_models:InvalidModel', 'Model must be an integer between 1 and 4.');
    end

    if ~(ischar(experiment) || isstring(experiment))
        error('simulate_models:InvalidExp', 'Experiment must be a string or character array.');
    end

    %% 2. Routing Logic
    % Experiments C1, C2, and C3 utilize the Chambon-style design 
    % (forced/free choice paradigms).
    
    is_chambon = any(strcmpi(experiment, {'C1', 'C2', 'C3'}));

    if is_chambon
        % Output: {sta, cho, out, cou, free}
        [sim_data] = simulate_models_C123(parameters, model, experiment, varargin{:});
    else
        % Output: {sta, cho, out, cou}
        [sim_data] = simulate_models_noC123(parameters, model, experiment, varargin{:});
    end

end