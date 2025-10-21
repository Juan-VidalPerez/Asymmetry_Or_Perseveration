function [sim_data] = simulate_models_Chambon(parameters, model, experiment, varargin)
% Simulates behavioral data for models designed for observational learning
% paradigms, such as those in Chambon et al. (2020).
%
% This function can simulate the default number of sessions for an experiment
% or a custom number of sessions if provided as an optional argument.
%
% USAGE:
%   sim_data = sim_pers_Chambon(params, 4, 'C1');      % Simulates the default 3 sessions
%   sim_data = sim_pers_Chambon(params, 4, 'C1', 5);   % Simulates 5 sessions
%
% INPUTS:
%   parameters  - [nSubs x nParams] matrix where each row contains the
%                 parameters for one simulated subject.
%   model       - An integer specifying the model to simulate:
%                 1: RW with symmetric learning for free trials and a separate
%                    learning rate (lr3) for forced trials.
%                 2: RW with asymmetric learning (lr1, lr2) for free trials
%                    and a separate learning rate (lr3) for forced trials.
%                 3: Model 1 with an added perseveration component (tau, phi).
%                 4: Full model (Model 2 with perseveration).
%   experiment  - A string specifying the experimental design: 'C1', 'C2', 'C3'.
%   varargin{1} (optional) - The number of sessions to simulate. If not
%                            provided, the default for the experiment is used.
%
% OUTPUT:
%   sim_data    - A 1x5 cell array containing the simulated data, structured
%                 as {states, choices, outcomes, counterfactuals, trial_type}.

%% 1. Configure Experiment Parameters
% Set up contingencies and trial structure based on the experiment string.
switch experiment
    case 'C1' % Observational rich/poor
        p_rew = [0.6 0.9; 0.6 0.9; 0.1 0.4; 0.1 0.4];
        % Block types: 0 = all free choice, 1 = intermixed free/forced
        observational_blocks = [0 1 0 1];
        default_n_sessions = 3;
        n_trials_per_block = 40;
        has_counterfactual = 0;
        
    case 'C2' % Observational complete
        p_rew = [0.6 0.9; 0.6 0.9; 0.1 0.4; 0.1 0.4];
        observational_blocks = [1 1 1 1]; % All blocks are intermixed
        default_n_sessions = 2;
        n_trials_per_block = 20;
        has_counterfactual = 1;
        
    case 'C3' % Observational symmetric/asymmetric
        p_rew = [0.5 0.5; 0.3 0.7; 0.5 0.5; 0.3 0.7];
        observational_blocks = [0 1 0 1];
        default_n_sessions = 3;
        n_trials_per_block = 20;
        has_counterfactual = 0;
end

% Set the number of sessions: use the default unless an override is provided.
if ~isempty(varargin)
    n_sessions = varargin{1};
else
    n_sessions = default_n_sessions;
end

n_blocks = size(p_rew, 1);

%% 2. Run Simulation Loop
% Pre-allocate output cell arrays
cho = cell(size(parameters, 1), 1); con = cell(size(parameters, 1), 1);
out = cell(size(parameters, 1), 1); cou = cell(size(parameters, 1), 1);
obs = cell(size(parameters, 1), 1);

% Loop over each subject (each row in the parameters matrix)
for ss = 1:size(parameters, 1)
    
    % --- Unpack parameters for the current subject and model ---
    params = squeeze(parameters(ss, :));
    if model == 1
        beta  = params(1); lr1 = params(2); lr2 = params(2); lr3 = params(3);
        tau = 0; phi = 0;
    elseif model == 2
        beta  = params(1); lr1 = params(2); lr2 = params(3); lr3 = params(4);
        tau = 0; phi = 0;
    elseif model == 3
        beta  = params(1); lr1 = params(2); lr2 = params(2); lr3 = params(3);
        tau = params(4); phi = params(5);
    elseif model == 4
        beta  = params(1); lr1 = params(2); lr2 = params(3); lr3 = params(4);
        tau = params(5); phi = params(6);
    end
    
    % Initialize temporary data arrays for the current subject
    s_sub = []; a_sub = []; r_sub = []; c_sub = []; o_sub = [];
    
    % --- Main simulation loops ---
    for session = 1:n_sessions
        for bb = 1:n_blocks
            Q = [0 0]; C = [0 0]; % Reset values for each block
            p_rew_block = p_rew(bb, :);
            
            % Determine trial types for the block (free, forced, or intermixed)
            if observational_blocks(bb) == 0 % All trials are free choice
                trial_types = ones(1, n_trials_per_block);
            else % Intermixed block: half free, half forced, randomized
                trial_types = repelem([1 0], n_trials_per_block / 2);
                trial_types = trial_types(randperm(length(trial_types)));
            end

            for tt = 1:n_trials_per_block
                trial_is_free = trial_types(tt); % 1 = free, 0 = forced
                
                % 1. Make a choice
                if trial_is_free
                    % Softmax for free-choice trials
                    prob_choice2 = 1 / (1 + exp(-beta * (Q(2) - Q(1)) - phi * (C(2) - C(1))));
                else
                    % Random choice for forced-choice trials
                    prob_choice2 = 0.5;
                end
                choice = 1 + (rand < prob_choice2);
                
                % 2. Get outcome
                outcome = 2 * (double(rand < p_rew_block(choice)) - 0.5);
                PE_chosen = outcome - Q(choice);
                
                % 3. Update values and traces based on trial type
                if trial_is_free
                    % Update Q-values with free-choice learning rates
                    Q(choice) = Q(choice) + lr1 * PE_chosen * (PE_chosen > 0) + lr2 * PE_chosen * (PE_chosen < 0);
                    % Update choice trace: strengthen chosen, weaken unchosen
                    C(choice)   = C(choice)   + tau * (1 - C(choice));
                    C(3-choice) = C(3-choice) + tau * (0 - C(3-choice));
                else % Forced-choice trial
                    % Update Q-values with the observational learning rate
                    Q(choice) = Q(choice) + lr3 * PE_chosen;
                    % Update choice trace: decay both as choice was not self-generated
                    C(choice)   = C(choice)   + tau * (0 - C(choice));
                    C(3-choice) = C(3-choice) + tau * (0 - C(3-choice));
                end
                
                % 4. Handle counterfactual outcome
                if has_counterfactual
                    counterfactual = 2 * (double(rand < p_rew_block(3 - choice)) - 0.5);
                    PE_unchosen = counterfactual - Q(3 - choice);
                    if trial_is_free
                        Q(3 - choice) = Q(3 - choice) + lr2 * PE_unchosen * (PE_unchosen > 0) + lr1 * PE_unchosen * (PE_unchosen < 0);
                    else
                        Q(3 - choice) = Q(3 - choice) + lr3 * PE_unchosen;
                    end
                else
                    counterfactual = 0;
                end
                
                % 5. Store trial data
                state = bb + n_blocks * (session - 1);
                s_sub = [s_sub; state]; a_sub = [a_sub; choice];
                r_sub = [r_sub; outcome]; c_sub = [c_sub; counterfactual];
                o_sub = [o_sub; trial_is_free];
            end
        end
    end
    
    % Store the complete data for the current subject
    con{ss} = s_sub; cho{ss} = a_sub;
    out{ss} = r_sub; cou{ss} = c_sub;
    obs{ss} = o_sub;
end

% Package the data for output
sim_data = {con, cho, out, cou, obs};
end