function [sim_data] = simulate_models(parameters, model, experiment, varargin)
% Simulates behavioral data for a given reinforcement learning model in a
% specified experimental paradigm.
%
% This function can simulate the default number of sessions for an experiment
% or a custom number of sessions if provided as an optional argument.
%
% USAGE:
%   sim_data = sim_pers(params, 4, 'P1');       % Simulates the default 2 sessions
%   sim_data = sim_pers(params, 4, 'P1', 10);  % Simulates 10 sessions
%
% INPUTS:
%   parameters  - [nSubs x nParams] matrix where each row contains the
%                 parameters for one simulated subject.
%   model       - An integer specifying the model to simulate:
%                 1: Simple Rescorla-Wagner (RW) model (symmetric learning).
%                 2: Contextual Bandit (CB) model (asymmetric learning).
%                 3: Perseveration (PSL) model (symmetric learning).
%                 4: Full (CBPERS) model (asymmetric learning + perseveration).
%   experiment  - A string specifying the experimental design to simulate:
%                 'P2', 'S1b', 'P1', 'S1a', 'L1', 'L2', 'C4'.
%   varargin{1} (optional) - The number of sessions to simulate. If not
%                            provided, the default for the experiment is used.
%
% OUTPUT:
%   sim_data    - A 1x4 cell array containing the simulated data, structured
%                 as {states, choices, outcomes, counterfactuals}.

%% 1. Configure Experiment Parameters
% Set up contingencies, trial structure, and default session number.
switch experiment
    case {'P2', 'S1b'} % Full feedback, reversal
        p_rew = [0.5 0.5; 0.25 0.75; 0.25 0.75; 0.17 0.83];
        default_n_sessions = 2;
        n_trials_per_block = 24;
        flip_block = 4;
        flip_trial = 13;
        has_counterfactual = 1;
        
    case {'P1', 'S1a'} % Partial feedback, reversal
        p_rew = [0.5 0.5; 0.25 0.75; 0.25 0.75; 0.17 0.83];
        default_n_sessions = 2;
        n_trials_per_block = 24;
        flip_block = 4;
        flip_trial = 13;
        has_counterfactual = 0;
        
    case {'L1', 'L2'} % Partial feedback, stable
        p_rew = [0.25 0.25; 0.25 0.75; 0.25 0.75; 0.75 0.75];
        default_n_sessions = 1;
        n_trials_per_block = 24;
        flip_block = NaN;
        flip_trial = NaN;
        has_counterfactual = 0;
        
    case 'C4' % Go/NoGo task
        p_rew = [0.5 0.5; 0.3 0.7];
        default_n_sessions = 3;
        n_trials_per_block = 100;
        flip_block = NaN;
        flip_trial = NaN;
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
% Pre-allocate output cell arrays for each subject's data
cho = cell(size(parameters, 1), 1);
con = cell(size(parameters, 1), 1);
out = cell(size(parameters, 1), 1);
cou = cell(size(parameters, 1), 1);

% Loop over each subject (each row in the parameters matrix)
for ss = 1:size(parameters, 1)
    
    % --- Unpack parameters for the current subject and model ---
    params = squeeze(parameters(ss, :));
    if model == 1 % RW model (symmetric learning, no perseveration)
        beta  = params(1);
        lr1   = params(2);
        lr2   = params(2); % Symmetric
        tau   = 0;
        phi   = 0;
    elseif model == 2 % CB model (asymmetric learning, no perseveration)
        beta  = params(1);
        lr1   = params(2);
        lr2   = params(3);
        tau   = 0;
        phi   = 0;
    elseif model == 3 % PSL model (symmetric learning, with perseveration)
        beta  = params(1);
        lr1   = params(2);
        lr2   = params(2); % Symmetric
        tau   = params(3);
        phi   = params(4);
    elseif model == 4 % Full CBPERS model (asymmetric learning, with perseveration)
        beta  = params(1);
        lr1   = params(2);
        lr2   = params(3);
        tau   = params(4);
        phi   = params(5);
    end
    
    % Initialize temporary data arrays for the current subject
    s_sub = []; a_sub = []; r_sub = []; c_sub = [];
    
    % --- Main simulation loops ---
    for session = 1:n_sessions
        for bb = 1:n_blocks
            % Initialize/reset Q-values and choice traces for each new block
            Q = [0 0];
            C = [0 0];
            p_rew_block = p_rew(bb, :);
            
            for tt = 1:n_trials_per_block
                % Check for contingency reversal within the block
                if tt == flip_trial && bb == flip_block
                    p_rew_block = fliplr(p_rew_block);
                end
                
                % 1. Make a choice using the softmax rule
                prob_choice2 = 1 / (1 + exp(-beta * (Q(2) - Q(1)) - phi * (C(2) - C(1))));
                choice = 1 + (rand < prob_choice2);
                
                % 2. Get outcome for the chosen option
                outcome = 2 * (double(rand < p_rew_block(choice)) - 0.5);
                
                % 3. Update Q-values for chosen option
                PE_chosen = outcome - Q(choice);
                Q(choice) = Q(choice) + lr1 * PE_chosen * (PE_chosen > 0) + lr2 * PE_chosen * (PE_chosen < 0);
                
                % 4. Handle counterfactual outcome and update if applicable
                if has_counterfactual
                    counterfactual = 2 * (double(rand < p_rew_block(3 - choice)) - 0.5);
                    PE_unchosen = counterfactual - Q(3 - choice);
                    Q(3 - choice) = Q(3 - choice) + lr2 * PE_unchosen * (PE_unchosen > 0) + lr1 * PE_unchosen * (PE_unchosen < 0);
                else
                    counterfactual = 0; % No feedback for unchosen
                end
                
                % 5. Update choice traces for perseveration
                C(choice)   = C(choice)   + tau * (1 - C(choice));
                C(3-choice) = C(3-choice) + tau * (0 - C(3-choice));
                
                % 6. Store trial data
                state = bb + n_blocks * (session - 1); % Unique state for each block
                s_sub = [s_sub; state];
                a_sub = [a_sub; choice];
                r_sub = [r_sub; outcome];
                c_sub = [c_sub; counterfactual];
            end
        end
    end
    
    % Store the complete data for the current subject
    con{ss} = s_sub;
    cho{ss} = a_sub;
    out{ss} = r_sub;
    cou{ss} = c_sub;
end

% Package the data for output
sim_data = {con, cho, out, cou};

end