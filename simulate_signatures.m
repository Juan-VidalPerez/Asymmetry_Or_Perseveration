function [sim_data] = simulate_signatures(parameters, model, experiment)
% Simulates behavioral data for specific paradigms used to generate characteristic
% learning signatures (e.g., win-stay/lose-shift, risk preference).
%
% INPUTS:
%   parameters  - [nSubs x nParams] matrix of parameters for simulation.
%                 Parameter order depends on the model:
%                 Model 1 (CB):   [beta, lr1, lr2]
%                 Model 2 (PSL):  [beta, lr1, tau, phi] (lr2=lr1 internally)
%                 Model 3 (CBPERS):[beta, lr1, lr2, tau, phi]
%   model       - Integer specifying the model: 1=CB, 2=PSL, 3=CBPERS.
%   experiment  - String specifying the paradigm:
%                 'Stable':     50/50 rewards, 24 trials.
%                 'Reversal':   80/20 rewards, reversal at trial 13/24.
%                 'Risk':       Sure reward (0) vs. Risky (50% +/-1), 24 trials.
%                 'Classic':    75/25 rewards, 50 trials.
%                 'Katahira':   70/30 rewards, multiple reversals, 200 trials.
%                 'Reversal50': 80/20 rewards, reversal at trial 26/50.
%                 'HighTask':   90/60 rewards, 50 trials.
%                 'LowTask':    40/10 rewards, 50 trials.
%
% OUTPUT:
%   sim_data    - A 1x4 cell array {states, choices, outcomes, counterfactuals}.

%% 1. Configure Experiment Parameters
counterfactual = 0; % All these paradigms are partial feedback
n_sessions = 1;     % All are single-session simulations

switch experiment
    case 'Stable'
        p_rew = [0.5 0.5];      mag_rew = [1 1];
        n_trials = 24;
        flip_block = NaN;       flip_trial = NaN;
    case 'Classic' % formerly Stable2
        p_rew = [0.25 0.75];    mag_rew = [1 1];
        n_trials = 50;
        flip_block = NaN;       flip_trial = NaN;
    case 'Reversal'
        p_rew = [0.2 0.8];      mag_rew = [1 1];
        n_trials = 24;
        flip_block = 1;         flip_trial = 13;
    case 'Katahira' % formerly Reversal2
        p_rew = [0.3 0.7];      mag_rew = [1 1];
        n_trials = 200;
        flip_block = 1;         flip_trial = [51, 101, 151];
     case 'Reversal50' % formerly Reversal3
        p_rew = [0.2 0.8];      mag_rew = [1 1];
        n_trials = 50;
        flip_block = 1;         flip_trial = [26];
    case 'Risk'
        p_rew = [1 0.5];        mag_rew = [0 1];
        n_trials = 24;
        flip_block = NaN;       flip_trial = NaN;
    case 'HighTask' % New
        p_rew = [0.6 0.9];      mag_rew = [1 1];
        n_trials = 50;
        flip_block = NaN;       flip_trial = NaN;
    case 'LowTask' % New
        p_rew = [0.1 0.4];      mag_rew = [1 1];
        n_trials = 50;
        flip_block = NaN;       flip_trial = NaN;
    otherwise
        error('Unknown experiment type: %s', experiment);
end

n_blocks = size(p_rew, 1); % Usually just 1 block for these signatures

%% 2. Run Simulation Loop
cho = cell(size(parameters, 1), 1); con = cell(size(parameters, 1), 1);
out = cell(size(parameters, 1), 1); cou = cell(size(parameters, 1), 1);

for ss = 1:size(parameters, 1) % Loop over subjects/parameter sets
    params = squeeze(parameters(ss, :));
    
    % --- Unpack parameters based on model ---
    if model == 1 % CB (asymmetric learning, no perseveration)
        beta  = params(1); lr1 = params(2); lr2 = params(3);
        tau   = 0;         phi = 0;
    elseif model == 2 % PSL (symmetric learning, with perseveration)
        beta  = params(1); lr1 = params(2); lr2 = params(2); % Symmetric LR
        tau   = params(3); phi = params(4);
    elseif model == 3 % CBPERS (asymmetric learning, with perseveration)
        beta  = params(1); lr1 = params(2); lr2 = params(3);
        tau   = params(4); phi = params(5);
    end
    
    s_sub = []; a_sub = []; r_sub = []; c_sub = []; % Temp arrays for subject
    
    for session = 1:n_sessions % This loop will only run once
        for bb = 1:n_blocks
            Q = [0 0]; C = [0 0]; % Reset for block
            p_rew_block = p_rew(bb, :);
            mag_rew_block = mag_rew;
            
            for tt = 1:n_trials
                % Check for reversal
                if ismember(tt, flip_trial) && bb == flip_block
                    p_rew_block = fliplr(p_rew_block);
                end
                
                % 1. Make choice
                prob_choice2 = 1 / (1 + exp(-beta * (Q(2) - Q(1)) - phi * (C(2) - C(1))));
                choice = 1 + (rand < prob_choice2);
                
                % 2. Get outcome
                rewarded = double(rand < p_rew_block(choice));
                outcome = (2 * rewarded - 1) * mag_rew_block(choice);
                
                % 3. Update Q-values
                PE_chosen = outcome - Q(choice);
                Q(choice) = Q(choice) + lr1 * PE_chosen * (PE_chosen > 0) + lr2 * PE_chosen * (PE_chosen < 0);
                
                % 4. Counterfactual (always 0 here)
                counterfactual = 0;
                
                % 5. Update choice traces
                C(choice)   = C(choice)   + tau * (1 - C(choice));
                C(3-choice) = C(3-choice) + tau * (0 - C(3-choice));
                
                % 6. Store trial data
                state = bb + n_blocks * (session - 1); % State is just the block number here
                s_sub = [s_sub; state]; a_sub = [a_sub; choice];
                r_sub = [r_sub; outcome]; c_sub = [c_sub; counterfactual];
            end % End trial loop
        end % End block loop
    end % End session loop
    
    con{ss} = s_sub; cho{ss} = a_sub;
    out{ss} = r_sub; cou{ss} = c_sub;
end % End subject loop

sim_data = {con, cho, out, cou};
end