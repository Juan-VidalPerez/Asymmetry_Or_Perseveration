function [sim_data] = simulate_newtask(parameters, model)
% Simulates behavioral data for a specific 4-condition task designed to
% elicit different choice patterns (e.g., LR vs WS).
%
% INPUTS:
%   parameters  - [nSubs x nParams] matrix of parameters for simulation.
%                 Parameter order depends on the model:
%                 Model 1 (CB):   [beta, lr1, lr2]
%                 Model 2 (PSL):  [beta, lr1, tau, phi] (lr2=lr1 internally)
%                 Model 3 (CBPERS):[beta, lr1, lr2, tau, phi]
%                 Model 4 (Q0):   [beta, lr1, q0] (lr2=lr1 internally)
%   model       - Integer specifying the model: 1=CB, 2=PSL, 3=CBPERS, 4=Q0.
%
% OUTPUT:
%   sim_data    - A 1x5 cell array {condition, choice, outcome, counterfactual, unchosen}.

%% 1. Configure Experiment Parameters
% This simulates a task with 4 distinct conditions (pairs of stimuli).
% Stimuli have different expected values (EV) and standard deviations (SD).
% Some conditions provide feedback, others don't.
ev = [1 1 2 2]; % Expected value for stimuli 1, 2, 3, 4
sd = [0.5 1 1 0.5]; % Standard deviation for stimuli 1, 2, 3, 4
n_trials_per_cond = 50;
pairs = [1 3; 2 4; 1 2; 3 4]; % Stimulus pairs presented in conditions 1, 2, 3, 4
feedback_given = [1 1 0 0]; % Whether feedback is provided in conditions 1, 2, 3, 4 (1=Yes, 0=No)
interleaved = 0; % 0 = Blocked presentation, 1 = Interleaved
counterfactual = 1; % Whether to calculate counterfactual outcomes (even if not shown)

n_conditions = size(pairs, 1);

%% 2. Run Simulation Loop
% Pre-allocate output cell arrays
ucho = cell(size(parameters, 1), 1); cho = cell(size(parameters, 1), 1);
con = cell(size(parameters, 1), 1); out = cell(size(parameters, 1), 1);
cou = cell(size(parameters, 1), 1);

for ss = 1:size(parameters, 1) % Loop over subjects
    params = squeeze(parameters(ss, :));
    
    % --- Unpack parameters based on model ---
    q0 = 0; % Default initial Q value
    if model == 1 % CB
        beta  = params(1); lr1 = params(2); lr2 = params(3);
        tau   = 0;         phi = 0;
    elseif model == 2 % PSL
        beta  = params(1); lr1 = params(2); lr2 = params(2); % Symmetric LR
        tau   = params(3); phi = params(4);
    elseif model == 3 % CBPERS
        beta  = params(1); lr1 = params(2); lr2 = params(3);
        tau   = params(4); phi = params(5);
    elseif model == 4 % Q0 model
        beta  = params(1); lr1 = params(2); lr2 = params(2); % Symmetric LR
        tau   = 0;         phi = 0;
        q0    = params(3); % Initial Q value is a parameter
    end
    
    % Initialize Q-values and choice traces
    Q = ones(1, 4) * q0;
    C = zeros(1, 4);
    
    % --- Determine trial order ---
    condition_seq = []; feedback_seq = [];
    if interleaved == 1
        condition_seq = repmat(1:n_conditions, [1, n_trials_per_cond]);
        feedback_seq = feedback_given(condition_seq);
        rand_order = randperm(length(condition_seq));
        condition_seq = condition_seq(rand_order);
        feedback_seq = feedback_seq(rand_order);
        n_total_trials = length(condition_seq);
    else % Blocked presentation
        for cond_idx = 1:n_conditions
            condition_seq = [condition_seq, repmat(cond_idx, [1, n_trials_per_cond])];
            feedback_seq = [feedback_seq, repmat(feedback_given(cond_idx), [1, n_trials_per_cond])];
        end
        n_total_trials = length(condition_seq);
        % Note: Original code randomized within feedback blocks, keeping it simple here.
        % If strict blocking by feedback type then randomization is needed, uncomment below:
        % condition_f1 = condition_seq(feedback_seq==1); rand_f1 = randperm(length(condition_f1));
        % condition_f0 = condition_seq(feedback_seq==0); rand_f0 = randperm(length(condition_f0));
        % condition_seq = [condition_f1(rand_f1), condition_f0(rand_f0)];
        % feedback_seq = [ones(1, length(condition_f1)), zeros(1, length(condition_f0))];
    end
    
    % Initialize arrays for subject data
    s_sub = NaN(n_total_trials, 1); a_sub = NaN(n_total_trials, 1);
    r_sub = NaN(n_total_trials, 1); c_sub = NaN(n_total_trials, 1);
    u_sub = NaN(n_total_trials, 1);
    
    % --- Trial Loop ---
    for tt = 1:n_total_trials
        current_condition = condition_seq(tt);
        provides_feedback = feedback_seq(tt);
        bandits = pairs(current_condition, :); % Get the two stimuli for this trial
        
        % 1. Make choice
        prob_choice2 = 1 / (1 + exp(-beta * (Q(bandits(2)) - Q(bandits(1))) - phi * (C(bandits(2)) - C(bandits(1)))));
        chosen_idx = 1 + (rand < prob_choice2); % Index within the pair (1 or 2)
        chosen_stim = bandits(chosen_idx);      % Actual stimulus ID (1-4)
        unchosen_stim = bandits(3 - chosen_idx);% Actual stimulus ID (1-4)
        
        % 2. Get outcome for chosen stimulus
        outcome = randn * sd(chosen_stim) + ev(chosen_stim);
        
        % 3. Update Q-values (only if feedback is provided for the condition)
        PE_chosen = outcome - Q(chosen_stim);
        if provides_feedback == 1
            Q(chosen_stim) = Q(chosen_stim) + lr1 * PE_chosen * (PE_chosen > 0) + lr2 * PE_chosen * (PE_chosen < 0);
            
            % 4. Calculate Counterfactual outcome and update if applicable
            if counterfactual == 1
                outcome_un = randn * sd(unchosen_stim) + ev(unchosen_stim);
                PE_unchosen = outcome_un - Q(unchosen_stim);
                Q(unchosen_stim) = Q(unchosen_stim) + lr2 * PE_unchosen * (PE_unchosen > 0) + lr1 * PE_unchosen * (PE_unchosen < 0);
                c_sub(tt) = outcome_un;
            else
                c_sub(tt) = 0;
            end
        else % No feedback provided for this condition
            c_sub(tt) = 0; % No counterfactual if no chosen feedback
        end
        
        % 5. Update choice traces (always updated)
        is_chosen = (1:4 == chosen_stim); % Logical vector: 1 for chosen, 0 otherwise
        C = C + tau * (is_chosen - C);
        
        % 6. Store trial data
        s_sub(tt) = current_condition;
        a_sub(tt) = chosen_stim;
        u_sub(tt) = unchosen_stim;
        r_sub(tt) = outcome;
        % c_sub is stored above
    end % End trial loop
    
    con{ss} = s_sub; cho{ss} = a_sub;
    out{ss} = r_sub; cou{ss} = c_sub;
    ucho{ss} = u_sub;
    
end % End subject loop

sim_data = {con, cho, out, cou, ucho};
end