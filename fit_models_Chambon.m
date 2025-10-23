function [parametersLPP,LPP] = fit_models_Chambon(data,model,fit_procedure)
% This function fits reinforcement learning models to behavioral data from
% experiments that include forced-choice (observational) trials, as in
% Chambon et al. (2020), Nat. Human Behavior.
%
% INPUTS:
%   data          - Can be a string to load a specific dataset or a cell
%                   array containing the behavioral data.
%                   * As a string, valid options are: 'C1', 'C2', 'C3'.
%                   * As a cell array, must be in the format {sta, cho, out, cou, obs}:
%                       - sta: {1 x nSubs} cell array of state vectors.
%                       - cho: {1 x nSubs} cell array of choice vectors (1 or 2).
%                       - out: {1 x nSubs} cell array of outcome vectors (-1 or 1).
%                       - cou: {1 x nSubs} cell array of counterfactual outcome vectors.
%                       - obs: {1 x nSubs} cell array of trial type vectors (1=free, 0=forced).
%
%   model         - An integer specifying the model to fit:
%                   1: RW model with separate learning rate for forced trials
%                      (beta, lr1, lr3).
%                   2: RW model with asymmetric learning for free choices and a
%                      separate learning rate for forced trials (beta, lr1, lr2, lr3).
%                   3: Model 1 with perseveration (beta, lr1, lr3, tau, phi).
%                   4: Full model with asymmetric learning, perseveration, and a
%                      separate learning rate for forced trials (beta, lr1, lr2, lr3, tau, phi).
%
%   fit_procedure - A string specifying the fitting method:
%                   'MAP': Maximum a Posteriori estimation (uses priors).
%                   'MLE': Maximum Likelihood Estimation (no priors).
%
% OUTPUTS:
%   parametersLPP - [nSubs x nParams] matrix of the fitted parameters for
%                   each subject. The column order depends on the model:
%                   * Model 1: [beta, lr1, lr3]
%                   * Model 2: [beta, lr1, lr2, lr3]
%                   * Model 3: [beta, lr1, lr3, tau, phi]
%                   * Model 4: [beta, lr1, lr2, lr3, tau, phi]
%
%   LPP           - [nSubs x 1] vector of the log posterior probability (for MAP)
%                   or log likelihood (for MLE), adjusted for Laplace approximation,
%                   for each subject's best fit.
%
% REQUIREMENTS:
%   - MATLAB Optimization Toolbox (for fmincon).
%   - MATLAB Parallel Computing Toolbox (for parfor loop).

%% INPUT VALIDATION
if ~(ischar(data) || iscell(data))
    error('fit_pers_Chambon:InvalidInput', "Input 'data' must be a string or a cell array.");
end
if iscell(data) && numel(data) ~= 5
    error('fit_pers_Chambon:InvalidInput', "If 'data' is a cell array, it must contain 5 elements: {states, choices, chosen_outcomes, unchosen_outcomes, trial_types}.");
end
if ~isnumeric(model) || ~isscalar(model) || ~ismember(model, 1:4)
    error('fit_pers_Chambon:InvalidInput', "Input 'model' must be an integer from 1 to 4.");
end
if ~ischar(fit_procedure) || ~ismember(fit_procedure, {'MAP', 'MLE'})
    error('fit_pers_Chambon:InvalidInput', "Input 'fit_procedure' must be either 'MAP' or 'MLE'.");
end
close all
% clear all
%% loading data
data_path = 'Data'; % Define the data subdirectory
if ischar(data)
    if strcmp(data,'C1')
        load(fullfile(data_path, 'C2020a.mat')); % Chambon 2020 (observational rich/poor) C1 in the paper
    elseif strcmp(data,'C2')
        load(fullfile(data_path, 'C2020b.mat')); % Chambon 2020 (observational complete) C2 in the paper
    elseif strcmp(data,'C3')
        load(fullfile(data_path, 'C2020c.mat')); % Chambon 2020 (observational symmetric/asymmetric) C3 in the paper
    else
        error('fit_pers_Chambon:InvalidDataString', 'Unknown data string provided.');
    end
elseif iscell(data) % Data provided directly
    sta=data{1};
    cho=data{2};
    out=data{3};
    cou=data{4};
    obs=data{5};
end
if ~exist('cou','var') || isempty(cou)
    for ii=1:length(out)
        cou{ii}=zeros(size(out{ii}));
    end
end
%% FIT MODEL
subjecttot=numel(sta);
nfpm=[3,4,5,6];
options = optimset('Algorithm', 'interior-point', 'Display', 'off', 'MaxIter', 10000,'MaxFunEval',10000);
parametersLPP = NaN(subjecttot, nfpm(model)); % Preallocate
LPP = NaN(subjecttot, 1); % Preallocate
check_conv = NaN(subjecttot, 1); % Preallocate

for k_sub = 1:subjecttot
        if model == 1
            lb = [0 0 0];       LB = [0 0 0];
            ub = [15 1 1];       UB = [Inf 1 1];
        elseif model == 2
            lb = [0 0 0 0];       LB = [0 0 0 0];
            ub = [15 1 1 1];       UB = [Inf 1 1 1];
        elseif model == 3
            lb = [0 0 0 0 -15];       LB = [0 0 0 0 -Inf];
            ub = [15 1 1 1 15];       UB = [Inf 1 1 1 Inf];
        elseif model == 4
            lb = [0 0 0 0 0 -15];       LB = [0 0 0 0 0 -Inf];
            ub = [15 1 1 1 1 15];       UB = [Inf 1 1 1 1 Inf];
        end

        if strcmp(fit_procedure,'MLE')
            LB=lb; UB=ub;
        end

        ddb = ub - lb;
        
        n_rep           = 5;
        parametersLPP_rep  = NaN(n_rep,nfpm(model));
        LPP_rep            = NaN(n_rep,1);
        FminHess        = NaN(n_rep,nfpm(model),nfpm(model));
       
        x0 = lb + rand(n_rep,length(lb)).*ddb;
        x0 = x0(:,1:nfpm(model));
        
        current_sta = sta{k_sub};
        current_cho = cho{k_sub};
        current_out = out{k_sub};
        current_cou = cou{k_sub};
        current_obs = obs{k_sub};

        parfor k_rep = 1:n_rep
            [parametersLPP_rep(k_rep,:),LPP_rep(k_rep),~,~,~,~,FminHess(k_rep,:,:)]=fmincon(@(x) Parameters_Priors_Complete_Final(x,current_sta,current_cho,current_out,current_cou,current_obs,model,fit_procedure),x0(k_rep,:),[],[],[],[],LB,UB,[],options);
        end
        
        [~,posLPP] = min(LPP_rep);
        parametersLPP(k_sub,1:nfpm(model)) = parametersLPP_rep(posLPP(1),1:nfpm(model));
        LPP(k_sub) = LPP_rep(posLPP(1),:) - nfpm(model)*log(2*pi)/2 + real(log(det(squeeze(FminHess(posLPP(1),:,:)))))/2;
        check_conv(k_sub) = ~any(eig(squeeze(FminHess(posLPP(1),:,:)))<0);
end
if any(~check_conv)
    warning('fit_pers_Chambon:Convergence', '%d subjects may not have converged properly (Hessian not positive definite).', sum(~check_conv));
end
end

function [post]=Parameters_Priors_Complete_Final(params,s,a,r,c,obs,model,fit)
if model == 1
    beta  = params(1); pbeta = log(gampdf(beta,1.2,5.0));
    lr1   = params(2); plr1  = log(betapdf(lr1,1.1,1.1));
    lr3   = params(3); plr3  = log(betapdf(lr3,1.1,1.1));
    p = [pbeta plr1 plr3];
elseif model == 2
    beta  = params(1); pbeta = log(gampdf(beta,1.2,5.0));
    lr1   = params(2); plr1  = log(betapdf(lr1,1.1,1.1));
    lr2   = params(3); plr2  = log(betapdf(lr2,1.1,1.1));
    lr3   = params(4); plr3  = log(betapdf(lr3,1.1,1.1));
    p = [pbeta plr1 plr2 plr3];
elseif model == 3
    beta  = params(1); pbeta = log(gampdf(beta,1.2,5.0));
    lr1   = params(2); plr1  = log(betapdf(lr1,1.1,1.1));
    lr3   = params(3); plr3  = log(betapdf(lr3,1.1,1.1));
    tau   = params(4); ptau  = log(betapdf(tau,1.1,1.1));
    phi   = params(5); pphi  = log(normpdf(phi,0,1));
    p = [pbeta plr1 plr3 ptau pphi];
elseif model == 4
    beta  = params(1); pbeta = log(gampdf(beta,1.2,5.0));
    lr1   = params(2); plr1  = log(betapdf(lr1,1.1,1.1));
    lr2   = params(3); plr2  = log(betapdf(lr2,1.1,1.1));
    lr3   = params(4); plr3  = log(betapdf(lr3,1.1,1.1));
    tau   = params(5); ptau  = log(betapdf(tau,1.1,1.1));
    phi   = params(6); pphi  = log(normpdf(phi,0,1));
    p = [pbeta plr1 plr2 plr3 ptau pphi];
end
p = -sum(p);
l=Computational_Models_Complete_Final(params,s,a,r,c,obs,model);
if strcmp(fit,'MAP')
    post = p + l;
elseif strcmp(fit,'MLE')
    post = l;
end
end

function lik = Computational_Models_Complete_Final(params,s,a,r,c,obs,model)
if isempty(s) || isempty(a) || isempty(obs)
    lik = Inf;
    return;
end
if model == 1
    beta  = params(1); lr1 = params(2); lr2 = params(2); lr3 = params(3); tau = 0; phi = 0;
elseif model == 2
    beta  = params(1); lr1 = params(2); lr2 = params(3); lr3 = params(4); tau = 0; phi = 0;
elseif model == 3
    beta  = params(1); lr1 = params(2); lr2 = params(2); lr3 = params(3); tau = params(4); phi = params(5);
elseif model == 4
    beta  = params(1); lr1 = params(2); lr2 = params(3); lr3 = params(4); tau = params(5); phi = params(6);
end

Q = zeros(max(s),2);
C = zeros(max(s),2);
lik=0;
for i = 1:length(a)
    if (a(i)==1 || a(i)==2) && ~isnan(s(i)) && s(i)>0
        if obs(i)==1
            lik = lik + log (1/(1+ exp(-beta*(Q(s(i),a(i))-Q(s(i),3-a(i))) - phi*(C(s(i),a(i))-C(s(i),3-a(i))))));
        end
        PEc =  r(i) - Q(s(i),a(i));
        if obs(i)==1
            Q(s(i),a(i)) = Q(s(i),a(i)) + lr1 * PEc * (PEc>0) +  lr2 * PEc * (PEc<0);
            C(s(i),a(i)) = C(s(i),a(i))     + tau * (1 - C(s(i),a(i)));
            C(s(i),3-a(i)) = C(s(i),3-a(i)) + tau * (0 - C(s(i),3-a(i)));
        elseif obs(i)==0
            Q(s(i),a(i)) = Q(s(i),a(i)) + lr3 * PEc;
            C(s(i),a(i)) = C(s(i),a(i))     + tau * (0 - C(s(i),a(i)));
            C(s(i),3-a(i)) = C(s(i),3-a(i)) + tau * (0 - C(s(i),3-a(i)));
        end
        if c(i)~=0
            PEu =  c(i) - Q(s(i),3-a(i));
            if obs(i)==1
                Q(s(i),3-a(i)) = Q(s(i),3-a(i)) + lr2 * PEu * (PEu>0) +  lr1 * PEu * (PEu<0);
            elseif obs(i)==0
                Q(s(i),3-a(i)) = Q(s(i),3-a(i)) + lr3 * PEu;
            end
        end
    end
end
lik = -lik;
end