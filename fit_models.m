function [parametersLPP,LPP] = fit_models(data,model,fit_procedure)
% This function fits reinforcement learning models to behavioral data using
% Maximum a Posteriori (MAP) or Maximum Likelihood Estimation (MLE).
% It is adapted from Palminteri et al., (2023), Nat Rev Neurosci.
%
% INPUTS:
%   data          - Can be a string to load a specific dataset or a cell 
%                   array containing the behavioral data.
%                   * As a string, valid options are: 'S1a', 'S1b', 'P2', 
%                     'P1', 'L1', 'L2', 'C4'.
%                   * As a cell array, must be in the format {sta, cho, out, cou}:
%                       - sta: {1 x nSubs} cell array of state vectors.
%                       - cho: {1 x nSubs} cell array of choice vectors (1 or 2).
%                       - out: {1 x nSubs} cell array of outcome vectors (-1 or 1).
%                       - cou: {1 x nSubs} cell array of counterfactual outcome vectors.
%
%   model         - An integer specifying the model to fit:
%                   1: Simple Rescorla-Wagner (RW) model (beta, lr1).
%                   2: RW model with asymmetric learning rates (beta, lr1, lr2).
%                   3: RW model with perseveration (beta, lr1, tau, phi).
%                   4: Full model with asymmetric learning and perseveration
%                      (beta, lr1, lr2, tau, phi).
%
%   fit_procedure - A string specifying the fitting method:
%                   'MAP': Maximum a Posteriori estimation (uses priors).
%                   'MLE': Maximum Likelihood Estimation (no priors).
%
% OUTPUTS:
%   parametersLPP - [nSubs x nParams] matrix of the fitted parameters for
%                   each subject. The column order depends on the model:
%                   * Model 1: [beta, lr1]
%                   * Model 2: [beta, lr1, lr2]
%                   * Model 3: [beta, lr1, tau, phi]
%                   * Model 4: [beta, lr1, lr2, tau, phi]
%
%   LPP           - [nSubs x 1] vector of the log posterior probability
%                   (or log likelihood adjusted for Laplace approximation)
%                   for each subject's best fit.
%
% REQUIREMENTS:
%   - MATLAB Optimization Toolbox (for fmincon).
%   - MATLAB Parallel Computing Toolbox (for parfor loop).

%% INPUT VALIDATION
if ~(ischar(data) || iscell(data))
    error('fit_pers:InvalidInput', "Input 'data' must be a string or a cell array.");
end
if iscell(data) && numel(data) ~= 4
    error('fit_pers:InvalidInput', "If 'data' is a cell array, it must contain 4 elements: {states, choices, chosen_outcomes, unchosen_outcomes}.");
end
if ~isnumeric(model) || ~isscalar(model) || ~ismember(model, 1:4)
    error('fit_pers:InvalidInput', "Input 'model' must be an integer from 1 to 4.");
end
if ~ischar(fit_procedure) || ~ismember(fit_procedure, {'MAP', 'MLE'})
    error('fit_pers:InvalidInput', "Input 'fit_procedure' must be either 'MAP' or 'MLE'.");
end
%% LOADING DATA
data_path = 'Data'; % Define the data subdirectory
if ischar(data)
    if strcmp(data,'S1a')
        load(fullfile(data_path, 'S2021.mat'));
        for ii=1:length(out)
                sta{ii}=sta{ii}(cou{ii}==0);
                cho{ii}=cho{ii}(cou{ii}==0);
                out{ii}=out{ii}(cou{ii}==0);
                cou{ii}=cou{ii}(cou{ii}==0);
        end
    elseif strcmp(data,'S1b')
        load(fullfile(data_path, 'S2021.mat'));
        for ii=1:length(out)
            sta{ii}=sta{ii}(cou{ii}~=0);
            cho{ii}=cho{ii}(cou{ii}~=0);
            out{ii}=out{ii}(cou{ii}~=0);
            cou{ii}=cou{ii}(cou{ii}~=0);
        end
    elseif strcmp(data,'P2')
        load(fullfile(data_path, 'P2017b.mat'));
    elseif strcmp(data,'P1')
        load(fullfile(data_path, 'P2017a.mat'));
    elseif strcmp(data,'L1')
        load(fullfile(data_path, 'L2017a.mat'));
    elseif strcmp(data,'L2')
        load(fullfile(data_path, 'L2017b.mat'));
    elseif strcmp(data,'C4')
        load(fullfile(data_path, 'C2020d.mat'));
    else
        error('fit_pers:InvalidDataString', 'Unknown data string provided.');
    end
elseif iscell(data) % Data provided directly
    sta=data{1};
    cho=data{2};
    out=data{3};
    cou=data{4};
end
if ~exist('cou','var') || isempty(cou)
    for ii=1:length(out)
        cou{ii}=zeros(size(out{ii}));
    end
end
%% FIT THE MODEL
subjecttot=numel(sta);
nfpm=[2,3,4 5];
options = optimset('Algorithm', 'interior-point', 'Display', 'off', 'MaxIter', 10000,'MaxFunEval',10000);
parametersLPP = NaN(subjecttot, nfpm(model));
LPP = NaN(subjecttot, 1);
check_conv = NaN(subjecttot, 1);

for k_sub = 1:subjecttot
        if model == 1
            lb = [0 0];       LB = [0 0];
            ub = [15 1];       UB = [15 1];
        elseif model == 2
            lb = [0 0 0];       LB = [0 0 0];
            ub = [15 1 1];       UB = [15 1 1];
        elseif model == 3
            lb = [0 0 0 -1];       LB = [0 0 0 -Inf];
            ub = [15 1 1 1];       UB = [Inf 1 1 Inf];
        elseif model == 4
            lb = [0 0 0 0 -1];       LB = [0 0 0 0 -Inf];
            ub = [15 1 1 1 1];       UB = [Inf 1 1 1 Inf];
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

        parfor k_rep = 1:n_rep
            [parametersLPP_rep(k_rep,:),LPP_rep(k_rep),~,~,~,~,FminHess(k_rep,:,:)]=fmincon(@(x) Parameters_Priors_Complete_Final(x,current_sta,current_cho,current_out,current_cou,model,fit_procedure),x0(k_rep,:),[],[],[],[],LB,UB,[],options);
        end
        
        [~,posLPP] = min(LPP_rep);
        parametersLPP(k_sub,1:nfpm(model)) = parametersLPP_rep(posLPP(1),1:nfpm(model));
        LPP(k_sub) = LPP_rep(posLPP(1),:) - nfpm(model)*log(2*pi)/2 + real(log(det(squeeze(FminHess(posLPP(1),:,:)))))/2;
        check_conv(k_sub) = ~any(eig(squeeze(FminHess(posLPP(1),:,:)))<0);
end
if any(~check_conv)
    warning('fit_pers:Convergence', '%d subjects may not have converged properly (Hessian not positive definite).', sum(~check_conv));
end
end

function [post]=Parameters_Priors_Complete_Final(params,s,a,r,c,model,fit)
if model == 1
    beta  = params(1); pbeta = log(gampdf(beta,1.2,5.0));
    lr1   = params(2); plr1  = log(betapdf(lr1,1.1,1.1));
    p = [pbeta plr1];
elseif model == 2
    beta  = params(1); pbeta = log(gampdf(beta,1.2,5.0));
    lr1   = params(2); plr1  = log(betapdf(lr1,1.1,1.1));
    lr2   = params(3); plr2  = log(betapdf(lr2,1.1,1.1));
    p = [pbeta plr1 plr2];
elseif model == 3
    beta  = params(1); pbeta = log(gampdf(beta,1.2,5.0));
    lr1   = params(2); plr1  = log(betapdf(lr1,1.1,1.1));
    tau   = params(3); ptau  = log(betapdf(tau,1.1,1.1));
    phi   = params(4); pphi  = log(normpdf(phi,0,1));
    p = [pbeta plr1 ptau pphi];
elseif model == 4
    beta  = params(1); pbeta = log(gampdf(beta,1.2,5.0));
    lr1   = params(2); plr1  = log(betapdf(lr1,1.1,1.1));
    lr2   = params(3); plr2  = log(betapdf(lr2,1.1,1.1));
    tau   = params(4); ptau  = log(betapdf(tau,1.1,1.1));
    phi   = params(5); pphi  = log(normpdf(phi,0,1));
    p = [pbeta plr1 plr2 ptau pphi];
end
p = -sum(p);
l=Computational_Models_Complete_Final(params,s,a,r,c,model);
if strcmp(fit,'MAP')
    post = p + l;
elseif strcmp(fit,'MLE')
    post = l;
end
end

function lik = Computational_Models_Complete_Final(params,s,a,r,c,model)
if isempty(s) || isempty(a)
    lik = Inf;
    return;
end
if model == 1
    beta  = params(1); lr1 = params(2); lr2 = params(2); tau = 0; phi = 0;
elseif model == 2
    beta  = params(1); lr1 = params(2); lr2 = params(3); tau = 0; phi = 0;
elseif model == 3
    beta  = params(1); lr1 = params(2); lr2 = params(2); tau = params(3); phi = params(4);
elseif model == 4
    beta  = params(1); lr1 = params(2); lr2 = params(3); tau = params(4); phi = params(5);
end

Q = zeros(max(s),2);
C = zeros(max(s),2);
lik=0;
for i = 1:length(a)
    if (a(i)==1 || a(i)==2) && ~isnan(s(i)) && s(i)>0 % Added check for s(i)>0
        lik = lik + log (1/(1+ exp(-beta*(Q(s(i),a(i))-Q(s(i),3-a(i))) - phi*(C(s(i),a(i))-C(s(i),3-a(i))))));
        PEc =  r(i) - Q(s(i),a(i));
        Q(s(i),a(i)) = Q(s(i),a(i)) + lr1 * PEc * (PEc>0) +  lr2 * PEc * (PEc<0);
        if c(i)~=0
            PEu =  c(i) - Q(s(i),3-a(i));
            Q(s(i),3-a(i)) = Q(s(i),3-a(i)) + lr2 * PEu * (PEu>0) +  lr1 * PEu * (PEu<0);
        end
        C(s(i),a(i)) = C(s(i),a(i))     + tau * (1 - C(s(i),a(i)));
        C(s(i),3-a(i)) = C(s(i),3-a(i)) + tau * (0 - C(s(i),3-a(i)));
    end
end
lik = -lik;
end