function logb = logb_term(th_cur,th_pr,sim_model)
% Computes the logb-term for the MCMC acceptance ratio. 
% Assumes symmetric proposal q.

logr = 0;
if nargin == 3 && isfield(sim_model,'logprior_eval')
    logr = sim_model.logprior_eval(th_pr) - sim_model.logprior_eval(th_cur);
elseif nargin == 3 && isfield(sim_model,'prior_eval')
    logr = log(sim_model.prior_eval(th_pr)./sim_model.prior_eval(th_cur));
end
logq = 0; % WE ASSUME SYMMETRIC Q HERE SO THAT THE PROPOSAL DENSITY CANCELS OUT!
logb = logr + logq;
end

