function [gp, gp_optim_opt] = init_gp(gp_opt, th_grid, y_tr, th_tr, sigma_tr)
% Sets up initial GP structure and its default init values. This function
% however does not optimise the GP hyperparameters. 
%
% We set up rather noninformative priors for the GP model. The lengthscale is set 
% according to the bounds of the parameter space. The signal variance is set acording to
% 'typical' magnitudes of loglik. The noise variance is given a very wide prior (needed 
% only when gp_opt.noise_model is not used). 

% TODO: better default init values?
th_range = th_grid.range(:,2)-th_grid.range(:,1);

%% set up GP model etc.
pl = prior_t('s2',(th_range(1)/2)^2); % we suppose parameters are scaled similarly
ydata_var = 1000^2;
%ydata_var = 100^2;
%pm = prior_unif();
pm = prior_sqrtt('s2',ydata_var);
%pn = prior_unif();
pn = prior_sqrtt('s2',50^2);

gpcf = gpcf_sexp('lengthScale', th_range(:)'/3, 'magnSigma2', ydata_var, ...
    'lengthScale_prior', pl, 'magnSigma2_prior', pm);

if gp_opt.noise_model
    n_tr = length(y_tr);
    lik = lik_gaussiansmt('ndata', n_tr, 'sigma2', sigma_tr.^2);
    % Note: there is no sigma_n parameter to infer!
else
    sigma_n = 10;
    lik = lik_gaussian('sigma2', sigma_n^2, 'sigma2_prior', pn);
    %lik = lik_gaussian('sigma2', sigma_n^2, 'sigma2_prior', prior_fixed()); % sigma fixed
end

if gp_opt.meanf
    % set different mean functions
    gpmf1 = gpmf_constant('prior_mean',0,'prior_cov',30^2);
    gpmf2 = gpmf_linear('prior_mean',0,'prior_cov',30^2);
    gpmf3 = gpmf_squared('prior_mean',0,'prior_cov',30^2);
    
    gp = gp_set('lik', lik, 'cf', gpcf, 'meanf', {gpmf1,gpmf2,gpmf3}, 'jitterSigma2', 1e-9);
    %gp = gp_set('lik', lik, 'cf', gpcf, 'meanf', {gpmf1}, 'jitterSigma2', 1e-9); % only const term
else
    % use the zero mean function
    gp = gp_set('lik', lik, 'cf', gpcf, 'jitterSigma2', 1e-9);
end
gp_optim_opt = optimset('TolFun',1e-3,'TolX',1e-3','display', 'off');

end

