function [post,mcmc_diag] = post_from_gp_mcmc(th_grid, sim_model, th_mcmc, opt)
% Determine the posterior estimate from the approx. MCMC samples in the 
% GP-MH case. The purpose of this function is analogous to 
% 'post_from_gp_surrogate' used in the MH-BLFI case.

% Leave out burn-in
% However, if the amount of samples is small (which happens when this 
% function is called during the early iterations), we do not ignore burnin.
if any(isnan(th_mcmc(:)))
    error('NaN values in approx. MCMC samples.');
end
samples = th_mcmc;
[ns,d] = size(samples);
if ns > 20
    ind1 = min(opt.n-20,max(1,floor(opt.burninfrac*ns)));
    samples = samples(ind1:end,:);
end

% Assess the convergence
mcmc_diag = basic_mcmc_diagnostics(samples, opt.misc.display_type);

% Thin to final length (to reduce the save size of the samples-matrix)
samples = thin_samples(samples,opt.nfinal);

% Compute marginals in the grid using KDE
if d <= 2
    post.epost = kde_for_abc(th_grid,samples,1);
    post.samples = samples;
else
    post.epost = kde_for_abc(th_grid,samples,1)';
    post.samples = samples;
end
end

