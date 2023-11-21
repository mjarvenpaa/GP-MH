function [gp,gp_optim_opt] = reinit_gp(gp, gp_optim_opt, gp_opt, th_grid, ...
    y_tr, th_tr, sigma_tr)
% Updates the GPstuff GP structure when a new data point is added but GP 
% hyperparameters are not to be optimised. This is actually only needed
% for GP with noise model; in the case of standard GP with constant noise 
% variance, this function does nothing.

if isempty(gp)
    error('Empty GP structure provided for reinit_gp.');
end
if gp_opt.noise_model
    %[gp, gp_optim_opt] = init_gp(gp_opt, th_grid, y_tr, th_tr, sigma_tr);
    gp.lik.r = zeros(size(sigma_tr));
    gp.lik.sigma2 = sigma_tr.^2;
    gp.lik.U = ones(size(sigma_tr));
    gp.lik.ndata = length(sigma_tr);
end
end

