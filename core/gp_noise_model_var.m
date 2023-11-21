function sigma2 = gp_noise_model_var(gp, th_tr, sigma_tr, th, gp_opt)
% Returns the noise variance of the GP model at input parameter 'th'. 

if gp_opt.noise_model %strcmp(gp.lik.type, 'Gaussian-smt')
    % Special noise model
    % Estimates for the noise variance 'sigma2' at \theta could be modelled 
    % and computed from a separate GP but here we either use a given 
    % ad-hoc estimate 'gp_opt.acq_tol2' (chosen e.g. based on the 
    % assumption that the noise level is negligible although this typically
    % does not quite hold in practice) or use the nearest neighbour 
    % interpolation to estimate it.
    n_th = size(th,1);
    if gp_opt.acq_noise_model
        if isempty(th_tr) || isempty(sigma_tr) || any(~isfinite(sigma_tr))
            error('Incorrect input to gp_noise_model_var.');
        end
        sigma2 = NaN(n_th,1);
        for i = 1:n_th
            sigma2(i) = variance_estim_nn(th_tr, sigma_tr.^2, th(i,:));
        end
    else
        sigma2 = gp_opt.acq_tol2 * ones(n_th,1); % use given estimate
    end
else
    % standard GP case: constant noise
    sigma2 = gp.lik.sigma2 * ones(size(th,1),1);
end
end


function tol2 = variance_estim_nn(th_tr,sigma2_tr,th)
% Uses 1-nearest neighbour to estimate the variance at 'th', not carefully
% tested.

dists = sum((th_tr - repmat(th,size(th_tr,1),1)).^2,2);
tol2 = mean(sigma2_tr(dists==min(dists))); % if multiple equal distances, take the mean
end

