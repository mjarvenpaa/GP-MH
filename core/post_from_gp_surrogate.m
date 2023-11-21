function [post, mcmc_diag] = post_from_gp_surrogate(th_grid, sim_model,...
    gp, gp_opt, log_lik_tr, th_tr, sigma_tr, P, mcmc_opt, full_output)
% Computes the estimate of the (unnormalised) posterior and various related
% quantities based on the fitted GP surrogate model in the MH-BLFI case.
%
% If dim <= 2, the computations are done in a grid, otherwise sampling is
% used.
% TODO: Could compute logprior directly instead of log(prior) which might 
% underflow to -Inf in some cases, though this does not matter in the 
% current test cases because uniform priors are always specified.

if nargin < 10
    full_output = 1;
end
mcmc_diag = [];
d = th_grid.dim;
if d <= 2
    if d == 1
        th_gr = th_grid.theta(:);
    else
        th_gr = th_grid.theta2d';
    end
    [eft,varft] = gp_pred_fast(gp,th_tr,log_lik_tr,th_gr,P);
    
    % loglik etc.
    if full_output
        post.eloglik = eft; % point estimate for loglik; GP mean
        post.varloglik = varft; % var of loglik
        post.varnloglik = varft + gp_noise_model_var(gp, th_tr, sigma_tr, th_gr, gp_opt);
        
        if any(post.varloglik < 0)
            error('Negative GP variance encountered.'); % should not happen anymore
        end
        
        % loglik CI (computed from  Gaussian densities)
        post.loglik_lb = eft - 1.96*sqrt(varft); % loglik (latent) uncertainty
        post.loglik_ub = eft + 1.96*sqrt(varft);
        post.nloglik_lb = eft - 1.96*sqrt(post.varnloglik); % new loglik measurement uncertainty
        post.nloglik_ub = eft + 1.96*sqrt(post.varnloglik);
    end
    
    % posterior
    pr_val = sim_model.prior_eval(th_gr);
    log_pr_val = log(pr_val);
    log_meanpost = log_pr_val + (eft + 0.5*varft);
    log_medpost = log_pr_val + (eft);
    log_modepost = log_pr_val + (eft - varft);
    % NOTE1: We scale these according to the main point estimator; this
    % ensures that all estimators are in the same scale and that the main 
    % estimator is computed robustly but this does not guarantee that 
    % under- or overflows were impossible for the other estimators.
    % NOTE2: Calling 'exponentiate_log_array' corresponds scaling by exp(-c)
    %[epost,c] = exponentiate_log_array(log_meanpost); % MEAN AS A POINT ESTIMATOR
    %[epost,c] = exponentiate_log_array(log_medpost); % MEDIAN AS A POINT ESTIMATOR
    [epost,c] = exponentiate_log_array(log_modepost); % MODE AS A POINT ESTIMATOR
    post.epost = epost;
    post.samples = [];
    
    % extra
    if full_output
        post.meanpost = exp(log_meanpost - c); % could still overflow
        post.medpost = exp(log_medpost - c);
        post.modepost = exp(log_modepost - c);
        
        %log_qpost = log_pr_val - sqrt(varft); % some lower quantile, TBD
        log_varpost = real(2*log_pr_val + 2*(eft + varft) + (log1p(-exp(-varft))));
        post.varpost = exp(log_varpost - 2*c);
        
        % posterior CI (computed from lognormal density)
        post.post_lb = exp(log_pr_val + post.loglik_lb - c);
        post.post_ub = exp(log_pr_val + post.loglik_ub - c);
    end
    
else % dim > 2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% compute 'slices' at the true parameter - these are used for visualisation only
    if isfield(sim_model,'true_theta') %&& full_output
        gl = size(th_grid.theta,2);
        post.slice.eloglik = NaN(d,gl);
        post.slice.varloglik = NaN(d,gl);
        post.slice.epost = NaN(d,gl);
        for i = 1:d
            th_gri = th_grid.theta;
            ii = setdiff(1:d,i);
            th_gri(ii,:) = repmat(sim_model.true_theta(ii),gl,1)';
            [efti,varfti] = gp_pred_fast(gp,th_tr,log_lik_tr,th_gri',P);
            
            % loglik (slice)
            post.slice.eloglik(i,:) = efti;
            post.slice.varloglik(i,:) = varfti; 
            
            % posterior (slice)
            pr_vali = sim_model.prior_eval(th_gri');
            log_pr_vali = log(pr_vali);
            %log_eposti = log_pr_vali + (efti + 0.5*varfti); % MEAN AS A POINT ESTIMATOR
            %log_eposti = log_pr_vali + (efti); % MEDIAN AS A POINT ESTIMATOR
            log_eposti = log_pr_vali + (efti - varfti); % MODE AS A POINT ESTIMATOR
            post.slice.epost(i,:) = exponentiate_log_array(log_eposti);
        end
    end
    
    %% sampling from the point estimate of posterior
    % If dim > 2, we sample from the estimated posterior density using MCMC and compute 
    % the marginals from the resulting samples using KDE.
    % Here we use the adaptive MCMC (DRAM) by Haario et al. 2006 but many other MCMC
    % algorithms could be alternatively used.
    
    display_type = mcmc_opt.display_type;
    
    % Set up variables and options for MCMC (some of these could be varied according to each particular problem)
    npar = sim_model.dim;
    model.ssfun = @(th,data) -2*log_post_point_estim(th,th_tr,log_lik_tr,gp,P,sim_model);
    data = [];
    model.N = 1;
    
    % Take the maximum posterior value evaluated at training data points as the initial point
    f_all = model.ssfun(th_tr,[]);
    [f_opt,opt_ind] = min(f_all); % min because -2*log(.)
    init = th_tr(opt_ind,:);
    
    if ~strcmp(display_type,'off')
        disp('Initial point for MCMC:');
        init
    end
    
    params = cell(npar,1);
    for i = 1:npar
        params{i} = {sprintf('\\theta_{%d}',i), init(i), th_grid.theta(i,1), th_grid.theta(i,end)};
    end
    
    % Additional MCMC settings
    options.nsimu = mcmc_opt.nsimu;
    options.qcov = 1/10^2*diag((th_grid.theta(:,1)-th_grid.theta(:,end)).^2);
    options.method = 'am';
    options.updatesigma = 0;
    options.verbosity = ~strcmp(display_type,'off'); % no printing from mcmc
    options.waitbar = 0;

    % Initialize results
    samples_all = NaN(mcmc_opt.nsimu,npar,mcmc_opt.nchains);
    results_all = cell(mcmc_opt.nchains,1);
    
    % Run MCMC chains!
    for i = 1:mcmc_opt.nchains
        if ~strcmp(display_type,'off')
            if i == 1
                disp('Running MCMC in the second stage of MH-BLFI...');
            end
            if mcmc_opt.nchains > 1
                disp(['Chain ', num2str(i), '/', num2str(mcmc_opt.nchains)]);
            end
        end
        [results,samples] = mcmcrun(model,data,params,options);
        results_all{i} = results;
        samples_all(:,:,i) = samples;
        if i == mcmc_opt.nchains && ~strcmp(display_type,'off')
            disp('Done.');
        end
    end
    
    % Leave out burn-in (the first half of each chain)
    cl = size(samples_all,1);
    samples_all = samples_all(ceil(cl/2:cl),:,:);
    
    % Assess the convergence
    mcmc_diag = basic_mcmc_diagnostics(samples,display_type);
    
    % Add the chains together
    cl = size(samples_all,1);
    samples = NaN(mcmc_opt.nchains*cl,npar);
    for i = 1:mcmc_opt.nchains
        samples(1 + (i-1)*cl:i*cl,:) = samples_all(:,:,i);
    end
    
    % Thin to final length (to reduce the save size of the samples-matrix)
    samples = thin_samples(samples,mcmc_opt.nfinal);
    
    % Compute marginals in the grid using KDE
    post.epost = kde_for_abc(th_grid,samples,1)';
    post.samples = samples;
end
end


function lp = log_post_point_estim(th_prop,th_tr,log_lik_tr,gp,P,sim_model)
% Wrapper for computing the log posterior value from the GP surrogate for the MCMC code.

[eft,varft] = gp_pred_fast(gp, th_tr, log_lik_tr, th_prop, P);
%ll = eft + 0.5*varft; % MEAN AS A POINT ESTIMATOR
%ll = eft; % MEDIAN AS A POINT ESTIMATOR
ll = eft - varft; % MODE AS POINT ESTIMATOR
lp = log(sim_model.prior_eval(th_prop)) + ll;
end



