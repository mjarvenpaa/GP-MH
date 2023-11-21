function [estim_density,init_samples,mcmc_diag] = standard_mcmc_inference(sim_model,...
    grid_th, opt, mcmc_opt, method, use_grid)
% Runs a fairly simple implementation of (exact/noisy/pseudo-marginal) MCMC
% or performs corresponding grid-based computations. Used in particular for
% SL-MCMC.
%
% NOTE: If dim <= 2 and 'use_grid' option is on, posterior values in a 
% grid are returned in 'estim_density', otherwise this variable contains 
% posterior samples generated using (exact/noisy/pseudo-marginal) MCMC. 
% NOTE2: Variable 'init_samples' contains the first 1000 samples generated
% after neglecting the burn-in (whenever applicable). 

if nargin < 6
    use_grid = 0;
end

if use_grid && sim_model.dim == 1
    % evaluate (SL) posterior in a 1d grid - for demonstration only
    
    grid_len = length(grid_th.theta);
    post_grid = NaN(grid_len,1);
    for i = 1:grid_len
        post_grid(i) = noisy_loglik_estim(sim_model,opt,grid_th.theta(i),method) ...
            + log(sim_model.prior_eval(theta(i)));
    end
    estim_density = post_grid;
    init_samples = [];
    mcmc_diag = [];
    
elseif use_grid && sim_model.dim == 2
    % evaluate (SL) posterior in a 2d grid
    
    error('Not implemented.');
    
else % grid not used and/or dim > 2
    
    % We use the adaptive MCMC (DRAM) by Haario et al. 2006 but many other MCMC
    % algorithms could be alternatively used.
    
    display_type = mcmc_opt.display_type;
    
    model.ssfun = @(x,data) -2*(noisy_loglik_estim(sim_model,opt,x,method) + log(sim_model.prior_eval(x)));
    data = [];
    model.N = 1;
    npar = sim_model.dim;
    params = cell(npar,1);
    for i = 1:npar
        params{i} = {sprintf('\\theta_{%d}',i), mcmc_opt.init(i), grid_th.theta(i,1), grid_th.theta(i,end)};
    end
    
    % Additional MCMC settings
    options.nsimu = mcmc_opt.nsimu;
    if isfield(mcmc_opt,'qcov') && ~isempty(mcmc_opt.qcov)
        options.qcov = mcmc_opt.qcov; % initial cov matrix provided, use it
    else
        options.qcov = 1/10^2*diag((grid_th.theta(:,1)-grid_th.theta(:,end)).^2); % use default
    end
    options.method = 'am';
    options.updatesigma = 0;
    options.verbosity = 0; % no printing from mcmc
    options.waitbar = 0;
    
    % Initialize results
    samples_all = NaN(mcmc_opt.nsimu,npar,mcmc_opt.nchains);
    results_all = cell(mcmc_opt.nchains,1);
    
    % Run MCMC chains!
    for i = 1:mcmc_opt.nchains
        if ~strcmp(display_type,'off')
            if i == 1
                disp('Running MCMC...');
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
    
    % Leave out burn-in
    cl = size(samples_all,1);
    if isfield(mcmc_opt,'burninfreq') && ~isempty(mcmc_opt.burninfreq)
        samples_all = samples_all(ceil(mcmc_opt.burninfreq*cl):cl,:,:);
    else
        % default option, leave out first half
        samples_all = samples_all(ceil(cl/2):cl,:,:);
    end
    
    % Assess the convergence
    mcmc_diag = basic_mcmc_diagnostics(samples_all,display_type);
    
    % Add the chains together
    cl = size(samples_all,1);
    samples = NaN(mcmc_opt.nchains*cl,npar);
    for i = 1:mcmc_opt.nchains
        samples(1 + (i-1)*cl:i*cl,:) = samples_all(:,:,i);
    end
    
    % Extract the first 1000 samples 
    % (If multiple chains, only 1000 samples from the first chain are returned)
    if nargout > 1
        init_n = min(1000,size(samples,1));
        init_samples = samples(1:init_n,:);
    end
    
    % Thin to final length (to reduce the save size of the samples-matrix)
    samples = thin_samples(samples,mcmc_opt.nfinal);
    
    % Print and plot results to visually examine whether convergence is reached
    if ~strcmp(display_type,'off')
        figure(50);
        clf;
        mcmcplot(samples,1:npar,results_all{1}.names,'chainpanel');
        set(gcf,'Position',[60 600 600 400]);
    end
    estim_density = samples;
end
end



