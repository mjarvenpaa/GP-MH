function res = GPMH_test_run()
% An example script for running an implementation of GP-MH (and simultaneously MH-BLFI).
% Use addpath(genpath('.')); to add files to MATLAB path.

clear; close('all'); format compact;

use_profiler = 0; % whether to use matlab profiler to analyze the code
seed = 123;
rng(seed);


%% some general settings
opt.init_n = 10; % number of initial points for GP fitting
opt.n = 0.5*10^5; % how many samples for GP-MH (NOTE: We use small value for demonstration, larger value may be necessary for real analysis)
opt.nbatch = 1; % how many evaluations acquired at a time (currently only 1 works)
opt.burninfrac = 0.25; % fraction of burn-in samples that are neglected (set to e.g. 0.25 or 0.5)
opt.nfinal = 10000; % maximum amount of samples after neglecting burn-in and thinning
opt.maxacq = 1000; % hard threshold for the maximum number of acquisitions
opt.maxacqiter = Inf; % hard threshold for the maximum number of acquisitions / iteration (Inf == no limit)
opt.viz = 1; % interactive plotting: plot only final result=1 / plot after each iter=2 / compute estimates but no plotting=-1 (cluster computation)
opt.ret_lvl = 1; % return basic output and debugging info=1 / return some more info=2


%% select test problem
%[th_grid,sim_model] = get_test_model('gaussian1d',[],100); 
%[th_grid,sim_model] = get_test_model('gaussian2d',[],100); 
%[th_grid,sim_model] = get_test_model('gaussian_6',[],100);
%[th_grid,sim_model] = get_test_model('ricker_1',[],100); 
%[th_grid,sim_model] = get_test_model('ricker_12',[],100); 
[th_grid,sim_model,samples] = get_test_model('ricker',[],100);
%[th_grid,sim_model,samples] = get_test_model('thetaricker',[],100); % 4-parameter version
%[th_grid,sim_model,samples] = get_test_model('thetaricker2',[],100); % 5-parameter version used in the paper
%[th_grid,sim_model] = get_test_model('simple2d',4,[]);
%[th_grid,sim_model] = get_test_model('banana2d',1,[]);
%[th_grid,sim_model] = get_test_model('bimodal2d',1,[]);
%[th_grid,sim_model,samples] = get_test_model('simple6d',4,[]);
%[th_grid,sim_model,samples] = get_test_model('banana6d',2,[]);
%[th_grid,sim_model,samples] = get_test_model('bimodal6d',2,[]);
%[th_grid,sim_model] = get_test_model('cell_model',[],2500);
%[th_grid,sim_model,samples] = get_test_model('gk_model',[],100);
%[th_grid,sim_model] = get_test_model('lorenz',[],100);
%[th_grid,sim_model] = get_test_model('ma2',[],100);

%[th_grid,sim_model,samples] = get_test_model('bacterial_w5',[],1);
%[th_grid,sim_model,samples] = get_test_model('bacterial_w10',[],1);
%[th_grid,sim_model,samples] = get_test_model('bacterial_w20',[],1);
%[th_grid,sim_model,samples] = get_test_model('bacterial_w40',[],1);


%% MH settings
opt.mcmc.th_init = th_grid.range(:,1) + 0.2*(th_grid.range(:,2)-th_grid.range(:,1)); % initial point
opt.mcmc.q_Sigma = diag((th_grid.range(:,1)-th_grid.range(:,2)).^2/10^2)/100; % initial cov matrix for Gaussian proposal
opt.mcmc.q_adapt = 1; % 1 if adapt covariance matrix as in Adaptive Metropolis method
opt.mcmc.q_i0 = 100; % when the cov matrix adaptation is started (set to at least 20)
opt.mcmc.extra_scaling = 1; % extra scaling for adaptive cov matrix (use 1)
opt.mcmc.ep = 1e-9; % for ensuring proposal cov is proper cov matrix


%% special initial locations etc.:
if 0
    opt.mcmc.th_init = sim_model.true_theta; % initialise at 'true' value
end
if strcmp(sim_model.name,'Cell')
    sc = 10^2;
    opt.mcmc.th_init = [0.5  sc*15e-4]; % true: [0.35, sc*1e-3]
    opt.mcmc.q_Sigma = diag([0.02 sc*2e-4].^2);
elseif strcmp(sim_model.name,'Theta-Ricker5')
    opt.mcmc.th_init = [3.4   0.9 3.0   8.0 0.30]; % true: [3.5, 1, 3.5, 10, 0.3]
    opt.mcmc.q_Sigma = diag([0.05   0.1 0.25   0.5 0.05].^2);
end


%% set the required accuracy of the accept/reject decision of MH
opt.mcmc.tol_crit = 'uncond'; % which MH accept/reject error criterion: 'cond' or 'uncond'
%opt.mcmc.tol_cond = 0.3; % maximum allowed error, \epsilon in the paper
opt.mcmc.tol_uncond = 0.3; % maximum allowed error, \epsilon in the paper


%% settings related to acquisition of new loglik evaluations
%opt.acq_opt.method = 'naive_random'; % selects current location or the proposal randomly, called just 'naive' in the paper
%opt.acq_opt.method = 'EPoE'; % minimise expected probability of error of the M-H acceptance/rejection decision
opt.acq_opt.method = 'EPoEr'; % as above but optimisation restricted to current and proposed point only
opt.acq_opt.optim_alg = 'fmincon'; % which optimisation algorithm for obtaining the design points in EPoE case (use 'fmincon', others not carefully tested)
%opt.acq_opt.optim_alg = 'cmaes';
%opt.acq_opt.optim_alg = 'grid';
%opt.acq_opt.optim_alg = 'rs'; % random search
%opt.acq_opt.optim_alg = 'direct';
opt.acq_opt.rs.evals = 1000; % number of evals in random search
opt.acq_opt.fmincon.nr_inits = 1000; % number of initial points for the multistart optimization
opt.acq_opt.fmincon.nr_inits_grad = 10; % number of the best initial points that are actually used for multistart optimization
opt.acq_opt.fmincon.tolf = 1e-5; % error tolerances for fmincon optimizer
opt.acq_opt.fmincon.tolx = 1e-5;
opt.acq_opt.direct.maxevals = 1000; % max. number of function evals  (default is 20)
opt.acq_opt.direct.maxits = 100; % max. number of iterations  (default is 10)
opt.acq_opt.direct.maxdeep = 100; % max. number of rect. divisions (default is 100)
opt.acq_opt.display_type = 'off'; % 'on' / 'off', whether to print some info on acq


%% GP settings
opt.gp_opt.noise_model = 1; % 0=constant GP noise term, 1=uses bootstrapped noise variance estimates in GP
opt.gp_opt.meanf = 1; % 0 is zero mean GP prior, 1 enables const/lin/quadratic terms
opt.gp_opt.acq_noise_model = 0; % 1 uses NN interp. to estim. variance at new \theta^*, otherwise acq_tol2 (used only when opt.gp_opt.noise_model=1)
opt.gp_opt.acq_tol2 = 1e-2;
opt.gp_opt.display_type = 'off'; % 'on' / 'off', whether to print some info on GP fitting


%% MCMC settings for the second stage of MH-BLFI
opt.mcmcgp.nchains = 5; % how many chains
opt.mcmcgp.nsimu = 20000; % how many samples for each chain (NOTE: We use small value for demonstration, larger value may be necessary for real analysis)
opt.mcmcgp.nfinal = 10000; % maximum amount of samples after concatenating the chains and thinning
opt.mcmcgp.display_type = 'on'; % 'on' / 'off', whether to print some info on MCMC


%% loglik computation
% which method to obtain loglik estimates:
if isfield(sim_model,'loglik_eval') % A function that outputs (noisy) log-likelihood evaluations is provided
    opt.lik_opt.method = 'exact';
    opt.gp_opt.noise_model = 0;
    % NOTE: loglik can in principle be replaced by any log-density of interest
else % Common LFI setting: A function that simulates (summarised) data from the model is only provided
    opt.lik_opt.method = 'sl'; % Synthetic likelihood is used in this case, other methods currently not implemented
end
opt.lik_opt.sl.estimator = 'sl'; % which SL estimator: 'sl', 'ubsl', 'ublogsl'
opt.lik_opt.sl.N = sim_model.N; % number of repeated samples computing SL at each evaluation location
opt.lik_opt.sl.robust_bootvar = 1; % if 1 uses robust variance estimator in bootstrap 
opt.lik_opt.sl.bootn = 2000; % how many bootstrap samples
opt.lik_opt.maxabsloglik = 1e6; % maximum absolute value of loglik that is considered valid output
opt.lik_opt.maxbootsigma = 1e3; % maximum loglik stdev that is considered valid output


%% other
opt.misc.res_ind = logspaced_inds(1,opt.n,25); % iterations when to compute TV/KL if cluster computation
opt.misc.display_type = 'on'; % 'on' / 'off', whether to print some info on the progression of the main algorithm


%% further adjustments for Cell model
if strcmp(sim_model.name,'Cell')
    opt.lik_opt.sl.bootn = 1000; % to speedup testing with Cell model here
    opt.gp_opt.acq_tol2 = 1^2;
end  

if exist('samples','var')
    sim_model.samples = samples; % usually not included to 'sim_model' to make its filesize smaller
end

% Note that some settings related to the GP prior and MCMC are currently hard-coded and not included here.
% However, it may be useful to adjust also these for particular inference problem at hand. 

if use_profiler
    profile on;
end

if 1
    % print initial point and the initial cov matrix:
    init = opt.mcmc.th_init
    init_cov_sqrt_diag = sqrt(diag(opt.mcmc.q_Sigma)')
end

%% run main algorithm
c = tic;
res = GPMH_main(sim_model, th_grid, opt);
elapsed_time_in_minutes = toc(c)/60

if use_profiler
    profile viewer;
end
end


