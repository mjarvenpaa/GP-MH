function [] = genpost_inference_baseline(cluster_run, run_id, m_id, M, plot_res)
% Computes the baseline (generalised, discrepancy loss-based) posterior 
% density for the real test problems whose exact likelihood is unavailable.
% Currently this means the bacterial infections model.
% The results are saved to file to be used for comparisons in the GP-MH 
% code. This same code can also be used for plotting the resulting 
% posterior once the MCMC has first been run and results saved.
%
% NOTE: The code structure is similar to 'sl_inference_baseline' which is 
% for the synthetic likelihood case.
%
% TODO: Refactor? - A lot of code is same or similar to 'sl_inference_baseline'.

% We first compute w10 case with M=1 and then with M=10. 
% Then w5,w20 with M=10.

TEST_MODE = 0;
SAVE_FIGS = 1;

if cluster_run
	warning off; 
    set_paths_for_cluster_run();
    plot_res = 0;
else
    disp('Local run...');
end

root = '../results/genpost_baseline/';
if TEST_MODE
    root = '../results/genpost_baseline_test/'; % For testing...
end

% bacterial model: 1 simul (M=1)  *1-2s*,  10^4 simul (M=1)  2.8-5.6h
%                  1 simul (M=10) 10-20s,  10^4 simul (M=10) 1d4h-2d8h
%                  1 simul (M=50) 50-100s, 10^4 simul (M=10) >5d

models = {'bacterial_w5','bacterial_w10','bacterial_w20','bacterial_w40'}; % 4 w's, widest gen post case first, 4th case added after running 3 first...
nmodels = length(models);
nruns = 50*ones(1,nmodels); % how many MCMC chains
Ms = M*ones(1,nmodels); % M for the Monte Carlo estimate of E(discrepancy)
ns = 1.25*1e4*ones(1,nmodels); % how long MCMC chain
qcovs = cell(nmodels,1); % kept the same for all M, based on pilot runs with GP-MH and EPoEr
qcovs{1} = 2.4^2/3*[5.57 0.25 -0.13; 0.25 0.15 0; -0.13 0 0.0094];
qcovs{2} = 2.4^2/3*diag([0.85 0.032 0.0014]);
qcovs{3} = 2.4^2/3*diag([0.32 0.012 0.00037]);
qcovs{4} = 2.4^2/3*diag([0.18 0.0046 0.00034]);

if TEST_MODE % For testing...
    nruns = 2*ones(1,nmodels);
    Ms = 3*ones(1,nmodels); % that is, M=3 case is forced
    ns = 100*ones(1,nmodels);
end
method = 'exact'; % that is, not SL

% increase the chain length if M=1 case (mixing likely worse then...)
if M==1 && ~TEST_MODE
    ns = 2*ns;
end

% results with different M are saved to different folders:
root = [root,'/','M',num2str(Ms(m_id)),'/'];

if ~plot_res
    if ~cluster_run
        % some debug printing:
        m_id, run_id, cur_model=models{m_id}, cur_M=Ms(m_id), cur_n = ns(m_id)
    end
    
    % get test model settings
    [grid_th,sim_model] = get_test_model(models{m_id},[],Ms(m_id));
    
    %% settings
    % mcmc related settings
    genpost_mcmc_opt.init = sim_model.true_theta; % USE THE 'TRUE VALUE' AS INITIAL POINT FOR MCMC
    genpost_mcmc_opt.nsimu = ns(m_id);   
    genpost_mcmc_opt.nchains = 1;
    genpost_mcmc_opt.nfinal = 10000;
    genpost_mcmc_opt.display_type = 'off';
    genpost_mcmc_opt.qcov = qcovs{m_id};
    genpost_mcmc_opt.burninfreq = 0.20; % default would be 0.5
    
    lik_opt = []; % M is the only option and is set via model specification in 'get_test_model(...)'
    
    %% run SL-MCMC inference
    genpost_seed = 42 + run_id; % different seed for each run
    rng(genpost_seed);
    if ~cluster_run
        t0=tic;
    end
    try
        [genpost_samples,init_genpost_samples,mcmc_diag] = standard_mcmc_inference(sim_model,...
            grid_th, lik_opt, genpost_mcmc_opt, method);
    catch err
        if cluster_run
            write_err(run_id, m_id, err);
            exit;
        else
            err.getReport('extended','hyperlinks','off')
            keyboard;
        end
    end
    if ~cluster_run
        toc(t0)
    end
    
    % compute marginals from the MCMC samples
    % NOTE: This is redone in 'get_test_model' so that small changes to 
    % e.g. parameter bounds can be done without rerunning MCMC.
    post_kde = kde_for_abc(grid_th,genpost_samples);
    if sim_model.dim == 2
        post_kde = vec_to_grid_matrix(post_kde, grid_th);
    end
    
    %% save thinned set of samples, initial 1000 samples and the final estimated density
    fn = [root,'/',models{m_id},'_run',num2str(run_id)];
    save(fn,'genpost_samples','init_genpost_samples','post_kde','mcmc_diag','grid_th',...
        'sim_model','genpost_mcmc_opt','lik_opt');
end


%% Plot results
if plot_res && ~cluster_run
    close all;
    % load already computed results from files...
    nru = nruns(m_id);
    %nru=5
    post = cell(nru,1);
    samples = cell(nru,1);
    for i = 1:nru
        fn = [root,'/',models{m_id},'_run',num2str(i)];
        load(fn,'genpost_samples','post_kde','mcmc_diag','grid_th','sim_model','genpost_mcmc_opt','lik_opt'); 
        post{i} = post_kde;
        samples{i} = genpost_samples;
        % print basic convergence info:
        R = mcmc_diag.R
        is_converged = mcmc_diag.is_converged
        
        % plot trajectory
        if i <= 10
            figure(10+i);
            d = size(genpost_samples,2);
            for j = 1:d
                subplot(d,1,j);
                plot(1:size(genpost_samples,1),genpost_samples(:,j)','-k');
            end
        end
    end
    d = sim_model.dim;
    
    % plot (generalised) posterior
    figure(1);
    cols = {'b','r','k','g','m','--b','--r','--k','--g','--m'};
    if d ~= 2
        % 1) plot all 1d marginals to one figure:
        set(gcf,'Position',[50 1000 400*d 450]);
        for i = 1:d
            subplot(1,d,i); hold on;
            for j = 1:min(length(cols),nru)
                plot(grid_th.theta(i,:),post{j}(:,i),[cols{j}]);
            end
            hold off;
            if isfield(sim_model,'theta_names')
                xlabel(sim_model.theta_names{i});
            end
            box on;
        end
        drawnow;
        if SAVE_FIGS
            fn = [root,'/sl_marg_post_',models{m_id}];
            my_export_fig(fn,'-png');
        end
        
        % 2) plot also all 2d marginals:
        if d > 2
            for j = 1:min(10,nru)
                figure(1+j);
                set(gcf,'Position',[50 50 1000 1000]);
                bivariate_marginals_plot_old(grid_th,sim_model,samples{j},post{j});
                drawnow;
                if SAVE_FIGS && j == 1
                    % only first figure saved to file
                    fn = [root,'/sl_post_',models{m_id}];
                    my_export_fig(fn,'-png');
                end
            end
        end
        
    else % d == 2
        % plot joint 2d posterior:
        nr_contour_lines = 25;
        thx = grid_th.theta(1,:);
        thy = grid_th.theta(2,:);
        set(gcf,'Position',[50 1000 400*nru 350]);
        for j = 1:min(nru,10)
            subplot(1,min(nru,10),j);
            contour(thx, thy, post{j}, nr_contour_lines);
            if isfield(sim_model,'theta_names')
                xlabel(sim_model.theta_names{1}); ylabel(sim_model.theta_names{2});
            end
        end
        drawnow;
        if SAVE_FIGS
            fn = [root,'/sl_post_',models{m_id}];
            my_export_fig(fn,'-png');
        end
    end
    
    % compute all pairwise TV/KL/L2 distances to get a (rough) estimate of variability
    if 1
        % KL is not symmetric so in KL case we could in fact compute all combinations
        % but this is not done here.
        ds = NaN(nru,nru);
        for i = 1:nru
            for j = (i+1):nru
                [kl,tv,l2] = compute_dist(grid_th, post{i}', post{j}');
                ds(i,j) = mean(tv);
            end
        end
        
        % print results
        ds
        meand = nanmean(ds(:))
        mediand = nanmedian(ds(:))
        maxd = nanmax(ds(:))
    end
    
    %% Combine the samples to create the final ground-truth posterior
    if m_id >= 1
        % Add the chains together
        cl = size(samples{1},1);
        genpost_samples = NaN(nru*cl,d);
        for i = 1:nru
            genpost_samples(1 + (i-1)*cl:i*cl,:) = samples{i};
        end
        
        % Thin to final length (to reduce the save size of the samples-matrix)
        % NOTE: WE USE ACTUALLY 20000 SAMPLES HERE
        genpost_samples = thin_samples(genpost_samples,20000);
        
        fn = [root,'/',models{m_id}];
        save(fn,'genpost_samples','grid_th','sim_model','genpost_mcmc_opt','lik_opt');
    end
end

if cluster_run
    exit;
end
end


