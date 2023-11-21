function [] = sl_inference_baseline(cluster_run, run_id, m_id, plot_res)
% Computes the baseline SL posterior density for the real test problems 
% whose exact likelihood is unavailable. The results are saved to file to 
% be used for comparisons in the GP-MH code. This same code can also be 
% used for plotting the resulting posterior once the MCMC has first been 
% run and results saved. 

TEST_MODE = 0;
SAVE_FIGS = 1;

if cluster_run
	warning off; 
    set_paths_for_cluster_run();
    plot_res = 0;
else
    disp('Local run...');
end

root = '../results/sl_baseline/';
if TEST_MODE
    root = '../results/sl_baseline_test/'; % For testing...
end

models = {'lorenz','ricker','ricker_12','gk_model','ma2','thetaricker','thetaricker2','cell_model'};
nmodels = length(models);
nruns = [5,5,5,5,5,5,5,20]; % 5 runs per model, 10 for slow cell_model
Ns = [100,100,100,100,100,100,100,2500]; % N for SL
ns = [1e6,1e6,1e6,1e6,1e6,1e6,2e6,2e4]; % how long MCMC chain, shorter for slow cell_model
qcovs = cell(nmodels,1);
qcovs{8} = diag([0.01 10^2*1e-4].^2); % a specific choice for cell_model, the default is ok for the rest

if TEST_MODE % For testing...
    %m_id = 8; run_id = 1;
    nruns = 2*ones(1,8);
    ns = [1000*ones(1,7),200];
end
method = 'sl';

if ~plot_res
    if ~cluster_run
        % some debug printing:
        m_id, run_id, cur_model=models{m_id}, cur_N=Ns(m_id), cur_n = ns(m_id)
    end
    
    % get test model settings
    [grid_th,sim_model] = get_test_model(models{m_id},[],Ns(m_id));
    
    %% SL settings
    % mcmc related settings
    sl_mcmc_opt.init = sim_model.true_theta; % USE THE 'TRUE VALUE' AS INITIAL POINT FOR MCMC
    sl_mcmc_opt.nsimu = ns(m_id);   
    sl_mcmc_opt.nchains = 1;
    sl_mcmc_opt.nfinal = 10000;
    sl_mcmc_opt.display_type = 'off';
    sl_mcmc_opt.qcov = qcovs{m_id};
    if m_id == 8
        %% NOTE: SOME SPECIFIC CHANGES FOR CELL MODEL
        sl_mcmc_opt.burninfreq = 0.25; % SHORTER BURNIN FOR CELL MODEL
        sl_mcmc_opt.nfinal = 15000; % return more samples (all when size is 2e4)
    end
    
    sl_opt.N = Ns(m_id);
    sl_opt.estimator = 'sl'; % 'sl','ubsl','ublogsl'
    
    %% run SL-MCMC inference
    sl_seed = 42 + run_id; % different seed for each run
    rng(sl_seed);
    if ~cluster_run
        t0=tic;
    end
    try
        [sl_samples,init_sl_samples,mcmc_diag] = standard_mcmc_inference(sim_model,...
            grid_th, sl_opt, sl_mcmc_opt, method);
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
    
    % compute marginals from the SL-MCMC samples
    % NOTE: This is redone in 'get_test_model' so that small changes to 
    % e.g. parameter bounds can be done without rerunning SL-MCMC.
    post_kde = kde_for_abc(grid_th,sl_samples);
    if sim_model.dim == 2
        post_kde = vec_to_grid_matrix(post_kde, grid_th);
    end
    
    %% save thinned set of samples, initial 1000 samples and the final estimated density
    fn = [root,'/',models{m_id},'_run',num2str(run_id)];
    save(fn,'sl_samples','init_sl_samples','post_kde','mcmc_diag','grid_th',...
        'sim_model','sl_mcmc_opt','sl_opt');
end


%% Plot results
if plot_res && ~cluster_run
    close all;
    % load already computed results from files...
    nru = nruns(m_id);
    post = cell(nru,1);
    samples = cell(nru,1);
    for i = 1:nru
        fn = [root,'/',models{m_id},'_run',num2str(i)];
        load(fn,'sl_samples','post_kde','mcmc_diag','grid_th','sim_model','sl_mcmc_opt','sl_opt'); 
        post{i} = post_kde;
        samples{i} = sl_samples;
        % print basic convergence info:
        R = mcmc_diag.R
        is_converged = mcmc_diag.is_converged
        
        % plot trajectory
        if i <= 10
            figure(10+i);
            d = size(sl_samples,2);
            for j = 1:d
                subplot(d,1,j);
                plot(1:size(sl_samples,1),sl_samples(:,j)','-k');
            end
        end
    end
    d = sim_model.dim;
    
    % plot SL-posterior
    figure(1);
    cols = {'b','r','k','g','m','--b','--r','--k','--g','--m'};
    if d ~= 2
        % 1) plot all 1d marginals to one figure:
        set(gcf,'Position',[50 1000 400*d 450]);
        for i = 1:d
            subplot(1,d,i); hold on;
            for j = 1:nru
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
            for j = 1:nru
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
    end
    
    %% Combine the samples to create the final ground-truth posterior
    if m_id <= 7
        % We just use the first MCMC run because the variability is small
        % enough and to keep compatibility with previous work. This is
        % handled manually. 
        
    elseif m_id == 8
        % cell_model case; we do combine all the fairly short chains to reduce MC error
        
        % Add the chains together
        cl = size(samples{1},1);
        sl_samples = NaN(nru*cl,d);
        for i = 1:nru
            sl_samples(1 + (i-1)*cl:i*cl,:) = samples{i};
        end
        
        % Thin to final length (to reduce the save size of the samples-matrix)
        % NOTE: WE USE ACTUALLY 20000 SAMPLES HERE
        sl_samples = thin_samples(sl_samples,20000);
        
        fn = [root,'/',models{m_id}];
        save(fn,'sl_samples','grid_th','sim_model','sl_mcmc_opt','sl_opt');
    end
end

if cluster_run
    exit;
end
end


