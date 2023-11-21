function h = plot_comparison_fig(th_grid, loglik_tr, th_tr, th_mcmc_all, ...
    estim_post_gp, estim_post_mcmc, sim_model, res, opt, fignr)
% General function for plotting the estimated posterior and comparing it 
% to the ground truth (if available).

% *red* line == exact or true (SL) posterior
% *black* line/points == estimated GP-based posterior (MH-BLFI)
% *blue* line/points == estimated approx. MCMC posterior (GP-MH)

nr_init = opt.init_n; % initial points for fitting GP
d = th_grid.dim;
% check labels for plotting
if ~isfield(sim_model,'theta_names')
    names = cell(1,d); 
    for i = 1:d
        names{i} = ['\theta_',num2str(i)];
    end
else
    names = sim_model.theta_names;
end

% plot at most 100 of the approx. MCMC samples
th_mcmc = thin_samples(estim_post_mcmc.samples,100);

h = figure(fignr);
clf;
if d == 1
    set(gcf,'Position',[50 1000 1400 400]);
    thx = th_grid.theta;
    nr_plots = 3;
    
    %% Fig 1:
    subplot(1,nr_plots,1); % log lik
    my_shadedplot(thx, estim_post_gp.nloglik_lb, estim_post_gp.nloglik_ub, ...
        0.8*[1,1,1], [1,1,1]); % new loglik uncertainty
    hold on;
    plot(thx, estim_post_gp.eloglik,'-k'); % log lik
    plot(thx, estim_post_gp.loglik_lb,'--k'); % log lik CI (latent)
    plot(thx, estim_post_gp.loglik_ub,'--k');
    
    plot(th_tr,loglik_tr,'*k'); % training data points
    
    xlabel(names{1});
    title('log likelihood');
    xlim([th_grid.range(1),th_grid.range(2)]);
    hold off; box on;

    %% Fig 2:
    subplot(1,nr_plots,2); % estimated posterior based on *GP surrogate* (MH-BLFI)
    
    % NOTE: WE RESCALE THESE AGAIN SO THAT THE POINT ESTIMATE OF THE DENSITY INTEGRATES TO 1
    epost_gp = estim_post_gp.epost; 
    [epost_gp,c] = normalise_pdf_in_grid(epost_gp, th_grid);
    
    my_shadedplot(thx, estim_post_gp.post_lb/c, estim_post_gp.post_ub/c, ...
        0.8*[1,1,1], [1,1,1]); % post uncertainty
    hold on;
    
    plot(thx, epost_gp,'-k'); % estimated post
    if isfield(sim_model,'true_post_pdf') && ~isempty(sim_model.true_post_pdf)
    	plot(thx, sim_model.true_post_pdf,'-r'); % true post
    else
        % plot only true parameter instead
        plot(sim_model.true_theta, 0, '+r','MarkerSize',14);
    end
    plot(th_tr,zeros(size(th_tr)),'*k'); % training data points
    
    xlabel(names{1});
    title('posterior (MH-BLFI)');
    xlim([th_grid.range(1),th_grid.range(2)]);
    hold off; box on;
    
    %% Fig 3
    subplot(1,nr_plots,3); % estimated posterior based on *approx. MCMC* (GP-MH)
    epost_mcmc = estim_post_mcmc.epost; 
    
    hold on;
    plot(thx, epost_mcmc,'-b'); % estimated post
    if isfield(sim_model,'true_post_pdf') && ~isempty(sim_model.true_post_pdf)
    	plot(thx, sim_model.true_post_pdf,'-r'); % true post
    else
        % plot only true parameter instead
        plot(sim_model.true_theta, 0, '+r','MarkerSize',14);
    end
    plot(th_mcmc,zeros(size(th_mcmc)),'*b'); % samples by approx. MCMC (GP-MH)
    if 1
        plot(opt.mcmc.th_init,0,'sm','MarkerSize',6,'MarkerFaceColor','m'); % initial point of MCMC
    end
    
    xlabel(names{1});
    title('posterior (GP-MH)');
    xlim([th_grid.range(1),th_grid.range(2)]);
    hold off; box on;
    
    add_general_title(res);
    
elseif d == 2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    set(gcf,'Position',[50 1000 1400 600]);
    nr_contour_lines = 25;
    thx = th_grid.theta(1,:);
    thy = th_grid.theta(2,:);
    is_truep = isfield(sim_model,'true_post_pdf2d') && ~isempty(sim_model.true_post_pdf2d);
    
    %% Fig 1: loglik - computed from the GP
    subplot(2,3,1);
    eloglik = vec_to_grid_matrix(estim_post_gp.eloglik, th_grid);
    hold on;
    contour(thx, thy, eloglik, nr_contour_lines); % mean
    colorbar;
    add_datapoints_2dplot(th_tr,nr_init,0,sim_model,opt);
    xlabel(names{1}); ylabel(names{2});
    hold off; box on;
    title('GP mean of loglik + eval. locations');
    
    %% Fig 2: stdev of loglik
    subplot(2,3,2);
    stdevloglik = vec_to_grid_matrix(sqrt(estim_post_gp.varloglik), th_grid);
    hold on;
    contour(thx, thy, stdevloglik, nr_contour_lines); % stdev
    colorbar;
    add_datapoints_2dplot(th_tr,nr_init,0,sim_model,opt);
    xlabel(names{1}); ylabel(names{2});
    hold off; box on;
    title('GP stdev of loglik + eval. locations');
    
    %% Fig 3: estimated posterior (MH-BLFI) + training data points
    % NOTE: WE RESCALE THESE AGAIN SO THAT THE POINT ESTIMATE OF THE DENSITY INTEGRATES TO 1
    subplot(2,3,3);
    epost = estim_post_gp.epost;
    epost = vec_to_grid_matrix(epost, th_grid);
    [epost,c] = normalise_pdf_in_grid(epost, th_grid);
    hold on;
    contour(thx, thy, epost, nr_contour_lines);
    add_datapoints_2dplot(th_tr,nr_init,0,sim_model,opt);
    xlabel(names{1}); ylabel(names{2});
    hold off; box on;
    title('estim. posterior (MH-BLFI) + eval. locations');
    
    %% Fig 4: estimated posterior (GP-MH) + training data points
    % NOTE: WE RESCALE THESE AGAIN SO THAT THE POINT ESTIMATE OF THE DENSITY INTEGRATES TO 1
    subplot(2,3,4);
    epost = estim_post_mcmc.epost;
    epost = vec_to_grid_matrix(epost, th_grid);
    [epost,c] = normalise_pdf_in_grid(epost, th_grid);
    hold on;
    contour(thx, thy, epost, nr_contour_lines);
    add_datapoints_2dplot(th_tr,nr_init,0,sim_model,opt);
    xlabel(names{1}); ylabel(names{2});
    hold off; box on;
    title('estim. posterior (GP-MH) + eval. locations');
    
    %% Fig 5&6: true posterior (plot only if available) + training/MCMC points
    if is_truep
        true_post = vec_to_grid_matrix(sim_model.true_post_pdf2d, th_grid);
    end
    for i = 5:6
        subplot(2,3,i);
        hold on;
        if is_truep
            contour(thx, thy, true_post, nr_contour_lines); % true post
        end
        if i==5
            add_datapoints_2dplot(th_tr,nr_init,0,sim_model,opt); % training data points
            title('exact posterior + eval. locations');
        else
            add_datapoints_2dplot(th_mcmc,0,1,sim_model,opt); % samples by approx. MCMC
            title('exact posterior + GP-MH samples');
        end
        xlabel(names{1}); ylabel(names{2});
        xlim([thx(1) thx(end)]); ylim([thy(1) thy(end)]);
        hold off; box on;
    end
    add_general_title(res);
    
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Plot slices
    set(gcf,'Position',[50 1000 1600 800]);
    pls = isfield(estim_post_gp.slice,'eloglik');
    k = 1; % height of the plot
    if pls
        k = k + 3;
        for i = 1:d
            limsi = [th_grid.range(i,1),th_grid.range(i,2)];
            % first plot row
            subplot(k,d,i);
            plot(th_grid.theta(i,:), estim_post_gp.slice.eloglik(i,:),'-k'); % slice of estimated loglik at true val
            xlim(limsi);
            xlabel(names{i});
            if i==1
                ylabel('GP mean loglik (slice)');
            end
            
            % second plot row
            subplot(k,d,d+i);
            plot(th_grid.theta(i,:), sqrt(max(estim_post_gp.slice.varloglik(i,:),0)),'-k'); % slice of stdev of loglik at true val
            xlim(limsi);
            xlabel(names{i});
            if i==1
                ylabel('GP stdev loglik (slice)');
            end
            
            % third plot row
            subplot(k,d,2*d+i);
            plot(th_grid.theta(i,:), estim_post_gp.slice.epost(i,:),'-k'); % slice of estimated post at true val
            xlim(limsi);
            xlabel(names{i});
            if i==1
                ylabel('estim post (MH-BLFI,slice)');
            end
        end
    end
    
    %% Plot posterior marginals (GP and approx. MCMC-based)
    for i = 1:d
        % either first or fourth plot row
        subplot(k,d,k*d-d+i);
        thi = th_grid.theta(i,:);
        hold on;
        h1=plot(thi, estim_post_gp.epost(i,:),'-k'); % point estimate of marginal post, MH-BLFI
        h2=plot(thi, estim_post_mcmc.epost(i,:),'-b'); % point estimate of marginal post, GP-MH
        xlim([th_grid.range(i,1),th_grid.range(i,2)]);
        if isfield(sim_model,'true_post_pdf') && ~isempty(sim_model.true_post_pdf)
            plot(thi, sim_model.true_post_pdf(i,:),'-r'); % true post
        else
            % true post unavailable, plot true parameter instead
            yl = ylim;
            plot(sim_model.true_theta(i)*[1,1],[0,yl(end)],'-r');
            ylim([0,yl(end)]);
        end
        add_datapoint_proj_ndplot(th_tr, nr_init, 0, i); % training data, projections
        add_datapoint_proj_ndplot(th_mcmc, 0, 1, i); % approx. MCMC samples, projections
        if 1
            plot(opt.mcmc.th_init(i),0,'sm','MarkerSize',6,'MarkerFaceColor','m'); % initial point of MCMC
        end
        xlabel(names{i});
        hold off; box on;
        if i==1
            ylabel('estim post');
        elseif i==d
            legend([h1, h2], {'MH-BLFI','GP-MH'}, 'Location' ,'Northeast');
        end
    end
    add_general_title(res);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting help functions:

function [] = add_general_title(res,gen_titl)
% Adds a general level title to the plot which shows TV/KL if these values are computed.

if nargin < 2
    gen_titl = 1;
end
if isfield(res.res_mcmc,'tv') && ~isempty(res.res_mcmc.tv)
    titl = ['TV(MH-BLFI/GP-MH)=',sprintf('%.2f',mean(res.res_gp.tv)),'/',sprintf('%.2f',mean(res.res_mcmc.tv)),...
        ', KL(MH-BLFI/GP-MH)=',sprintf('%.2f',mean(res.res_gp.kl)),'/',sprintf('%.2f',mean(res.res_mcmc.kl))];
    if gen_titl
        suptitle(titl);
    else
        title(titl);
    end
end
end


function [] = add_datapoints_2dplot(ths, nr_init, is_amcmc, sim_model, opt)
% Adds the training data points/MCMC samples to the figure. Also adds the
% true point if it is available. 

if ~isempty(ths)
    col = 'k';
    if is_amcmc
        col = 'b';
    end
    nr_pts = size(ths,1);
    plot(ths(1:nr_init,1),ths(1:nr_init,2),['x',col]); % init/burnin points
    if nr_pts > nr_init
        plot(ths(nr_init+1:end,1),ths(nr_init+1:end,2),['*',col]); 
    end
end
if isfield(sim_model,'true_theta') && ~isempty(sim_model.true_theta)
    plot(sim_model.true_theta(1),sim_model.true_theta(2),'dr','MarkerSize',6,'MarkerFaceColor','r'); % true param
end
if 1
    plot(opt.mcmc.th_init(1),opt.mcmc.th_init(2),'sm','MarkerSize',6,'MarkerFaceColor','m'); % initial point of MCMC
end
end


function [] = add_datapoint_proj_ndplot(ths, nr_init, is_amcmc, i)
% Adds the datapoints projected to coordinate i.

if ~isempty(ths)
    col = 'k';
    if is_amcmc
        col = 'b';
    end
    nr_pts = size(ths,1);
    plot(ths(1:nr_init,i),zeros(nr_init,1),['x',col]); % init/burnin points
    if nr_pts > nr_init % at least one batch
        plot(ths(nr_init+1:end,i),zeros(size(ths(nr_init+1:end,1))),['*',col]); 
    end
end
end



