function [] = demo_loglik_gp_modelling()
% Plots a 1x2 figure that demonstrates the robust likelihood estimator
% that downweights the estimate in regions where computations have not been
% performed and where the uncertainty is thus especially large. 
% ***Used to create figures for v1/v2 of the paper (Figure 2, B.1, B.2 in v1).***

subplot = @(m,n,p) subtightplot (m, n, p, [0.07 0.09], [0.2 0.1], [0.1 0.1]); % better subplots

close all;
rng(1234);

%% set up model etc.
save_fig = 1;
load_from_file = 0;

%root = '../results/demo_loglik/'; % for v1 of the paper
root = '../results/demo_loglik_v2/'; % for v2  of the paper
%model = 'ricker_1'; N = 100;
model = 'gaussian1d'; N = 50;
[th_grid,sim_model] = get_test_model(model,[],N);
if 1 && strcmp(model,'gaussian1d')
    % slightly adjust the parameter space of the default test case for a better demonstration
    th_grid.range(2) = th_grid.range(2)+10;
    th_grid.theta = linspace(th_grid.range(1),th_grid.range(2),1000);
    mo = 'g';
elseif 1 && strcmp(model,'ricker_1')
    th_grid.range(2) = th_grid.range(2)+0;
    th_grid.theta = linspace(th_grid.range(1),th_grid.range(2),500);
    mo = 'r';
end

%% Adjust these two settings for different plots:
extra_eval = 0;
gp_opt.meanf = 0; % 0 is zero mean GP, 1 enables const/lin/quadratic terms

%% other fixed settings
gp_opt.noise_model = 1; % if 1, then use bootstrapped noise variance estimates in GP
gp_opt.display_type = 'off';
gp_opt.acq_noise_model = 0;
gp_opt.acq_tol2 = 1e-4;

lik_opt.method = 'sl'; 
lik_opt.sl.estimator = 'sl'; % 'sl', 'ubsl', 'ublogsl'
lik_opt.sl.N = N;

%% generate some logSL values (here eval points set manually for the demo)
if strcmp(model,'gaussian1d')
    th_tr = [-18 -15 -13 -11.5 -9 -2 0.5 1 7]';
    if extra_eval
        th_tr = [th_tr; 25]; % add an extra point...
    end
else
    th_tr = linspace(th_grid.range(1),th_grid.range(2)-1,20)';
end
t = length(th_tr);
fn_simul = [root,'loglik_cached_simulations'];
if load_from_file
    % load cached simulations
    load(fn_simul);
else
    % simulate now and save the results
    loglik_tr = NaN(t,1);
    bootvar = NaN(t,1);
    for j = 1:t
        [loglik_tr(j),bootvar(j)] = noisy_loglik_estim(sim_model,lik_opt.sl,...
            th_tr(j,:),lik_opt.method,gp_opt.noise_model);
    end
    save(fn_simul,'loglik_tr','bootvar');
end
sigma_tr = sqrt(bootvar);

%% fit GP and get the loglik estimates from it etc. 
gp = []; gp_optim_opt = [];
[gp,gp_optim_opt] = fit_gp_model(gp, gp_optim_opt, gp_opt,th_grid, loglik_tr,th_tr,sigma_tr);
if 1
    % decrease the lengthscale from the optimal value to make a better demonstration
    gp.cf{1}.lengthScale = gp.cf{1}.lengthScale/2;
end
P = precompute_gp_pred(gp, th_tr, loglik_tr, gp_opt);
estim_post = post_from_gp_surrogate(th_grid,sim_model,gp,gp_opt,loglik_tr,th_tr,sigma_tr,P,[]);

%% plot everything
figure(1);
set(gcf,'Position',[40 40 900 300]);

lw = 0.95; % linewidth, matlab default 0.5

% 1/2: GP on loglik
subplot(1,2,1);
thx = th_grid.theta;
my_shadedErrorBar(thx, estim_post.eloglik, estim_post.loglik_lb, estim_post.loglik_ub,'r',lw); % loglik + its uncertainty
hold on;
plot(th_tr,loglik_tr,'k.','MarkerSize',10); % loglik observations
if 1
    % plot 'errorbars' of the loglik observations
    for i = 1:t
        plot(th_tr(i)*[1 1],loglik_tr(i)+[1 -1]*1.96*sigma_tr(i),'-k','Linewidth',0.7);
    end
end
hold off;
box on;
xlim([thx(1),thx(end)]);
ylim([min(estim_post.eloglik)-0.1*abs(min(estim_post.eloglik)),max(estim_post.eloglik)+1.3*abs(max(estim_post.eloglik))]);
xlabel('\theta');
ylabel('Log-synthetic likelihood');
title('(a)');


% 2/2: uncertainty of lik
subplot(1,2,2);
smed=my_shadedErrorBar(thx, estim_post.medpost, estim_post.post_lb, estim_post.post_ub,'r',lw); % med + uncertainty
hold on; 
if 1
    hmode=plot(thx,estim_post.modepost,'--b','Linewidth',lw); % mode as the more sensible estimator in our case
end
if 0 % not plotted for v2 anymore; was redundant 
    % plot true value
    plot(sim_model.true_theta,0,'dr','MarkerSize',6,'MarkerFaceColor','r');
end
hold off;
box on;
xlim([thx(1),thx(end)]);
ylim([0,1.3*max(estim_post.epost)]);
xlabel('\theta');
ylabel('Synthetic likelihood (rescaled)');
title('(b)');
if 1
    %legend([smed.mainLine,hmode],{'med (35)','mode (37)'},'Location','NorthWest','FontSize',9);
    legend([smed.mainLine,hmode],{'median','mode'},'Location','NorthWest','FontSize',9);
end

% save to file
if save_fig
    fn = [root,'fig_loglik_gp_', mo, '_ep',num2str(extra_eval),'_mf',num2str(gp_opt.meanf)];
    my_export_fig(fn,'-transparent','-painters','-pdf');
end
end


