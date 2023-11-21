function [ths,f_ths] = acquire(i,u,th_grid,th_cur,th_pr,th_tr,y_tr,sigma_tr,gp,gpo_opt,gp_opt,P,sim_model,opt,fignr)
% Acquire next evaluation location(s) for improving the accuracy of the 
% MCMC acceptance check.

if i < 1 || isempty(gp)
    error('Wrong iteration.');
end
if opt.nbatch ~= 1
    error('Only sequential case b=1 currently implemented.');
end
interp_acq = 1;
rep_th = 1; % no joint batch methods here
acq_opt = opt.acq_opt;
viza = 0; f_ths = NaN;
f_alt = []; % alternative approach to compute acq or some baseline acq
inits = []; th_range_loc = [];

% record time used to acquire the next point(s) in the current iteration
start = tic;

%% Select the design criterion (aka acquisition function)
if strcmp(acq_opt.method,'naive_random')
    % Select randomly either current or proposed point.
    % Called just 'naive' in the paper.
    ind = (rand(1) > 0.5);
    ths = ind*th_pr(:)' + (~ind)*th_cur(:)';
    
elseif strcmp(acq_opt.method,'naive_prop')
    % Select proposed point.
    % NOTE: This method may lead to the algorithm getting stuck. 
    ths = th_pr(:)';
    
elseif strcmp(acq_opt.method,'EValpha')
    % extra acq method not carefully tested
    f = @(th_new)EV_acq_numint(1,u,th_new,th_cur,th_pr,th_grid,gp,gpo_opt,gp_opt,th_tr,y_tr,sigma_tr,sim_model,P);
    [ths,f_ths,inits,th_range_loc] = optimize_acq(f, th_grid, th_tr, th_cur, th_pr, gp,gpo_opt,gp_opt,acq_opt, rep_th);
    viza = 1;
    
elseif strcmp(acq_opt.method,'EVp')
    % extra acq method not carefully tested
    f = @(th_new)EV_acq_numint(2,u,th_new,th_cur,th_pr,th_grid,gp,gpo_opt,gp_opt,th_tr,y_tr,sigma_tr,sim_model,P);
    [ths,f_ths,inits,th_range_loc] = optimize_acq(f, th_grid, th_tr, th_cur, th_pr, gp,gpo_opt,gp_opt,acq_opt, rep_th);
    viza = 1;
    
elseif strcmp(acq_opt.method,'EPoE')
    %f_alt = @(th_new)EV_acq_numint(3,u,th_new,th_cur,th_pr,th_grid,gp,gpo_opt,gp_opt,th_tr,y_tr,sigma_tr,sim_model,P);
    f = @(th_new)EPoE_acq(th_new,th_cur,th_pr,th_grid,gp,gpo_opt,gp_opt,th_tr,y_tr,sigma_tr,sim_model,P);
    [ths,f_ths,inits,th_range_loc] = optimize_acq(f, th_grid, th_tr, th_cur, th_pr, gp,gpo_opt,gp_opt,acq_opt, rep_th);
    viza = 1;
    
elseif strcmp(acq_opt.method,'EPoEr') || strcmp(acq_opt.method,'EPoE_restricted')
    % As EPoE but optimisation is over only current and proposed point!
    f = @(th_new)EPoE_acq(th_new,th_cur,th_pr,th_grid,gp,gpo_opt,gp_opt,th_tr,y_tr,sigma_tr,sim_model,P);
    fs = f([th_cur(:)';th_pr(:)']);
    if fs(1) < fs(2)
        ths = th_cur(:)'; f_ths = fs(1);
    else
        ths = th_pr(:)'; f_ths = fs(2);
    end
else
    error('Incorrect acq method.');
end

tm = toc(start);
if strcmp(acq_opt.display_type,'on')
    %disp(' ');
    disp(['Time used for acquiring the next point(s) = ', num2str(tm), ' seconds.']);
    disp(' ');
end

%% Plot the surface of the design criterion (aka acquisition function)
if viza && opt.viz >= 2
    viz_acq(f,f_alt,th_grid,f_ths,ths,inits,th_range_loc,th_cur,th_pr,th_tr,y_tr,opt,interp_acq,fignr);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function acq = EPoE_acq(th_new,th_cur,th_pr,th_grid,gp,gpo_opt,gp_opt,th_tr,y_tr,sigma_tr,sim_model,P)
% Computes the (negative logarithm of) xi-function whose minimum value is 
% the optimal design point for the three design criteria.  

nnew = size(th_new,1);
acq = NaN(nnew,1);
for i = 1:nnew
    exi = gp_lookahead_cov_fast(gp,gp_opt,th_tr,y_tr,sigma_tr,[th_cur;th_pr],th_new(i,:),P,0);
    acq(i) = exi(1,1) + exi(2,2) - 2*exi(1,2);
end
acq = -log(max(eps,acq));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extra methods/code:

function acq = EV_acq_numint(method,u,th_new,th_cur,th_pr,th_grid,gp,gpo_opt,gp_opt,th_tr,y_tr,sigma_tr,sim_model,P)
% Computes the acquisition functions using raw numerical simulation. This
% is slow and intended for testing purposes only. 

li = Inf;
nnew = size(th_new,1);
acq = NaN(nnew,1);
for i = 1:nnew
    f = @(ep)EV_obj(ep,method,u,th_new(i,:),th_cur,th_pr,th_tr,y_tr,sigma_tr,th_grid,gp,gpo_opt,gp_opt,P,sim_model);
    acq(i) = integral(f,-li,li,'AbsTol',1e-4,'RelTol',1e-4);
end
%acq = log(max(eps,acq));
end


function in = EV_obj(ep,method,u,th_new,th_cur,th_pr,th_tr,y_tr,sigma_tr,th_grid,gp,gpo_opt,gp_opt,P,sim_model)

[mt,s2t] = gp_pred_fast(gp,th_tr,y_tr,th_new,P);
sigma2_new = gp_noise_model_var(gp,th_tr,sigma_tr,th_new,gp_opt);
y_new = mt + sqrt(s2t + sigma2_new).*ep;

logb = logb_term(th_cur,th_pr,sim_model);

% new point (y_new,th_new)
th_trnew = [th_tr;th_new];
sigma_trnew = [sigma_tr;sqrt(sigma2_new)];
sep = numel(ep);
in = NaN(size(ep));
for i = 1:sep
    y_trnew = [y_tr;y_new(i)];
    gpnew = reinit_gp(gp, gpo_opt, gp_opt, th_grid, y_trnew, th_trnew, sigma_trnew);
    Pnew = precompute_gp_pred(gpnew, th_trnew, y_trnew, gp_opt);
    
    [mtnew,s2tnew] = gp_pred_fast(gpnew,th_trnew,y_trnew,[th_cur;th_pr],Pnew);
    ctnew = gp_pred_cov_fast(gpnew,th_trnew,y_trnew,th_cur,th_pr,Pnew,0);
    
    mutnew = mtnew(2) - mtnew(1) + logb;
    sigmatnew = sqrt(max(0,s2tnew(1) + s2tnew(2) - 2*ctnew));
    if method == 1 % expected variance of \alpha
        %ini = varf_alpha_numint(mutnew,sigmatnew);
        [~,~,ini] = alpha_stats(mtnew(1),mtnew(2),s2tnew(1),s2tnew(2),ctnew,logb,[]);
    elseif method == 2 % expected variance of p
        ini = varf_p(mutnew,sigmatnew,u);
    elseif method == 3
        ini = p1mp(mutnew,sigmatnew,u);
    else
        error('Incorrect method.');
    end
    in(i) = 1/sqrt(2*pi)*exp(log(ini) - 0.5*ep(i)^2);
end
end


function v = varf_alpha_numint(mut,sigmat)
% Computes the variance of \alpha using numerical integration.

% This is fast and we use default values for MATLAB integral
f2 = @(alpha)alpha.*normcdf_fast((mut-log(alpha))./sigmat);
m2 = integral(f2,0,1);
f1 = @(alpha)normcdf_fast((mut-log(alpha))./sigmat);
m1 = integral(f1,0,1);
v = max(0, 2*m2 - m1.^2);
end


function v = varf_p(mut,sigmat,u)
% Computes the variance of p (indicator rv of whether \alpha<u or not)

logu = log(u);
ps(1,:) = normcdf_fast((logu-mut)./sigmat);
ps(2,:) = 1 - ps(1,:);
v = max(0,ps(1,:).*ps(2,:));
end


function minp1p = p1mp(mut,sigmat,u)
% Computes min{p,(1-p)} (where p is the indicator rv of whether \alpha<u or not)
% which equals the conditional error. 

logu = log(u);
minp1p = max(0,normcdf_fast(-abs(logu-mut)./sigmat));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h = viz_acq(f,f_alt,th_grid,f_ths,ths,inits,th_range_loc,th_cur,th_pr,th_tr,y_tr,opt,interp_acq,fignr)
% Wrapper for plotting the acq function surface

grid_size2d = 70;

d = th_grid.dim;
acq_grid.opt_f = f_ths;
acq_grid.opt_th = ths;
acq_grid.inits = inits;
acq_grid.th_range_loc = th_range_loc;

% acq surface is only to be plotted if dim <= 2
if d == 1
    acq_grid.acq = f(th_grid.theta(:));
    if ~isempty(f_alt)
        acq_grid.acq_alt = f_alt(th_grid.theta(:));
    end
elseif d == 2
    % Plotting in full 2d grid may be costly for some acq functions so evaluate it in
    % a smaller 2d grid and then interpolate to the larger grid. Computing exactly
    % here is unnecessary because this computation is only for visualisation.
    
    if ~interp_acq
        acq_grid.acq = vec_to_grid_matrix(f(th_grid.theta2d'), th_grid);
        if ~isempty(f_alt)
            acq_grid.acq_alt = vec_to_grid_matrix(f_alt(th_grid.theta2d'), th_grid);
        end
    else
        interp_method = 'cubic';
        th_s = theta_grid(th_grid.range, grid_size2d);
        [th1_s,th2_s] = meshgrid(th_s.theta(1,:),th_s.theta(2,:));
        [th1,th2] = meshgrid(th_grid.theta(1,:),th_grid.theta(2,:));
        acq_s = vec_to_grid_matrix(f(th_s.theta2d'), th_s);
        acq_grid.acq = interp2(th1_s, th2_s, acq_s, th1, th2, interp_method);
        if ~isempty(f_alt)
            acq_s = vec_to_grid_matrix(f_alt(th_s.theta2d'), th_s);
            acq_grid.acq_alt = interp2(th1_s, th2_s, acq_s, th1, th2, interp_method);
        end
    end
else
    return;
end
h = figure(fignr);
clf;
plot_acq(th_grid, acq_grid, th_cur, th_pr, th_tr, y_tr, opt.acq_opt);
drawnow;
end


