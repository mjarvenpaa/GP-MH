function results = GPMH_main(sim_model, th_grid, opt)
% Runs the main GP-MH algorithm, visualises its output and compares the 
% result to the ground truth posterior if available. Extra GP-based MCMC 
% sampling is also run at the end or during specified iterations so that 
% corresponding MH-BLFI results are also obtained simultaneously. 

% TODO: Is there a better way to deal with occasional failed/nonsensical loglik outputs
% TODO: Take advantage of repeated evals at the same points in SL case

%% initialize
d = sim_model.dim;
n = opt.n; n0 = opt.init_n;
th_mcmc0 = opt.mcmc.th_init(:)'; % initial point
th_mcmc = NaN(n,d); % n x dim, samples by GP-MH
th_tr = []; % t x dim, training points, obtained using design criterion (aka acquisition function)
loglik_tr = []; % t x 1, log-likelihood evaluations at training points
sigma_tr = []; % t x 1
results.iter = cell(n,1); % results are saved here
evals = zeros(1,n+1); evals_invalid = zeros(1,n+1); 
errs = Inf(1,n);
acc_nr = 0;
n_tr_old = 0;

%% some additional settings
%max_invalid_logliks = 20;
n0_max_tries = 5*n0;

%% gather evaluations for initial GP fitting
q_Sigma0 = opt.mcmc.q_Sigma;
if isfield(opt.mcmc,'q_Sigma0') && ~isempty(opt.mcmc.q_Sigma0)
    q_Sigma0 = opt.mcmc.q_Sigma0;
end
for i0 = 1:n0_max_tries
    th0 = q_init(th_mcmc0, q_Sigma0, 1, th_grid);
    [loglik0,bootvar0] = noisy_loglik_estim(sim_model, opt.lik_opt.sl,...
        th0, opt.lik_opt.method, opt.gp_opt.noise_model);
    if check_loglik(loglik0,opt,bootvar0)
        % valid loglik output, add it to initial training data
        th_tr = [th_tr;th0(:)'];
        loglik_tr = [loglik_tr;loglik0];
        sigma_tr = [sigma_tr;sqrt(bootvar0)];
    else
        evals_invalid(1) = evals_invalid(1) + 1;
        if strcmp(opt.misc.display_type,'on')
            disp(['!!! Invalid loglik: ',num2str(loglik0),', bootstdev: ',num2str(sqrt(bootvar0))]);
        end
    end
    if length(loglik_tr) == n0
        break;
    elseif i0 == n0_max_tries
        error(['Could not obtain ',num2str(n0),' valid initial evaluations.'...
            ' Better initialisation might be needed.']);
        % In principle we could e.g. try to decrease opt.mcmc.q_Sigma but
        % this is not implemented here. 
    end
end
evals(1) = i0;

%% fit initial GP
[gp,gpo_opt] = fit_gp_model([],[], opt.gp_opt, th_grid,loglik_tr, th_tr, sigma_tr);
P = precompute_gp_pred(gp, th_tr, loglik_tr, opt.gp_opt);

%% main iteration of GP-MH:
for i = 1:n
    pr_iter = print_iter_info(i,opt,0);
    
    %% Propose new point \theta'
    % We here use Gaussian proposal whose covariance matrix is adaptively updated
    if i == 1
        th_cur = th_mcmc0;
        [th_pr,C,mc] = q_gen([], [], th_mcmc0, 1, opt.mcmc);
    else
        th_cur = th_mcmc(i-1,:);
        [th_pr,C,mc] = q_gen(C, mc, th_mcmc(1:(i-1),:), 1, opt.mcmc);
    end
    u = rand(1);
    % In our implementation here we require the prior is some density 
    % truncated to a bounded rectangle specified in th_grid (uniform 
    % density is used in all of our examples in get_test_model.m) and we 
    % check here whether the proposed point th_pr is inside these prior 
    % bounds:
    ll_ok = inside_prior_bounds(th_grid, th_pr);
    
    erri = NaN;
    j = 1;
    while 1
        if ~ll_ok
            break;
        end
        
        %% compute the error in accept/reject decision based on current GP
        [mt,s2t] = gp_pred_fast(gp,th_tr,loglik_tr,th_cur,P);
        [mt_pr,s2t_pr] = gp_pred_fast(gp,th_tr,loglik_tr,th_pr,P);
        ct = gp_pred_cov_fast(gp,th_tr,loglik_tr,th_cur,th_pr,P,0);
        logb = logb_term(th_cur,th_pr,sim_model);
        
        % Note: Many of the computations below are actually not needed in GP-MH:
        [a_me,a_med,a_va,a_stdev] = alpha_stats(mt,mt_pr,s2t,s2t_pr,ct,logb,[]);
        [ps,p_me,p_va,p_stdev,ce,uce] = p_stats(u,mt,mt_pr,s2t,s2t_pr,ct,logb);
        % tol_eps is \epsilon in the paper
        if strcmp(opt.mcmc.tol_crit,'cond')
            erri = ce; tol_eps = opt.mcmc.tol_cond;
        elseif strcmp(opt.mcmc.tol_crit,'uncond')
            erri = uce; tol_eps = opt.mcmc.tol_uncond;
        else
            error('Incorrect criterion.');
        end
        
        %% check if error smaller than the given tolerance (\epsilon)
        if erri <= tol_eps
            break;
        end
        
        %% check if maximum evaluation budget is already used
        if sum(evals)+j-1 >= opt.maxacq || j-1 >= opt.maxacqiter
            break;
        end
        
%         %% some debug printing
%         if j == 1 && strcmp(opt.misc.display_type,'on')
%             if ~pr_iter
%                 pr_iter = print_iter_info(i,opt,1);
%             end
%             disp(['th_cur: ', num2str(th_cur(:)')]);
%             disp(['th_pr: ', num2str(th_pr(:)')]);
%             disp(['mt_cur, mt_pr, u: ',num2str([mt, mt_pr, u])]);
%         end
        
        %% too much uncertainty and budget not used yet -> acquire new loglik evals
        fignr = 1;
        th = acquire(i,u,th_grid,th_cur,th_pr,th_tr,loglik_tr,sigma_tr,...
            gp,gpo_opt,opt.gp_opt,P,sim_model,opt,fignr);
        [loglik,bootvar] = noisy_loglik_estim(sim_model, opt.lik_opt.sl, th,...
            opt.lik_opt.method, opt.gp_opt.noise_model);
        
        if ~check_loglik(loglik,opt,bootvar)
            % If the acquired point 'th' resulted invalid loglik, we check
            % if the point was current or proposed point. If it was the 
            % current point we stop. If it was the proposed point we reject
            % it and move on. In the tricky case where neither holds, we 
            % randomly choose either the current or proposed point, perform
            % a new loglik evaluation there and then proceed as above to 
            % (hopefully!) avoid the difficulty in the future.
            % TODO: Evaluate at both points yet proceed as above??
            evals_invalid(i+1) = evals_invalid(i+1) + 1;
            if strcmp(opt.misc.display_type,'on')
                disp(['!!! Invalid loglik: ',num2str(loglik),', bootstdev: ',num2str(sqrt(bootvar))]);
            end
            if eq_theta(th,th_cur)
                error(['Got stuck in invalid region, i=',num2str(i),', evals=',num2str(sum(evals)+j),...
                    ', loc: ',num2str(th(:)'),'. Better initialisation might be needed.']);
            elseif eq_theta(th,th_pr)
                ll_ok = 0;
                erri = NaN;
                break;
            end
            ind = (rand(1) > 0.5);
            th = ind*th_cur + (~ind)*th_pr;
            [loglik,bootvar] = noisy_loglik_estim(sim_model, opt.lik_opt.sl, th,...
                opt.lik_opt.method, opt.gp_opt.noise_model);
            j = j + 1;
            if ~check_loglik(loglik,opt,bootvar)
                evals_invalid(i+1) = evals_invalid(i+1) + 1;
                if strcmp(opt.misc.display_type,'on')
                    disp(['!!! Invalid loglik: ',num2str(loglik),', bootstdev: ',num2str(sqrt(bootvar))]);
                end
                if eq_theta(th,th_cur)
                    error(['Got stuck in invalid region, i=',num2str(i),', evals=',num2str(sum(evals)+j),...
                        ', loc: ',num2str(th(:)'),'. Better initialisation might be needed.']);
                end
                % Invalid loglik at the proposed point
                ll_ok = 0;
                erri = NaN;
                break;
            end
        end
        
        %% update GP training data and GP fit
        th_tr = [th_tr;th(:)'];
        loglik_tr = [loglik_tr;loglik];
        sigma_tr = [sigma_tr;sqrt(bootvar)];
        if sum(evals)+j <= 300 || rem(sum(evals)+j,10) == 1
            % GP hypers always updated until 300th iteration and every 10th
            % iteration thereafter
            [gp,gpo_opt] = fit_gp_model(gp, gpo_opt, opt.gp_opt, th_grid,...
                loglik_tr, th_tr, sigma_tr);
        else
            [gp,gpo_opt] = reinit_gp(gp, gpo_opt, opt.gp_opt, th_grid, ...
                loglik_tr, th_tr, sigma_tr); % does not update GP hypers
        end
        P = precompute_gp_pred(gp, th_tr, loglik_tr, opt.gp_opt);
        j = j + 1;
    end
    evals(i+1) = j-1; 
    errs(i) = erri;
    
%     %% Stop if too many subsequent invalid loglik evaluations
%     if i > max_invalid_logliks && all(evals_invalid(i+1-max_invalid_logliks:i+1)>0)
%         error(['Last ',num2str(max_invalid_logliks),' iterations all produced invalid'...
%             'loglik evaluation. Better initialisation might be needed.']);
%     end

    %% Acceptance test is now computed accurately enough (or budget is used)
    acc = (ll_ok && p_me >= 0.5); % whether proposed point is accepted
    if acc % accept
        th_mcmc(i,:) = th_pr;
        acc_nr = acc_nr + 1;
    else % reject
        th_mcmc(i,:) = th_cur;
    end
    
    %% some more debug printing
    if strcmp(opt.misc.display_type,'on') && pr_iter
        if acc; disp('ACCEPTED');
        elseif ~ll_ok; disp('REJECTED (out of prior bounds or evaluated loglik is invalid)');
        else; disp('REJECTED');
        end
        if ll_ok
            disp(['Estimated \alpha: ',num2str(a_med),', estimated error: ', num2str(erri)]);
        end
        disp(['Acceptance ratio: ', num2str(acc_nr/i)]);
        disp(['Loglik evals (this iteration/cumulative): ', num2str(j-1), '/', num2str(sum(evals))]);
    end
    
    %% compute current posterior estimates and compare to the ground truth
    if i >= 5 && (opt.viz >= 2 || (opt.viz == 1 && i == n) || (opt.viz == -1 && any(opt.misc.res_ind==i)))
        % Posterior estimate based on the approx. MCMC samples (GP-MH):
        % This quantity changes on each iteration i.
        [epost_mcmc,mcmc_diag_mcmc] = post_from_gp_mcmc(th_grid, sim_model, th_mcmc(1:i,:), opt);
        
        % Posterior estimate based on the current GP surrogate (MH-BLFI):
        % This quantity changed only if new loglik evaluation was collected.
        if length(loglik_tr) > n_tr_old
            [epost_gp,mcmc_diag_gp] = post_from_gp_surrogate(th_grid, sim_model, gp, opt.gp_opt, ...
                loglik_tr, th_tr, sigma_tr, P, opt.mcmcgp, opt.viz~=-1);
            n_tr_old = length(loglik_tr);
        elseif length(loglik_tr) == n_tr_old
            % GP fit has not changed; use the old values of 'epost_gp', 'mcmc_diag_gp' to avoid new MCMC run
        else
            error('Problem with GP training data size.'); % should not happen but check just in case
        end
        
        % assess convergence of both methods visually
        if opt.viz >= 1 && i == n
            fignr = 2;
            mcmc_chain_plot(th_mcmc(1:i,:),sim_model,1,opt,fignr);
            if ~isempty(epost_gp.samples) % MCMC used only for dim>2 case
                mcmc_chain_plot(epost_gp.samples,sim_model,0,opt,fignr+1);
            end
        end
        
        % compare to the ground truth (if available)
        res_i.res_mcmc = compare_to_true_baseline(th_grid, epost_mcmc, sim_model);
        res_i.res_gp = compare_to_true_baseline(th_grid, epost_gp, sim_model);
        res_i.mcmc_diag_mcmc = mcmc_diag_mcmc;
        res_i.mcmc_diag_gp = mcmc_diag_gp;
        res_i.evals = sum(evals); % saved also here for convenience
        if opt.ret_lvl >= 2 && i == opt.misc.res_ind(end)
            res_i.epost_mcmc = epost_mcmc;
            res_i.epost_gp = epost_gp;
        end
        results.iter{i} = res_i;
    end

    %% compare to the ground truth visually
    if i >= 5 && (opt.viz >= 2 || (opt.viz == 1 && i == n))
        fignr = 4;
        plot_comparison_fig(th_grid, loglik_tr, th_tr, th_mcmc(1:i,:), epost_gp, epost_mcmc, ...
            sim_model, res_i, opt, fignr);
        % plot marginals in a separate plot(s)
        bivariate_marginals_plot(th_grid, sim_model, epost_gp, epost_mcmc, fignr+1);
        drawnow;
    end
end

%% collect statistics and handle output
fignr = 6;
results.stats = gather_and_plot_stats(evals,evals_invalid,errs,acc_nr,opt,fignr);
results.opt = opt;
if opt.ret_lvl >= 2
    % additional output:
    results.loglik_tr = loglik_tr;
    results.th_tr = th_tr;
    results.sigma_tr = sigma_tr;
    %results.th_mcmc = th_mcmc; % all approx. MCMC samples by GP-MH
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [is_ok,ok_ind] = check_loglik(logliks,opt,vars)
% Returns 1 if the loglik values are ok (finite, not NaN, real-valued, not
% too small/large). Returns also the indices which show which evals were ok
% an which were not. If the corresponding noise variance estimates are also 
% provided, their values are also similarly checked.

ok_ind = (isfinite(logliks) & isreal(logliks));
if isfield(opt.lik_opt,'maxabsloglik')
    ok_ind = ok_ind & (abs(logliks)<=opt.lik_opt.maxabsloglik);
end
if nargin == 3 && opt.gp_opt.noise_model
    ok_ind = ok_ind & (isfinite(vars) & isreal(vars) & vars>=0);
    if isfield(opt.lik_opt,'maxbootsigma')
        ok_ind = ok_ind & (vars<=opt.lik_opt.maxbootsigma^2);
    end
end
is_ok = all(ok_ind(:));
end


function th = q_init(th_cur,q_Sigma,n,th_grid)
% Generate n values from the proposal q but so that the points are inside 
% prior bounds. This function is for obtaining evaluations for initial GP 
% fitting.

maxtries = 10000;
d = length(th_cur);
th = NaN(n,d);
for i = 1:n
    for j = 1:maxtries
        th_try = mvnrnd(th_cur(:)',q_Sigma,1);
        if inside_prior_bounds(th_grid,th_try)
            th(i,:) = th_try;
            break;
        end   
        if j == maxtries % should not happen unless q_Sigma is poorly chosen
            error('Could not generate from q so that prior bounds are satisfied.');
        end
    end
end
end


function [th,C,mc,ep] = q_gen(C,mc,ths,n,optmcmc)
% Generate n values from the (adaptive) Gaussian proposal q. 

[i,d] = size(ths);
th_cur = ths(end,:);
if optmcmc.q_adapt && i >= max(optmcmc.q_i0, 2*d) % adapt the covariance
    %C = cov(ths);
    if isempty(C)
        [C,mc] = covupd(ths,1); % taken from MCMC toolbox
    else
        [C,mc] = covupd(th_cur,1,C,mc,i-1); % taken from MCMC toolbox
    end
    % We first try the specified 'opt.ep' value. If this does not produce a valid 
    % cov matrix, we try increasing 'ep' (temporarily) so that we can still proceed. 
    sd = optmcmc.extra_scaling*2.4^2/d;
    epl = [optmcmc.ep,10.^(-10:-2)];
    for ep = epl
        q_Sigma1 = sd*(C + ep*eye(size(C)));
        [~,p] = chol(q_Sigma1,'lower');
        if p == 0
            q_Sigma = q_Sigma1;
            break;
        elseif ep == epl(end)
            error('New adaptive proposal cov could not be computed.');
        end
    end
else % fixed proposal; no adaptation
    ep = NaN;
    q_Sigma = optmcmc.q_Sigma;
end
th = mvnrnd(th_cur,q_Sigma,n);
end


function is_inside = inside_prior_bounds(th_grid,th)
% Returns 1 if the point 'th' is inside the specified prior bounds and 0
% otherwise.

is_inside = all(th_grid.range(:,1) <= th(:) & th(:) <= th_grid.range(:,2));
end


function pr = print_iter_info(i,opt,force_pr)
% Prints basic information about the progress of the algorithm on the screen.

freq = max(100,floor(opt.n/20));
pr = false;
% print for 10 first iterations and then after each 'freq' iterations:
if strcmp(opt.misc.display_type,'on') && (force_pr || i<=10 || i==opt.n || mod(i,freq)==0)
    if i > 1; disp(' '); end
    disp(['### iteration: ', num2str(i), '/', num2str(opt.n)]);
    pr = true;
end
end


function s = eq_theta(th1,th2)
% Test if two theta parameters are (essentially) the same.

s = max(abs(th1(:)-th2(:))) <= 1e-12;
end


