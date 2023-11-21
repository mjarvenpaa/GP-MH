function [llik,bootvar] = noisy_loglik_estim(sim_model, opt, theta, method, use_boot)
% Outputs a (noisy) log-likelihood evaluation at parameter 'theta' for a 
% given model specified as 'sim_model'. 
% Either a given function handle in 'sim_model' is directly used for this 
% or N datasets are simulated and LFI approximation 'method' (currently 
% only SL is implemented though) is then used. In the latter case an 
% estimate of the variance of the loglik value is also computed using the 
% bootstrap.

if nargin < 5
    use_boot = 0;
end
bootvar = NaN;

%% (noisy) loglik is computed without simulation runs:
if strcmp(method,'exact')
    llik = sim_model.loglik_eval(theta);
    return;
end

%% LFI case:
N = opt.N;
summaries_th = NaN(sim_model.summary_dim,N); % each column one summary vector

% repeated sampling for SL (this could be parallelized in practice)
for i = 1:N
    data_i = sim_model.gen(theta, sim_model.n_data);
    summaries_th(:,i) = sim_model.comp_summaries(data_i,sim_model.data);
end

if strcmpi(method,'sl')
    llik = eval_sl_loglik(sim_model.summary_true, summaries_th, opt.estimator);
    if ~isfinite(llik)
        msg = ['SL computation resulted non-finite value at ', num2str(theta(:)')];
        error(msg);
    end
    
    if use_boot && nargout > 1
        % uses bootstrap to compute estimate of the variance of the loglik estimate
        if isfield(opt,'bootn')
            bootn = opt.bootn;
        else
            bootn = 2000;
        end
        llik_boot = NaN(bootn,1);
        for i = 1:bootn
            boot_inds = randsample(N,N,1); % sample with replacement
            s_thi = summaries_th(:,boot_inds);
            llik_boot(i) = eval_sl_loglik(sim_model.summary_true, s_thi, opt.estimator);
        end
        llik_boot = llik_boot(~isnan(llik_boot));
        if isfield(opt,'robust_bootvar') && opt.robust_bootvar == 1
            % use 1.48*MAD as a robust estimator for standard deviation:
            % (multiplication with 1.48 comes from a Gaussian assumption)
            bootvar = (1.48*median(abs(llik_boot-median(llik_boot))))^2;
        else
            bootvar = var(llik_boot);
        end
        if ~isfinite(bootvar)
            msg = ['Computing bootstrapped variance of SL failed at ', num2str(theta(:)')];
            error(msg);
        end
        
        if 0
            % show bootstrap samples
            figure(25);
            histogram(llik_boot,30);
            title('Bootstrap samples of SL loglik');
            xlabel('SL loglik values');
            xlim([min(llik_boot) max(llik_boot)]);
            pause;
        end
    end

else % lfire or other methods
    error('Not implemented.');
end

% %% The magnitude of the log-likelihood near boundaries can be large for some models which
% %% can make fitting the GP problematic. We deal with such cases by lower bounding the loglik
% %trunclik = -Inf;
% trunclik = -10^5;
% llik = max(trunclik,llik);
% bootvar(llik<=trunclik) = 100^2;
end


