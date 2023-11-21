function [a_me,a_med,a_va,a_stdev,a_cdf,a_grid] = alpha_stats(mt,mt_pr,s2t,s2t_pr,ct,logb,opt)
% Computes the mean, median, variance, stdev of the distribution of M-H acceptance
% check value \alpha = min{ratio,1}.
% Both exact and simulation based computation are included. 
%
% mt: GP mean at current \theta
% mt_pr: GP mean at proposed \theta'
% s2t: GP variance at current \theta
% s2t_pr: GP variance at proposed \theta'
% ct: GP covariance between \theta and \theta'
% logb: log ratio term involving prior and proposal q densities

mut = mt_pr - mt + logb;
sigma2t = s2t_pr + s2t - 2*ct;
sigmat = sqrt(max(0,sigma2t));

if isfield(opt,'method') && strcmp(opt.method,'sim')
    % compute using simulation
    [a_me,a_med,a_va] = alpha_stats_simul(mt,mt_pr,s2t,s2t_pr,ct,logb,opt);
else
    % use exact formulas (slightly reformulated to avoid numerical issues)
    a_me = normcdf_fast(mut./sigmat) + ...
        exp(mut+sigma2t/2 + log_normcdf_fast(-(mut+sigma2t)./sigmat));
    a_va = normcdf_fast(mut./sigmat) + ...
        exp(2*(mut+sigma2t) + log_normcdf_fast(-(mut+2*sigma2t)./sigmat))...
        - a_me.^2;
    a_va = max(0,a_va);
    a_med = min(exp(mut),1);
    
%     % use numerical integration:
%     alpha_grid = linspace(1/opt.grid_pts,1,opt.grid_pts);
%     v1_grid = normcdf_fast((mut - log(alpha_grid))./sigmat);
%     a_me = trapz(alpha_grid,v1_grid);
%     a_me = max(a_me,0); % as \alpha >= 0
%     a_med = []; % TBD!
%     
%     % variance
%     v2_grid = alpha_grid.*v1_grid;
%     a_va = 2*trapz(alpha_grid,v2_grid) - a_me.^2;
%     a_va = max(a_va,0); 
end
a_stdev = sqrt(a_va);

% compute the cdf of \alpha
if nargout > 4
    a_grid = linspace(1/opt.grid_pts,1,opt.grid_pts);
    a_cdf = normcdf_fast((log(a_grid) - mut)./sigmat);
    a_cdf(a_grid<0) = 0;
    a_cdf(a_grid>=1) = 1; % discontinuity point at \alpha == 1
end
end


function [a_me,a_med,a_va] = alpha_stats_simul(mt,mt_pr,s2t,s2t_pr,ct,logb,opt)
% Computes the stats by raw simulation.

logr = mvnrnd([mt mt_pr],[s2t ct;ct s2t_pr],opt.n_simul);
logrd = logr(:,2)-logr(:,1);
alpha_sim = exp(min(0,logb + logrd)); % simulated \alphas

a_me = mean(alpha_sim);
a_med = median(alpha_sim);
a_va = var(alpha_sim);
end


