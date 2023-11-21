function [ps,p_me,p_va,p_stdev,ce,uce] = p_stats(u,mt,mt_pr,s2t,s2t_pr,ct,logb)
% Computes the statistics of the indicator variable that is 0 if \alpha < u
% and 1 otherwise. It follows Bernoulli distribution given u and supposing
% the loglik follows the GP. Computes also the conditional and
% unconditional error which are denoted with 'ce' and 'uce', respectively.

mut = mt_pr - mt + logb;
sigma2t = s2t_pr + s2t - 2*ct;
sigmat = sqrt(max(0,sigma2t));

% pmf
ps = NaN(2,length(mut));
ps(1,:) = normcdf_fast((log(u)-mut)./sigmat);
ps(2,:) = 1 - ps(1,:);

% mean & var
p_me = ps(2,:);
p_va = max(0,ps(1,:).*ps(2,:));
p_stdev = sqrt(p_va);

% conditional error (which is the same as min of p and 1-p)
% Alternatively, we could compute this using 'conderr'-function.
ce = min(ps(1,:),ps(2,:));

% unconditional error
uce = unconderr(sigmat, mut);
% computation might fail due to numerics, check this:
if isempty(uce) || any(~isfinite(uce)) || any(~isreal(uce))
    error('Uncond.error comput. failed.');
end
end

