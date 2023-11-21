function res = basic_post_sample_stats(samples,d)
% Computes sample mean and covariance of posterior samples.

res.post_mean = mean(samples);
res.post_cov = cov(samples);

% just in case, check dimensions
if nargin > 1 && (numel(res.post_mean) ~= d || numel(res.post_cov) ~= d^2)
    error('Incorrect dimensions with post summaries.');
end
end