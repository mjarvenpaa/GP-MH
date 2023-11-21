function samples_thinned = thin_samples(samples,m)
% A simple function that thins the given sample set 'samples' to length 'm'. 
% 'Samples' must be n x dim matrix, where n is the amount of samples.

% some very basic input checkings to avoid common errors
if isempty(samples) || numel(m) ~= 1 || ~isfinite(m) || m <= 0 
    error('Incorrect input to thin_samples.');
end
[ns,d] = size(samples);
if d > 100
    error('Probably incorrect input to thin_samples.');
end

% thin to desired length
n_final = min(m,ns);
samples_thinned = samples(floor(linspace(1,ns,n_final)),:);
end

