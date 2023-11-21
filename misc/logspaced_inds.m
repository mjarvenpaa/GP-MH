function v = logspaced_inds(a,b,n,n0)
% Returns (at most) n distinct integers between a and b which are 
% approximately on log scale and so that the points 1:(n0-1) have been 
% neglected. This function is intended for determining those iterations
% when the approx. MCMC samples etc. are to be saved for which case an
% approximate solution suffices.
% NOTE: May not work as intended if e.g. b is very large. If n is large 
% enough, the output cannot contain n distinct integers in which case a 
% smaller number of points (max(n0,a):b) is returned instead. 

if nargin < 4
    n0 = 5; % no results are saved until 5th iteration anyway
elseif nargin == 4 && (numel(n0) ~= 1 || n0 <= 0)
    error('Incorrect input to logspaced_pts.');
end
if numel(a) ~= 1 || a <= 0 || numel(b) ~= 1 || max(n0,a) >= b || numel(n) ~= 1 || n <= 0
    error('Incorrect input to logspaced_pts.');
end
a = max(n0,a);
v = [];
m = n;
while ~(length(v) >= n || isequal(v,a:b) || m-n >= 10^3)
    % Main idea: we increase m from n until we either get n points after 
    % rounding or after some other inevitable stopping condition is met.
    v = unique(floor(logspace(log10(a),log10(b),m)));
    m = m + 1;
end
%disp(num2str(m-1));
if length(v) > n % if too many points gets included
    v = [v(1:(n-1)), v(end)];
end
% ensure that the last point is exactly b (and e.g. not b-1)
if length(v) > 1
    v(end) = b;
end
end
