function ce = conderr(sigmat, mut, u, logoutput)
% Computes the conditional error (or logarithm of conditional error if
% logoutput == 1). sigmat and mut need to be arrays with the same size or
% either one needs to be a scalar. 

if nargin < 4
    logoutput = 0;
end
if ~(numel(sigmat) == 1 || numel(mut) == 1 || all(size(sigmat)==size(mut)))
    error('Incorrect inputs for conderr.');
end

if logoutput
    ce = log_normcdf_fast(-abs(mut-log(u))./sigmat);
else
    ce = normcdf_fast(-abs(mut-log(u))./sigmat);
end
end
