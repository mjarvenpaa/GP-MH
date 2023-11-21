function p = log_normcdf_fast(z)
% This code is much faster than log(normcdf) if called repeatedly

p = -log(2) + log(erfc(-z ./ sqrt(2)));

% If z << 0 so that p becomes -Inf due to underflow we use the following approximation to 
% avoid -Inf cases.
% There is a little bit of discontinuity near those z, where the approximation below is 
% first used but this shouldn't cause any issues in practice.
% The approximation below is likely accurate also for larger z than z<-38.4
ind = (z<-38.4) | ((z<0) & ~isfinite(p));
p(ind) = -0.5*log(2*pi) - 0.5*z(ind).^2 - log(-z(ind));
end
