function mu_t = rand_mut(m,k,s)
% Generate m random variables for the density of \mu_t with dimensionality 
% k when proposal scaling is s 
% (i.e. q(\theta'|\theta)=Normal(\theta'|\theta,s^2\Sigma))

mu_t = zeros(1,m);
for i = 1:k
    r = randn(2,m);
    r = chol([1 1; 1 1 + s^2],'lower')*r;
    mu_t = mu_t + 0.5*(r(1,:).^2 - r(2,:).^2);
end
end

