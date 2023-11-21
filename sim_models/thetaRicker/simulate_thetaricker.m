function y = simulate_thetaricker(param,N,T)
% Simulates one data set with length T from the theta-Ricker model.
% If theta==1 and K==r0, then the special case of (scaled) Ricker model, 
% which is often used as an ABC benchmark problem, is obtained.
%
% INPUT:
% N - the starting population (equal to 1 in our application)
% T - the length of the data set
% param - parameters: log(r) i.e. r0, theta, K, phi, sigma_e
%                 OR: log(r) i.e. r0, theta, phi, sigma_e
%
% OUTPUT:
% y - the simulated data set

if numel(param) == 4 % K==r
    r0 = param(1);
    theta = param(2);
    phi = param(3);
    sigma_e = param(4);
    K = r0;
elseif numel(param) == 5 % K is given as an additional parameter
    r0 = param(1);
    theta = param(2);
    K = param(3);
    phi = param(4);
    sigma_e = param(5);
else
    error('Incorrect parameters for theta-Ricker.');
end

Ns = [N; zeros(T,1)];
srn = sigma_e*randn(T,1);
for t = 1:T
    g_t = r0*(1-(Ns(t)/K)^theta);
    Ns(t+1) = Ns(t)*exp(g_t + srn(t)); % population size
end
y = poissrnd_fast(phi*Ns(2:end)); % the actual observed random variable

end



