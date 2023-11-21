function mcmc_diag = basic_mcmc_diagnostics(samples,display_type)
% Wrapper for some basic assessment of the convergence of MCMC.

mcmc_diag = [];
[ns,npar,nchains] = size(samples);
if ns <= 1 || nchains > 100
    error('Incorrect dimensions for samples matrix.');
elseif ns <= 20
    return; % too few samples to even run the diagnostics
end

% psrf is taken from GPstuff/diag
[R,neff,Vh,W,B,tau,thin] = psrf(samples);
mcmc_diag.R = R;
mcmc_diag.neff = neff;
mcmc_diag.is_converged = (max(abs(R-1)) < 0.1); % one number summary of convergence assessment
if 0 && ~strcmp(display_type,'off')
    disp(' ');
    disp('***');
    disp(['nr chains = ', num2str(nchains)]);
    disp(['R = ',num2str(R)]);
    disp(['"is_converged" = ', num2str(mcmc_diag.is_converged)]);
    disp(['neff = ', num2str(neff)]);
end
if mcmc_diag.is_converged ~= 1
    warning('Convergence not reached.');
end
end

