function [gp,gp_optim_opt] = fit_gp_model(gp, gp_optim_opt, gp_opt, th_grid, ...
    y_tr, th_tr, sigma_tr)
% Optimizes the GP hyperparameters using GPstuff 'gp_optim'.
%
% NOTE: GP_OPTIM_OPT are the settings for GP hyperparam optimisation, 
% GP_OPT are the general settings structure used for specifying GP model settings etc.

if nargin < 7
    sigma_tr = [];
end

%% Set up the GP model priors, initial settings etc.
init = isempty(gp);
if init
    [gp,gp_optim_opt] = init_gp(gp_opt, th_grid, y_tr, th_tr, sigma_tr);
else
    % This is needed for noise model GP as the GP structure needs to be
    % updated to account for the new data points with sigma_tr values;
    % in standard GP case adding new data does not affect the information
    % in GP structure.
    
    % If "data dependent priors" were used, then this would need to be
    % changed. 
    [gp,gp_optim_opt] = reinit_gp(gp, gp_optim_opt, gp_opt, th_grid, ...
        y_tr, th_tr, sigma_tr);
end

%% fit GP hyperparameters
try
    [gp,fval,exitflag] = gp_optim(gp, th_tr, y_tr(:), 'opt', gp_optim_opt);
catch
    % Optimisation failed, try again by starting from default settings
    [gp,gp_optim_opt] = init_gp(gp_opt, th_grid, y_tr, th_tr, sigma_tr);
    [gp,fval,exitflag] = gp_optim(gp, th_tr, y_tr(:), 'opt', gp_optim_opt);
    % If this still fails, let the error happen 
    
    %error('Optimisation of GP hypers failed due to code error.');
    %warning('Optimisation of GP hypers failed due to code error.');
    %keyboard;
end

%% debugging
if ~strcmp(gp_opt.display_type,'off')
    [w,s] = gp_pak(gp);
    disp('Hyperparameters:');
    disp(s);
    %disp('Values of log(hyperparameters):')
    %disp(w);
    disp('Values of hyperparameters:');
    disp(exp(w(1))); 
    if gp_opt.noise_model
        disp(exp(w(2:end)));
    else
        disp(exp(w(2:end-1)));
        disp(exp(w(end)));
    end
    disp(' ');
end
end



