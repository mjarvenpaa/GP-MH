function h = mcmc_chain_plot(samples,sim_model,is_amcmc,opt,fignr)
% Plots the (combined) MCMC chain for e.g. visual convergence assessment.

plot_burnin_sep = 1; % whether to plot burnin samples with separate colour

[ns,d] = size(samples);
if is_amcmc
    %titl = 'Approx. MCMC';
    titl = 'GP-MH';
else
    %titl = 'MCMC for GP-based approx.';
    %titl = 'MCMC (GP)';
    titl = 'MH-BLFI';
end

% check labels for plotting
if ~isfield(sim_model,'theta_names')
    names = cell(1,d); 
    for i = 1:d
        names{i} = ['\theta_',num2str(i)];
    end
else
    names = sim_model.theta_names;
end

h = figure(fignr);
clf;
if 0
    % use code from MCMC toolbox:
    mcmcplot(samples,1:d,names,'chainpanel'); 
    suptitle(titl);
else
    % my own simple plottings:
    for i = 1:min(d,10)
        subplot(min(d,10),1,i);
        hold on;
        if plot_burnin_sep && is_amcmc
            % plot burnin and the actual approx. MCMC samples separately
            bi = floor(opt.burninfrac*ns);
            plot(1:bi,samples(1:bi,i),'.b');
            plot(bi+1:ns,samples(bi+1:end,i),'.k');
        else
            % just plot all samples; here burnin is already neglected or is
            % not to be separated from the rest of the samples anyway
            plot(1:ns,samples(:,i),'.k');
        end
        hold off;
        box on;
        ylabel(names{i});
        %xlabel('iteration');
        xlim([1,ns]);
        if i == 1
            title(titl);
        end
    end
end
set(gcf,'Position',[50 600 500 min(d*200,800)]);
end

