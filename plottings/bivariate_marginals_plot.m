function h = bivariate_marginals_plot(th_grid, sim_model, epost_gp, epost_mcmc, fignr, th_init)
% Plots all bivariate marginal densities (GP,MCMC-based and the true one)
% nicely to the same figure.

FULLRANGE = 1; % if 1, plots full range of the parameters in the figures
TRUE_PARAM = 1; % if 1, plots the true parameter to the figures
MAX_S = 500; % maximum number of samples to plot
MARSIZE = 3.5; % size of marker of true value
PTSIZE = 4; % size of the sample points

if nargin < 6
    th_init = [];
end
h = [];
d = sim_model.dim;
if d <= 2
    return;
end
h = figure(fignr);
set(gcf,'Position',[50 1000 700 700]);
subplot = @(m,n,p)subtightplot(m, n, p, [0.055 0.055], [0.2 0.1], [0.1 0.1]);

samples_gp = thin_for_plotting(epost_gp.samples,MAX_S);
samples_mcmc = thin_for_plotting(epost_mcmc.samples,MAX_S);
plot_true = isfield(sim_model,'samples') && ~isempty(sim_model.samples);
if plot_true
    samples_true = thin_for_plotting(sim_model.samples,MAX_S);
end

for i = 1:d
    for j = 1:d
        if j == i % diagonal -> draw marginal densities
            subplot(d,d,(i-1)*d+j);
            hold on;
            plot(th_grid.theta(i,:),epost_gp.epost(i,:),'-k');
            plot(th_grid.theta(i,:),epost_mcmc.epost(i,:),'-b');
            if plot_true
                plot(th_grid.theta(i,:),sim_model.true_post_pdf(i,:),'-r');
            end  
            if TRUE_PARAM && isfield(sim_model,'true_theta')
                plot(sim_model.true_theta(i),0,'dr','MarkerSize',MARSIZE,'MarkerFaceColor','r');
            end
            hold off;
            if isfield(sim_model,'theta_names')
                xlabel(sim_model.theta_names{i}); 
            end
            if FULLRANGE
                xlim([th_grid.range(i,1),th_grid.range(i,2)]);
            end
            set(gca,'ytick',[]);
        elseif j > i % upper diagonal -> approx bivariate marginal posteriors
            subplot(d,d,(i-1)*d+j);
            hold on;
            
            if 1
                hgp = plot(samples_gp(:,j),samples_gp(:,i),'.k','MarkerSize',PTSIZE);
                ham = plot(samples_mcmc(:,j),samples_mcmc(:,i),'.b','MarkerSize',PTSIZE);
                %hgp = scatter(samples_gp(:,j),samples_gp(:,i),'.k');
                %ham = scatter(samples_mcmc(:,j),samples_mcmc(:,i),PTSIZE,'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
            else
                matlabkde = 1; % default bandwidth does not work well...
                pind = [j i];
                th_grid_ij = theta_grid(th_grid.range(pind,:),size(th_grid.theta,2));
                gp_post = vec_to_grid_matrix(kde_for_abc(th_grid_ij,samples_gp(:,pind),matlabkde),th_grid_ij);
                contour(th_grid_ij.theta(1,:), th_grid_ij.theta(2,:), gp_post, 1, '-b');
                % TBD...
            end
            
            if TRUE_PARAM && isfield(sim_model,'true_theta')
                plot(sim_model.true_theta(j),sim_model.true_theta(i),'dr',...
                    'MarkerSize',MARSIZE,'MarkerFaceColor','r');
            end
            if ~isempty(th_init)
                plot(th_init(j),th_init(i),'sm','MarkerSize',MARSIZE,'MarkerFaceColor','m');
            end
            hold off;
            if isfield(sim_model,'theta_names')
                xlabel(sim_model.theta_names{j}); 
                ylabel(sim_model.theta_names{i}); 
            end
            if FULLRANGE
                xlim([th_grid.range(j,1),th_grid.range(j,2)]);
                ylim([th_grid.range(i,1),th_grid.range(i,2)]);
            end
        elseif plot_true % lower diagonal -> true bivariate marginal posteriors
            subplot(d,d,(i-1)*d+j);
            hold on;
            plot(samples_true(:,i),samples_true(:,j),'.r','MarkerSize',PTSIZE);
            if TRUE_PARAM && isfield(sim_model,'true_theta')
                plot(sim_model.true_theta(i),sim_model.true_theta(j),'dk',...
                    'MarkerSize',MARSIZE,'MarkerFaceColor','r');
            end
            if ~isempty(th_init)
                plot(th_init(i),th_init(j),'sm','MarkerSize',MARSIZE,'MarkerFaceColor','m');
            end
            hold off;
            if isfield(sim_model,'theta_names')
                xlabel(sim_model.theta_names{i}); 
                ylabel(sim_model.theta_names{j}); 
            end
            if FULLRANGE
                xlim([th_grid.range(i,1),th_grid.range(i,2)]);
                ylim([th_grid.range(j,1),th_grid.range(j,2)]);
            end
        end
        box on;
    end
end
end

function samples1 = thin_for_plotting(samples,MAX_S)
MAX_S = min(size(samples,1),MAX_S);
samples1 = samples(floor(linspace(1,size(samples,1),MAX_S)),:);
end

