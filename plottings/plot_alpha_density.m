function [] = plot_alpha_density(a_grid,a_cdf)
% Plots the cdf of the \alpha for GP-MH. 

plot(a_grid,a_cdf,'-k');
xlabel('\alpha');
ylabel('F(\alpha)');
xlim([0,1]); ylim([0,1]);
end