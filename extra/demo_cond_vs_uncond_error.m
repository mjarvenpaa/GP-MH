function [] = demo_cond_vs_uncond_error()
% Plots conditional error (with some u) and unconditional error both with
% fixed \sigma_t and as \mu_t is varied. 
% Note: This plots the errors, not their upper bounds. 
% ***Used to create one figure for v1/v2 of the paper (Figure 1 v1/v2).***

close all;

save_fig = 1;
%root = '../results/bounds/'; % for v1 of the paper
root = '../results/bounds_v2/'; % for v2  of the paper
u = [0.01,0.25,0.5,0.75,1];

sigma_ts = [10 5 1 0.1];
for sigma_t = sigma_ts
    mu_t = linspace(-30,20,5000);
    % some adjustments for nicer plottings:
    if sigma_t < 0.5
        mu_t = linspace(-6,1,5000);
    elseif sigma_t < 2
        mu_t = linspace(-9,4,5000);
    end
    
    % conditional error, depends on u:
    nu = numel(u);
    cei = cell(nu,1);
    for i = 1:nu
        cei{i} = conderr(sigma_t, mu_t, u(i));
    end
    
    % unconditional error:
    uce = unconderr(sigma_t,mu_t);
    
    figure(1); 
    clf;
    set(gcf,'Position',[50 1000 400 300]);
    hold on;
    leuc = cell(nu,1);
    lw = 1.2;
    for i = 1:nu
        plot(mu_t,cei{i},'Color',col_get(i),'Linewidth',lw);
        %leuc{i} = ['Cond. error, $u=',num2str(u(i)),'$'];
        % legend to show only u-values:
        leuc{i} = ['$u=',num2str(u(i)),'$'];
    end
    lec = 'Uncond. error';
    plot(mu_t,uce,'k--','Linewidth',lw);
    hold off;
    box on; grid on;
    xlabel('$\mu_t$','interpreter','latex');
    if sigma_t == sigma_ts(1)
        %ylabel('error');
        ylabel('Conditional/unconditional error');
    else
        %ylabel('Conditional/unconditional error','color','white');
        ylabel('Conditional/unconditional error');
    end
    xlim([mu_t(1) mu_t(end)]);
    if sigma_t == sigma_ts(1)
        %legend([leuc;lec],'interpreter','latex','FontSize',7.5,... %'color','none',...
        %    'Location','Northwest');
        legend(leuc,'interpreter','latex','FontSize',9,'Location','Northwest');
    end
    if 1
        title(['$\sigma_t=',num2str(sigma_t),'$'],'interpreter','latex');
    end
    
    if save_fig
        fn = [root,'fig_cond_vs_uncond_error','_sigma_t',num2str(sigma_t)];
        my_export_fig(fn,'-transparent','-pdf');
    end
end
end

