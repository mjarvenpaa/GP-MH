function [] = demo_bounds_uncond_error()
% Analysis of the *conditional error*. These are deterministic bounds. 
% ***Used to create figures for v1/v2 of the paper (Figure 3b,c in v1).***

close all;

%% 1/3: check the derived exact formula for the unconditional error
if 0
    sigma_t = 1;
    mu_t = linspace(-10,10,10000);
    n_u = 1000;
    ugrid = linspace(1/n_u,1,n_u);
    
    uce_int = NaN(length(mu_t),1);
    for i = 1:length(mu_t)
        % integral in the error formula computed numerically
        uce_int(i) = trapz(ugrid,normcdf_fast(abs(mu_t(i)-log(ugrid))./sigma_t));
    end
    uce_int = 1 - uce_int;
    
    uce_exact = unconderr(sigma_t,mu_t); % exact formula
    
    figure(1);
    hold on;
    plot(mu_t,uce_int,'-r');
    plot(mu_t,uce_exact,'-b');
    hold off;
    box on;
    xlabel('\mu_t');
    ylabel('unconderr');
    legend(['max \mu_t:', num2str(mu_t(uce_exact==max(uce_exact)))]);
    
    % Numerically computed seems to match the analytical formula as it should!
end


%root = '../results/bounds/'; % for v1 of the paper
root = '../results/bounds_v2/'; % for v2  of the paper

%% 2/3: plots *worst case unconditional error* **figure to the paper**
%% Plotted using dashed lines now in v2 and sigma_n also in reserved order. 
if 1
    save_fig = 1;
    
    %n = 1:1:200;
    n = 1:1:100;
    sigma_n = flip([0.5 1 2]); % Note: Changed the order for v2
    % deterministic bound, no epsilon (e). Also no u as we integrate over it
    
    ns = numel(sigma_n);
    nn = numel(n);
    mu_t_gr = linspace(-0.7,0,10000); % we simpy maximise in a dense grid
    mu_t_gr = repmat(mu_t_gr(:),1,nn);
    res_wc = NaN(nn,ns);
    
    % worst case error
    for i = 1:ns
        ci = 2*sqrt(2)*sigma_n(i)./repmat(sqrt(n),length(mu_t_gr),1);
        ri = normcdf_fast(mu_t_gr./ci) + exp(mu_t_gr + 0.5*ci.^2).*...
            (normcdf_fast(-(mu_t_gr+ci.^2)./ci) - 2*normcdf_fast(-ci));
        res_wc(:,i) = max(ri);
    end
    
    figure(2);
    set(gcf,'Position',[50 1000 400 300]);
    lw = 1.3;
    hold on;
    for i = 1:ns
        h(i) = plot(n,res_wc(:,i),'Color',col_get(i),'LineStyle','--','Linewidth',lw);
        le{i} = ['$\bar{\sigma}_n=',num2str(sigma_n(i)),'$'];
    end
    hold off;
    box on; grid on;
    legend(h,le,'interpreter','latex')
    stry = 'Uncond. error';
    ylabel(stry,'interpreter','latex');
    %xlabel('n');
    xlabel('b');
    title('(b)');
    
    %% save figure
    if save_fig
        fn = [root,'/fig_wc_uncond_error'];
        my_export_fig(fn,'-transparent','-pdf');
    end
end

%% 3/3: plots unconditional error **figure to the paper**
if 1
    save_fig = 1;
    
    %n = 1:1:200;
    n = 1:1:100;
    e = [0.1:0.1:0.4];
    sigma_n = 1;
    k = 5; % dim of target
    %s = 0.2; % proposal scaling
    s = sqrt(2.4^2/k);
    % No u as we integrate over it
    
    nn = numel(n);
    ne = numel(e);
    res = NaN(nn,ne);
    
    m = 200000;
    mu_t = rand_mut(m,k,s);
    mu_t = repmat(mu_t(:)',nn,1);
    cn = 2*sqrt(2)*sigma_n./repmat(sqrt(n(:)),1,size(mu_t,2));
    r = unconderr(cn,mu_t); % note that cn in place of sigma_t -> upper bound
    for i = 1:ne
        res(:,i) = mean(r>=e(i),2);
    end
    
    figure(3);
    set(gcf,'Position',[500 1000 400 300]);
    lw = 1.3;
    hold on;
    for i = 1:ne
        h(i) = plot(n,res(:,i),'Color',col_get(i),'LineStyle','-','Linewidth',lw);
        le{i} = ['$\varepsilon=',num2str(e(i)),'$'];
    end
    hold off;
    box on; grid on;
    ylim([0,0.7]);
    legend(h,le,'interpreter','latex')
    stry = '$\bf P \rm (\textnormal{Uncond. error} \geq \varepsilon)$';
    %stry = 'Prob. uncond. error larger than \epsilon';
    ylabel(stry,'interpreter','latex');
    %xlabel('n');
    xlabel('b');
    title('(c)');
    
    %% save figure
    if save_fig
        fn = [root,'/fig_uncond_error'];
        my_export_fig(fn,'-transparent','-pdf');
    end
end
end


