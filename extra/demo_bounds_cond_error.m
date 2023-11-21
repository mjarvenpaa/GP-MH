function [] = demo_bounds_cond_error()
% Error bounds for the *conditional error*. These are probabilistic bounds
% wrt. u and possibly also \mu_t
% ***Used to create figures for v1/v2 of the paper (Figure 3a in v1).***

close all;

if 0
    %% WORST CASE ERROR --- TABLE
    %% Note: This analysis is not in the paper and is kind of obsolete. 
    
    n = [1 10 50 100 500];
    e = 0.4:-0.1:0.1;
    sigman = 1;
    
    % n values on different rows, e values on different columns
    nn = numel(n);
    ne = numel(e);
    sc = 4*sqrt(2);
    %sc = 2;
    p_errs = 1-exp(sc*repmat(norminv(e(:)'),nn,1).*sigman./repmat(sqrt(n(:)),1,ne));
    
    print_perr_table(p_errs, n, e);
end

if 0
    %% GAUSSIAN TARGET ASSUMPTION TO ACCOUNT DIFFERENT \MU_T VALUES --- TABLE
    %% Note: This analysis is not in the paper and is kind of obsolete. 
    
    n = [1 10 50 100 500];
    e = 0.4:-0.1:0.1;
    sigman = 1;
    k = 10; % dim of target
    s = sqrt(2.4^2/k); % proposal scaling
    
    % n values on different rows, e values on different columns
    nn = numel(n);
    ne = numel(e);
    
    m = -0.5*k*s^2;
    sigma2 = 0.5*k*s^2*(s^2+2);
    sigma = sqrt(sigma2);
    sc = 2*sqrt(2);
    %sc = 2;
    sqrtw = -sc*repmat(norminv(e(:)'),nn,1).*sigman./repmat(sqrt(n(:)),1,ne);
    p_errs = ...
        exp(m+0.5*sigma2)*(exp(sqrtw).*normcdf_fast(-(m+sqrtw+sigma2)./sigma) ...
        - exp(-sqrtw).*normcdf_fast(-(m-sqrtw+sigma2)./sigma)) ...
        + normcdf_fast((m+sqrtw)./sigma) ...
        - normcdf_fast((m-sqrtw)./sigma); % this is approximate; mu_t approx. with Gaussian
    
    print_perr_table(p_errs, n, e);
end

if 1
    %% BOTH CONDITIONAL ERRORS --- This plots **the figure to the paper** 
    
    save_fig = 1;
    %root = '../results/bounds/'; % for v1 of the paper
    root = '../results/bounds_v2/'; % for v2  of the paper
    
    n = 1:1:200; % Note: Renamed to b in v2 to surely distinguish it from n as for 'noise'
    %n = 1:1:100;
    e = [0.1:0.1:0.4];
    sigman = 1;
    k = 5; % dim of target
    s = sqrt(2.4^2/k); % proposal scaling
    
    % n values on different rows, e values on different columns
    nn = numel(n);
    ne = numel(e);
    
    %% worst case error
    sc_w = 4*sqrt(2);
    %sc_w = 2;
    p_errs_w = 1-exp(sc_w*repmat(norminv(e(:)'),nn,1).*sigman./repmat(sqrt(n(:)),1,ne));
    
    %% error where we average over typical mu_t values
    if 0
        % this is approximate; mu_t approx. with Gaussian:
        m = -0.5*k*s^2;
        sigma2 = 0.5*k*s^2*(s^2+2);
        sigma = sqrt(sigma2);
        sc = 2*sqrt(2);
        %sc = 2;
        sqrtw = -sc*repmat(norminv(e(:)'),nn,1).*sigman./repmat(sqrt(n(:)),1,ne);
        p_errs = ...
            exp(m+0.5*sigma2)*(exp(sqrtw).*normcdf_fast(-(m+sqrtw+sigma2)./sigma) ...
            - exp(-sqrtw).*normcdf_fast(-(m-sqrtw+sigma2)./sigma)) ...
            + normcdf_fast((m+sqrtw)./sigma) ...
            - normcdf_fast((m-sqrtw)./sigma);
    else
        p_errs = cond_bound_simul(k,e,ne,n,nn,s,sigman); % exact, computed using simulation
    end
    
%     %% error - simpler approx. upper bound (testing)
%     extra_method = 0;
%     if extra_method
%         p_errs_ub = cond_bound_simul(k,e,ne,n,nn,s,sigman); % exact, computed using simulation
%         if 0
%         p_errs_ub = ...
%             exp(m+0.5*sigma2)*(exp(sqrtw).*normcdf_fast(-(m+sqrtw+sigma2)./sigma) ...
%             - exp(-sqrtw).*normcdf_fast(-(m-sqrtw+sigma2)./sigma)) ...
%             + sqrt(2)*sqrtw./(sigma.*sqrt(pi));
%         end
%     end

    figure(1);
    set(gcf,'Position',[50 1000 400 300]);
    lw = 1.3;
    hold on;
    le = []; h = [];
    for i = 1:ne
        hi = plot(n,p_errs_w(:,i)','Color',col_get(i),'LineStyle','--','Linewidth',lw);
        h = [h,hi];
        le = [le,{['$\varepsilon=',num2str(e(i)),'$']}];
    end
    for i = 1:ne
        hi = plot(n,p_errs(:,i)','Color',col_get(i),'LineStyle','-','Linewidth',lw);
        h = [h,hi];
        le = [le,{['$\varepsilon=',num2str(e(i)),'$']}];
    end
%     if extra_method
%         for i = 1:ne
%             hi = plot(n,p_errs_ub(:,i)','Color',col_get(i),'LineStyle','-.','Linewidth',lw);
%             h = [h,hi];
%             le = [le,{['$\varepsilon=',num2str(e(i)),'$']}];
%         end
%     end
    hold off;
    box on; grid on;
    ylim([0,0.8]);
    legend(h,le,'interpreter','latex','NumColumns',2);
    stry = '$\bf P \rm (\textnormal{Cond. error} \geq \varepsilon)$';
    %stry = 'Prob. cond. error larger than \epsilon';
    ylabel(stry,'interpreter','latex');
    %xlabel('n');
    xlabel('b');
    title('(a)');
    
    %% save figure
    if save_fig
        fn = [root,'/fig_cond_error'];
        my_export_fig(fn,'-transparent','-pdf');
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function p = print_perr_table(p, m, e)
% plot numbers in a table format
p = [NaN,e(:)';m(:),p]
disp(' ');
end


function p = cond_bound_simul(k,e,ne,n,nn,s,sigman)
% Computes an exact bound but uses simulation instead of the Gaussian
% approximation.

m = 100000; % samples from the target
mu_t = rand_mut(m,k,s);
mu_t = repmat(mu_t(:)',nn,1);
p = NaN(nn,ne); % result as nn x ne matrix
for i = 1:ne
    lambda = -2*sqrt(2)*norminv(e(i)).*sigman./repmat(sqrt(n(:)),1,m);
    p1 = max(1-exp(mu_t-lambda),0) + min(exp(mu_t+lambda)-1,0);
    p(:,i) = mean(p1,2);
end
end


