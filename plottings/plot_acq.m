function [] = plot_acq(th_grid, acq_grid, th_cur, th_pr, th_tr, y_tr, acq_opt)
% Plots the results of the acq optimization. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Additional settings
plot_acq_alt = 1;
plot_opt_region = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d = th_grid.dim;
% check the input struct
opt_f = acq_grid.opt_f;
opt_th = acq_grid.opt_th;
acq = []; acq_alt = []; 
if isfield(acq_grid,'acq')
    acq = acq_grid.acq;
end
if isfield(acq_grid,'acq_alt')
    acq_alt = acq_grid.acq_alt;
end

plot_acq_alt = (plot_acq_alt && ~isempty(acq_alt));
nr_subplots = 1 + plot_acq_alt;
titl = [acq_opt.method, ', t = ', num2str(length(y_tr)+1)];

%""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""  
%% plot acq. function values in the grid for closer inspection
if d == 1
    max_a = max(acq(isfinite(acq)));
    min_a = min(acq(isfinite(acq)));
    if isempty(max_a) || isempty(min_a)
        max_a = 1; min_a = 0;
    end
    if max_a <= min_a
        max_a = min_a+1; % for avoiding rare issue
    end
    
    if plot_acq_alt
        max_var = max(acq_alt);
        min_var = min(acq_alt);
        if max_var <= min_var
            max_var = min_var+1; % for avoiding rare issue
        end
        
        subplot(1,nr_subplots,2);
        hold on;
        plot(th_grid.theta,acq_alt,'-k');
        if isfield(acq_grid,'cur_acq')
            plot([th_grid.range(1), th_grid.range(2)],acq_grid.cur_acq*[1,1],'--r'); % current acq value
            max_var = max(max_var,acq_grid.cur_acq);
            min_var = min(min_var,acq_grid.cur_acq);
        end
        for i = 1:length(opt_th)
            ha=plot(opt_th(i)*[1,1],[min_var,max_var],'b-'); % plot optimum (i.e. the latest eval point)
        end
        if 1
            hc=plot(th_cur*[1,1],[min_a,max_a],'r-'); % current and proposed point
            hp=plot(th_pr*[1,1],[min_a,max_a],'r--');
        end
        hold off;
        box on;
        xlim([th_grid.range(1), th_grid.range(2)]);
        ylim([min_var,max_var]);
        ylabel('Acquisition value');
        title('Alt acq');
        set_legend(hc,hp,ha,d,0);
    end
    
    subplot(1,nr_subplots,1);
    hold on;
    plot(th_grid.theta,acq,'-k'); % plot acq-value
    if isfield(acq_grid,'cur_acq')
        plot([th_grid.range(1), th_grid.range(2)],acq_grid.cur_acq*[1,1],'--r'); % current uncertainty value
        max_a = max(max_a,acq_grid.cur_acq);
        min_a = min(min_a,acq_grid.cur_acq);
    end
    for i = 1:length(opt_th)
        ha=plot(opt_th(i)*[1,1],[min_a,max_a],'b-'); % plot optimum (i.e. the latest eval point)
    end
    if 1
        hc=plot(th_cur*[1,1],[min_a,max_a],'r-'); % current and proposed point
        hp=plot(th_pr*[1,1],[min_a,max_a],'r--');
    end
    plot_loc_rectangle(acq_grid,d,~plot_opt_region);
    hold off;
    box on;
    % data points not plotted here...
    xlim([th_grid.range(1), th_grid.range(2)]);
    ylim([min_a,max_a]);
    set(gcf,'Position',[50 50 nr_subplots*500 380]);
    ylabel('Acquisition value');
    title(titl,'Interpreter','none');
    set_legend(hc,hp,ha,d,0);

%""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""      
elseif d == 2
    
    if plot_acq_alt
        acq_alt = vec_to_grid_matrix(acq_alt, th_grid);
    end
    
    thx = th_grid.theta(1,:);
    thy = th_grid.theta(2,:);
    nr_cont_lines = 50;
    if plot_acq_alt
        subplot(1,nr_subplots,2);
        hold on;
        contour(thx, thy, acq_alt, nr_cont_lines);
        ha=plot(opt_th(:,1),opt_th(:,2),'k*','MarkerSize',8,'LineWidth',1.5); % plot optimum (i.e. the latest eval point)
        plot(th_tr(:,1),th_tr(:,2),'k.','MarkerSize',10); % plot also datapoints
        if 1
            hc=plot(th_cur(1),th_cur(2),'rx','MarkerSize',12,'LineWidth',1.5); % current and proposed point
            hp=plot(th_pr(1),th_pr(2),'r+','MarkerSize',12,'LineWidth',1.5);
        end
        colorbar;
        hold off;
        box on;
        title('Alt acq');
        set_legend(hc,hp,ha,d,0);
    end
    
    subplot(1,nr_subplots,1);
    hold on;
    contour(thx, thy, acq, nr_cont_lines);
    ha=plot(opt_th(:,1),opt_th(:,2),'k*','MarkerSize',8,'LineWidth',1.5); % plot optimum(s) i.e. the latest eval point(s)
    plot(th_tr(:,1),th_tr(:,2),'k.','MarkerSize',10); % plot also datapoints
    if 1
        hc=plot(th_cur(1),th_cur(2),'rx','MarkerSize',12,'LineWidth',1.5); % current and proposed point
        hp=plot(th_pr(1),th_pr(2),'r+','MarkerSize',12,'LineWidth',1.5);
    end
    plot_loc_rectangle(acq_grid,d,~plot_opt_region);
    hold off;
    colorbar;
    xlim([th_grid.range(1,1), th_grid.range(1,2)]);
    ylim([th_grid.range(2,1), th_grid.range(2,2)]);
    box on;
    set(gcf,'Position',[50 10 nr_subplots*450 350]);
    title(titl,'Interpreter','none');
    set_legend(hc,hp,ha,d,0);
    
%""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""    
else % dim > 2

    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    batch_size = size(opt_th,1);
    for i = 1:d
        for j = 1:batch_size
            loc = [0.05 + 0.2*(j-1), (d-i+1)/d/1.15];
            str = num2str(opt_th(j,i));
            text(loc(1),loc(2),str);
        end
    end
    box on;
    title('Acquired location(s)');
    set(gcf,'Position',[50, 50, 100 + 100*batch_size, 250]);
    
%     disp(' ');
%     if min(size(opt_th)) == 1
%         disp(['The latest acquired point = ', num2str(opt_th(:)')]);
%     else
%         disp('The latest batch of acquired points:');
%         opt_th
%     end
%     disp(' ');
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = set_legend(hc,hp,ha,d,skip)
if skip
    return;
end
if d==1
    legend([hc,hp,ha],{'current','proposed','acquired'});
elseif d==2
    %legend([hc,hp,ha],{'current','proposed','acquired'},'Location','Southeast');
    legend([hc,hp],{'current','proposed'},'Location','Southeast');
end
end


function plot_loc_rectangle(acq_grid,d,skip)
% Plots the initial points and/or the local optimisation rectangle
if skip
    return;
end
% points:
if d==1
    if ~isempty(acq_grid.inits)
        plot(acq_grid.inits,zeros(size(acq_grid.inits)),'m.'); % initial points used in optimisation
    end
elseif d==2
    if ~isempty(acq_grid.inits)
        plot(acq_grid.inits(:,1),acq_grid.inits(:,2),'md','MarkerSize',6,'LineWidth',1.4);
    end
end

% rectangle:
th_range = acq_grid.th_range_loc;
if d==1
    if ~isempty(th_range)
        plot(th_range(1),0,'m>'); plot(th_range(2),0,'m<'); % TODO: y-value needs to be adjusted
    end
elseif d==2
    if ~isempty(th_range)
        plot([th_range(1,1) th_range(1,2) th_range(1,2) th_range(1,1) th_range(1,1)],...
            [th_range(2,1) th_range(2,1) th_range(2,2) th_range(2,2) th_range(2,1)],'-m');
    end
end
end




