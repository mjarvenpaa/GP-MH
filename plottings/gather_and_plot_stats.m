function [stats,h] = gather_and_plot_stats(evals,evals_invalid,errs,acc_nr,opt,fignr)
% Computes and plots some statistics on the progression of GP-MH. 

bbs = min(evals):max(evals);
stats.counts_evals = histcounts(evals(:)',[bbs bbs(end)+1]);
%stats.cum_counts_evals = cumsum(bbs.*stats.counts_evals);
stats.acc_nr = acc_nr;
stats.acc_freq = acc_nr/opt.n;
stats.sevals = basic_stats(evals);
stats.serrs = basic_stats(errs);
stats.evals = sparse(evals); % no evals at most iterations -> sparse useful
stats.evals_invalid = sparse(evals_invalid);
if opt.ret_lvl >= 2
    stats.errs = errs;
end
h = [];

if opt.viz >= 1
    h = figure(fignr);
    set(gcf,'Position',[50 1000 1600 600]);
    nx = 2; ny = 2;
    
    % errors associated with MH accept/reject decision as a function of MCMC iteration
    subplot(ny,nx,1);
    errsn = errs;
    errsn(isnan(errsn)) = 0; % out of bounds values are NaNs -> set them 0
    x = 1:length(errs);
    plot(x,errsn,'-k');
    hold on;
    plot(x(isnan(errs)),zeros(size(x(isnan(errs)))),'*b'); % NaN cases
    hold off;
    xlim([x(1),x(end)]);
    title('Remaining error');
    xlabel('iteration');
    
    % errors associated with MH accept/reject decision as histogram
    subplot(ny,nx,3);
    h1 = histogram(errsn);
    ylim([0,max(1.1*h1.Values(2),1)]); % truncate first large bin
    title('Remaining error');
    
    % number of loglik evaluations as a function of MCMC iteration 
    subplot(ny,nx,2);
    plot(1:length(evals),evals,'-k');
    xlim([1,length(errs)]);
    ylim([0,max(evals)]);
    title('Evaluations');
    xlabel('iteration');
    
    % number of loglik evaluations as histogram
    subplot(ny,nx,4);
    h2 = histogram(evals);
    ylim([0,max(1.1*h2.Values(2),1)]); % truncate first large bin
    title('Evaluations/iteration');

    % some statistics as text boxes:
    str = {'Errors:',['mean=',num2str(nanmean(errs))],['med=',num2str(nanmedian(errs))],...
        ['90%q=',num2str(quantile(errs,0.9))],['max=',num2str(nanmax(errs))]};
    annotation('textbox',[0.35 0.15 0.25 0.25],'String',str,'FitBoxToText','on');

    str = {'Evals:',['mean=',num2str(nanmean(evals))],['med=',num2str(nanmedian(evals))],...
        ['90%q=',num2str(quantile(evals,0.9))],['max=',num2str(nanmax(evals))],...
        ['nr. invalid evals=',num2str(sum(evals_invalid))]};
    annotation('textbox',[0.8 0.15 0.25 0.25],'String',str,'FitBoxToText','on');
end
end

function s = basic_stats(v)
% Some basic statistics, computed from vector 'v' and so that possible NaN's are ignored.  
s.mean = nanmean(v);
s.med = nanmedian(v);
s.q90 = quantile(v,0.9);
s.max = nanmax(v);
end

