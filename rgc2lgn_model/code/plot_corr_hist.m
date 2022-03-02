function [fig, ax] = plot_corr_hist(model_corr_mean, axPars, axParsAux)
    %[fig, ax] = plot_corr(model_corr_mean, axPars, axParsAux);
    % Plots histogram of correlations between cell responses and model
    % INPUT: 
    % model_corr_mean : mean correlation across cross-validation repeats 
    % axPars : axis parameters 
    % axParsAux : auxilliary axis pars that could not be passed to set()
    % OUPTUT:
    % fig, ax : figure and axis handles
    % TODO: implement for varying thresholds
    
    fig = figure; hold on;

    % Get histogram counts
    edges = linspace(min(model_corr_mean),1,15);
    n_all = histcounts(model_corr_mean,edges);
    n_all = n_all/sum(n_all)*100; % get percent cells
    edges(1) = [];

    % Plot hist
    bh = bar(edges,n_all,'histc');
    set(bh, 'FaceColor', axParsAux.barcolor, 'EdgeColor', [1 1 1], 'LineWidth',1.5);

    % Plot average correlation marker line
    x = zeros(100,1);
    x(:,1) = nanmedian(model_corr_mean);
    plot(x,linspace(0,ceil(max(n_all)),100)','--k','LineWidth',1.5)
    
    % Adjust plot
    ax = gca();
    set(gca, axPars);

    set(ax,'TickDir','out','XLim', [min(model_corr_mean) 1],'XTick',round(min(model_corr_mean)*10)/10:0.2:1,'FontSize',axParsAux.AxisFontSize)
    % set(newT1,'YLim',[0 ceil(max(n_all))])
    xlabel(ax,'Correlation','FontSize',axParsAux.LabelFontSize);
    ylabel(ax,'% cells','FontSize',axParsAux.LabelFontSize);
    
    t = text(ax.XLim(1)+0.05, ax.YLim(2), ...
        sprintf('median = %.2f \n(n = %d)',round(nanmedian(model_corr_mean),2),nnz(model_corr_mean)),...
         'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 12);
    
%     lg = legend(ax, sprintf('median = %.2f \n(n = %d)',round(median(model_corr_mean),2),nnz(model_corr_mean)));
%     lg.Box = 'off';
%     lg.Location = 'northwest';
%     lg.FontSize = 12;
%     title(ax, 'Correlation');
    
    % Display results
    fprintf('Population correlation: \n\tmean = %.2f,  median = %.2f\n', nanmean(model_corr_mean), nanmedian(model_corr_mean));

end