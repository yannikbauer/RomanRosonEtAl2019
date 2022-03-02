function [fig, ax] = plot_rmse_hist(model_rmse_mean, axPars, axParsAux)
    % Plots histogram of RMSE for model

    fig = figure; hold on;
    
    % Remove some extreme values from failed model convergence
    model_rmse_mean = model_rmse_mean(model_rmse_mean < 10);
    % Get histogram counts
    minval = min(model_rmse_mean);
    maxval = max(model_rmse_mean);
    edges = linspace(-0.01,maxval,25);
    n_all = histcounts(model_rmse_mean,edges);
    n_all = n_all/sum(n_all)*100; % get percent cells
    edges(1) = [];

    % Plot histogram
    bh = bar(edges,n_all,'histc');
    set(bh, 'FaceColor', axParsAux.barcolor, 'EdgeColor', [1 1 1], 'LineWidth',1.5);

    % Plot average RMSE marker line
    x = zeros(100,1);
    x(:,1) = nanmedian(model_rmse_mean);
    plot(x,linspace(0,ceil(max(n_all)),100)','--k','LineWidth',1.5)
    
    % Adjust plot
    ax = gca();
    set(ax, axPars)
    set(ax,'TickDir','out','XLim', [minval maxval], 'FontSize', axParsAux.AxisFontSize)
    xlabel(ax,'RMSE','FontSize',axParsAux.LabelFontSize)
    ylabel(ax,'% cells','FontSize',axParsAux.LabelFontSize)
    
    t = text(ax.XLim(1)+0.005, ax.YLim(2), ...
        sprintf('median = %.2f',round(nanmedian(model_rmse_mean),2)),...
         'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 12);
    
%     lg = legend(ax, sprintf('median = %.2f',round(median(model_rmse_mean),2)));
%     lg.Box = 'off';
%     lg.Location = 'northwest';
%     lg.FontSize = 12;
%     title(ax,'RMSE')

    % Display results
    fprintf('Population RMSE: \n\tmean = %.2f,  median = %.2f\n', nanmean(model_rmse_mean), nanmedian(model_rmse_mean));

end