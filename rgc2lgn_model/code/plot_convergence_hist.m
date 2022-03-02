function [fig, ax] = plot_convergence_hist(model_w_norm, weight_threshold, axPars, axParsAux)
    % Plots RGC>dLGN convergence = # of RGC clusters needed for dLGN cell
    % reconstruction

    % Get mean # of RGC clusters for each cell
    mnr = zeros(size(model_w_norm,1),1);
    for i = 1:size(model_w_norm,1)
        mnr(i) = mean(sum(model_w_norm(i,:,:) > weight_threshold | model_w_norm(i,:,:) < -weight_threshold,2));
    %     counts(:, ithreshold) = histcounts(mean(sum(norm_w > threshold,2),3), edges);
    end

    % Get histogram counts
    edges = min(round(mnr))-1:1:max(round(mnr))+1;
    edges = edges+0.5;
    % n_all = histcounts(round(mnr),edges);
    n_all = histcounts(mnr,edges);
    n_all = n_all/sum(n_all)*100; % get percent cells
    edges(end) = [];
    
    % Range of mean # of RGCs needed for reconstruction
    minval = min(round(mnr));
    maxval = max(round(mnr));    
    prctile_5 = round(prctile(mnr, 5));
    prctile_95 = round(prctile(mnr, 95));

    % Plot histogram
    fig = figure; hold on;
    
    bh = bar(edges,n_all,'histc');
    ax = gca();
%     set(bh, 'FaceColor', axParsAux.barcolor, 'EdgeColor', [1 1 1], 'LineWidth',1.5);
%     set(ax,'Xtick',minval:2:maxval,'XLim',[minval-1 maxval+1])

    % Plot average number marker line
    x = zeros(100,1);
    x(:,1) = median(mnr);
    plot(x,linspace(0,ceil(max(n_all)),100)','--k','LineWidth',1.5)
    
    % Adjust plot
    set(ax, axPars)
    set(bh, 'FaceColor', axParsAux.barcolor, 'EdgeColor', [1 1 1], 'LineWidth',1.5);
    set(ax,'Xtick',minval:2:maxval,'XLim',[minval-1 maxval+1])
    set(ax,'TickDir','out','FontSize',axParsAux.AxisFontSize)
    xlabel(ax,'# RGC types','FontSize',axParsAux.LabelFontSize)
    ylabel(ax,'% of cells','FontSize',axParsAux.LabelFontSize)
    title(ax, sprintf('threshold = %.3f', weight_threshold));
    t = text(0, ax.YLim(2), ...
         sprintf('median = %.2f',round(median(mnr),2)),...
         'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 12);
     
%     lg = legend(ax,sprintf('median = %.2f',round(median(mnr),2)));
%     lg.Box = 'off';
%     lg.Location = 'northwest';
%     lg.FontSize = 10;
    
    % Display results
    fprintf('Population convergence (threshold = %.3f): \n', weight_threshold)
    fprintf('\tmean = %.2f,  median = %.2f, min = %i, max = %i, 5th %%ile = %i, 95th %%ile = %i\n',...
        mean(mnr), median(mnr), minval, maxval, prctile_5, prctile_95);
end