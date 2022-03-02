function [fig, ax, cbar] = plot_weights_per_threshold(model_w, model_w_norm, axPars)
    
    % Set pars
    [nneurons, ntypes, nrepeats] = size(model_w); % dimensions
    thresholds = 0:0.01:1;     

    % Go through the thresholds to get # of clusters as a function of threshold
    edges = 0:ntypes+1;
%     edges = edges + 0.5;
    edges = edges - 0.5;
    counts = nan(numel(edges)-1, numel(thresholds));
    for ithreshold = 1 : numel(thresholds)
        threshold = thresholds(ithreshold);
        % sum: sum over the RGC types for counting, mean over repeats
        counts(:, ithreshold) = histcounts(mean(sum(model_w_norm > threshold | model_w_norm < -threshold,2),3), edges);
    end
    counts(counts == 0) = NaN;

    % Plot figure with number of clusters as a function of threshold
    fig = figure;
%     imagesc(1:ntypes, thresholds, counts');
    im = imagesc(0:ntypes, thresholds, counts');
    
    % Adjust plot
    ax = gca();
    ylabel('Weight threshold');
    xlabel('Number of RGC types');
    axis square;
    box off;
    set(ax, axPars);
    set(ax, 'YDir', 'normal');    

%     colormap(flipud(gray));
    cmap = [1 1 1 ; flipud(bone(max(max(im.CData))))]; % add ones to represent NaNs as white
    cmap = [cmap(1:2, :); cmap(0.1*max(max(im.CData)):end,:)]; % clip cmap to show values close to white as more grey
    colormap(cmap);

    cbar = colorbar;
    title(cbar, '# dLGN cells')
    cbar.FontSize = axPars.FontSize;
    cbar.TickDirection = 'out';
    
    ylim = max(max(mean(model_w_norm,3)));
    ax.YLim = [ax.YLim(1) ylim];
    
%     view([90 -90]); % flip x-y axis view
end
