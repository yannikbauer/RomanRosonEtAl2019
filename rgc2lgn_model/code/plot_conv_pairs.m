function [fig, ax] = plot_conv_pairs(convMat_norm, model_clu_idx, maxRGC, axPars, axParsAux)
    %% Plot convergence pairs per cluster   
    
    nRGC = length(convMat_norm);
    color_ind = round(linspace(1,size(axParsAux.cMap_igor,1),maxRGC));

    fig = figure; hold on;
    
    vals = sum(convMat_norm);
    vals = vals/max(vals);
    for i = 1:length(vals)
        tmp = zeros(1,length(vals));
        tmp(i) = vals(i);

        bh = bar(tmp);
        set(bh, 'FaceColor', axParsAux.cMap_igor(color_ind(model_clu_idx(i)),:),'LineWidth',1.5);
    end

    % Adjust plot
    ax = gca();
    set(ax, axPars)%, 'Position', [26 h-26 10.5 2]);
    set(ax,'FontSize',axParsAux.AxisFontSize,'TickDir','out','Xlim',[0 nRGC + 0.5])
    set(ax,'XTick',[],'XTickLabel',[],'YTick',[0 1],'YTickLabel',[0 1])
    ylabel(ax,'Mean weights','FontSize',axParsAux.LabelFontSize)
end