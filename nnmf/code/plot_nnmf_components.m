function [fig, axes, OUTPERM] = plot_nnmf_components(A, psth, Y_norm)
    % Plots NNMF components, arranged in a dendrogram, and a bar plot showing
    % the % of dLGN cells using a particular component
    % TODO: recode, so that analysis is separate from plotting, e.g. by
    % calling [outperm] = get_nnmf_comps, and plot_nnmf_comps
    do_one_plot = 1;
    if do_one_plot
        disp('NNMF plot: dendrogram, components and weights');
        k_best = size(A,1);

        fig=figure;

        % Plot dendrogram
        subplot(1,20,1:3);
        d = pdist(A, 'euclidean');
        tree = linkage(d,'ward');
        leafOrder = optimalleaforder(tree,d, 'CRITERIA', 'adjacent');
        [~,~,OUTPERM] = dendrogram(tree, 0, 'Reorder', leafOrder, 'Orientation', 'left', ... 
                                    'ColorThreshold','default');
        title('Component #');
        set(gca, 'Color', 'none');
        set(gca,'XColor','none');
        set(gca,'YTickLabel', {fliplr([1:size(A,1)])});
        axes(1) = gca();

        % Plot components
        subplot(1,20,5:19);
        PlotDisplaced(psth.ts, A(fliplr(OUTPERM),:), 10, 'std', 'k');
        set(gca,'TickDir','out','LineWidth',1);
        title('chirp components');
    %     axis off
    %     axes('Color','none','XColor','none');
    %     set(gca(),'Visible','off');
        set(gca,'color', [0.975,0.975,0.975]); % axis background color
        % Plot stimulus trigger lines
        [~, ~, stim_cumTs] = plotChirpStim();    
        nLines = length(stim_cumTs);
        vals = linspace(min(get(gca,'ylim')),max(get(gca,'ylim')),100);
        ycor = zeros(100,nLines);
        ycor(:,1:nLines)=repmat(vals',1,nLines);
        xcor = ycor;
        xcor(1:100,:)=repmat(stim_cumTs,100,1);
        plot(xcor,ycor,':', 'color', [0.7 0.7 0.7],'LineWidth',1.5);
        axes(2) = gca();

        % Plot weights
        nCompSum = sum(Y_norm,1);
        subplot(1,20,20);
        b(1) =  barh(nCompSum(OUTPERM)/max(nCompSum), ...
            'FaceColor', [0.75,0.75,0.75], 'EdgeColor', [0.5,0.5,0.5]);
        set(gca,'YLim',[0 k_best+1],'YTick',1:k_best,'TickDir','out','LineWidth',1);
        set(gca, 'YTickLabel', OUTPERM);
        axes(3) = gca();
        box off;
        axis off;

        num = round(max(nCompSum(OUTPERM))*100/size(Y_norm,1),1); % percentage of the maximal component
        title({'weights', 'max feat:', num2str(num), '%'});

        tightfig;
        set(gcf, 'Position', [1 1 570 754]);
    end
end
