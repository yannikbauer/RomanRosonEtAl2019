%% Plot example cell reconstructions from NNMF components
function plot_nnmf_example_cells(cell_idx, psth, A, X, Y, Y_norm, OUTPERM, save_fig, data_type)
    % Plots example cell responses and their NNMF reconstructions.
    % inputs: cell_idx, A, X, Y, Y_norm
    % outputs: 
    % TODO: make documentation, plot whole figure not multiple subfigures,
    % outsource save command

    for ic = 1:length(cell_idx)
        
        fprintf('plot_nnmf_example_cells: plotting ... \n')

        % Plot NNMF cell reconstruction
        figure;
        %set(gcf,'Position',[245 572 1528 275])
        subplot(2,1,1)
        hold on
        plot(psth.ts,X(:,cell_idx(ic)), 'color', 'k')
        plot(psth.ts,Y(cell_idx(ic),:)*A)
    %     plot(psth.ts,Y(cells(ic),:)*A, 'color', [0, 0.4470, 0.7410])    
        set(gca,'FontSize',11,'TickDir','out','XLim', [0 32], 'xtick',[], 'Color', 'none')
    %     xlabel('time (s)','FontSize',13);
        ylabel('Firing rate','FontSize',13)
        title(['cell # ' num2str(cell_idx(ic))]);
        legend('dLGN cell','nnmf reconstruction')
        box off;    
        % remove all tick labels except first and last - find more efficient way
        tmp = get(gca,'yticklabel');
        for i = 1:length(tmp)
            if (i > 1 & i < length(tmp)), tmp(i) = {' '}; end
        end
        set(gca,'yticklabel', tmp);

        % Plot reconstruction component weights as bar plot
        subplot(2,5,6:9)
        dist = Y_norm(cell_idx(ic),:);
        bar(dist(fliplr(OUTPERM)), 'facecolor', [0, 0.4470, 0.7410]);
    %     set(gca,'XTick',1:1:17,'XTickLabel',{fliplr(OUTPERM)})
        set(gca,'FontSize',11,'TickDir','out', 'Color', 'none')%,'YLim',[0 1])
        xlabel('Features','FontSize',13);
        ylabel('Weight','FontSize',13)
        box off;
        % remove all tick labels except first and last - find more efficient way
        tmp = get(gca,'yticklabel');
        for i = 1:length(tmp)
            if (i > 1 & i < length(tmp)), tmp(i) = {' '}; end
        end
        set(gca,'yticklabel', tmp);

        % Plot reconstruction component weights as sorted line plot
        subplot(2,5,10)
        plot(fliplr(sort(Y_norm(cell_idx(ic),:))), 'k')
        set(gca, 'XLim',[0 size(A,1)], 'TickDir','out', 'Color', 'none')    
        box off;
        % remove all tick labels except first and last - find more efficient way
        tmp = get(gca,'yticklabel');
        for i = 1:length(tmp)
            if (i > 1 & i < length(tmp)), tmp(i) = {' '}; end
        end
        set(gca,'yticklabel', tmp);

        set(gcf, 'Position', [72, 218, 1277, 468])

        % save fig
%         save_fig = false;
        if save_fig
            % save figure
            fname = ['nnmf_', data_type, '_cell_', num2str(cell_idx(ic))];
            ffname = fullfile('../results', 'cell_reconstructions', fname);
            fprintf('Saving figure to %s ...\n', ffname);            
            hgexport(gcf, ffname, hgexport('factorystyle'),'Format','eps');
            saveas(gcf, ffname, 'png');
        %     close(gcf)
        end
    end
end