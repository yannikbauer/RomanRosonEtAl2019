%% Correlate NNMF features with dLGN-p RGC type mean responses
function [max_corrs, fig] = corr_nnmf_rgc(OUTPERM, psth, A)
    
    % Correlates NNMF and RGCs after getting dLGN-p RGC type response traces 
    % and convolving LGN components with ogb1 kernel
    % inputs: OUTPERM, psth, A

    % Get the RGC cluster means
    load('data/o1rgc.mat')
    load('data/g6rgc.mat')
    load('data/kerns.mat')

    % Order components by dendrogram
    order = fliplr(OUTPERM);
    cluIdx = unique(g6rgc.cluIdx);
    nClust = length(cluIdx);

    % Get RGC distribution
    rgc_dist = zeros(nClust,1);
    for iu=1:nClust
        rgc_dist(iu) = nnz(g6rgc.cluIdx == cluIdx(iu));
    end

    % Get cluster means / get it from gcamp later
    rgc_mean = zeros(nClust,length(o1rgc.chirp_time));
    for iclu = 1:nClust
        rgc_mean(iclu,:) = mean(o1rgc.chirp(:,o1rgc.cluIdx==cluIdx(iclu)),2);
    end

    % Downsample components
    comp_mean = zeros(size(A,1),size(rgc_mean,2));
    for icomp = 1:size(A,1)
        trace = A(icomp,:);
        comp_mean(icomp,:) = interp1(psth.ts,trace,o1rgc.chirp_time','linear');
    end

    % Convolve components
    comp_conv = zeros(size(comp_mean));
    for icomp = 1:size(A,1)
        trace = comp_mean(icomp,:);
        tmp = conv(trace, o1Kern);
        comp_conv(icomp,:) = tmp(1:size(comp_mean,2));
    end

    % correct for convolution artifact % ?
    comp_conv(3,1) = 0.1598;
    comp_conv(3,2) = 0.1678;
    comp_conv(3,3) = 0.1725;
    comp_conv(3,4) = 0.1735;
    comp_conv(16,1) = 0.4831;
    comp_conv(16,2) = 0.4745;
    comp_conv(16,3) = 0.4678;
    comp_conv(16,4) = 0.4591;
    comp_conv(16,5) = 0.4516;

    % Normalize components
    comp_norm = zeros(size(comp_conv));
    for icomp = 1:size(A,1)
        trace = comp_conv(icomp,:);

%         % mean normalization
%         comp_norm(icomp,:) = trace - mean(trace(1:8));
%         comp_norm(icomp,:) = comp_norm(icomp,:) / max(abs(comp_norm(icomp,:)));
        
        % range normalization
        comp_norm(icomp,:) = trace - min(trace);
        comp_norm(icomp,:) = comp_norm(icomp,:) / max(comp_norm(icomp,:));
    end

    % Normalize RGC means
    rgc_mean_norm = zeros(size(rgc_mean));
    for irgc = 1:size(rgc_mean,1)
        trace = rgc_mean(irgc,:);

%         % mean normalization
%         rgc_mean_norm(irgc,:) = trace - mean(trace(1:8));
%         rgc_mean_norm(irgc,:) = rgc_mean_norm(irgc,:) / max(abs(rgc_mean_norm(irgc,:)));
            
        % range normalization
        rgc_mean_norm(irgc,:) = trace - min(trace);
        rgc_mean_norm(irgc,:) = rgc_mean_norm(irgc,:) / max(rgc_mean_norm(irgc,:));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Correlate components with RGC Cluster means
    corrMat = zeros(size(A,1),size(rgc_mean_norm,1));
    for i = 1:size(A,1)
        for j = 1:size(rgc_mean_norm,1)
            corrMat(i,j) = corr(comp_norm(i,:)',rgc_mean_norm(j,:)');
        end
    end

    %  PLOT STUFF
    h = 45;
    [figPars, axPars] = aux_setPlotPars();
    fh = figure(figPars, 'Position', [32 14 21 h]);


    % Plot Correlation matrix
    f1 = figure;

%     corrMat(corrMat < 0) = 0; % ?? How is this legitimized?
    imagesc(corrMat(order,:));
    figHandles = findall(f1, 'Type', 'axes');
    newT1 = copyobj(figHandles(1), fh);
    set(newT1, axPars, 'Position', [3 h-9 16 8]);
    close(f1)

    set(gca,'TickDir','out','YTick',1:2:size(corrMat,1),'XTick',1:1:length(cluIdx),'FontSize',10)
    set(gca,'XTickLabel',cluIdx,'XTickLabelRotation',90)
    set(gca, 'CLim', [-1 1]);
    ylabel(gca,'NNMF components','FontSize',14)
    xlabel(gca,'RGC clusters','FontSize',14)
    t = title(gca,'Correlation of NNFM components & RGC clusters');
    t.FontSize = 14;

    % daspect([1 1 1])
    
    colorbar(gca)

    % Plot best corr
    f2 = figure; hold on

    [corrs,ind_rgc] = max(corrMat,[],2);
    labels = cell(1,length(ind_rgc));
    for i = 1:length(ind_rgc)
        labels(i) = {num2str(cluIdx(ind_rgc(i)))};
    end

    b1 = bar(corrs(order), 'FaceColor', [0.75,0.75,0.75], 'EdgeColor', [0.5,0.5,0.5]);
%     set(b1, 'FaceColor', [0.00 0.44 0.74], 'LineWidth',1.0);
    set(gca,'XTick',1:1:size(corrMat, 1))
    xt = get(gca, 'XTick');
    text(xt, corrs(order), labels(order), 'HorizontalAlignment','center', 'VerticalAlignment','bottom')

    figHandles = findall(f2, 'Type', 'axes');
    newT2 = copyobj(figHandles(1), fh);
    set(newT2, axPars, 'Position', [3 h-14 11.7 2.5]);
    close(f2)

    set(newT2,'Xlim',[0 size(corrMat, 1)+1],'YLim',[0 1])
    set(newT2,'TickDir','out','XTick',1:2:size(corrMat, 1),'YTick',[0 0.5 1],'FontSize',12)
    ylabel(newT2,'Correlation','FontSize',14)
    xlabel(newT2,'NNMF components','FontSize',14)

    t = title(newT2,'Best NNMF-RGC Cluster Correlation');
    t.FontSize = 14;

    
    % Plot max correlation histogram
%     figure; hist(vals)
%     figure; histogram(vals, 15)
    fig = figure; 
    hist = histfit(corrs, 20, 'kernel');
%     xlabel('Correlation')
    ylabel('Count')
%     title('NNMF-RGC cluster correlation: distribution of max correlations')
    set(gca, 'view', [90 -90]);
    set(hist(1), 'FaceColor', [0.75,0.75,0.75], 'EdgeColor', [0.75,0.75,0.75]);
    set(hist(2), 'Color', 'k');
    
    figHandles = findall(fig, 'Type', 'axes');
    newT = copyobj(figHandles(1), fh);
    set(newT, axPars, 'Position', [15.3373   31.0162    2.8928    2.5047]);
    set(newT,'Xlim',[0 1])%,'YLim',[0 1])
    set(newT,'TickDir','out','FontSize',12)
    close(fig)
    
    % Plot example components
%     [corrs, ind_rgc] = max(corrMat,[],2);
    [~,ind_clu] = sort(corrs,'descend');
    
    num_exs = 6; % set desired number of example components to plot
    if num_exs > length(ind_clu) % check that desired n does not exceed n_components, else decrease it
        num_exs = length(ind_clu);
    end
    for i = 1:1:num_exs
        f3 = figure; hold on
        plot(o1rgc.chirp_time,comp_norm(ind_clu(i),:), 'LineWidth',1.5, 'Color', 'k');
        ph = plot(o1rgc.chirp_time,rgc_mean_norm(ind_rgc(ind_clu(i)),:), 'LineWidth',1.5, 'Color', [0    0.4470    0.7410]);

        minval = min([rgc_mean_norm(ind_rgc(ind_clu(i)),:) comp_norm(ind_clu(i),:)]);
        maxval = max([rgc_mean_norm(ind_rgc(ind_clu(i)),:) comp_norm(ind_clu(i),:)]);

        figHandles = findall(f3, 'Type', 'axes');
        newT3 = copyobj(figHandles(1), fh);
        set(newT3, axPars, 'Position', [3 h-16-i*4 14.2 2.5]);
        close(f3)

        set(newT3,'TickDir','out','YLim',[minval maxval])
        set(newT3,'XLim', [0 32],'XTick',0:5:32,'FontSize',12)
%             xlabel(newT3,'Time (s)','FontSize',14)

        % Adjust subplot descriptions differently for 1st and last plot
        if i == 1
            lg = legend(newT3, sprintf('NNMF comp. # %d',find(ind_clu(i) == order)),...
                sprintf('RGC cluster # %d (corr = %.2f)', cluIdx(ind_rgc(ind_clu(i))), corrs(ind_clu(i))),...
                'Orientation', 'horizontal');
            lg.Position = lg.Position + [0 0.026 0 0];

            t = title(newT3,sprintf('Example NNMFs and & best-correlated RGC clusters')); %corr = %.2f',corrs(ind_clu(i)))            
            t.Position = t.Position + [0 0.3 0];
            newT3.YLabel.String = 'A.U.';
        elseif i ~= 1
            [lg, ~, plots] = legend(newT3, sprintf('%d',find(ind_clu(i) == order)),...
                sprintf('%d (%.2f)', cluIdx(ind_rgc(ind_clu(i))), corrs(ind_clu(i))),...
                'Orientation', 'horizontal');
            shift_lg = [0 0.0255 0 0];
            
            % Change legend text color
            for ilg = 1:length(lg.String)
                lg.String{ilg} = ['\color[rgb]{' num2str(plots(ilg).Color) '} ' lg.String{ilg}];
            end
            
            lg.LineWidth = 0.0001;
            lg.Position = lg.Position + shift_lg;
        end

        if (i ~= num_exs) % remove xticklabels for all but last plot
            set(newT3, 'XTickLabel', []);            
        end
        if (i ~= 1) % remove yticklabels for all but first plot
            set(newT3, 'YTickLabel', []);            
        end
        
        set(lg, 'EdgeColor', 'none')
        set(lg.BoxFace, 'ColorType', 'truecoloralpha', 'ColorData', uint8([255;255;255;0.5*255]));
        t.FontSize = 14;
    end
    
    
    xlabel(newT3,'Time (s)','FontSize',14); % set xlabel only for last plot
    fig=fh;
    
    % Save output var
    max_corrs = corrs(order);
end