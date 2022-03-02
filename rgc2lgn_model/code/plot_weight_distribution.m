function [fig, ax] = plot_weight_distribution(model_w_norm, model_w_norm_mean,...
    model_corr_mean, weight_threshold, rgc, model_clu_idx, clu_vs_group, mean_vs_percent, axPars, axParsAux)
    % Plots bar chart showing weight distribution over RGCs for (1) clusters or types and (2)
    % as mean or percent.
    
    fprintf('Plotting %s weights per %s (weight threshold = %.3f).\n', mean_vs_percent, clu_vs_group, weight_threshold);
    
    if strcmp(clu_vs_group, 'clu') % Plot clusters
        if strcmp(mean_vs_percent, 'mean') % Plot mean
            %% Plot mean weights per cluster

            % Set weights = NaN if they are between 0 and abs(threshold),
            % i.e. in the range [-threshold, treshold]
            % i.e. if they are positive and below pos threshold, or negative
            % and larger than the negative treshold
            vals = nan(size(model_w_norm_mean));
%             vals(model_w_norm_mean > weight_threshold) = model_w_norm_mean(model_w_norm_mean > weight_threshold);
            vals(model_w_norm_mean > weight_threshold | model_w_norm_mean < -weight_threshold) = ...
                model_w_norm_mean(model_w_norm_mean > weight_threshold | model_w_norm_mean < -weight_threshold);
            vals = nanmean(vals);

            % extend matrix to size of clus
            maxRGC = max(rgc.cluIdx);
            vals_ext = zeros(1,maxRGC);
            vals_ext(model_clu_idx) = vals;
            color_ind = round(linspace(1,size(axParsAux.cMap_igor,1),maxRGC));

            % Plot bar chart
            fig = figure('Position', [1000, 1146, 686, 192]); 
            hold on;
            % vals_ext = vals;
            for i = 1:maxRGC
                tmp = zeros(size(vals_ext));
                tmp(i) = vals_ext(i);

                bh = bar(tmp);
                set(bh, 'FaceColor', axParsAux.cMap_igor(color_ind(i),:),'LineWidth',1.5);
            end
            
            % Adjust plot
            ax = gca();            
            set(ax, axPars, 'Position', [2.568, 1.63, 19.551, 3.944]);

            set(ax,'FontSize',axParsAux.AxisFontSize,'TickDir','out','Xlim',[0 maxRGC+0.5])
            set(ax,'XTick',1:1:maxRGC,'XTicklabel',rgc.names_clu,'XTickLabelRotation',90)
            
%             ax.XTickLabel{1} = ['color/{gray}', ax.XTickLabel{1}];
            
            set(ax,'YAxisLocation','left')
            xlabel(ax,'RGC types','FontSize',axParsAux.LabelFontSize);
            ylabel(ax,'Mean weights','FontSize',axParsAux.LabelFontSize)

            title(ax,['Weight threshold = ' num2str(weight_threshold)]);
        elseif strcmp(mean_vs_percent, 'percent') % Plot percent
            %% Plot percent per cluster

            fig = figure; hold on;
            
            % Get # cells w mean weight (across cross-validation repeats) above threshold
            % or below negative threshold
            vals = sum(model_w_norm_mean > weight_threshold | model_w_norm_mean < -weight_threshold);
            % Convert to percentage
            vals = round(vals/size(model_w_norm_mean,1)*100, 1);

            % extend matrix to size of clus
            maxRGC = max(rgc.cluIdx);
            vals_ext = zeros(1,maxRGC);
            vals_ext(model_clu_idx) = vals;
            color_ind = round(linspace(1,size(axParsAux.cMap_igor,1),maxRGC));

            for i = 1:maxRGC
                tmp = zeros(size(vals_ext));
                tmp(i) = vals_ext(i);

                bh = bar(tmp);
                set(bh, 'FaceColor', axParsAux.cMap_igor(color_ind(i),:),'LineWidth',1.5);
            end

            % Adjust plot
            ax = gca();
            set(ax, axPars);
            set(ax,'FontSize',axParsAux.AxisFontSize,'TickDir','out','Xlim',[0 maxRGC+0.5])
            set(ax,'XTick',1:1:maxRGC,'XTicklabel',rgc.names_clu,'XTickLabelRotation',90)
            set(ax,'YAxisLocation','left')
            xlabel(ax,'RGC types','FontSize',axParsAux.LabelFontSize);
            ylabel(ax,'% of cells','FontSize',axParsAux.LabelFontSize)
            t = title(ax,['Weight threshold = ' num2str(weight_threshold)]);
            t.FontSize = 14;
        end
    elseif strcmp(clu_vs_group, 'group') % Plot groups
        if strcmp(mean_vs_percent, 'mean') % Plot mean
            %% Plot group mean weights          

            % Set weights = NaN if they are between 0 and abs(threshold),
            % i.e. in the range [-threshold, treshold]
            % i.e. if they are positive and below pos threshold, or negative
            % and larger than the negative treshold
            vals = nan(size(model_w_norm_mean));
%             vals(model_w_norm_mean > weight_threshold) = model_w_norm_mean(model_w_norm_mean > weight_threshold);
            vals(model_w_norm_mean > weight_threshold | model_w_norm_mean < -weight_threshold) = ...
                model_w_norm_mean(model_w_norm_mean > weight_threshold | model_w_norm_mean < -weight_threshold);
            vals = nanmean(vals);

            % Convert from clusters to groups
            % TODO: improve code: get clu-to-group correspondence list
            vals_group = zeros(1,32);
            vals_group(1)  = mean(vals(model_clu_idx == 1));
            vals_group(2)  = mean(vals(model_clu_idx == 2));
            vals_group(3)  = mean(vals(model_clu_idx == 3));
            vals_group(4)  = mean([vals(model_clu_idx == 4) vals(model_clu_idx == 5)]);
            vals_group(5)  = mean([vals(model_clu_idx == 6) 0 vals(model_clu_idx == 8)]);
            vals_group(6)  = mean(0);
            vals_group(7)  = mean(vals(model_clu_idx == 10));
            vals_group(8)  = mean([vals(model_clu_idx == 11) vals(model_clu_idx == 12)]);
            vals_group(9)  = mean(vals(model_clu_idx == 13));
            vals_group(10) = mean(0);
            vals_group(11) = mean([vals(model_clu_idx == 15) vals(model_clu_idx == 16)]);
            vals_group(12) = mean([vals(model_clu_idx == 17) vals(model_clu_idx == 18)]);
            vals_group(13) = mean(0);
            vals_group(14) = mean(0);
            vals_group(15) = mean(vals(model_clu_idx == 21));
            vals_group(16) = mean(vals(model_clu_idx == 22));
            vals_group(17) = mean([vals(model_clu_idx == 23) vals(model_clu_idx == 24) 0]);
            vals_group(18) = mean([vals(model_clu_idx == 26) 0]);
            vals_group(19) = mean(vals(model_clu_idx == 28));
            vals_group(20) = mean(vals(model_clu_idx == 29));
            vals_group(21) = mean(vals(model_clu_idx == 30));
            vals_group(22) = mean([0 vals(model_clu_idx == 32)]);
            vals_group(23) = mean(vals(model_clu_idx == 33));
            vals_group(24) = mean(vals(model_clu_idx == 34));
            vals_group(25) = mean(0);
            vals_group(26) = mean(0);
            vals_group(27) = mean(vals(model_clu_idx == 37));
            vals_group(28) = mean([vals(model_clu_idx == 38) vals(model_clu_idx == 39)]);
            vals_group(29) = mean(0);
            vals_group(30) = mean(0);
            vals_group(31) = mean([0 vals(model_clu_idx == 43) 0 0 0]);
            vals_group(32) = mean([vals(model_clu_idx == 47) vals(model_clu_idx == 48) vals(model_clu_idx == 49)]);

            fig = figure; hold on;
            color_ind = round(linspace(1,size(axParsAux.cMap_igor,1),32));

            for i = 1:32
                tmp = zeros(size(vals_group));
                tmp(i) = vals_group(i);

                bh = bar(tmp);
                set(bh, 'FaceColor', axParsAux.cMap_igor(color_ind(i),:),'LineWidth',1.5);
            end
            
            % Adjust plot
            ax = gca();            
            set(ax, axPars);
            set(ax,'FontSize',axParsAux.AxisFontSize,'TickDir','out','Xlim',[0 32+0.5])
            set(ax,'XTick',1:1:32,'XTicklabel',rgc.names_group,'XTickLabelRotation',90)
            set(ax,'YAxisLocation','left')
            xlabel(ax,'RGC types','FontSize',axParsAux.LabelFontSize);
            ylabel(ax,'Mean weights','FontSize',axParsAux.LabelFontSize)
            t = title(ax,['Weight threshold = ' num2str(weight_threshold)]);
            t.FontSize = 14;
        elseif strcmp(mean_vs_percent, 'percent') % Plot percent
            %% Plot group percent           

            fig = figure; hold on;
            color_ind = round(linspace(1,size(axParsAux.cMap_igor,1),32));

            % Get # cells w mean weight (across cross-validation repeats) above threshold
            % or below negative threshold
            vals = sum(model_w_norm_mean > weight_threshold | model_w_norm_mean < -weight_threshold);
            % Convert to percentage
            vals = round(vals/size(model_w_norm_mean,1)*100, 1);

            % vals = round(sum(mw > weight_threshold)*100/size(mw,1),1); % Same in one line

            vals_group = zeros(1,32);
            vals_group(1)  = mean(vals(model_clu_idx == 1));
            vals_group(2)  = mean(vals(model_clu_idx == 2));
            vals_group(3)  = mean(vals(model_clu_idx == 3));
            vals_group(4)  = mean([vals(model_clu_idx == 4) vals(model_clu_idx == 5)]);
            vals_group(5)  = mean([vals(model_clu_idx == 6) 0 vals(model_clu_idx == 8)]);
            vals_group(6)  = mean(0);
            vals_group(7)  = mean(vals(model_clu_idx == 10));
            vals_group(8)  = mean([vals(model_clu_idx == 11) vals(model_clu_idx == 12)]);
            vals_group(9)  = mean(vals(model_clu_idx == 13));
            vals_group(10) = mean(0);
            vals_group(11) = mean([vals(model_clu_idx == 15) vals(model_clu_idx == 16)]);
            vals_group(12) = mean([vals(model_clu_idx == 17) vals(model_clu_idx == 18)]);
            vals_group(13) = mean(0);
            vals_group(14) = mean(0);
            vals_group(15) = mean(vals(model_clu_idx == 21));
            vals_group(16) = mean(vals(model_clu_idx == 22));
            vals_group(17) = mean([vals(model_clu_idx == 23) vals(model_clu_idx == 24) 0]);
            vals_group(18) = mean([vals(model_clu_idx == 26) 0]);
            vals_group(19) = mean(vals(model_clu_idx == 28));
            vals_group(20) = mean(vals(model_clu_idx == 29));
            vals_group(21) = mean(vals(model_clu_idx == 30));
            vals_group(22) = mean([0 vals(model_clu_idx == 32)]);
            vals_group(23) = mean(vals(model_clu_idx == 33));
            vals_group(24) = mean(vals(model_clu_idx == 34));
            vals_group(25) = mean(0);
            vals_group(26) = mean(0);
            vals_group(27) = mean(vals(model_clu_idx == 37));
            vals_group(28) = mean([vals(model_clu_idx == 38) vals(model_clu_idx == 39)]);
            vals_group(29) = mean(0);
            vals_group(30) = mean(0);
            vals_group(31) = mean([0 vals(model_clu_idx == 43) 0 0 0]);
            vals_group(32) = mean([vals(model_clu_idx == 47) vals(model_clu_idx == 48) vals(model_clu_idx == 49)]);

            for i = 1:32
                tmp = zeros(size(vals_group));
                tmp(i) = vals_group(i);

                b9 = bar(tmp);
                set(b9, 'FaceColor', axParsAux.cMap_igor(color_ind(i),:),'LineWidth',1.5);
            end
            
            % Adjust plot
            ax = gca();
            set(ax, axPars);
            
            set(ax,'FontSize',axParsAux.AxisFontSize,'TickDir','out','Xlim',[0 32+0.5])
            set(ax,'XTick',1:1:32,'XTicklabel',rgc.names_group,'XTickLabelRotation',90)
%             set(ax,'YAxisLocation','right')
            xlabel(ax,'RGC types','FontSize',axParsAux.LabelFontSize);
            ylabel(ax,'% of cells','FontSize',axParsAux.LabelFontSize)
            t = title(ax,['Weight threshold = ' num2str(weight_threshold)]);
            t.FontSize = 14;

        end
    end
end