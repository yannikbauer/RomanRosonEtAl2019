function [fig, lg] = plot_units_model(cell_idx, rgc, model_clu_idx, model_units_mean, model_units_y_hat_mean,...
    model_corr_mean, model_w_norm_mean, model_w_norm_sd, weight_threshold,...
    model_type, figPars, axPars, axParsAux, h, color_ind)
    %% Plot example cells
    % TODO: This function has not yet ben modified at all, so change
    % everything
    % This is how it should be: write function that makes single cell trace and
    % weights plot using subplot(1,2,1or2), call that in a loop, and position the output fig into
    % fig_all using copyobj()
    % NOTE: could even position all axes in one go by passing them all together
    % into copyobj([ax(1), ax(2)..., ax(8)]) or better: copyobj(ax(:)), if we
    % output ax as an array from the function
    % INPUT: cell_idx, model_units_mean, model_units_y_hat_mean,
    % model_corr_mean, model_w_norm_mean, model_w_norm_sd, weight_threshold,
    % model_type, axPars, axParsAux
    
    % Display progress
    fprintf('Plotting example cells.\n')
    
    % Make figure for all subplots
    fig = figure('Position', [1040, 812, 520, 526]);

%     h = 24;
    pos = linspace(12,22,length(cell_idx)); % subplot positions; TODO: replace w actual subplot call
    
    % Get overall highest n of bars amongst subplots to normalize their widths in summary plot
    max_n_bars = 0; % initialize var
    
    % Get overall highest weight amongst subplots to normalize their widths in summary plot
    max_bar_height = 0; % initialize var
    min_bar_height = 0; % initialize var
    
    % TODO: call function
    % inputs: model_units_mean, model_units_y_hat_mean, model_corr_mean,
    % weight_threshold, axPars, model_type
    for i = 1:length(cell_idx)        

        iunit = cell_idx(i);    
        
        % Make individual unit traces figure
        fig_trace = figure; hold on

        % Plot cell response vs model time series    
        plot(rgc.chirp_time, model_units_mean(iunit,:),'k','LineWidth',1.5)   % data
        plot(rgc.chirp_time, model_units_y_hat_mean(iunit,:),'LineWidth',1.5) % model

        % TODO: outsource this later
        figHandles = findall(fig_trace, 'Type', 'axes');
        copy = copyobj(figHandles(1), fig);
        set(copy, axPars, 'Position', [2 h-pos(i) 8 2]);
        close(fig_trace)

        % Adjust plot
        minval = min([min(model_units_mean(iunit,:)) min(model_units_y_hat_mean(iunit,:))]); % min value of cell response or model 
        set(copy,'TickDir','out','Xlim',[0 32],'YLim',[minval 1],'FontSize',axPars.FontSize);
        set(copy,'YTick',[0 1], 'XTick', [0 2]);
        % Info only in first subplot
        if i == 1
            lg = legend(copy,'data',sprintf('%s, p=%.3f', model_type, model_corr_mean(iunit)));
            set(lg, 'EdgeColor', 'none')
            set(lg.BoxFace, 'ColorType', 'truecoloralpha', 'ColorData', uint8([255;255;255;0.5*255]));
            ylabel('Norm. \DeltaF/F')
        else
            th = text(copy, copy.XLim(2), copy.YLim(2), sprintf('p=%.3f', model_corr_mean(iunit)), ...
                'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', 12);
        end

        % Plot cell model weights barplot
        % Make individual unit traces figure
        fig_weights = figure; hold on
        if i == 4
            a = 1;
        end
        % Find weights above threshold
        vals = model_w_norm_mean(iunit,:);
        ind = find(vals > weight_threshold | vals < -weight_threshold);        
        
        for ic = 1:length(ind)
            tmp = zeros(1,length(ind));
            tmp_mw = zeros(1,length(ind));
            tmp_sd = zeros(1,length(ind));

            tmp(ic) = vals(ind(ic));
            tmp_mw(ic) = model_w_norm_mean(iunit,ind(ic));
            tmp_sd(ic) = model_w_norm_sd(iunit,ind(ic))/2;

            bh = bar(tmp);
            set(bh, 'FaceColor', axParsAux.cMap_igor(color_ind(model_clu_idx(ind(ic))),:),'LineWidth',0.5, 'BarWidth', 0.8);
            errorbar(tmp_mw,tmp_sd,'k','marker','none','LineStyle','none')
        end
        
        % Get overall highest n of bars amongst subplots to normalize their widths in summary plot
        max_n_bars_tmp = length(ind);
        if max_n_bars_tmp > max_n_bars; max_n_bars=max_n_bars_tmp; end
        
        % Get overall highest & lowest weight amongst subplots to normalize their widths in summary plot
        max_bar_height_tmp = max(vals);
        if max_bar_height_tmp > max_bar_height; max_bar_height = max_bar_height_tmp; end
        min_bar_height_tmp = min(vals);
        if min_bar_height_tmp < min_bar_height; min_bar_height = min_bar_height_tmp; end

        figHandles = findall(fig_weights, 'Type', 'axes');
        copy = copyobj(figHandles(1), fig);
        set(copy, axPars, 'Position', [11 h-pos(i) 3 2]);
        close(fig_weights)

        set(copy,'FontSize',axParsAux.AxisFontSize,'TickDir','out','XLim',[0.5 length(ind)+0.5]) % 'YLim',[0 1],'YTick',[0 1],
        set(copy,'XTick',1:1:length(ind),'XTicklabel',model_clu_idx(ind),'XTickLabelRotation',90)
        if i == length(cell_idx)
            xlabel(copy,'RGC clusters')%,'FontSize',axParsAux.LabelFontSize);
        end

    end
    
    % Summary plot adjustments
    figHandles = findall(fig, 'Type', 'axes');
    for i = 1:2:length(figHandles)
        bar_plot_ax = figHandles(i);
        bar_plot_child = bar_plot_ax.Children;
        for j = 2:2:length(bar_plot_child)
            % Adjust bar widths for summary plot
            bh = bar_plot_child(j);
            bh.BarWidth = bh.BarWidth * (length(bar_plot_child)/2) / max_n_bars;
            
        end
        % Adjust bar plot heights 
        if isempty(regexp(model_type, 'nonneg')) % adjust for models without non-negativity constraint
            bar_plot_ax.YLim = [round(min_bar_height,1), round(max_bar_height,1)];
            bar_plot_ax.YTick = [round(min_bar_height,1), round(max_bar_height,1)];
            bar_plot_ax.YTickLabel = [round(min_bar_height,1), round(max_bar_height,1)];
        else % set to [0,1] for models with non-neg
            bar_plot_ax.YLim = [0, 1];
            bar_plot_ax.YTick = [0, 1];
            bar_plot_ax.YTickLabel = [0, 1];
        end

    end
    


end