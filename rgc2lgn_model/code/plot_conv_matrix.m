function [fig, ax, convMat_norm] = plot_conv_matrix(model_w_norm, model_w_norm_mean, weight_threshold, model_clu_idx, axPars, axParsAux)
    %% PLOT CONVERGENCE MATRIX   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % z.B. Type A kommt insgesamt 3x vor, Type B kommt 2x vor, beide zusammen
    % kommen 2x vor, dann h?tten wir 2/((3+2)/2) -> 80%
    
    nRGC = size(model_w_norm_mean, 2);

    fig = figure;

    convMat = zeros(size(model_w_norm_mean,2));
    for ic = 1:size(model_w_norm_mean,1)

        var   = model_w_norm_mean(ic,:);
        nums  = find(var > weight_threshold);    
        if length(nums) < 2 %|| length(nums) > 2
            continue
        end

        tmp = nchoosek(1:length(nums),2);
        pairs = nums(tmp);

        for ip = 1:size(pairs,1)
            convMat(pairs(ip,1),pairs(ip,2)) = convMat(pairs(ip,1),pairs(ip,2)) + 1;
            convMat(pairs(ip,2),pairs(ip,1)) = convMat(pairs(ip,2),pairs(ip,1)) + 1;
        end
    end

    % Laura
    convMat_norm = zeros(size(model_w_norm,2));
    for i = 1:size(model_w_norm,2)
        for j = 1:size(model_w_norm,2)
            a = convMat(i,j);
            b = sum(model_w_norm(:,i) > 0.001);
            c = sum(model_w_norm(:,j) > 0.001);        
            res = a/((b+c)/2);
            convMat_norm(i,j) = res;
        end
    end

    % % Phil
    % convMat_norm = zeros(size(norm_w,2));
    % for i = 1:size(norm_w,2)
    %     for j = 1:size(norm_w,2)
    %         
    %         a = convMat(i,j);
    %         var = sum(norm_w(:,j) > 0.001);
    %         res = a/var;
    %         
    %         if i==j
    %             convMat_norm(i,j) = 0;
    %         elseif i<j                  
    %             convMat_norm(i,j) = res;
    %         elseif i>j            
    %             convMat_norm(i,j) = res;
    %         end
    %     end
    % end

    convMat_norm(isnan(convMat_norm)) = 0;
    convMat_norm=round(convMat_norm*100,1);
    imagesc(convMat_norm)
    
    % Adjust plot
    ax = gca();    
    set(ax, axPars)%, 'Position', [26 h-40 12 15]);   

    set(ax,'XLim',[0.5 nRGC+0.5],'YLim',[0.5 nRGC+0.5]);
    set(ax,'FontSize',axParsAux.AxisFontSize,'LineWidth',1.5,'TickDir','out')
    set(ax,'XTick',1:1:nRGC,'XTicklabel',model_clu_idx,'XTickLabelRotation',90)
    set(ax,'YTick',flipud(1:1:nRGC),'YTicklabel',model_clu_idx,'XTickLabelRotation',90)
    title(ax,'RGC cluster co-occurrence matrix')

    daspect(ax,[1 1 1]);
    cl = colorbar(ax);
    set(cl,'TickDir','out','XTick',[0 (round(max(max(convMat_norm))*10))*10]);
    box off;
end