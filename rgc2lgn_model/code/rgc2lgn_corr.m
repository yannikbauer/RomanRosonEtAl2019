% One more thought for the modeling: once we have, for each single LGN
% neuron, the combination of weights for the RGC clusters, we could try
% to see if similarly responding neurons have similar weights. Maybe we
% can try to use the NNMF results to determine "similarity" and then
% compare the weight distributions?

% There are cells, whose NNMF weights can be correlated with a lot of
% cells. In contrast, there are cells, whose NNMF weights do not correlate
% well with others. Use Gaussian Mixture Models to determine the groups

%% Load data

load('data/nnmf_best.mat')
load('workspace/lin_range_o1g6.mat')


%% SET PARAMETERS

saveFigs = 1;                          % 1 = save figures; 0 figs won't be saved
fname    = 'figs/lin_range_corr';  % figure folder name
varname  = 'lin_range_corr';

nnmf_thr = 0.3;
lgn_norm = psth.psth;

%% COMPUTE AND PLOT 

nCells = 2;

% specify colors
c1 = [0.85    0.32    0.09];
c2 = [0.35    0.62    0.29];
c3 = [0.00    0.44    0.74];

% compute correlations between NNMF weights
[R_nmf,P_nmf] = corrcoef(Y,'rows','pairwise');

h = 37;
[figPars, axPars] = aux_setPlotPars();

for iunit = 1:1:nLGN
    
    %%%%%%%%
    % NNMF %
    %%%%%%%%
    
    % most similarly responding cells based on NNMF weights
    [corr_nnmf,ind] = sort(R_nmf(iunit,:),'descend');
    [ind_nnmf] = ind(1:nCells+1);
    var_corr_nnmf(iunit,:) = corr_nnmf(2:3);
    
    % compute correlations between convolved LGN cells
    corr_lgn = corr(lgn_norm(ind_nnmf,:)');
    var_corr_lgn(iunit,:) = corr_lgn(1,2:3);
    
    % # of nnmf components
    nComp = size(Y,1);
    
    
    %%%%%%%%%%%%%%
    % PREDICTION %
    %%%%%%%%%%%%%%
    
    % compute correlations between MODEL weights
    corr_w = corrcoef(var_w(ind_nnmf,:)');
    var_corr_w(iunit,:) = corr_w(1,2:3);
    
    % compute correlations between LGN predictions
    corr_mod = corr(var_yhat2(ind_nnmf,:)');
    var_corr_mod(iunit,:) = corr_mod(1,2:3);
 
    
    %%%%%%%%
    % PLOT %
    %%%%%%%%
    
    fh = figure(figPars, 'Position', [33 7 20 h]);
    
    % PLOT LGN WEIGHTS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    f2 = figure; hold on    
    
    dist = Y(:,ind_nnmf(1)) <= nnmf_thr;    
    b11 = bar(Y(:,ind_nnmf(1)),'FaceColor',c1);
    b12 = bar(Y(:,ind_nnmf(1)) .* dist,'FaceColor',[0.8 0.8 0.8]);    
    figHandles = findall(f2, 'Type', 'axes');
    newT2 = copyobj(figHandles(1), fh);
    set(newT2, axPars, 'Position', [2 h-5 15 2.5]);
    close(f2)
    
    set(newT2,'FontSize',10,'TickDir','out','YLim',[0 1],'YTick',[0 1])
    set(newT2,'XTick',1:1:nComp,'XTicklabel',1:1:25,'XTickLabelRotation',90)
    set(newT2,'XLim', [0 nComp+1])
    title(newT2,'NNMF Weights');
    lg = legend(newT2, 'Cell 1');
    lg.Box = 'off';
    lg.Location = 'eastoutside';
    
    %%%%%%%%%%%%%
    
    f3 = figure; hold on
    
    dist = Y(:,ind_nnmf(2)) <= nnmf_thr;    
    b21 = bar(Y(:,ind_nnmf(2)),'FaceColor',c2);
    b22 = bar(Y(:,ind_nnmf(2)) .* dist,'FaceColor',[0.8 0.8 0.8]);
    figHandles = findall(f3, 'Type', 'axes');
    newT3 = copyobj(figHandles(1), fh);
    set(newT3, axPars, 'Position', [2 h-8.5 16.5 2.5]);
    close(f3)
    
    set(newT3,'FontSize',10,'TickDir','out','YLim',[0 1],'YTick',[0 1])
    set(newT3,'XTick',1:1:nComp,'XTicklabel',1:1:25,'XTickLabelRotation',90)
    set(newT3,'XLim', [0 nComp+1])    
    lg = legend(newT3, sprintf('Cell 2, corr: %.2f',corr_nnmf(2)));
    lg.Box = 'off';
    lg.Location = 'eastoutside';
    
    %%%%%%%%%%%%%
    
    f4 = figure; hold on
    
    dist = Y(:,ind_nnmf(3)) <= nnmf_thr;    
    b11 = bar(Y(:,ind_nnmf(3)),'FaceColor',c3);
    b12 = bar(Y(:,ind_nnmf(3)) .* dist,'FaceColor',[0.8 0.8 0.8]);
    figHandles = findall(f4, 'Type', 'axes');
    newT4 = copyobj(figHandles(1), fh);
    set(newT4, axPars, 'Position', [2 h-12 16.5 2.5]);
    close(f4)
    
    set(newT4,'FontSize',10,'TickDir','out','YLim',[0 1],'YTick',[0 1])
    set(newT4,'XTick',1:1:nComp,'XTicklabel',1:1:25,'XTickLabelRotation',90)
    set(newT4,'XLim', [0 nComp+1])
    xlabel(newT4,'RGC Types','FontSize',11);
    lg = legend(newT4, sprintf('Cell 3, corr: %.2f',corr_nnmf(3)));
    lg.Box = 'off';
    lg.Location = 'eastoutside';
    
    
    % PLOT LGN TRACES  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    f1 = figure; hold on
    plot(psth.ts, lgn_norm(ind_nnmf(1),:),'Color',c1,'LineWidth',1.5)
    plot(psth.ts, lgn_norm(ind_nnmf(2),:),'Color',c2,'LineWidth',1.5)
    plot(psth.ts, lgn_norm(ind_nnmf(3),:),'Color',c3,'LineWidth',1.5)
    
    figHandles = findall(f1, 'Type', 'axes');
    newT1 = copyobj(figHandles(1), fh);
    set(newT1, axPars, 'Position', [2 h-17.5 16.7 2.5]);
    close(f1)
    
    set(newT1,'FontSize',10,'TickDir','out','YLim',[0 1],'YTick',[0 1],'XLim',[0 32])
    lg = legend(newT1,'Cell 1',sprintf('Cell 2, corr: %.2f',corr_lgn(2)),sprintf('Cell 3, corr: %.2f',corr_lgn(3)));
    lg.Location = 'eastoutside';
    lg.Box = 'off';
    title(newT1,'LGN Cells');
    
    
    
    %%%%%%%%%%%%%%%%%%%%
    % PLOT PREDICTIONS %
    %%%%%%%%%%%%%%%%%%%%
    
    % PLOT MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    f5 = figure; hold on
    plot(rgc.chirp_time, var_yhat2(ind_nnmf(1),:),'Color',c1,'LineWidth',1.5)
    plot(rgc.chirp_time, var_yhat2(ind_nnmf(2),:),'Color',c2,'LineWidth',1.5)
    plot(rgc.chirp_time, var_yhat2(ind_nnmf(3),:),'Color',c3,'LineWidth',1.5)
    
    figHandles = findall(f5, 'Type', 'axes');
    newT5 = copyobj(figHandles(1), fh);
    set(newT5, axPars, 'Position', [2 h-22.5 16.5 2.5]);
    close(f5)
    
    set(newT5,'FontSize',10,'TickDir','out','YLim',[0 1],'YTick',[0 1],'XLim',[0 32])
    
    lg = legend(newT5,'Cell 1',sprintf('Cell 2, corr: %.2f',corr_mod(2)),sprintf('Cell 3, corr: %.2f',corr_mod(3)));
    lg.Location = 'eastoutside';
    lg.Box = 'off';  
    
    title(newT5,'LGN Predictions');
    
    % PLOT MODEL WEIGHTS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    f6 = figure; hold on
    
    b6 = bar(var_w(ind_nnmf(1),:),'FaceColor',c1);
    figHandles = findall(f6, 'Type', 'axes');
    newT6 = copyobj(figHandles(1), fh);
    set(newT6, axPars, 'Position', [2 h-27.5 15 2.5]);
    close(f6)
    
    set(newT6,'FontSize',10,'TickDir','out','YLim',[0 1],'YTick',[0 1])
    set(newT6,'XTick',1:1:nRGC,'XTicklabel',cluIdx,'XTickLabelRotation',90)
    title(newT6, 'Model Weights');
    lg = legend(newT6, 'Cell 1');
    lg.Box = 'off';
    lg.Location = 'eastoutside';
    
    %%%%%%%%%%%%%
    
    f7 = figure; hold on
    
    b7 = bar(var_w(ind_nnmf(2),:),'FaceColor',c2);
    figHandles = findall(f7, 'Type', 'axes');
    newT7 = copyobj(figHandles(1), fh);
    set(newT7, axPars, 'Position', [2 h-31 16.5 2.5]);
    close(f7)
    
    set(newT7,'FontSize',10,'TickDir','out','YLim',[0 1],'YTick',[0 1])
    set(newT7,'XTick',1:1:nRGC,'XTicklabel',cluIdx,'XTickLabelRotation',90)
    lg = legend(newT7, sprintf('Cell 2, corr: %.2f',corr_w(1,2)));
    lg.Box = 'off';
    lg.Location = 'eastoutside';
    
    %%%%%%%%%%%%%
    
    f8 = figure; hold on
    
    b8 = bar(var_w(ind_nnmf(3),:),'FaceColor',c3);
    figHandles = findall(f8, 'Type', 'axes');
    newT8 = copyobj(figHandles(1), fh);
    set(newT8, axPars, 'Position', [2 h-34.5 16.5 2.5]);
    close(f8)
    
    set(newT8,'FontSize',10,'TickDir','out','YLim',[0 1],'YTick',[0 1])
    set(newT8,'XTick',1:1:nRGC,'XTicklabel',cluIdx,'XTickLabelRotation',90)
    xlabel(newT8,'RGC Types','FontSize',11);    
    lg = legend(newT8, sprintf('Cell 3, corr: %.2f',corr_w(1,3)));
    lg.Box = 'off';
    lg.Location = 'eastoutside';
    
    st = supertitle(sprintf('Cell %d; Qi: %.2f; ranksum: %.3f',iunit,lgn.qi(iunit),lgn.corr_p(iunit)),0.99);
    st.FontSize = 14;
    st.FontWeight = 'Bold';
    
    % save figs
    if saveFigs == 1
        if(exist(fname, 'dir') ~= 7)
            mkdir(fname)
        end
        titlestr = sprintf('M%d',iunit);
        hgexport(fh,fullfile(fname,titlestr),hgexport('factorystyle'),'Format','png');
        close(fh)
    end
end

if saveFigs == 1
    save(fullfile('workspace',varname));
end


return

%% PLOT SOME STUFF

% create mask with failed fits
fit_mask = ones(nLGN,1);
fit_mask(isnan(var_corr_exp2)) = 0;
fit_mask(var_corr_exp2 <= 0) = 0;

% apply mask
var_corr_mod(fit_mask == 0,:)  = [];
var_corr_lgn(fit_mask == 0,:)  = [];
var_corr_nnmf(fit_mask == 0,:) = [];
var_corr_w(fit_mask == 0,:)    = [];


% model correlations   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
edges = linspace(0,1,50);

n_c1 = histcounts(var_corr_mod(:,1),edges);
n_c2 = histcounts(var_corr_mod(:,2),edges);
n_c1 = n_c1/max(n_c1);
n_c2 = n_c2/max(n_c2);
edges(end) = [];

figure; hold on;
b11 = bar(edges,n_c1,'histc');
b12 = bar(edges,n_c2,'histc','FaceColor',c2);

set(b11, 'FaceColor', [0.8 0.8 0.8], 'LineWidth',1.5);
set(b12, 'FaceColor', [0.7 0 0],'FaceAlpha',0.4,'LineWidth',1.5);
set(gca,'TickDir','out','XLim', [0 1.1],'XTick',0:0.2:1.1,'FontSize',13)

xlabel(gca,'Correlation','FontSize',14);
ylabel(gca,'Norm counts','FontSize',14);
title(gca, 'LGN Model Correlation');
box off

lg = legend('Cell2','Cell3','Location','NorthWest');
lg.Box = 'off';


% model correlations   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
edges = linspace(0,1,50);

n_c1 = histcounts(var_corr_lgn(:,1),edges);
n_c2 = histcounts(var_corr_lgn(:,2),edges);
n_c1 = n_c1/max(n_c1);
n_c2 = n_c2/max(n_c2);
edges(end) = [];

figure; hold on;
b11 = bar(edges,n_c1,'histc');
b12 = bar(edges,n_c2,'histc','FaceColor',c2);

set(b11, 'FaceColor', [0.8 0.8 0.8], 'LineWidth',1.5);
set(b12, 'FaceColor', [0.7 0 0],'FaceAlpha',0.4,'LineWidth',1.5);
set(gca,'TickDir','out','XLim', [0 1.1],'XTick',0:0.2:1.1,'FontSize',13)

xlabel(gca,'Correlation','FontSize',14);
ylabel(gca,'Norm counts','FontSize',14);
title(gca, 'LGN Traces Correlation');
box off

lg = legend('Cell2','Cell3','Location','NorthWest');
lg.Box = 'off';


% model correlations   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
edges = linspace(0,1,50);

n_c1 = histcounts(var_corr_nnmf(:,1),edges);
n_c2 = histcounts(var_corr_nnmf(:,2),edges);
n_c1 = n_c1/max(n_c1);
n_c2 = n_c2/max(n_c2);
edges(end) = [];

figure; hold on;
b11 = bar(edges,n_c1,'histc');
b12 = bar(edges,n_c2,'histc','FaceColor',c2);

set(b11, 'FaceColor', [0.8 0.8 0.8], 'LineWidth',1.5);
set(b12, 'FaceColor', [0.7 0 0],'FaceAlpha',0.4,'LineWidth',1.5);
set(gca,'TickDir','out','XLim', [0 1.1],'XTick',0:0.2:1.1,'FontSize',13)

xlabel(gca,'Correlation','FontSize',14);
ylabel(gca,'Norm counts','FontSize',14);
title(gca, 'NNMF Weights Correlation');
box off

lg = legend('Cell2','Cell3','Location','NorthWest');
lg.Box = 'off';


% model correlations   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
edges = linspace(0,1,50);

n_c1 = histcounts(var_corr_w(:,1),edges);
n_c2 = histcounts(var_corr_w(:,2),edges);
n_c1 = n_c1/max(n_c1);
n_c2 = n_c2/max(n_c2);
edges(end) = [];

figure; hold on;
b11 = bar(edges,n_c1,'histc');
b12 = bar(edges,n_c2,'histc','FaceColor',c2);

set(b11, 'FaceColor', [0.8 0.8 0.8], 'LineWidth',1.5);
set(b12, 'FaceColor', [0.7 0 0],'FaceAlpha',0.4,'LineWidth',1.5);
set(gca,'TickDir','out','XLim', [0 1.1],'XTick',0:0.2:1.1,'FontSize',13)

xlabel(gca,'Correlation','FontSize',14);
ylabel(gca,'Norm counts','FontSize',14);
title(gca, 'Model Weights Correlation');
box off

lg = legend('Cell2','Cell3','Location','NorthWest');
lg.Box = 'off';


% plot both together   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zero_mask = ones(nnz(fit_mask),1);
zero_mask(var_corr_mod(:,1) <= 0 | isnan(var_corr_mod(:,1))) = 0;
zero_mask(var_corr_mod(:,2) <= 0 | isnan(var_corr_mod(:,2))) = 0;

v1 = var_corr_nnmf(logical(zero_mask),:);
v2 = var_corr_lgn(logical(zero_mask),:);

figure; hold on;
scatter(mean(v1,2),mean(v2,2))
plot(0:0.1:1,0:0.1:1)

set(gca,'TickDir','out','XTick',0:0.25:1,'YTick',0:0.25:1)
xlabel(gca,'mean NNMF correlation','FontSize',14);
ylabel(gca,'mean Model correlation','FontSize',14);
title(gca, 'LGN Traces vs Model Correlation (mean from 2 cells)');







return 

%% SORT CORRELATIONS
cut_off = 0.5;

R_nmd_sorted = sort(R_nmf,2,'descend');
[~,ind] = sort(mean(R_nmd_sorted > 0.5,2));
figure; imagesc(R_nmd_sorted(ind,:)>0.5)


%%
for ip = 30%10:10:100
    
    ydata = tsne(R_nmd_sorted, [], 2, [], ip);
    F = ksdensity(ydata,'Bandwidth',4);
    
    figure
    %plot(ydata(:,1),ydata(:,2),'r.','MarkerSize',5), hold on
    imagesc(reshape(F,30,30));
    
    set(gca,'FontSize',11,'LineWidth',1.5,'TickDir','out')
    %set(gca,'XLim',[-60 60], 'YLim', [-60 60])
    box off; xlabel('dim 1','FontSize',13); ylabel('dim 2','FontSize',13)
    
    title(sprintf('tSNE, Perplexity = %d',ip))
end


















