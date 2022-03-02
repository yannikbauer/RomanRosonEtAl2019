
%% Load data
% NOTES: YB: seems redundant w make_fig_x.m

%load('workspace/lin_range_o1g6');


%% Compute stuff

mcorr = mean(var_corr_lin,2);
mrmse = mean(var_rmse_lin,2);
mw    = mean(mean(var_w,3));
msd   = mean(std(var_w,[],3));

% mean # of RGC clusters
mnr = zeros(size(var_w,1),1);
for i = 1:size(var_w,1)
    mnr(i) = mean(sum(var_w(i,:,:) > 0.01,2));
end


%% sort

% no model
ind = find(mnr == 0);
mcorr(ind) = [];
mrmse(ind) = [];
mnr(ind,:) = [];


%% POPULATION PLOTS

h = 29;
[figPars, axPars] = aux_setPlotPars();
fh2 = figure(figPars, 'Position', [24 17 29 h]);


% PLOT RMSE   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f1 = figure; hold on;

minval = min(mrmse);
maxval = max(mrmse);

edges = linspace(-0.01,maxval,25);
n_all = histcounts(mrmse,edges);
n_all = n_all/max(n_all);
edges(1) = [];

b1 = bar(edges,n_all,'histc');
set(b1, 'FaceColor', [0.8 0.8 0.8],'LineWidth',1.5);

figHandles = findall(f1, 'Type', 'axes');
newT1 = copyobj(figHandles(1), fh2);
set(newT1, axPars, 'Position', [2 h-8 7 6]);
close(f1)

set(newT1,'TickDir','out','XLim', [minval maxval],'FontSize',12)
xlabel(newT1,'RMSE','FontSize',15)
ylabel(newT1,'Norm counts','FontSize',15)
title(newT1,'RMSE')
lg = legend(newT1,sprintf('median = %.2f, n = %d',round(median(mrmse),2), length(mcorr)));
lg.Box = 'off';
lg.Location = 'northwest';


% PLOT CORRELATION   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f2 = figure; hold on;

minval = min(mcorr);

edges = linspace(minval,1,20);
n_all = histcounts(mcorr,edges);
n_all = n_all/max(n_all);
edges(1) = [];

b2 = bar(edges,n_all,'histc');
set(b2, 'FaceColor', [0.8 0.8 0.8],'LineWidth',1.5);

figHandles = findall(f2, 'Type', 'axes');
newT1 = copyobj(figHandles(1), fh2);
set(newT1, axPars, 'Position', [11 h-8 7 6]);
close(f2)

set(newT1,'TickDir','out','XLim', [minval 1],'XTick',round(minval*10)/10:0.2:1,'FontSize',13)
xlabel(newT1,'Correlation','FontSize',14)
ylabel(newT1,'Norm counts','FontSize',14)
title(newT1, 'Correlation')
lg = legend(newT1,sprintf('median = %.2f, n = %d', round(median(mcorr),2), length(mcorr)));
lg.Box = 'off';
lg.Location = 'northwest';


% PLOT # OF RGCs NEEDED FOR RECONSTRUCTION   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f3 = figure; hold on;

edges = min(round(mnr))-1:1:max(round(mnr))+1;
n_all = histcounts(round(mnr),edges);
n_all = n_all/max(n_all);
edges(end) = [];

minval = min(round(mnr));
maxval = max(round(mnr));

b3 = bar(edges,n_all);
set(b3, 'FaceColor', [0.8 0.8 0.8],'LineWidth',1.5);
set(gca,'Xtick',0:2:maxval,'XLim',[0 maxval+1])

figHandles = findall(f3, 'Type', 'axes');
newT6 = copyobj(figHandles(1), fh2);
set(newT6, axPars, 'Position', [20 h-8 7 6]);
close(f3)

set(newT6,'TickDir','out','FontSize',13)
xlabel(newT6,'# of RGCs','FontSize',14)
ylabel(newT6,'Norm counts','FontSize',14)
title(newT6, 'RGC-LGN Convergence')
lg = legend(newT6,sprintf('median %.2f, n = %d',round(median(mnr),2), length(mcorr)));
lg.Box = 'off';


% PLOT MEAN WEIGHTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f4 = figure; hold on;

b4 = bar(mw);
set(b4, 'FaceColor', [0.8 0.8 0.8], 'LineWidth',1.5);
errorbar(mw,msd/2,'k','marker','o','MarkerSize',6,'LineStyle','none')

figHandles = findall(f4, 'Type', 'axes');
newT4 = copyobj(figHandles(1), fh2);
set(newT4, axPars, 'Position', [2 h-17 20 6]);
close(f4)

set(newT4,'FontSize',8,'TickDir','out')
set(newT4,'XTick',1:1:nRGC,'XTicklabel',cluIdx,'XTickLabelRotation',90)
xlabel(newT4,'RGC types','FontSize',14);
ylabel(newT4,'Mean weights','FontSize',14)

t = title(newT4,['Mean Weight Distribution, # RGC types = ' num2str(nRGC)]);
t.FontSize = 14;


% PLOT WEIGHT/PERCENTIGE OF CELLS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f5 = figure; hold on;

wpc = round(sum(mean(var_w,3) > 0.001)*100/size(mcorr,1),1);

b5 = bar(wpc);
set(b5, 'FaceColor', [0.8 0.8 0.8],'LineWidth',1.5);

figHandles = findall(f5, 'Type', 'axes');
newT5 = copyobj(figHandles(1), fh2);
set(newT5, axPars, 'Position', [2 h-26 20 6]);
close(f5)

set(newT5,'FontSize',8,'TickDir','out','XTickLabelRotation',90)
set(newT5,'XTick',1:1:nRGC,'XTicklabel',cluIdx)
xlabel(newT5,'RGC types','FontSize',14);
ylabel(newT5,'# of cells (%)','FontSize',14)

t = title(newT5,'Cell/Weight Distribution');
t.FontSize = 14;





