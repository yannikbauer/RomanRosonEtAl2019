%% EVERYTHING TOGETHER
% TODO: 
% - improve description :-D
% - improve variable saving, e.g. using structure, to distinguish diff models

clear all; close all; clc;

%% Parameters

% Model pars
% Model types: 
% - 'lin': linear
% - 'lin_nonneg': linear w non-neg constraint
% - 'lin_lasso': linear lasso
model_type = 'lin'; % OPTIONS: 'lin', 'lin_nonneg', 'lin_lasso'

% Directory pars
data_dir = '../data/';   % figure folder name
filename = sprintf('model_%s', model_type);


%% Load data

% load('../data/model_lin_nonneg.mat'); % 6.5 GB, former file name: lin_range_o1g6
load(fullfile(data_dir, filename)); % Currently loads a whole workspace - CHANGE

%% Set parameters

weight_threshold = 0.2; % Set weight


str      = 'lin_range';    % text for the legend
titlestr = 'Linear Model';

afs = 11;  % axis font size
lfs = 15;  % lable font size
barsize = 1;

barcolor = [0.75 0.75 0.75];
maxRGC = max(o1rgc.cluIdx);

rgcstrc = {'1 OFF local, OS    1'; '2 OFF DS    2'; '3 OFF step    3'; '4a OFF slow    4';...
    '4b OFF slow    5';'5a OFF alpha sust.    6';'5b OFF alpha sust.    7';...
    '5c OFF alpha sust.    8';'6 (ON-)OFF "JAM-B" mix    9';'7 OFF sust.   10';...
    '8a OFF alpha trans.   11'; '8b OFF alpha trans.   12'; '9 OFF "mini" alpha trans.   13';...
    '10 ON-OFF local-edge "W3"   14';'11a ON-OFF local   15'; '11b ON-OFF local   16';...
    '12a ON-OFF DS 1   17'; '12b ON-OFF DS 1   18'; '13 ON-OFF DS 2   19'; ...
    '14 (ON-)OFF local, OS   20'; '15 ON step   21'; '16 ON DS trans.   22';...
    '17a ON local trans., OS   23'; '17b ON local trans., OS   24'; ...
    '17c ON local trans., OS   25'; '18a ON trans.   26'; '18b ON trans.   27';...
    '19 ON trans., large   28'; '20 ON high freq.   29'; '21 ON low freq.   30';...
    '22a ON sust.   31'; '22b ON sust   32'; '23 ON "mini" alpha   33'; '24 ON alpha   34'; ...
    '25 ON DS sust. 1   35'; '26 ON DS sust. 2   36'; '27 ON slow   37';...
    '28a ON contrast suppr.   38'; '28b ON contrast suppr.   39'; '29 ON DS sust. 3   40';...
    '30 ON local sust., OS   41'; '31a OFF suppr. 1   42'; '31b OFF suppr. 1   43';...
    '31c OFF suppr. 1   44'; '31d OFF suppr. 1   45'; '31e OFF suppr. 1   46';...
    '32a OFF suppr. 2   47'; '32b OFF suppr. 2   48'; '32c OFF suppr. 2   49'};

rgcstrg = {'OFF local, OS    1'; 'OFF DS    2'; 'OFF step    3'; ' OFF slow     4';...
    'OFF alpha sust.    5';'(ON-)OFF "JAM-B" mix    6';'OFF sust.   7';...
    'OFF alpha trans.   8'; 'OFF "mini" alpha trans.    9'; ' ON-OFF local-edge "W3"   10';...
    'ON-OFF local   11'; 'ON-OFF DS 1   12'; 'ON-OFF DS 2   13'; ...
    '(ON-)OFF local, OS   14'; 'ON step   15'; 'ON DS trans.   16';...
    'ON local trans., OS   17'; 'ON trans.   18'; 'ON trans., large   19';...
    'ON high freq.   20'; 'ON low freq.   21'; 'ON sust.   22'; 'ON "mini" alpha   23';...
    'ON alpha   24'; 'ON DS sust. 1   25'; 'ON DS sust. 2   26'; 'ON slow   27';...
    'ON contrast suppr.   28'; 'ON DS sust. 3   29'; 'ON local sust., OS   30';...
    'OFF suppr. 1   31'; 'OFF suppr. 2   32'};


load('cMaps/cMap_igor.mat')
color_ind = round(linspace(1,size(cMap_igor,1),maxRGC));
[figPars, axPars] = aux_setPlotPars();


%% normalize weights

[nneurons, ntypes, nrepeats] = size(var_w);
norm_w = var_w;
for ineuron = 1 : nneurons
    for irepeat = 1 : nrepeats
%         this_max = max(squeeze(var_w(ineuron, :,irepeat))); % normalize to max weight
        this_sum = sum(squeeze(var_w(ineuron, :,irepeat))); % normalize to sum of weights
%         norm_w(ineuron,:,irepeat) = var_w(ineuron, :,irepeat)/this_max;
        norm_w(ineuron,:,irepeat) = var_w(ineuron, :,irepeat)/this_sum;
    end
end

%% COMPUTE MEANS
mvu = mean(var_units,3);
mvy = mean(var_yhat,3);
mvc = mean(var_corr_lin,2);

mw  = mean(norm_w,3);
sd = std(norm_w,[],3);


%% PLOT EXAMPLE CELLS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = 29.7;
% fh = figure(figPars, 'Position', [17 5 54 h]);
fh = figure(figPars, 'Position', [2 2 21 h]);


% [~,corr_ind] = sort(var_corr_lin,'descend');
best_corr_ind = zeros(4,1);

% manual or automatic assigment
best_corr_ind(1) = 18; % old: 37
best_corr_ind(2) = 133; % 153
best_corr_ind(3) = 379-1; % 280
best_corr_ind(4) = 810-1; % 132
best_corr_ind(5) = 452; % 452

pos = linspace(4,22,5);
barwidth = 2.5;

time  = rgc.chirp_time;
for i = 1:5
    
    f5 = figure; hold on
    
    iunit = best_corr_ind(i);
    minval = min([min(mvu(iunit,:)) min(mvy(iunit,:))]);

    % CELL RESPONSE VS MODEL  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    plot(time,mvu(iunit,:),'k','LineWidth',1.5)   % data
    plot(time,mvy(iunit,:),'LineWidth',1.5)       % lin
    
    figHandles = findall(f5, 'Type', 'axes');
    newT5 = copyobj(figHandles(1), fh);
    set(newT5, axPars, 'Position', [2 h-pos(i) 9 2]);
    close(f5)
    
    set(newT5,'TickDir','out','Xlim',[0 32],'YLim',[minval 1],'FontSize',afs);
    set(newT5,'YTick',[0 1], 'XTick', [0 2]);
    lg = legend(newT5,'data',sprintf('LIN,p=%.3f',mvc(iunit)));
    lg.EdgeColor = 'w';
    
    
    % CELL WEIGHTS BAR PLOT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    f6 = figure; hold on
    if i == 4
        a = 1;
    end
%     vals = round(mw(iunit,:),2);
    vals = mw(iunit,:);
    ind = find(vals > weight_threshold);
    for ic = 1:length(ind)
        tmp = zeros(1,length(ind));
        tmp_mw = zeros(1,length(ind));
        tmp_sd = zeros(1,length(ind));
        
        tmp(ic) = vals(ind(ic));
        tmp_mw(ic) = mw(iunit,ind(ic));
        tmp_sd(ic) = sd(iunit,ind(ic))/2;
        
        b7 = bar(tmp);
        set(b7, 'FaceColor', cMap_igor(color_ind(cluIdx(ind(ic))),:),'LineWidth',1.5);
        errorbar(tmp_mw,tmp_sd,'k','marker','none','LineStyle','none')
        
    end
    
    figHandles = findall(f6, 'Type', 'axes');
    newT6 = copyobj(figHandles(1), fh);
    set(newT6, axPars, 'Position', [12 h-pos(i) 4 2]);
    close(f6)
    
    set(newT6,'FontSize',afs,'TickDir','out','YLim',[0 1],'YTick',[0 1],'XLim',[0.5 length(ind)+0.5])
    set(newT6,'XTick',1:1:length(ind),'XTicklabel',cluIdx(ind),'XTickLabelRotation',90)
    xlabel(newT6,'RGC Types','FontSize',lfs);
    
%     for ibar = 2:2:numel(newT6.Children)
%         newT6.Children(ibar).BarWidth = 1;%barwidth/mnc;
%     end
end

%% Plot mean # of clusters, correlation and RMSE

% CLEAN DATASETS FOR POPULATION PLOTS
mcorr = mean(var_corr_lin,2);
mrmse = mean(var_rmse_lin,2);
mmw   = mean(mean(norm_w,3));
msd   = mean(std(norm_w,[],3));


%% PLOT CORRELATION   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO: implement for varying thresholds
f1 = figure; hold on;

edges = linspace(minval,1,15);
n_all = histcounts(mcorr,edges);
n_all = n_all/sum(n_all)*100; % get percent cells
edges(1) = [];

minval = min(mcorr);

% plot hist
b1 = bar(edges,n_all,'histc');
%set(b1, 'FaceColor', [1 1 1], 'EdgeColor', [0 0 0], 'LineWidth',1.5);
set(b1, 'FaceColor', barcolor, 'EdgeColor', [1 1 1], 'LineWidth',1.5);

% plot mean
x = zeros(100,1);
x(:,1) = median(mcorr);
plot(x,linspace(0,ceil(max(n_all)),100)','--k','LineWidth',1.5)

% display
disp('Correlation')
fprintf('mean = %.2f\n',mean(mcorr));
fprintf('median = %.2f\n',median(mcorr));

set(gca, axPars)

figHandles = findall(f1, 'Type', 'axes');
newT1 = figHandles(1);

% YB: commented out for standalone figure
% newT1 = copyobj(figHandles(1), fh);
% set(newT1, axPars, 'Position', [18.5 h-6 6 5]);
% close(f1)

set(newT1,'TickDir','out','XLim', [minval 1],'XTick',round(minval*10)/10:0.2:1,'FontSize',afs)
%set(newT1,'YLim',[0 ceil(maxval)])
xlabel(newT1,'Correlation','FontSize',lfs)
ylabel(newT1,'% cells','FontSize',lfs)
title(newT1, 'Correlation')
lg = legend(newT1,sprintf('%.2f (n = %d)',round(median(mcorr),2),nnz(mcorr)));
lg.Box = 'off';
lg.Location = 'northwest';


%% PLOT RMSE   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO: implement for varying thresholds

f2 = figure; hold on;

minval = min(mrmse);
maxval = max(mrmse);

edges = linspace(-0.01,maxval,25);
n_all = histcounts(mrmse,edges);
n_all = n_all/sum(n_all)*100; % get percent cells
edges(1) = [];


% plot hist
b2 = bar(edges,n_all,'histc');
%set(b2, 'FaceColor', [1 1 1], 'EdgeColor', [0 0 0], 'LineWidth',1.5);
set(b2, 'FaceColor', barcolor, 'EdgeColor', [1 1 1], 'LineWidth',1.5);

% plot mean
x = zeros(100,1);
x(:,1) = median(mrmse);
plot(x,linspace(0,ceil(max(n_all)),100)','--k','LineWidth',1.5)

% display
disp('RMSE')
fprintf('mean = %.2f\n',mean(mrmse));
fprintf('median = %.2f\n',median(mrmse));

figHandles = findall(f2, 'Type', 'axes');
newT2 = figHandles(1);
set(gca, axPars)


% YB: commented out for standalone figure
% newT2 = copyobj(figHandles(1), fh);
% set(newT2, axPars, 'Position', [18.5 h-14 6 5]);
% close(f2)

set(newT2,'TickDir','out','XLim', [minval maxval],'FontSize',afs)
xlabel(newT2,'RMSE','FontSize',lfs)
ylabel(newT2,'# of cells (%)','FontSize',lfs)
title(newT2,'RMSE')
lg = legend(newT2,sprintf('%.2f',round(median(mrmse),2)));
lg.Box = 'off';
lg.Location = 'northwest';


%% PLOT # OF RGCs NEEDED FOR RECONSTRUCTION   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% mean # of RGC clusters
fprintf('Using weight treshold: %.3f\n', weight_threshold)
mnr = zeros(size(norm_w,1),1);
for i = 1:size(norm_w,1)
    mnr(i) = mean(sum(norm_w(i,:,:) > weight_threshold,2));
%     counts(:, ithreshold) = histcounts(mean(sum(norm_w > threshold,2),3), edges);
end


edges = min(round(mnr))-1:1:max(round(mnr))+1;
edges = edges+0.5;
% n_all = histcounts(round(mnr),edges);
n_all = histcounts(mnr,edges);

n_all = n_all/sum(n_all)*100; % get percent cells
edges(end) = [];

minval = min(round(mnr));
maxval = max(round(mnr));

f3 = figure; hold on;
% plot histo
b3 = bar(edges,n_all,'histc');
%set(b3, 'FaceColor', [1 1 1], 'EdgeColor', [0 0 0], 'LineWidth',1.5);
set(b3, 'FaceColor', barcolor, 'EdgeColor', [1 1 1], 'LineWidth',1.5);
set(gca,'Xtick',minval:2:maxval,'XLim',[minval-1 maxval+1])

% plot mean
x = zeros(100,1);
x(:,1) = median(mnr);
plot(x,linspace(0,ceil(max(n_all)),100)','--k','LineWidth',1.5)

% display
disp('Convergence')
fprintf('mean = %.2f\n',mean(mnr));
fprintf('median = %.2f\n',median(mnr));

figHandles = findall(f3, 'Type', 'axes');
newT3 = figHandles(1);
set(gca, axPars)

% YB: commented out for standalone figure
% newT3 = copyobj(figHandles(1), fh);
% set(newT3, axPars, 'Position', [18.5 h-22 6 5]);
% close(f3)

set(newT3,'TickDir','out','FontSize',afs)
xlabel(newT3,'# of RGCs','FontSize',lfs)
ylabel(newT3,'% of cells','FontSize',lfs)
title(newT3, 'RGC-LGN Convergence')
lg = legend(newT3,sprintf('%.2f',round(median(mnr),2)));
lg.Box = 'off';



%% PLOT MEAN WEIGHTS (CLUSTERS)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USE only for clusters; mean weights per group below

fprintf('Using weight treshold: %.3f\n', weight_threshold)

% vals = zeros(size(mw));
% vals(mw > weight_threshold) = mw(mw > weight_threshold);
% vals = mean(vals);

% Set weights below threshold = NaN
vals = nan(size(mw));
vals(mw > weight_threshold) = mw(mw > weight_threshold);
vals = nanmean(vals);

% extend matrix to size of clus
maxRGC = 49;
vals_ext = zeros(1,maxRGC);
vals_ext(cluIdx) = vals;
color_ind = round(linspace(1,size(cMap_igor,1),maxRGC));


f7 = figure; hold on;
% vals_ext = vals;
for i = 1:maxRGC
    tmp = zeros(size(vals_ext));
    tmp(i) = vals_ext(i);
    
    b7 = bar(tmp);
    set(b7, 'FaceColor', cMap_igor(color_ind(i),:),'LineWidth',1.5);
end

figHandles = findall(f7, 'Type', 'axes');
newT7 = figHandles(1);
set(newT7, axPars);

% newT7 = copyobj(figHandles(1), fh);
% set(newT7, axPars, 'Position', [2 h-28 20 3]);
% close(f7)

set(newT7,'FontSize',afs,'TickDir','out','Xlim',[0 maxRGC+0.5])
set(newT7,'XTick',1:1:maxRGC,'XTicklabel',1:1:49,'XTickLabelRotation',90)
set(newT7,'YAxisLocation','left')
xlabel(newT7,'RGC types','FontSize',lfs);
ylabel(newT7,'Mean weights','FontSize',lfs)

title(newT7,['Mean Weight Distribution, # RGC types = ' num2str(nRGC)]);


%% PLOT PERCENTAGE OF CELLS (CLUSTERS) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USE ONLY FOR CLUSTERS, percentage cells per group below

f8 = figure; hold on;

vals = round(sum(mean(norm_w,3) > weight_threshold)/size(mcorr,1)*100,1);
% Do for either clusters or groups

% extend matrix to size of clus
maxRGC = 49;
vals_ext = zeros(1,maxRGC);
vals_ext(cluIdx) = vals;
color_ind = round(linspace(1,size(cMap_igor,1),maxRGC));


for i = 1:maxRGC
    tmp = zeros(size(vals_ext));
    tmp(i) = vals_ext(i);
    
    b8 = bar(tmp);
    set(b8, 'FaceColor', cMap_igor(color_ind(i),:),'LineWidth',1.5);
end


figHandles = findall(f8, 'Type', 'axes');
newT8 = figHandles(1);
set(newT8, axPars);
% newT8 = copyobj(figHandles(1), fh);
% set(newT8, axPars, 'Position', [2 h-33.5 20 3]);
% close(f8)

set(newT8,'FontSize',afs,'TickDir','out','Xlim',[0 maxRGC+0.5])
set(newT8,'XTick',1:1:maxRGC,'XTicklabel',rgcstrc,'XTickLabelRotation',90)
set(newT8,'YAxisLocation','left')
xlabel(newT8,'RGC types','FontSize',lfs);
ylabel(newT8,'% of cells','FontSize',lfs)
t = title(newT8,'Cell/Weight Distribution');
t.FontSize = 14;


%% PLOT MEAN WEIGHTS (GROUPS)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For groups only; mean weights per cluster above

fprintf('Using weight treshold: %.3f\n', weight_threshold);

% % Set weights below threshold = 0
% vals = zeros(size(mw));
% vals(mw > weight_threshold) = mw(mw > weight_threshold);
% vals = mean(vals);

% Set weights below threshold = NaN
vals = nan(size(mw));
vals(mw > weight_threshold) = mw(mw > weight_threshold);
vals = nanmean(vals);

% Convert from clusters to groups
% TODO: find more elegant coding solution
vals_group = zeros(1,32);
vals_group(1)  = mean(vals(cluIdx == 1));
vals_group(2)  = mean(vals(cluIdx == 2));
vals_group(3)  = mean(vals(cluIdx == 3));
vals_group(4)  = mean([vals(cluIdx == 4) vals(cluIdx == 5)]);
vals_group(5)  = mean([vals(cluIdx == 6) 0 vals(cluIdx == 8)]);
vals_group(6)  = mean(0);
vals_group(7)  = mean(vals(cluIdx == 10));
vals_group(8)  = mean([vals(cluIdx == 11) vals(cluIdx == 12)]);
vals_group(9)  = mean(vals(cluIdx == 13));
vals_group(10) = mean(0);
vals_group(11) = mean([vals(cluIdx == 15) vals(cluIdx == 16)]);
vals_group(12) = mean([vals(cluIdx == 17) vals(cluIdx == 18)]);
vals_group(13) = mean(0);
vals_group(14) = mean(0);
vals_group(15) = mean(vals(cluIdx == 21));
vals_group(16) = mean(vals(cluIdx == 22));
vals_group(17) = mean([vals(cluIdx == 23) vals(cluIdx == 24) 0]);
vals_group(18) = mean([vals(cluIdx == 26) 0]);
vals_group(19) = mean(vals(cluIdx == 28));
vals_group(20) = mean(vals(cluIdx == 29));
vals_group(21) = mean(vals(cluIdx == 30));
vals_group(22) = mean([0 vals(cluIdx == 32)]);
vals_group(23) = mean(vals(cluIdx == 33));
vals_group(24) = mean(vals(cluIdx == 34));
vals_group(25) = mean(0);
vals_group(26) = mean(0);
vals_group(27) = mean(vals(cluIdx == 37));
vals_group(28) = mean([vals(cluIdx == 38) vals(cluIdx == 39)]);
vals_group(29) = mean(0);
vals_group(30) = mean(0);
vals_group(31) = mean([0 vals(cluIdx == 43) 0 0 0]);
vals_group(32) = mean([vals(cluIdx == 47) vals(cluIdx == 48) vals(cluIdx == 49)]);

f10 = figure; hold on;
color_ind = round(linspace(1,size(cMap_igor,1),32));

for i = 1:32
    tmp = zeros(size(vals_group));
    tmp(i) = vals_group(i);
    
    b9 = bar(tmp);
    set(b9, 'FaceColor', cMap_igor(color_ind(i),:),'LineWidth',1.5);
end

figHandles = findall(f10, 'Type', 'axes');
newT10 = figHandles(1);
set(newT10, axPars);
% newT10 = copyobj(figHandles(1), fh);
% set(newT10, axPars, 'Position', [25 h-15 18 3]);
% close(f10)

set(newT10,'FontSize',afs,'TickDir','out','Xlim',[0 32+0.5])
set(newT10,'XTick',1:1:32,'XTicklabel',rgcstrc,'XTickLabelRotation',90)
set(newT10,'YAxisLocation','left')
xlabel(newT10,'RGC types','FontSize',lfs);
ylabel(newT10,'Mean weights','FontSize',lfs)
t = title(newT10,'Cell/Weight Distribution');
t.FontSize = 14;

%% PLOT PERCENTAGE OF CELLS (GROUPED)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For groups only; percentage cells per cluster above

f10 = figure; hold on;
color_ind = round(linspace(1,size(cMap_igor,1),32));

% Get # cells w mean weight (across cross-validation repeats) above threshold
vals = sum(mw > weight_threshold);
% Convert to percentage
vals = round(vals/size(mw,1)*100, 1);

% vals = round(sum(mw > weight_threshold)*100/size(mw,1),1); % Same in one line

vals_group = zeros(1,32);
vals_group(1)  = mean(vals(cluIdx == 1));
vals_group(2)  = mean(vals(cluIdx == 2));
vals_group(3)  = mean(vals(cluIdx == 3));
vals_group(4)  = mean([vals(cluIdx == 4) vals(cluIdx == 5)]);
vals_group(5)  = mean([vals(cluIdx == 6) 0 vals(cluIdx == 8)]);
vals_group(6)  = mean(0);
vals_group(7)  = mean(vals(cluIdx == 10));
vals_group(8)  = mean([vals(cluIdx == 11) vals(cluIdx == 12)]);
vals_group(9)  = mean(vals(cluIdx == 13));
vals_group(10) = mean(0);
vals_group(11) = mean([vals(cluIdx == 15) vals(cluIdx == 16)]);
vals_group(12) = mean([vals(cluIdx == 17) vals(cluIdx == 18)]);
vals_group(13) = mean(0);
vals_group(14) = mean(0);
vals_group(15) = mean(vals(cluIdx == 21));
vals_group(16) = mean(vals(cluIdx == 22));
vals_group(17) = mean([vals(cluIdx == 23) vals(cluIdx == 24) 0]);
vals_group(18) = mean([vals(cluIdx == 26) 0]);
vals_group(19) = mean(vals(cluIdx == 28));
vals_group(20) = mean(vals(cluIdx == 29));
vals_group(21) = mean(vals(cluIdx == 30));
vals_group(22) = mean([0 vals(cluIdx == 32)]);
vals_group(23) = mean(vals(cluIdx == 33));
vals_group(24) = mean(vals(cluIdx == 34));
vals_group(25) = mean(0);
vals_group(26) = mean(0);
vals_group(27) = mean(vals(cluIdx == 37));
vals_group(28) = mean([vals(cluIdx == 38) vals(cluIdx == 39)]);
vals_group(29) = mean(0);
vals_group(30) = mean(0);
vals_group(31) = mean([0 vals(cluIdx == 43) 0 0 0]);
vals_group(32) = mean([vals(cluIdx == 47) vals(cluIdx == 48) vals(cluIdx == 49)]);

for i = 1:32
    tmp = zeros(size(vals_group));
    tmp(i) = vals_group(i);
    
    b9 = bar(tmp);
    set(b9, 'FaceColor', cMap_igor(color_ind(i),:),'LineWidth',1.5);
end

figHandles = findall(f10, 'Type', 'axes');
newT10 = figHandles(1);
set(newT10, axPars);
% newT10 = copyobj(figHandles(1), fh);
% set(newT10, axPars, 'Position', [25 h-15 18 3]);
% close(f10)

set(newT10,'FontSize',afs,'TickDir','out','Xlim',[0 32+0.5])
set(newT10,'XTick',1:1:32,'XTicklabel',rgcstrc,'XTickLabelRotation',90)
set(newT10,'YAxisLocation','right')
xlabel(newT10,'RGC types','FontSize',lfs);
ylabel(newT10,'% of cells','FontSize',lfs)
t = title(newT10,'Cell/Weight Distribution');
t.FontSize = 14;


% Extra:
% 1.choose best cells for the supplementary figure (relay mode)
% rmc_ind = find(sum(norm_w > 0.001,2) == 2 & norm_w(:,26) > 0.01 & norm_w(:,32) > 0.01);


%% PLOT CONVERGENCE MATRIX   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% z.B. Type A kommt insgesamt 3x vor, Type B kommt 2x vor, beide zusammen
% kommen 2x vor, dann h?tten wir 2/((3+2)/2) -> 80%

f4 = figure;

convMat = zeros(size(mw,2));
for ic = 1:size(mw,1)
    
    var   = mw(ic,:);
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
convMat_norm = zeros(size(norm_w,2));
for i = 1:size(norm_w,2)
    for j = 1:size(norm_w,2)
        a = convMat(i,j);
        b = sum(norm_w(:,i) > 0.001);
        c = sum(norm_w(:,j) > 0.001);        
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
figHandles = findall(f4, 'Type', 'axes');
newT4 = copyobj(figHandles(1), fh);
set(newT4, axPars, 'Position', [26 h-40 12 15]);
close(f4)

set(newT4,'XLim',[0.5 nRGC+0.5],'YLim',[0.5 nRGC+0.5]);
set(newT4,'FontSize',afs,'LineWidth',1.5,'TickDir','out')
set(newT4,'XTick',1:1:nRGC,'XTicklabel',cluIdx,'XTickLabelRotation',90)
set(newT4,'YTick',flipud(1:1:nRGC),'YTicklabel',cluIdx,'XTickLabelRotation',90)
title(newT4,'RGC cluster co-occurrence matrix')

daspect(newT4,[1 1 1]);
cl = colorbar(newT4);
set(cl,'TickDir','out','XTick',[0 (round(max(max(convMat_norm))*10))*10]);
box off;


%% PLOT CONVERGENCE PAIRS PER CLUSTER   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f9 = figure; hold on;
color_ind = round(linspace(1,size(cMap_igor,1),maxRGC));

vals = sum(convMat_norm);
vals = vals/max(vals);
for i = 1:length(vals)
    tmp = zeros(1,length(vals));
    tmp(i) = vals(i);
    
    b7 = bar(tmp);
    set(b7, 'FaceColor', cMap_igor(color_ind(cluIdx(i)),:),'LineWidth',1.5);
end

figHandles = findall(f9, 'Type', 'axes');
newT9 = copyobj(figHandles(1), fh);
set(newT9, axPars, 'Position', [26 h-26 10.5 2]);
close(f9)

set(newT9,'FontSize',afs,'TickDir','out','Xlim',[0 nRGC + 0.5])
set(newT9,'XTick',[],'XTickLabel',[],'YTick',[0 1],'YTickLabel',[0 1])
ylabel(newT9,'Mean weights','FontSize',lfs)

