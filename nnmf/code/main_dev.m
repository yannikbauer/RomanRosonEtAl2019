
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Choose data set to load

data_type = 'crossval'; % OPTIONS: 'crossval', 'original'

if data_type == 'original'
    % original components
    load('data/nnmf_best.mat') % Miro's 25 components
    A = A'; % A' = h | FxT (components x time)
    Y = Y'; % Y' = w | NxF (cells x components)
    
elseif data_type == 'crossval'
    % cross-validated components (comment out if nec)
    load('nnmf_cv_model.mat'); % cross-validated components
    A = Vt;
    Y = U;
    k_best = size(U,2);
end
%
%% COMPUTE FEATURES AND WEIGHTS
disp('... features and weights')

% normalize the weights to 1
Y_norm = zeros(size(Y));
for iunit = 1:size(Y,1)
    Y_norm(iunit,:) = (Y(iunit,:)-min(Y(iunit,:)))/(max(Y(iunit,:))-min(Y(iunit,:)));
end


%% WEIGHT DISTRIBUTION
% disp('... weight distribution')
% 
% edges   = -0.1:0.1:1;
% n = histcounts(Y_norm,edges);
% 
% n = n/max(n);
% edges(end) = [];
% 
% b = bar(edges,n,'histc');
% set(b, 'FaceColor', [0 0 1],'FaceAlpha',0.5,'LineWidth',1.5);
% set(gca,'TickDir','out','XLim', [0 1],'FontSize',13)
% box off;
% 
% xlabel('Weight','FontSize',14)
% ylabel('Distribution','FontSize',14)
% title('Weight Distribution')
% 
% 
% %% COMPLETE NNMF PLOT (DISPLACED, DENDRO, DISTRIBUTION)
% disp('... complete nnmf plot')
% 
% d = pdist(A, 'euclidean');
% tree = linkage(d,'ward');
% leafOrder = optimalleaforder(tree,d, 'CRITERIA', 'adjacent');
% 
% figure
% [~,~,OUTPERM] = dendrogram(tree, 0, 'Reorder', leafOrder, 'Orientation', 'left');
% set(gca,'TickDir','out','LineWidth',1);
% title('sparsennmf');
% 
% figure; PlotDisplaced(A(fliplr(OUTPERM),:));
% set(gca,'TickDir','out','LineWidth',1);
% title('sparsennmf');
% 
% nCompSum = sum(Y_norm,1);
% 
% figure
% b(1) =  barh(nCompSum(OUTPERM)/max(nCompSum));
% set(gca,'YLim',[0 k_best+1],'YTick',1:k_best,'TickDir','out','LineWidth',1);
% set(gca, 'YTickLabel', OUTPERM)
% box off;
% 
% num = round(max(nCompSum(OUTPERM))*100/size(Y,1),1); % percentage of the maximal component
% title(['sparsennmf, max feat: ' num2str(num) '%']);

%% COMPLETE NNMF PLOT (DISPLACED, DENDRO, DISTRIBUTION) AS ONE PLOT

do_one_plot = 1;
if do_one_plot
    disp('... complete nnmf plot')

    figure;
    subplot(1,9,1:2);
    d = pdist(A, 'euclidean');
    tree = linkage(d,'ward');
    leafOrder = optimalleaforder(tree,d, 'CRITERIA', 'adjacent');

    [~,~,OUTPERM] = dendrogram(tree, 0, 'Reorder', leafOrder, 'Orientation', 'left');
    set(gca,'TickDir','out','LineWidth',1);
    title('sparsennmf');

    subplot(1,9,3:8);
    PlotDisplaced(A(fliplr(OUTPERM),:));
    set(gca,'TickDir','out','LineWidth',1);
    title('sparsennmf');

    nCompSum = sum(Y_norm,1);

    subplot(1,9,9);
    b(1) =  barh(nCompSum(OUTPERM)/max(nCompSum));
    set(gca,'YLim',[0 k_best+1],'YTick',1:k_best,'TickDir','out','LineWidth',1);
    set(gca, 'YTickLabel', OUTPERM)
    box off;

    num = round(max(nCompSum(OUTPERM))*100/size(Y,1),1); % percentage of the maximal component
    title(['sparsennmf, max feat: ' num2str(num) '%']);
end
%% save fig
do_save = 1;
if do_save
    % save figure
    str = ['nnmf_cv_' num2str(size(U,2))];
    hgexport(gcf,fullfile('figures',str),hgexport('factorystyle'),'Format','eps');
%     close(gcf)
end

%% COMPUTE FEATURE PAIRS
% the matrix disps how often pairs of features occur together. # out of
% 815 (# of cells)
disp('... feature pairs')

mat = Y_norm';

cellpairs = zeros(size(mat,1));
pairnums = nchoosek(1:size(cellpairs,1),2);

for ipair = 1:size(pairnums,1)
    
    pair = pairnums(ipair,:);
    
    var1 = pair(1);
    var2 = pair(2);
    
    % compares the occurrence of the same unit pairs in one feature space
    x = nnz(mat(var1,:) & mat(var2,:));
    
    cellpairs(var1,var2) = x;
    cellpairs(var2,var1) = x;
end

figure; imagesc(cellpairs)
set(gca,'FontSize',11,'LineWidth',1.5,'TickDir','out')
title('Feature co-ccurrence matrix')
box off;

cl = colorbar;
nums = round(linspace(0,max(max(cellpairs)),5));
set(cl,'Ticks',nums,'TickLabels',nums,'TickDir','out')


%% COMPUTE CORRELATIONS MATRIX

figure
[R,P]=corrcoef(Y); imagesc(R);

set(gca,'FontSize',11,'LineWidth',1.5,'TickDir','out')
title('Feature correlation matrix')
box off;

cl = colorbar;
set(cl,'TickDir','out')


%% GM ALGORITHM
disp('... GM algorithm')

% first order sorting
Y_sorted = sort(Y_norm,2,'descend');
Y_sorted = Y_sorted'; % W * N

% normalize
Y_ns = zeros(size(Y_sorted)); % ns = norm & sorted
for iunit = 1:size(Y_sorted,2)
    Y_ns(:,iunit) = (Y_sorted(:,iunit)-min(Y_sorted(:,iunit)))/(max(Y_sorted(:,iunit))-min(Y_sorted(:,iunit)));
end

mod = performClustering(Y_ns');


%% Plot BIC
disp('... plot BIC')

figure; plot(mod.bic,'LineWidth',2)
set(gca,'FontSize',11,'LineWidth',1.5,'TickDir','out')
box off; xlabel('# of clusters','FontSize',13); ylabel('BIC','FontSize',13)
title('BIC')


%% PLOT WEIGHT-CLUSTERS - SORTED BY CLUSTERS
disp('... weight clusters - sorted by clusters')

% sort by idx
idx = mod.idx;
num = max(unique(idx));
order = zeros(1,num);
for imean = 1:num
    order(imean) = mean(sum(Y_ns(2:end,idx==imean)));
end

[~,ind_order] = sort(order); % cluster order

% create cluster matrix
mat = zeros(size(Y_ns,1),size(Y_ns,2));
vec_line = zeros(1,num-1);
ss = 1;
for iclu = 1:num
    tmp = Y_ns(:,idx == ind_order(iclu));
    si = ss+size(tmp,2)-1;
    mat(:,ss:si) = tmp;
    ss = ss+size(tmp,2);
    vec_line(iclu) = ss;
end

figure; imagesc(mat); colorbar; hold on
set(gca,'FontSize',11,'LineWidth',1.5,'TickDir','out')
box off; xlabel('# of units','FontSize',13); ylabel('# of Comp','FontSize',13)
title(sprintf('GMM Clustering of weights, k = %d',mod.K))

% draw line
nLines = length(vec_line);
vals = linspace(min(get(gca,'ylim')),max(get(gca,'ylim')),100);
ycor = zeros(100,nLines);
ycor(:,1:nLines)=repmat(vals',1,nLines);
xcor = ycor;
xcor(1:100,:)=repmat(vec_line,100,1);
plot(xcor,ycor,'-', 'color', [1 1 1],'LineWidth',1.5)


%% PLOT WEIGHT-CLUSTERS - SORTED BY SIZE

disp('... weight clusters - sorted by size')

[~,ind] = sort(mean(Y_ns),2,'descend');
figure; imagesc(Y_ns(:,ind)); colorbar;
set(gca,'FontSize',11,'LineWidth',1.5,'TickDir','out')
box off; xlabel('# of units','FontSize',13); ylabel('# of Comp','FontSize',13)
title('Sorted weights by size')


%% PLOT MEAN
disp('... plot mean')

mat = mat';
figure
sem = std(mat)/sqrt(length(mat)); % standard error of the mean
ftr = shadedErrorBar(1:1:size(mat,2),mean(mat),sem,'r','transparent');
set(ftr.mainLine,'LineWidth',1,'Color','k')
set(gca,'FontSize',11,'LineWidth',1.5,'TickDir','out')
set(gca,'xlim',[1 size(mat,2)])
box off; xlabel('# of comp','FontSize',13); ylabel('mean','FontSize',13)
title('averaged weights')


%% plot non-normalized weights

h = 10;
[figPars, axPars] = aux_setPlotPars();
fh = figure(figPars, 'Position', [24 18 20 h]);

% Plot weights
[~,ind] = sort(mean(sort(Y','descend')),'ascend');
tmp = sort(Y(ind,:)','descend');
f1 = figure; imagesc(tmp)

figHandles = findall(f1, 'Type', 'axes');
newT6 = copyobj(figHandles(1), fh);
set(newT6, axPars, 'Position', [2 h-8 14 7]);
close(f1)

title(newT6,'NNMF weight Distribution')

% Plot distribution of weights
dist = sum(tmp,2);
dist = dist/max(dist);

f2 = figure; barh(flipud(dist),'FaceColor', [1 1 1])
set(gca,'Ylim',[0 25.5])
axis off

figHandles = findall(f2, 'Type', 'axes');
newT6 = copyobj(figHandles(1), fh);
set(newT6, axPars, 'Position', [16 h-8 2.5 7]);

title(newT6,'Sum of weights')

close(f2)

%% FIT SUM OF WEIGHTS

x = 1:1:25;
y = dist';

f = @(b,x) b(1).*exp(-b(2).*x) + b(3);
nrmrsd = @(b) norm(y - f(b,x));
B0 = rand(3,1);
[B,rnrm] = fminsearch(nrmrsd, B0);
x_plot = linspace(min(x), max(x));
figure(1)
plot(x, y, 'pg')
hold on
plot(x_plot, f(B,x_plot), '-r')
hold off
grid
xlabel('Components')
ylabel('Weight')
legend('Data', 'Fit', 'Location','NE')


%% PLOT NUMBER OF WEIGHTS WITH THE MEAN
figure; hold on
dist = sum(tmp>0,2);
dist_n = dist/max(dist);
bar(dist_n)

x = zeros(100,1);
x(:,1) = mean(sum(tmp>0,1));
plot(x,linspace(0,1,100),'r','LineWidth',1.5)

set(gca,'TickDir','out', 'FontSize',12,'Xlim',[0 26],'YTick',[0 0.5 1])
set(gca,'Xtick',1:2:25,'XTickLabel',1:2:25)

ylabel('Components','FontSize',14)
xlabel('Counts','FontSize',14)
t = title('Number of weights > 0');
t.FontSize = 14;


%% cluster means
disp('... cluster means')
load(fullfile('cMaps', 'cMap_igor.mat'));
cInd = round(linspace(1,size(cMap_igor,1), num));
colorCode = cMap_igor(cInd,:);
CM = zeros(num,1);
lgd = cell(num,1);
figure; hold on
for iclu = 1:num
    tmp = Y_ns(:,idx == ind_order(iclu));
    tmp = tmp';
    
    sem = std(tmp)/sqrt(length(tmp)); % standard error of the mean
    ftr = shadedErrorBar(1:1:size(tmp,2),mean(tmp),sem,{'Color',colorCode(iclu,:)},'transparent');
    
    set(ftr.mainLine,'LineWidth',1)
    set(gca,'FontSize',11,'LineWidth',1.5,'TickDir','out')
    box off; xlabel('# of comp','FontSize',13); ylabel('mean','FontSize',13)
    title('averaged weights/group')
    
    CM(iclu) = plot(NaN,NaN,'Color',colorCode(iclu,:)); % legend
    lgd(iclu,1) = {sprintf('group %d',iclu)};
end
legend(CM, lgd, 'Location', 'NorthEast');


%% Convolve LGN components with ogb1 kernel and correlate features with RGC types

% order components
order = fliplr(OUTPERM);

% get the RGC cluster means
load('data/o1rgc.mat')
load('data/g6rgc.mat')
load('data/kerns.mat')

cluIdx = unique(g6rgc.cluIdx);
nClust = length(cluIdx);

% get RGC distribution
rgc_dist = zeros(nClust,1);
for iu=1:nClust
    rgc_dist(iu) = nnz(g6rgc.cluIdx == cluIdx(iu));
end

% get cluster means / get it from gcamp later
rgc_mean = zeros(nClust,length(o1rgc.chirp_time));
for iclu = 1:nClust
    rgc_mean(iclu,:) = mean(o1rgc.chirp(:,o1rgc.cluIdx==cluIdx(iclu)),2);
end

% downsample components
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

% correct for convolution artifact
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
    
    % mean normalization
    comp_norm(icomp,:) = trace - mean(trace(1:8));
    comp_norm(icomp,:) = comp_norm(icomp,:) / max(abs(comp_norm(icomp,:)));
end

% Correlate components with RGC Cluster means
corrMat = zeros(size(A,1),size(rgc_mean,1));
for i = 1:size(A,1)
    for j = 1:size(rgc_mean,1)
        corrMat(i,j) = corr(comp_norm(i,:)',rgc_mean(j,:)');
    end
end

%  PLOT STUFF
h = 45;
[figPars, axPars] = aux_setPlotPars();
fh = figure(figPars, 'Position', [32 14 21 h]);


% Plot Correlation matrix
f1 = figure;

corrMat(corrMat < 0) = 0;
imagesc(corrMat(order,:))
figHandles = findall(f1, 'Type', 'axes');
newT1 = copyobj(figHandles(1), fh);
set(newT1, axPars, 'Position', [3 h-9 16 8]);
close(f1)

set(gca,'TickDir','out','YTick',1:2:25,'XTick',1:1:length(cluIdx),'FontSize',10)
set(gca,'XTickLabel',cluIdx,'XTickLabelRotation',90)
ylabel(gca,'NNMF components','FontSize',14)
xlabel(gca,'RGC clusters','FontSize',14)
t = title(gca,'Correlation of NNFM components & RGC clusters');
t.FontSize = 14;

daspect([1 1 1])
colorbar(gca)

% plot best corr
f2 = figure; hold on

[vals,ind_rgc] = max(corrMat,[],2);
labels = cell(1,length(ind_rgc));
for i = 1:length(ind_rgc)
    labels(i) = {num2str(cluIdx(ind_rgc(i)))};
end

b1 = bar(vals(order));
set(b1, 'FaceColor', [0.00 0.44 0.74], 'LineWidth',1.5);
set(gca,'XTick',1:1:25)
xt = get(gca, 'XTick');
text(xt, vals(order), labels(order), 'HorizontalAlignment','center', 'VerticalAlignment','bottom')

figHandles = findall(f2, 'Type', 'axes');
newT2 = copyobj(figHandles(1), fh);
set(newT2, axPars, 'Position', [3 h-14 14.2 2.5]);
close(f2)

set(newT2,'Xlim',[0 26],'YLim',[0 1])
set(newT2,'TickDir','out','XTick',1:2:25,'YTick',[0 0.5 1],'FontSize',12)
ylabel(newT2,'Correlation','FontSize',14)
xlabel(newT2,'NNMF components','FontSize',14)

t = title(newT2,'Best NNMF-RGC Cluster Correlation');
t.FontSize = 14;

% plot example componnets
[vals, ind_rgc] = max(corrMat,[],2);
[~,ind_clu] = sort(vals,'descend');

for i = 1:1:5
    f3 = figure; hold on
    plot(o1rgc.chirp_time,comp_norm(ind_clu(i),:), 'LineWidth',1.5);
    plot(o1rgc.chirp_time,rgc_mean(ind_rgc(ind_clu(i)),:), 'LineWidth',1.5);
    
    minval = min([rgc_mean(ind_rgc(ind_clu(i)),:) comp_norm(ind_clu(i),:)]);
    maxval = max([rgc_mean(ind_rgc(ind_clu(i)),:) comp_norm(ind_clu(i),:)]);
    
    figHandles = findall(f3, 'Type', 'axes');
    newT3 = copyobj(figHandles(1), fh);
    set(newT3, axPars, 'Position', [3 h-15-i*5 14.2 2.5]);
    close(f3)
    
    set(newT3,'TickDir','out','YLim',[minval maxval])
    set(newT3,'XLim', [0 32],'XTick',1:4:32,'FontSize',12)
    xlabel(newT3,'Time (s)','FontSize',14)
    lg = legend(newT3,sprintf('nnmf comp # %d',find(ind_clu(i) == order)),sprintf('RGC Cluster # %d',cluIdx(ind_rgc(ind_clu(i)))));
    
    t = title(newT3,sprintf('Best RGC Cluster/NNMF Component, corr = %.2f',vals(ind_clu(i))));
    t.FontSize = 14;
end

return


%% example reconstructions

% pool of example cells
% cells = [153 179 217 252 279 315 316 317 318 319 349 352 372 420 439 441 ...
%     446 456 462 467 471 475 484 489 493 494 496 508 509 517 522 526 527 ...
%     528 530 543 546 549 552 553 554 561 575 576 578 579 581 583 589 620 ...
%     623 637 667 672 696 700 702 726 728 732 733 734 735 760 762 772 785 ...
%     790 791 793 795 801 805 806 808 811 812 813 814];

cells = [226 507 537];

for ic = 1:length(cells)
    
    % plot reconstruction
    figure;
    %set(gcf,'Position',[245 572 1528 275])
    %subplot(1,2,1)
    hold on
    plot(psth.ts,X(:,cells(ic)))
    plot(psth.ts,Y(cells(ic),:)*A)
    
    set(gca,'FontSize',11,'TickDir','out','XLim', [0 32])
    xlabel('time (s)','FontSize',13);
    ylabel('normalized','FontSize',13)
    title(['cell # ' num2str(cells(ic))]);
    legend('dLGN cell','nnmf reconstruction')
    box off;
    
    figure;
    %subplot(1,2,2)
    dist = Y(cells(ic),:);
    bar(dist(fliplr(OUTPERM)))
    set(gca,'FontSize',11,'TickDir','out','XLim',[0 18],'YLim',[0 1])
    set(gca,'XTick',1:1:17,'XTickLabel',{fliplr(OUTPERM)})
    xlabel('Features','FontSize',13);
    ylabel('Weight','FontSize',13)
    title(['cell # ' num2str(cells(ic))]);
    box off;
end




% ADDITIONAL COMP STUFF

return
%% plot more cells
cells = [12 18 22 37 74 133 140 146 147 150 151 152 153 154 169 190 193 226 ...
    227 229 237 240 250 252 254 256 259 280 281 297 302 332 346 350 358 366 ...
    379 398 403 407 408 409 411 449 453 460 480 483 486 488 489 490 491 495 ...
    507 517 529 532 534 537 543 549 553 583 616 638 641 646 658 666 673 674 ...
    678 679 680 714 717 720 722 725 726 727 759 760 778 781];

for ic = 1:length(cells)
    figure;
    plot(psth.ts,X(:,cells(ic)))
    set(gca,'FontSize',11,'TickDir','out','XLim', [0 32])
    xlabel('time (s)','FontSize',13);
    ylabel('normalized','FontSize',13)
    title(['cell # ' num2str(cells(ic))]);
    box off;
end


%% plot an overview of all cells
cc = 1;
for i=1:8
    figure
    set(gcf,'Position',[16 765 2518 575])
    for j = 1:98
        subplot(7,14,j,'align')
        plot(psth.ts,X(:,cc))
        set(gca,'XTick',[],'YTick',[],'XLim', [0 32])
        title(['cell # ' num2str(cc)]);
        cc = cc + 1;
    end
end


%% PERFORM tSNE ANALYSIS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('... perform ')
for ip = 30%20:4:60
    
    ydata = tsne(psth.psth, [], 2, [], ip);
    F = ksdensity(ydata,'Bandwidth',4);
    
    figure
    %plot(ydata(:,1),ydata(:,2),'r.','MarkerSize',5), hold on
    imagesc(reshape(F,30,30));
    
    set(gca,'FontSize',11,'LineWidth',1.5,'TickDir','out')
    %set(gca,'XLim',[-60 60], 'YLim', [-60 60])
    box off; xlabel('dim 1','FontSize',13); ylabel('dim 2','FontSize',13)
    
    title(sprintf('tSNE, Perplexity = %d',ip))
end


%% PLOT MULTIPLE NNMF FEATURES
warning off MATLAB:nargchk:deprecated

savepath = '/Users/miro/Google Drive/Lab_Work/Matlab/LGN_Cluster/Dynamic/figures';
[~, ~, stim_cumTs] = plotChirpStim();

% use seed from the sequence
S = RandStream('mt19937ar','Seed',seed);
RandStream.setGlobalStream(S);

X = psth.psth';

for ik = 2:30
    
    % A: matrix, the basis matrix
    % Y: matrix, the coefficient matrix
    % Sparse-NMF optimized by NNLS
    [A, Y] = sparsenmfnnls(X,ik); txt = 'sparsenmfnnls';
    
    A = A'; % A' = h
    Y = Y'; % Y' = w;
    
    d = pdist(A, 'euclidean');
    tree = linkage(d,'ward');
    leafOrder = optimalleaforder(tree,d, 'CRITERIA', 'adjacent');
    
    [~,~,OUTPERM] = dendrogram(tree, 0, 'Reorder', leafOrder, 'Orientation', 'left');
    close(gcf)
    
    
    figure; PlotDisplaced(psth.psth_ts,A(OUTPERM,:));
    
    set(gca,'FontSize',11,'LineWidth',1.5,'TickDir','out')
    box off; xlabel('Time (s)','FontSize',14);
    set(gca,'YColor',[1 1 1],'YTick',[])
    title(['NNMF, k = ', num2str(ik)])
    
    % plot vertical stimulus lines
    nLines = length(stim_cumTs);
    vals = linspace(min(get(gca,'ylim')),max(get(gca,'ylim')),100);
    
    ycor = zeros(100,nLines);
    ycor(:,1:nLines)=repmat(vals',1,nLines);
    xcor = ycor;
    xcor(1:100,:)=repmat(stim_cumTs,100,1);
    plot(xcor,ycor,':', 'color', [0.7 0.7 0.7],'LineWidth',1.5)
    
    % save figure
    str = ['nnmf_' num2str(ik)];
    hgexport(gcf,fullfile('figures',str),hgexport('factorystyle'),'Format','eps');
    close(gcf)
end































