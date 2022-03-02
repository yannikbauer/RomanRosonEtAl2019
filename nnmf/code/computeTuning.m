

% load unit keys
load('data/nnmf_4tuning.mat')

% restricted table from SFTFOriTuning
sftfoTable = data.SFTFOriTuning(unitList) & food.SFTFOriUnitStability('r>=0.65');
list = fetch(data.ClusterInfo & sftfoTable);


%% choose experiment with better rsq for ori tuning
warning('off','DataJoint:longCondition')

oInfo = [];
cc = 10;
for iunit = 1 : numel(list)
    if mod(iunit,floor(numel(list)/10)) == 0 % prints the progress
        fprintf('... %d%% completed\n',cc)
        cc = cc+10;
    end
    sftfoInfo = fetch(sftfoTable & list(iunit), '*');
    oInfo = cat(2, oInfo, sftfoInfo([sftfoInfo.sftfo_ori_rsq] == max([sftfoInfo.sftfo_ori_rsq])));
end
[list.exp_num] = deal(oInfo.exp_num);


%% DIVIDE CELLS BASED ON PERMUTATION TEST   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pval_os = zeros(numel(list),1);
pval_ds = zeros(numel(list),1);

cc = 10;
for iunit = 1:numel(list)
    
    if mod(iunit,floor(numel(list)/10)) == 0 % prints the progress
        fprintf('... %d%% completed\n',cc)
        cc = cc+10;
    end
    
    myUnit = list(iunit);
    
    % Get Stim info
    prefSF = oInfo(iunit).sftfo_best_sf;
    prefTF = oInfo(iunit).sftfo_best_tf;
    
    % Get Data
    kk = fetch(data.GratingConditions(myUnit) & sprintf('grat_drift_rate LIKE %d',prefTF)...
        & sprintf('grat_x_frequency LIKE %d',prefSF) & 'grat_contrast LIKE 1', 'grat_num');
    [kk.unit_id] = deal(myUnit.unit_id);
    
    % get orientation and spike times
    [ori,spikeTimes] = fetchn(data.GratingConditions(kk) * data.TrialSpikes(kk),'grat_orientation','spike_times');
    
    %  directions in radiants (#directions x 1)
    dirs = unique(deg2rad(ori));
    
    % get spike counts
    counts = cellfun(@length,spikeTimes);
    cpc = reshape(counts,20,8);  % counts per condition
    
    % test for the orientation variance
    p = testTuning(dirs, cpc, 2);
    pval_os(iunit) = p;
    
    % cell ds selective - test for the ds
    p = testTuning(dirs, cpc, 1);
    pval_ds(iunit) = p;
end


%% TREATMENT OF OS & DS POSITIVE CELLS
% Take the higher dsi or osi value as an indicator for better tuning

p = 0.001; % set threshold 

dsi   = [oInfo.sftfo_ori_dsi]';
osi   = [oInfo.sftfo_ori_osi]';

ind = find(pval_ds<p & pval_os<p);
pval_ds(ind(osi(ind) > dsi(ind)))  = 1;
pval_os(ind(dsi(ind) >= osi(ind))) = 1;


%% COMPARE PERMUTATION TEST WITH DSI
% the cells in the last bar are cells which were manually dicarded, because
% they had both ds and os significant. Their assigment was placed based on
% higher dsi or osi

thr   = 0.33;
edges = [-5 0.001:0.03:0.5 5];

% DSI
[n_pval_ds, bin] = histc(pval_ds, edges);

edges(1) = -0.03;
edges(end) = edges(end-1)+ 0.03;

cellsInGroup = zeros(length(unique(bin))+1,1);
for iclu = 1:max(unique(bin)+1)
    cellsInGroup(iclu,1) = nnz(dsi(bin == iclu) > thr);
end

figure; hold on
b(1) = bar(edges,n_pval_ds,'histc');
b(2) = bar(edges,cellsInGroup,'histc');

set(gca,'LineWidth',1.5,'FontSize',12,'XLim',[-0.03 0.51])
set(gca,'TickDir','out','XTick',edges(2:2:end))

set(b(1), 'FaceColor', [0.8 0.8 0.8]);
set(b(2), 'FaceColor', [0 0 1],'FaceAlpha',0.4);

xlabel('# of cells','FontSize',14);
ylabel('p values','FontSize',14);
title('P-Value & DSi distribution')


%% COMPARE PERMUTATION TEST WITH OSI

thr   = 0.33;
edges = [-5 0.001:0.03:0.5 5];

[n_pval_os, bin] = histc(pval_os, edges);

edges(1) = -0.03;
edges(end) = edges(end-1)+ 0.03;

cellsInGroup = zeros(length(unique(bin))+1,1);
for iclu = 1:max(unique(bin)+1)
    cellsInGroup(iclu,1) = nnz(osi(bin == iclu) > thr);
end

figure; hold on
b(1) = bar(edges,n_pval_os,'histc');
b(2) = bar(edges,cellsInGroup,'histc');

set(gca,'LineWidth',1.5,'FontSize',12,'XLim',[-0.03 0.51])
set(gca,'TickDir','out','XTick',edges(2:2:end))

set(b(1), 'FaceColor', [0.8 0.8 0.8]);
set(b(2), 'FaceColor', [0 0 1],'FaceAlpha',0.4);

ylabel('# of cells','FontSize',14);
xlabel('p values','FontSize',14);
title('P-Value & OSi distribution')



%% PLOTS

ori_dp = [oInfo.sftfo_ori_dp]';

% ALL
figure;
nCells = length(ori_dp);
circ_plot(circ_ang2rad(ori_dp),'hist',[],(0:pi/4:(2*pi - pi/4)),true,false,'linewidth',2,'color','r');
title(['All cells;  Cells: ' num2str(nCells)]);
figure;circ_plot(circ_ang2rad(ori_dp),'pretty',[],(0:pi/4:(2*pi - pi/4)),true,false,'linewidth',2,'color','r');
title(['All cells;  Cells: ' num2str(nCells)]);

% TUNED
figure;
nCells = length(ori_dp(pval_ds<p | pval_os<p));
circ_plot(circ_ang2rad(ori_dp(pval_ds<p | pval_os<p)),'hist',[],(0:pi/4:(2*pi - pi/4)),true,false,'linewidth',2,'color','r');
title(['OS & DS cells;  Cells: ' num2str(nCells)]);
figure;circ_plot(circ_ang2rad(ori_dp(pval_ds<p | pval_os<p)),'pretty',[],(0:pi/4:(2*pi - pi/4)),true,false,'linewidth',2,'color','r');
title(['OS & DS cells;  Cells: ' num2str(nCells)]);

% DS
figure;
nCells = length(ori_dp(pval_ds<p));
circ_plot(circ_ang2rad(ori_dp(pval_ds<p)),'hist',[],(0:pi/4:(2*pi - pi/4)),true,false,'linewidth',2,'color','r');
title(['DS cells;  Cells: ' num2str(nCells)]);
figure;circ_plot(circ_ang2rad(ori_dp(pval_ds<p)),'pretty',[],(0:pi/4:(2*pi - pi/4)),true,false,'linewidth',2,'color','r');
title(['DS cells;  Cells: ' num2str(nCells)]);

% OS
figure;
nCells = length(ori_dp(pval_os<p));
circ_plot(circ_ang2rad(ori_dp(pval_os<p)),'hist',[],(0:pi/4:(2*pi - pi/4)),true,false,'linewidth',2,'color','r');
title(['OS cells;  Cells: ' num2str(nCells)]);
figure;circ_plot(circ_ang2rad(ori_dp(pval_os<p)),'pretty',[],(0:pi/4:(2*pi - pi/4)),true,false,'linewidth',2,'color','r');
title(['OS cells;  Cells: ' num2str(nCells)]);


%% ANALYZE THE 45deg (315deg) DIRECTION

% ALL 45deg
vals_all = ori_dp(ori_dp<67.5 & ori_dp>22.5);
figure; circ_plot(circ_ang2rad(vals_all),'hist',[],(0:pi/4:(2*pi - pi/4)),true,false,'linewidth',2,'color','-r');
title(['All cells; n = ' num2str(length(vals_all))])
figure; circ_plot(circ_ang2rad(vals_all),'pretty',[],(0:pi/4:(2*pi - pi/4)),true,false,'linewidth',2,'color','-r');
title(['All cells; n = ' num2str(length(vals_all))])

% DS 45deg
vals_ds = ori_dp(ori_dp<67.5 & ori_dp>22.5 & pval_ds<p);
figure; circ_plot(circ_ang2rad(vals_ds),'hist',[],(0:pi/4:(2*pi - pi/4)),true,false,'linewidth',2,'color','-r');
title(['DS cells; n = ' num2str(length(vals_ds))])
figure; circ_plot(circ_ang2rad(vals_ds),'pretty',[],(0:pi/4:(2*pi - pi/4)),true,false,'linewidth',2,'color','-r');
title(['DS cells; n = ' num2str(length(vals_ds))])

% OS 45deg
vals_os = ori_dp(ori_dp<67.5 & ori_dp>22.5 & pval_os<p);
figure; circ_plot(circ_ang2rad(vals_os),'hist',[],(0:pi/4:(2*pi - pi/4)),true,false,'linewidth',2,'color','-r');
title(['OS cells; n = ' num2str(length(vals_os))])
figure; circ_plot(circ_ang2rad(vals_os),'pretty',[],(0:pi/4:(2*pi - pi/4)),true,false,'linewidth',2,'color','-r');
title(['OS cells; n = ' num2str(length(vals_os))])


% get keys & chirp exp_num
ds45keys = list(pval_ds<p & (ori_dp<67.5 & ori_dp>22.5));
os45keys = list(pval_os<p & (ori_dp<67.5 & ori_dp>22.5));

ds45keys = fetch(miro.ChirpQuality(stripKey(ds45keys,'data.ClusterInfo')),'corr_p','berens_qi');
os45keys = fetch(miro.ChirpQuality(stripKey(os45keys,'data.ClusterInfo')),'corr_p','berens_qi');

corr_p = 0.001;
qi = 0.05;

mask_chirp_ds = [ds45keys.corr_p] <= corr_p & [ds45keys.berens_qi] >= qi; % cell responding / not responding to chirp
mask_chirp_os = [os45keys.corr_p] <= corr_p & [os45keys.berens_qi] >= qi; % cell responding / not responding to chirp

disp('STATS:')
fprintf('\nAll Cells (tuned and not tuned): %d\n',length(vals_all))
fprintf('DS Cells: %d\n',length(vals_ds))
fprintf('OS Cells: %d\n',length(vals_os))

fprintf('%d/%d ds cells respond to chirp\n',nnz(mask_chirp_ds), length(mask_chirp_ds))
fprintf('%d/%d os cells respond to chirp\n',nnz(mask_chirp_os), length(mask_chirp_os))


%% COMPUTE DISTRIBUTION FOR FEATURES COMPUTATION

disp('compute feature distribution')
list = stripKey(list,'data.ClusterInfo');

% ALL CELLS
disp('... processing all cells')
mask_all = zeros(numel(unitList),1);
for i = 1:numel(unitList)
    for j = 1:numel(list)
        if isequal(unitList(i),list(j))
            mask_all(i) = 1;
            continue;
        end
    end
end

% NOT TUNED
disp('... processing not tuned cells')
mask_nt = zeros(numel(unitList),1);
list_tmp = list(pval_ds>=p & pval_os>=p);
for i = 1:numel(unitList)
    for j = 1:numel(list_tmp)
        if isequal(unitList(i),list_tmp(j))
            mask_nt(i) = 1;
            continue;
        end
    end
end

% OS + DS
disp('... processing OS & DS')
mask_osds = zeros(numel(unitList),1);
list_tmp = list(pval_ds<p | pval_os<p);
for i = 1:numel(unitList)
    for j = 1:numel(list_tmp)
        if isequal(unitList(i),list_tmp(j))
            mask_osds(i) = 1;
            continue;
        end
    end
end

% OS
disp('... processing OS')
mask_os = zeros(numel(unitList),1);
list_tmp = list(pval_os<p);
for i = 1:numel(unitList)
    for j = 1:numel(list_tmp)
        if isequal(unitList(i),list_tmp(j))
            mask_os(i) = 1;
            continue;
        end
    end
end

% DS
disp('... processing DS')
mask_ds = zeros(numel(unitList),1);
list_tmp = list(pval_ds<p);
for i = 1:numel(unitList)
    for j = 1:numel(list_tmp)
        if isequal(unitList(i),list_tmp(j))
            mask_ds(i) = 1;
            continue;
        end
    end
end


%% PLOT TUNING

nCompALL  = sum(Y_norm(mask_all==1,:),1);
nCompDS   = sum(Y_norm(mask_ds==1,:),1);
nCompOS   = sum(Y_norm(mask_os==1,:),1);
nCompOSDS = sum(Y_norm(mask_osds==1,:),1);
nCompNT   = sum(Y_norm(mask_nt==1,:),1);

figure; hold on
b(2) =  bar(nCompALL(OUTPERM)/max(nCompALL));
b(3) =  bar(nCompDS(OUTPERM)/max(nCompDS));

set(gca,'XLim',[0 k_best+1],'XTick',1:k_best,'TickDir','out','LineWidth',1);
set(gca, 'XTickLabel', OUTPERM,'FontSize',12)

set(b(2), 'FaceColor',  [0.7 0.7 0.7], 'LineWidth',1.5); % all data
set(b(3), 'FaceColor',[0 0 1], 'LineWidth',1.5);
set(b(3),'FaceAlpha',0.4)

box off;
xlabel('Features','FontSize',14)
ylabel('norm sum of weights','FontSize',14)
title('ALL vs DS')

%%%%%%%%%%%%%%%%%%%

figure; hold on
b(2) =  bar(nCompALL(OUTPERM)/max(nCompALL));
b(3) =  bar(nCompOS(OUTPERM)/max(nCompOS));

set(gca,'XLim',[0 k_best+1],'XTick',1:k_best,'TickDir','out','LineWidth',1);
set(gca, 'XTickLabel', OUTPERM,'FontSize',12)

set(b(2), 'FaceColor',  [0.7 0.7 0.7], 'LineWidth',1.5); % all data
set(b(3), 'FaceColor',[0 0 1], 'LineWidth',1.5);
set(b(3),'FaceAlpha',0.4)

box off;
xlabel('Features','FontSize',14)
ylabel('norm sum of weights','FontSize',14)
title('ALL vs OS')

%%%%%%%%%%%%%%%%%%%

figure; hold on
b(2) =  bar(nCompALL(OUTPERM)/max(nCompALL));
b(3) =  bar(nCompOSDS(OUTPERM)/max(nCompOSDS));

set(gca,'XLim',[0 k_best+1],'XTick',1:k_best,'TickDir','out','LineWidth',1);
set(gca, 'XTickLabel', OUTPERM,'FontSize',12)

set(b(2), 'FaceColor',  [0.7 0.7 0.7], 'LineWidth',1.5); % all data
set(b(3), 'FaceColor',[0 0 1], 'LineWidth',1.5);
set(b(3),'FaceAlpha',0.4)

box off;
xlabel('Features','FontSize',14)
ylabel('norm sum of weights','FontSize',14)
title('ALL vs OSDS')

%%%%%%%%%%%%%%%%%%%

figure; hold on
b(2) =  bar(nCompOS(OUTPERM)/max(nCompOS));
b(3) =  bar(nCompDS(OUTPERM)/max(nCompDS));

set(gca,'XLim',[0 k_best+1],'XTick',1:k_best,'TickDir','out','LineWidth',1);
set(gca, 'XTickLabel', OUTPERM,'FontSize',12)

set(b(2), 'FaceColor',  [0.7 0.7 0.7], 'LineWidth',1.5); % all data
set(b(3), 'FaceColor',[0 0 1], 'LineWidth',1.5);
set(b(3),'FaceAlpha',0.4)

box off;
xlabel('Features','FontSize',14)
ylabel('norm sum of weights','FontSize',14)
title('OS vs DS')

%%%%%%%%%%%%%%%%%%%

figure; hold on
b(2) =  bar(nCompALL(OUTPERM)/max(nCompALL));
b(3) =  bar(nCompNT(OUTPERM)/max(nCompNT));

set(gca,'XLim',[0 k_best+1],'XTick',1:k_best,'TickDir','out','LineWidth',1);
set(gca, 'XTickLabel', OUTPERM,'FontSize',12)

set(b(2), 'FaceColor', [0.7 0.7 0.7], 'LineWidth',1.5); % all data
set(b(3), 'FaceColor',[0 0 1], 'LineWidth',1.5);
set(b(3),'FaceAlpha',0.4)

box off;
xlabel('Features','FontSize',14)
ylabel('norm sum of weights','FontSize',14)
title('ALL vs NT')






























