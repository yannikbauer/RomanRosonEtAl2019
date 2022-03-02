

%% PARAMETERS
corr_p      = 0.001;         % 0.001; corr_p represents the wilcoxon ransum correlation for between/withing segments
qi          = 0.05;          % 0.05; Philipps quality index
cluster_qi  = 3;             % 3; cluster quality index of likely single- vs multiunit
r           = 0.65;          % 0.7; regression index

dt          = 0.001;         % sampling
timeLims    = [0 32.05];     % time limits of the stimulus

fprintf('\nFUNCTION: GET UNITS\n')


%%   GET SFTFORI KEYS

% COMPUTE VISUAL RESPONSIVENESS %
fprintf('... Getting SFTFOri units\n');

% x Number of points have to be over the threshold
nVis = 10; % 10% -> 10
allUnits = fetch(food.SFTFOriUnitStability(sprintf('r>=%f',r)) & data.ClusterInfo(sprintf('quality <= %f',cluster_qi)));
unit_mask = zeros(numel(allUnits), 1);

for iunit = 1:numel(allUnits)
    
    stimInfo = fetch(data.StimInfo(allUnits(iunit)) & fetch(data.SFTFOriTuning(allUnits(iunit))), '*'); % contains both SFTFORIexps
    
    % check for 2 experiments
    if numel(stimInfo) == 1
        continue
    end
    
    % Go through each experiment separately
    visResExp = zeros(2,1);
    for iexp = 1:2
        
        key = allUnits(iunit);
        key.exp_num = stimInfo(iexp).exp_num;
        
        c = fetch(data.ConditionSpikes(key), 'grat_num', 'cond_rate', 'cond_sem');
        activeLog = ismember([c.grat_num], [stimInfo.num_active_grats]);
        blankLog  = ismember([c.grat_num], [stimInfo.num_blank_grats]);
        
        % positive response || negative response
        if nnz([c(activeLog).cond_rate] - 2.58*[c(activeLog).cond_sem] > mean([c(blankLog).cond_rate])) > nVis || ...
                nnz([c(activeLog).cond_rate] + 2.58*[c(activeLog).cond_sem] < mean([c(blankLog).cond_rate])) > nVis
            visResExp(iexp) = true;
        end
    end
    
    if(visResExp(1) == 1 && visResExp(2) == 1)
        unit_mask(iunit) = true;
    end
end

allUnits(unit_mask == 0) = [];


%% GET NON-RESPONSIVE CHIRP KEYS

fprintf('... Getting chirp units\n');

% Create a table of units that respond satisfactorily to chirp stimulus
chirpUnits = fetch(miro.ChirpQuality(allUnits) ...
    & data.ClusterInfo(sprintf('quality <= %d',cluster_qi)) ...
    & miro.ChirpQuality(sprintf('corr_p > %f',corr_p)) ...
    & miro.ChirpQuality(sprintf('berens_qi < %f',qi)), ...
    'corr_p', 'berens_qi');


%% GET BETTER RSQ FOR SFTFORI STIM
sftfoTable = data.SFTFOriTuning(stripKey(chirpUnits,'data.ClusterInfo')) & food.SFTFOriUnitStability('r>=0.65');
unitList = fetch(data.ClusterInfo & sftfoTable);

oInfo = [];
for iunit = 1 : numel(unitList)
    sftfoInfo = fetch(sftfoTable & unitList(iunit), '*');
    oInfo = cat(2, oInfo, sftfoInfo([sftfoInfo.sftfo_ori_rsq] == max([sftfoInfo.sftfo_ori_rsq])));
end
[unitList.exp_num] = deal(oInfo.exp_num);

[ori_dp, dsi, osi, tf] = fetchn(data.SFTFOriTuning(unitList), 'sftfo_ori_dp', 'sftfo_ori_dsi', 'sftfo_ori_osi', 'sftfo_best_tf');


%% PLOT DS DISTRIBUTION

ds_keys = fetch(data.SFTFOriTuning(unitList), 'sftfo_ori_dp', 'sftfo_ori_dsi');

% ALL
figure;
circ_plot(circ_ang2rad(ori_dp),'hist',[],(0:pi/4:(2*pi - pi/4)),true,false,'linewidth',2,'color','r');
figure;circ_plot(circ_ang2rad(ori_dp),'pretty',[],(0:pi/4:(2*pi - pi/4)),true,false,'linewidth',2,'color','r');
nCells = length(ori_dp);
title(['ALL;  Cells: ' num2str(nCells)]); % 8 directions

figure
circ_plot(circ_ang2rad(ori_dp(dsi > 0.33)),'hist',[],(0:pi/4:(2*pi - pi/4)),true,false,'linewidth',2,'color','r');
figure;circ_plot(circ_ang2rad(ori_dp(dsi > 0.33)),'pretty',[],(0:pi/4:(2*pi - pi/4)),true,false,'linewidth',2,'color','r');
nCells = length(ori_dp(dsi > 0.33));
title(['DSI > 0.33; Cells: ' num2str(nCells)]);


%% EXTRACT 45-DEG DS CELLS

ind = ori_dp<67.5 & ori_dp>22.5 & dsi > 0.33;

% DS 45deg
figure;
ds45_keys = ds_keys(ind);
circ_plot(circ_ang2rad([ds45_keys.sftfo_ori_dp]),'hist',[],(0:pi/4:(2*pi - pi/4)),true,false,'linewidth',2,'color','-r');
hold on
circ_plot(circ_ang2rad([ds45_keys.sftfo_ori_dp]),'pretty',[],(0:pi/4:(2*pi - pi/4)),true,false,'linewidth',2,'color','-r');
title(['DS > 0.33; n = ' num2str(nnz(ind))])

% PLOT LINES INDICATING THE STRENGTH OF THE FIRING RATES


%% COMPUTER DS FR & SEM RATES FOR EACH CELL

sftfoTable45 = fetch(sftfoTable & stripKey(ds45_keys,'data.Units'), '*');

% get chirp psth for single traces
psth = parseChirp(stripKey(ds45_keys,'data.ClusterInfo'), dt, timeLims);
foo.chirp    = psth.psth;
foo.chirp_ts = psth.ts;

% get cell info
ori_frRate  = zeros(numel(ds45_keys), 8);
ori_semRate = zeros(numel(ds45_keys), 8);
ori_fit  = zeros(numel(ds45_keys), size(oInfo(1).sftfo_ori_fit,2));
tf_fit = zeros(numel(ds45_keys), size(oInfo(1).sftfo_tf_fit,2));

for iu = 1:numel(ds45_keys)
    
    sftfo_best_tf = sftfoTable45(iu).sftfo_best_tf;
    sftfo_best_sf = sftfoTable45(iu).sftfo_best_sf;
    
    % get firing & sem rates
    frInfo = fetch(data.ConditionSpikes(ds45_keys(iu)) * data.GratingConditions('grat_contrast > 0'), '*');
    blankRate = sftfoTable45(iu).sftfo_rblank;
    
    oLog = [frInfo.grat_y_frequency] == sftfo_best_sf & [frInfo.grat_drift_rate] == sftfo_best_tf;
    ori_frRate(iu,:) = [frInfo(oLog).cond_rate];
    ori_semRate(iu,:) = [frInfo(oLog).cond_sem];
    
    % get fits
    ori_fit(iu,:) = sftfoTable45(iu).sftfo_ori_fit;
    tf_fit(iu,:) = sftfoTable45(iu).sftfo_tf_fit;
end

foo.ori_frRate  = ori_frRate;
foo.ori_semRate = ori_semRate;
foo.ori_fit     = ori_fit;
foo.tf_fit      = tf_fit;


%% MAKE DS PLOT
ori_fit_norm = zeros(size(ori_fit));
for i = 1:size(ori_fit,1)
    ori_fit_norm(i,:) = (ori_fit(i,:)-min(ori_fit(i,:)))/(max(ori_fit(i,:))-min(ori_fit(i,:)));
end

figure
ftr = shadedErrorBar([],mean(ori_fit_norm),std(ori_fit_norm),'k');
set(ftr.mainLine,'LineWidth',2)
set(gca,'FontSize',11,'LineWidth',1.5,'TickDir','out','Xlim', [0 360])

ylabel('normalized sp/s','FontSize',14)
xlabel('Orientation (deg)','FontSize',14)
title('DS Plot')
box off;


%% PLOT CHIRP CLUSTER

warning off MATLAB:nargchk:deprecated

psth = parseChirp(stripKey(ds45_keys,'data.ClusterInfo'), dt, timeLims);
myData = psth.psth;

% computer chirp stimulus trace
[stim_t, stim_y, stim_cumTs] = plot_chirpStim();
stim_cumTs = stim_cumTs(1:end-1);

figure;
imagesc(psth.ts,1:size(myData,1),myData)
set(gca,'Visible','off')

hold on

% plot vertical stimulus lines
nLines = length(stim_cumTs);
vals = linspace(min(get(gca,'ylim')),max(get(gca,'ylim')),100);

ycor = zeros(100,nLines);
ycor(:,1:nLines)=repmat(vals',1,nLines);
xcor = ycor;
xcor(1:100,:)=repmat(stim_cumTs,100,1);
plot(xcor,ycor,':', 'color', [1 0 0],'LineWidth',1.5)

% divide space
orgPosition = get(gca, 'Position');
thisPosition = orgPosition;
thisPosition(4) = orgPosition(4)/2;
thisPosition(2) = orgPosition(2) + orgPosition(4)/2;
set(gca, 'Position', thisPosition);
axes('Position',[orgPosition(1),orgPosition(2),orgPosition(3),thisPosition(4)*1.3],'xlim',psth.tl);

%%%%%%%%%%%%%%%%%%%%%%

ftr = shadedErrorBar(psth.ts,median(myData),std(myData),'k');
set(ftr.mainLine,'LineWidth',2)
set(gca,'xlim',psth.tl,'ylim',[-0.4 1.3])
set(gca,'XColor','w','YColor','w')
box off

hold on

% plot vertical stimulus lines
nLines = length(stim_cumTs);
vals = linspace(min(get(gca,'ylim')),max(get(gca,'ylim')),100);

ycor = zeros(100,nLines);
ycor(:,1:nLines)=repmat(vals',1,nLines);
xcor = ycor;
xcor(1:100,:)=repmat(stim_cumTs,100,1);
plot(xcor,ycor,':', 'color', [1 0 0],'LineWidth',1.5)

st = supertitle(sprintf('# of cells: %d; qi: %.3f; corr_p: %.3f',size(myData,1),qi,corr_p),0.98);
set(st,'FontWeight','Bold','FontSize',12);



return



%% TF PLOTS (BEST TF AT THEIR PRFERRED ORI)

edges = [0 0.75 1.5 3 6 12 24];

% ALL
figure;
subplot(2,2,1)
bar(histcounts(tf, edges))
nCells = length(ori_dp);
title(['TF for ALL;  Cells: ' num2str(nCells)]); % 8 directions


% DSI > 0.33
subplot(2,2,2)
bar(histcounts(tf(dsi > 0.33), edges))
nCells = length(ori_dp(dsi > 0.33));
title(['TF for DSI > 0.33; Cells: ' num2str(nCells)]);

% OSI > 0.33
subplot(2,2,3)
bar(histcounts(tf(osi > 0.33), edges))
nCells = length(ori_dp(osi > 0.33));
title(['TF for OSI > 0.33; Cells: ' num2str(nCells)]);

% OSI <= 0.33 & DSI <= 0.33
subplot(2,2,4)
bar(histcounts(tf(osi <= 0.33 & dsi <= 0.33), edges))
nCells = length(ori_dp(osi <= 0.33 & dsi <= 0.33));
title(['TF for OSI <= 0.33 & DSI <= 0.33; Cells: ' num2str(nCells)]);


%% DS PLOTS

% ALL
figure;
subplot(2,2,1);
circ_plot(circ_ang2rad(ori_dp),'hist',[],(0:pi/4:(2*pi - pi/4)),true,false,'linewidth',2,'color','r');
nCells = length(ori_dp);
title(['ALL;  Cells: ' num2str(nCells)]); % 8 directions

% DSI > 0.33
subplot(2,2,2);
circ_plot(circ_ang2rad(ori_dp(dsi > 0.33)),'hist',[],(0:pi/4:(2*pi - pi/4)),true,false,'linewidth',2,'color','r');
nCells = length(ori_dp(dsi > 0.33));
title(['DSI > 0.33; Cells: ' num2str(nCells)]);

% OSI > 0.33
subplot(2,2,3);
circ_plot(circ_ang2rad(ori_dp(osi > 0.33)),'hist',[],(0:pi/4:(2*pi - pi/4)),true,false,'linewidth',2,'color','r');
nCells = length(ori_dp(osi > 0.33));
title(['OSI > 0.33; Cells: ' num2str(nCells)]);

% OSI <= 0.33 & DSI <= 0.33
subplot(2,2,4);
circ_plot(circ_ang2rad(ori_dp(osi <= 0.33 & dsi <= 0.33)),'hist',[],(0:pi/4:(2*pi - pi/4)),true,false,'linewidth',2,'color','r');
nCells = length(ori_dp(osi <= 0.33 & dsi <= 0.33));
title(['OSI <= 0.33 & DSI <= 0.33; Cells: ' num2str(nCells)]);

st = supertitle('Based on OSI & DSI',1);
set(st,'FontWeight','Bold','FontSize',12);

%
% edges = [0 22.5 67.5 112.5 157.5 202.5 247.5 292.5 337.5 360];
% figure;
% counts = histcounts(ori_dp(osi > 0.33), edges);
% counts(1) = counts(1) + counts(end);
% counts(end) = [];
% bar(counts)
% set(gca, 'xticklabels', [0 45 90 135 180 225 270 315])



















