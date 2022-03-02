% LOCOMOTION ANALYSIS

% speed                        # speed during trial (m/s)
% trial_direction=null         # direction during trial (rad)
% loco_ts                      # timestamps of speed measurements (s)
% trial_percent_moving         # percentage speed running during trial
% trial_movement_threshold     # running threshold
% trial_mean_speed             # mean speed of trial
% trial_max_speed              # mean direction of trial


%% GET KEYS FOR LOCOMOTION
corr_p      = 0.001;         % 0.001; corr_p represents the wilcoxon ransum correlation for between/withing segments
qi          = 0.05;          % 0.05; Philipps quality index
cluster_qi  = 3;             % 3; cluster quality index of likely single- vs multiunit
r           = 0.65;          % 0.7; regression index
type        = 3;

unitList = getAllUnits(corr_p, qi, cluster_qi, r, type);
skeys = fetch(data.Locomotion(unitList) & data.Locomotion('sampling_rate > 85'));
trial_keys = fetch(data.ExtraTrials(skeys),'trial_num');
disp('... trial_keys done')


%% GET CHIRP LOCOMOTION

nTrials = numel(trial_keys);

% initialize
loco_speed = [];
loco_ts = [];
loco_percent_moving = zeros(nTrials,1);
loco_mean_speed = zeros(nTrials,1);
loco_max_speed = zeros(nTrials,1);

% get sampling rate
loco_samp_rate = mean(fetchn(data.Locomotion(skeys),'sampling_rate'));

% compute trial info
cc = 10;
for itrial = 1:nTrials
    
    if mod(itrial,floor(nTrials/10)) == 0 % prints the progress
        fprintf('... %d%% completed\n',cc)
        cc = cc+10;
    end
    
    [speed, ts, mean_speed, max_speed, percent_moving] = p_trialLocomotion(trial_keys(itrial));
    
    if isempty(speed)
        warning('Locomotion computation failed for trial #%d',itrial)
        loco_percent_moving(itrial) = nan;
        loco_mean_speed(itrial) = nan;
        loco_max_speed(itrial) = nan;
        continue
    end
    
    % locomotion info
    loco_percent_moving(itrial) = percent_moving;
    loco_max_speed(itrial) = max_speed;
    loco_mean_speed(itrial) = mean_speed;
    
    if itrial > 1
        
        % speed
        if size(loco_speed(end,:),2) ~= size(speed,2)
            df = abs(size(loco_speed(end,:),2) - size(speed,2)); % compute difference
            if size(loco_speed(end,:),2) < size(speed,2)
                loco_speed = [loco_speed, zeros(size(loco_speed,1),df)];   %#ok<AGROW>
            else
                speed = [speed, zeros(1,df)];                              %#ok<AGROW>
            end
        end
        
        % timestamp
        if size(loco_ts(end,:),2) ~= size(ts,2)
            df = abs(size(loco_ts(end,:),2) - size(ts,2));
            if size(loco_ts(end,:),2) < size(ts,2)
                loco_ts = [loco_ts, zeros(size(loco_ts,1),df)];            %#ok<AGROW>
            else
                ts = [ts, zeros(1,df)];                                    %#ok<AGROW>
            end
        end
    end
    
    loco_speed = [loco_speed; speed];                                      %#ok<AGROW>
    loco_ts = [loco_ts; ts];                                               %#ok<AGROW>
end

loco_percent_moving(isnan(loco_percent_moving)) = [];
loco_mean_speed(isnan(loco_mean_speed)) = [];
loco_max_speed(isnan(loco_max_speed)) = [];


%% PLOT MEAN OF TRIALS

tl = [min(loco_ts(1,:)) max(loco_ts(1,1:end - 1))];

figure;
imagesc(loco_speed)
set(gca,'TickDir','out','YTick',[250 500 750 1000],'XTick',1:896:size(loco_ts,2),'XtickLabel',round(loco_ts(1,1:896:end)))
box off

ylabel(gca, '# of trials','FontSize', 14);
xlabel(gca,'Time (s)','FontSize', 14);
title('Trial locomotion')

colorbar

%% PLOT LOCOMOTION

[t, y, cumTs] = plotChirpStim();
stim_cumTs = cumTs(1:end-1);

h = 10.5;
[figPars, axPars] = aux_setPlotPars();

fh = figure(figPars, 'Position', [32 30 21 h]);

% CHIRP LOCOMOTION
floco = figure; hold on
sem = std(loco_speed(:,1:end - 1))/sqrt(length(loco_speed(:,1:end - 1))); % standard error of the mean
ftr = shadedErrorBar(loco_ts(1,1:end - 1),mean(loco_speed(:,1:end - 1)),sem,'k');
set(ftr.mainLine,'LineWidth',1,'Color','k')

% plot vertical stimulus lines
nLines = length(stim_cumTs);
vals = linspace(0,30,100);

ycor = zeros(100,nLines);
ycor(:,1:nLines)=repmat(vals',1,nLines);
xcor = ycor;
xcor(1:100,:)=repmat(stim_cumTs,100,1);
plot(xcor,ycor,':', 'color', [0.4 0.4 0.4],'LineWidth',1.5)

figHandles = findall(floco, 'Type', 'axes');
newT = copyobj(figHandles(1), fh);
set(newT,'Ylim',[0 30],'YScale', 'log', 'YTick', [0 1:10 20 30],'YMinorTick', 'off')
set(newT,'YTickLabel', {'0' , '1', '', '', '', '', '', '', '', '', '10', '', '30'})
set(newT, axPars, 'Position', [2 h-6 18 5]);

close(floco)

ylabel(newT,'Speed (cm/s)','FontSize', 14);
set(newT, 'XLim', tl, 'XColor','w','XTick',[]);
box off;

st = supertitle('Chirp Locomotion',0.99);
st.FontWeight = 'Bold';
st.FontSize = 15;

% HEAT MAP
% fhm = figure;
% imagesc(mean(loco_speed(:,1:end - 1)))
% 
% figHandles = findall(fhm, 'Type', 'axes');
% newM = copyobj(figHandles(1), fh);
% set(newM, axPars, 'Position', [2 h-7 18 0.5]);
% 
% close(fhm)
% 
% ylabel(newM,'Speed (cm/s)','FontSize', 14);
% newM.YAxis.Visible = 'off';
% newM.XAxis.Visible = 'off';


% STIMULUS

fch = figure; 
plot(t,y,'-k');
set(gca,'XLim', tl)

figHandles = findall(fch, 'Type', 'axes');
newB = copyobj(figHandles(1), fh);
set(newB, axPars, 'Position', [2 h-7.5 18 1], 'YLim', [-2 1]);

close(fch) 
hold on

newB.YAxis.Visible = 'off';
xlabel(newB,'Time (s)','FontSize', 14);

% mean trace correlation
ls_down = interp1(loco_ts(1,1:end - 1)',mean(loco_speed(:,1:end - 1)),t');

% remove nans
y(1) = [];
ls_down(1) = [];

[r,p] = corr(y',ls_down);

%==========================================================================

%% GET KEYS FOR STATISTICS

unit_keys = fetch(data.ClusterInfo(unitList) * data.Locomotion(unitList));
unit_keys = rmfield(unit_keys,'exp_num');

% remove sampling_rate error
key_mask = zeros(numel(unit_keys),1);
for ikey = 1:numel(unit_keys)
    
    myKey = unitList(ikey);
    try
        smr = fetch1(data.Locomotion(myKey),'sampling_rate');
    catch
        continue
    end
    if smr > 85
        key_mask(ikey) = 1;
    end
end
unit_keys(key_mask==0) = [];

% get Aspontaneous units (contains series with multiple Aspontaneous)
tmp = fetch(data.Experiments('exp_name LIKE "%spontaneous"') * data.ClusterInfo(unit_keys));

% single units for Aspontaneous with sinificat firing rate
unit_keys = fetch(food.RunSpeedTuningFit(tmp)); 

% save info in other keys
rstf_keys = fetch(food.RunSpeedTuningFit(unit_keys),'run_speedmodulated','run_rsq','run_tune_type','run_sp');
rst_keys = fetch(food.RunSpeedTuning(unit_keys),'run_meanfr');

disp('... unit_keys done')


%% SPEED MODULATED PLOT
modulated = [rstf_keys.run_speedmodulated] == 1;   % units modulated by speed
not_modulated = [rstf_keys.run_speedmodulated] == 0;  % not modulated units

figure;
b1 = bar([nnz(modulated) nnz(not_modulated)],'w');

set(gca,'YLim',[0 nnz(not_modulated)],'YTick',linspace(0,nnz(not_modulated),5))
set(gca,'TickDir','out','Box','off','LineWidth',2,'FontSize',12)
set(gca,'XtickLabel',{'modulated' 'not modulated'},'YTickLabel',0:25:100)
ylabel('# of Neurons in %','FontSize',14)
title('Neurons modulated by running')
b1.LineWidth = 2;


%% MODULATED UNIT EXAMPLE TUNING TYPES WITH GOOD RSQ

rsq_threshold = 0.65;
mod_keys = rstf_keys(modulated);

mod_rsq = mod_keys([mod_keys.run_rsq] > rsq_threshold);  % modulated units with a good fit

% 1 = monotone increase
mi_keys = mod_rsq([mod_rsq.run_tune_type] == 1);
[~, ind] = sort([mi_keys.run_rsq]);
plot(food.RunSpeedTuningFit(stripKey(mi_keys(ind(end)),'data.ClusterInfo')));

% 2 = monoton decrease
md_keys = mod_rsq([mod_rsq.run_tune_type] == 2);
[~, ind] = sort([md_keys.run_rsq]);
%plot(food.RunSpeedTuningFit(stripKey(md_keys(ind(end)),'data.ClusterInfo')));
plot(food.RunSpeedTuningFit(stripKey(md_keys,'data.ClusterInfo')));

% 3 = broadband
bd_keys = mod_rsq([mod_rsq.run_tune_type] == 3);
[~, ind] = sort([bd_keys.run_rsq]);
%plot(food.RunSpeedTuningFit(stripKey(bd_keys(ind(end)),'data.ClusterInfo')));
plot(food.RunSpeedTuningFit(stripKey(bd_keys,'data.ClusterInfo')));


%% add some stats
disp('STATISTICS')
disp('==========')
fprintf('# significantly modulated neurons: %d out of %d\n', nnz(modulated), nnz(not_modulated))
fprintf('# monotone increase: %d out of %d\n', numel(mi_keys), numel(mod_rsq))
fprintf('# monotone descrease: %d out of %d\n', numel(md_keys), numel(mod_rsq))
fprintf('# broad band: %d out of %d\n', numel(bd_keys), numel(mod_rsq))

%##########################################################################
% unitList: 815
% wrong sampling_rate -> 482 units
% Aspontaneous, repeated -> 398 units
% Aspontaneous, RunSpeedTuning -> 323 units
% Aspontaneous, RunSpeedTuningFit -> 302 units (low spike rate?)
%##########################################################################

%%

% plot responses at preferred speed
edges = [-inf 2 5 10 15 20 25 inf];
run_sp = sort([mod_rsq.run_sp],'descend');

n = histcounts(run_sp,edges);

edges(1) = 0;
edges(end) = [];

figure;
b = bar(edges,n,'k');
set(b, 'LineWidth',1.5);
set(gca,'TickDir','out','FontSize',13)
box off;

ylabel('# of neurons','FontSize',14)
xlabel('Preferred run speed (cm/s)','FontSize',14)
title('Correct x-labels in illustrator due to bins')


return

%% POPULATE

keys = x;

for i = 1:numel(keys)
    fprintf('... processing %d out of %d\n', i,numel(keys))
    try
        populate(food.RunSpeedTuningFit,keys(i))
    catch
        warning('Mouse_counter: %d, series_num: %d\n',keys.mouse_counter, keys.series_num)
    end
end

disp('... done')


% data.RunSpeedInfo
% data.RunSpeedTuning
% food.RunSpeedTuning
% food.RunSpeedTuningFit

% xx_spon = fetch(data.Experiments('exp_name LIKE "%spontaneous"') & stripKey(unitList,'data.ClusterInfo'));

%%

%################
% Suff to do
%################

% RF ANALYSIS
%ns6 files oder was man haben will
%alteste oben, neueste unten

% SPEED ANALYSIS
%nev files nachschauen data.Locomotion 69 try to reach ballmovement from nev files
%load nevfiled
%170 zeile



























































