
%% Load data

load('workspace/lin_range_o1g6.mat')


%% SET PARAMETERS

saveFigs   = 1;                      % 1 = save figures; 0 figs won't be saved
fname      = 'figs/lin_range_o1g6';  % figure folder name

% mean # of RGC clusters
mnr = zeros(size(var_w,1),1);
for i = 1:size(var_w,1)
    mnr(i) = mean(sum(var_w(i,:,:) > 0.01,2));
end


%% PLOT SINGLE UNITS

h = 18;
[figPars, axPars] = aux_setPlotPars();

time  = rgc.chirp_time;

for iunit = 1:1:nLGN
    
    if mnr(iunit) == 0
        continue
    end
    
    fh = figure(figPars, 'Position', [14 7 22 h]);
    
    % NON-LINEAR MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    f1 = figure; hold on
    
    mvu = mean(var_units(iunit,:,:),3);
    mvy = mean(var_yhat(iunit,:,:),3);
    mvc = mean(var_corr_lin(iunit,:));
    
    minval = min([min(mvu) min(mvy)]);
    plot(time,mvu,'k','LineWidth',1) % data
    plot(time,mvy,'LineWidth',1)  % lin
    
    figHandles = findall(f1, 'Type', 'axes');
    newT1 = copyobj(figHandles(1), fh);
    set(newT1, axPars, 'Position', [2 h-8 17 5]);
    close(f1)
    
    set(newT1,'TickDir','out','Xlim',[0 32],'YLim',[minval 1],'FontSize',12);
    set(newT1,'YTick',[0 1], 'XTick', [0 2:5:32]);
    t = title(newT1,sprintf('LIN: %.2f',round(mvc*100)/100));
    t.FontSize = 14;
    xlabel(newT1,'Time (s)','FontSize', 15)
    lg = legend(newT1,'data','prediction');
    lg.EdgeColor = 'w';
    
    st1 = supertitle(sprintf('Cell %d; Qi: %.2f; ranksum: %.3f',iunit,lgn.qi(iunit),lgn.corr_p(iunit)),0.99);
    st1.FontSize = 14;
    st1.FontWeight = 'Bold';
    
    
    % WEIGHT DISTRIBUTION  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    f2 = figure; hold on
    
    sd = std(var_w(iunit,:,:),[],3);
    mw = mean(var_w(iunit,:,:),3);
    
    b1 = bar(mw); hold on
    set(b1, 'FaceColor', [0 0.45 0.75], 'LineWidth',1.5); % all data
    errorbar(mw,sd/2,'k','marker','none','LineStyle','none')
    
    figHandles = findall(f2, 'Type', 'axes');
    newT2 = copyobj(figHandles(1), fh);
    set(newT2, axPars, 'Position', [2 h-15 17 5]);
    close(f2)
    
    set(newT2,'FontSize',12,'TickDir','out','YLim',[0 1],'YTick',[0 1])
    set(newT2,'XTick',1:1:nRGC,'XTicklabel',cluIdx,'XTickLabelRotation',90)
    xlabel(newT2,'RGC Types','FontSize',15);
    ylabel(newT2,'Norm weight','FontSize',15)
    
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















































