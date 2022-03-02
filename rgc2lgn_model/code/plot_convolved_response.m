%% Plot convolved vs non-convolved dLGN responses

clear all
close all

%% Load data
% Load files saved in rgc2dlgn_model.m
load('../data/lgn_f7_mean_trials.mat', 'lgn_mean_trials'); % non-convolved traces
load('../data/lgn_f7_convolved.mat', 'lgn_conv_trials'); % convolved traces

% Load dLGN PSTHs saved in convolve_dlgn_traces.m
% save(fullfile(save_dir, 'lgn_psth.mat'), 'PSTH', 't');
load('../data/lgn_psth_s50_f1000.mat', 'PSTH_smooth', 't');
load('../data/lgn_psth_convolved.mat', 'PSTH_conv', 't');
load('../data/lgn_psth_convolved_resampled.mat', 'PSTH_conv_res', 't_interp');

%% Set parameters
doSave = false; % Save figure or not

% Choose example cells: cells used in make_fig_v2.m
cell_idx = [18, 133, 379, 450, 453, 456, 459, 460]; % [18, 133, 379-1, 810-1, 452];

%% Plot
% [fig, ax] = plot_conv_response(lgn_mean_trials, lgn_conv_trials, cell_idx); % function defined below
[fig, ax] = plot_conv_response_comparison(PSTH_smooth, t, lgn_conv_trials, t_interp, PSTH_conv_res, t_interp, cell_idx)


%% Save figure
if doSave
    fig_dir = '../figs/dlgn_convolved/';
    if (~exist(fig_dir, 'dir')); mkdir(fig_dir); end
    saveas(fig, fullfile(fig_dir, 'dlgn_convolved.eps'), 'epsc');
    saveas(fig, fullfile(fig_dir, 'dlgn_convolved.png'));
end

%%
function [fig, ax] = plot_conv_response(resp_orig, resp_conv, cell_idx)
    % function plot_conv_response(resp_orig, resp_conv, cell_idx)
    % Plots convolved (e.g. normalized DeltaF/F) against original cell response
    % (e.g. normalized firing rate) for selected example cells.
    %
    % Arguments:
    %   resp_orig : array [cell, time, trial] of original cell response time series
    %   resp_conv : array [cell, time, trial] of convolved cell response time series
    %   cell_idx :  array of cell indices for desired cells.
    %               DEFAULT: [1]
    % Returns:
    %   Plots time-series of the two trace types for specified example cells.
    %   fig, ax : optionally returns figure and axes handle
    % TODO: 
    % - include option to do single trial vs mean plotting
    % - insert hyphenated lines for chirp triggers?
    
    % Check for example cells to plot
    if nargin < 3 || isempty(cell_idx)
        cell_idx = [1]; % Default to first cell
    end
    
    % Call default figure parameters
    [figPars, axPars] = aux_setPlotPars();
    fig = figure(figPars, 'Position', [8.89,2.1519,30.6917,24.6239]); %, 'Position', [2 2 21 h]);


    % Loop through selected cells
    ax = gobjects(length(cell_idx)+1,1);
    for i = 1:length(cell_idx)
        
        % Plot chirp stimulus on top
        if i == 1
            ax(i) = subplot(length(cell_idx+1), 1, i);
            addpath(genpath('../../Functions/'));
            [t, y, cumTs] = plotChirpStim();
            plot(t,y, '-k', 'LineWidth', 1);
            set(ax(i), 'XLim', [0 cumTs(end)]);
            set(ax(i), axPars, 'YLim', [-1.5 1.5], 'YColor', 'w');
            axis off;
        end
                
        % Plot mean unconvolved response
        ax(i+1) = subplot(length(cell_idx)+1, 1, i+1);
        plot(mean(resp_orig(cell_idx(i),:,:),3), ...
            'Color', 'k', ...
            'LineWidth',1, ...
            'DisplayName','dLGN PSTH'); 
        hold on;
        % Plot mean convolved response
        plot(mean(resp_conv(cell_idx(i),:,:),3), ...
            'Color', [0.00,0.45,0.74], ...            
            'LineWidth',1, ...
            'DisplayName','dLGN PSTH convolved'); 
        
        % Set axis parameters
        set(ax(i+1), axPars);
        axis off;
        
        % Cell number text
        ylabel(ax(i+1),string(i)); 
        ylabh = get(ax(i+1),'ylabel');
        pos = ylabh.Position;
        text(pos(1), pos(2), pos(3), num2str(i), 'FontSize', axPars.FontSize);
        
        % Text and legend for second axis only (to avoid repetition)
        if i == 1 
            ylabel(ax(i+1),'A.U.');            
            text(pos(1)+8, pos(2)-2, pos(3), 'A.U.', 'Rotation', 90, ...
                'FontSize', axPars.FontSize);
            text(pos(1)-6, pos(2)+10, pos(3), 'cell #', ...
                'FontSize', axPars.FontSize);
            lg = legend(ax(i+1));
            legend('boxoff');
        end
%         set(ax,'TickDir','out','Xlim',[0 32],'YLim',[minval 1],'FontSize',afs);
%         set(ax,'YTick',[0 1], 'XTick', [0 2]);

        % Plot time scale line at the bottom
        if i == length(cell_idx)
            ax(i+1).XLim(1)
            line([ax(i+1).XLim(1),16], [-2,-2], 'Color', 'k', 'LineWidth', 3);
            text(4, -6, 0, '2 s',  'FontSize', axPars.FontSize);
        end
        
        % Plot legend
        if i == length(cell_idx)
            lg = legend('PSTH', 'PSTH convolved');
            lg.Box = 'off';
            lg.Orientation = 'horizontal';
%             lg.Location = 'southeast';
            lg.Position = lg.Position + [0 -0.03 0 0];
        end

    end
    
    % Add info to last plot axis
%     xlabel(ax,'RMSE','FontSize',lfs)
%     ylabel(ax,'# of cells (%)','FontSize',lfs)
%     title(newT2,'RMSE')
    % Add supertitle
    h = supertitle('dLGN cell electrophysiology responses (trial-mean) vs. Ca++-kernel-convolved responses (trial-mean)');
%     h.FontSize = axPars.FontSize
%     % plot mean or trials? >< do trials for now
%     figure;
%     plot(resp_unconv(ex_cell_idx(1),:,1));
%     hold on;
%     plot(resp_conv(ex_cell_idx(1),:,1));

    
end

%%
function [fig, ax] = plot_conv_response_comparison(resp_orig, t1, resp_conv, t2, resp_conv_new, t3, cell_idx)
    % function plot_conv_response(resp_orig, resp_conv, cell_idx)
    % Plots convolved (e.g. normalized DeltaF/F) against original cell response
    % (e.g. normalized firing rate) for selected example cells.
    %
    % Arguments:
    %   resp_orig : array [cell, time, trial] of original cell response time series
    %   resp_conv : array [cell, time, trial] of convolved cell response time series
    %   cell_idx :  array of cell indices for desired cells.
    %               DEFAULT: [1]
    % Returns:
    %   Plots time-series of the two trace types for specified example cells.
    %   fig, ax : optionally returns figure and axes handle
    % TODO: 
    % - include option to do single trial vs mean plotting
    % - insert hyphenated lines for chirp triggers?
    
    % Check for example cells to plot
    if nargin < 3 || isempty(cell_idx)
        cell_idx = [1]; % Default to first cell
    end
    
    % Call default figure parameters
    [figPars, axPars] = aux_setPlotPars();
    fig = figure(figPars, 'Position', [8.89,2.1519,30.6917,24.6239]); %, 'Position', [2 2 21 h]);


    % Loop through selected cells
    ax = gobjects(length(cell_idx)+1,1);
    for i = 1:length(cell_idx)
        
        % Plot chirp stimulus on top
        if i == 1
            ax(i) = subplot(length(cell_idx+1), 1, i);
            addpath(genpath('../../Functions/'));
            [t, y, cumTs] = plotChirpStim();
            plot(t,y, '-k', 'LineWidth', 1);
            set(ax(i), 'XLim', [0 cumTs(end)]);
            set(ax(i), axPars, 'YLim', [-1.5 1.5], 'YColor', 'w');
            axis off;
        end
                
        % Plot unconvolved response
        ax(i+1) = subplot(length(cell_idx)+1, 1, i+1);
        plot(t1, resp_orig(cell_idx(i),:), ...
            'Color', 'k', ...
            'LineWidth',1, ...
            'DisplayName','dLGN PSTH'); 
        hold on;
        
        %%% Plot comparison of convolved traces
        % Plot old mean convolved response        
        resp_conv_mean = mean(resp_conv(cell_idx(i),:,:),3);
        resp_conv_mean_norm = (resp_conv_mean - min(resp_conv_mean)) / (max(resp_conv_mean) - min(resp_conv_mean));
        plot(t2(1:249), resp_conv_mean_norm, ...
            'Color', [0.00,0.45,0.74], ...            
            'LineWidth',1, ...
            'DisplayName','dLGN PSTH convolved');
        
        % Plot new convolved response
        plot(t3(1:249), resp_conv_new(cell_idx(i),2:250), ...
            'Color', [1,0,0], ...            
            'LineWidth',1, ...
            'DisplayName','dLGN PSTH convolved');
        
        % Set axis parameters
        set(ax(i+1), axPars);
        axis off;
        
        % Cell number text
        ylabel(ax(i+1),string(i)); 
        ylabh = get(ax(i+1),'ylabel');
        pos = ylabh.Position;
        text(pos(1), pos(2), pos(3), num2str(i), 'FontSize', axPars.FontSize);
        
        % Text and legend for second axis only (to avoid repetition)
        if i == 1 
            ylabel(ax(i+1),'A.U.');            
            text(pos(1)+8, pos(2)-2, pos(3), 'A.U.', 'Rotation', 90, ...
                'FontSize', axPars.FontSize);
            text(pos(1)-6, pos(2)+10, pos(3), 'cell #', ...
                'FontSize', axPars.FontSize);
            lg = legend(ax(i+1));
            legend('boxoff');
        end
%         set(ax,'TickDir','out','Xlim',[0 32],'YLim',[minval 1],'FontSize',afs);
%         set(ax,'YTick',[0 1], 'XTick', [0 2]);

        % Plot time scale line at the bottom
        if i == length(cell_idx)
%             ax(i+1).XLim(1)
            line([ax(i+1).XLim(1),16], [-2,-2], 'Color', 'k', 'LineWidth', 3);
            text(4, -6, 0, '2 s',  'FontSize', axPars.FontSize);
        end
        
        % Plot legend
        if i == length(cell_idx)
            lg = legend('PSTH', 'PSTH convolved old', 'PSTH convolved new');
            lg.Box = 'off';
            lg.Orientation = 'horizontal';
%             lg.Location = 'southeast';
            lg.Position = lg.Position + [0 -0.03 0 0];
        end

    end
    
    % Add info to last plot axis
%     xlabel(ax,'RMSE','FontSize',lfs)
%     ylabel(ax,'# of cells (%)','FontSize',lfs)
%     title(newT2,'RMSE')
    % Add supertitle
    h = supertitle('dLGN cell electrophysiology responses (trial-mean) vs. Ca++-kernel-convolved responses (trial-mean)');
%     h.FontSize = axPars.FontSize
%     % plot mean or trials? >< do trials for now
%     figure;
%     plot(resp_unconv(ex_cell_idx(1),:,1));
%     hold on;
%     plot(resp_conv(ex_cell_idx(1),:,1));

    
end