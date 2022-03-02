%% Make summary figure for RGC-dLGN model
% Includes: 
% (1) example cell responses against model prediction + weights bar plot
% (2) Population correlation histogram
% (3) Population RMSE histogram
% (4) RGC-dLGN convergence histograms for weak and strong threshold
% (5) Histograms for number of RGCs used for reconstruction as a function of (continuous) threshold
% (6) Bar plots for population mean weight / percentage per cluster/group 
% TODO: 
% - improve dealing with cells w no model


clear all;
close all;
%% Setup parameters

% Model pars
% Model type OPTIONS:
% - 'lin': linear
% - 'lin_nonneg': linear w non-neg constraint
% - 'lin_nonneg_exp2': linear with non-neg constraint + two-term exponential
% - 'lasso': linear model w lasso regularization
% - 'elnet': linear model w elastic net (50/50) regularization
% - 'glm': Generalized linear model with exponential non-linearity
% - 'lassoglm': Generalized linear model with exponential non-linearity & lasso regularization
% - 'elnetglm': Generalized linear model with exponential non-linearity & elastic net (50/50) regularization
model_type = 'lin_nonneg';

% Choose whether to get subpopulation of dLGN cells (top or bottom of LGN)
get_subpop_ = 'top'; % OPTIONS: 'top' 'bottom', []

% Directory pars
data_dir = '../data/';   % figure folder name
filename = sprintf('model_%s', model_type);

% Load data
if ~exist('model_clu_idx')
    fprintf('Loading data ...\n');
    load(fullfile(data_dir, filename)); % 6.5 GB, former file name: lin_range_o1g6
end

% Plot pars
load(fullfile(data_dir, 'cMaps/cMap_igor.mat'))
[figPars, axPars] = aux_setPlotPars();
axPars.LabelFontSizeMultiplier = 1.2;

% Auxilliary axis pars that cannot be passed to set()
axParsAux.AxisFontSize = 11; 
axParsAux.LabelFontSize = 14;
axParsAux.barcolor = [0.75 0.75 0.75];
axParsAux.cMap_igor = cMap_igor;
color_ind = round(linspace(1,size(axParsAux.cMap_igor,1), max(model_clu_idx)));

% Get RGC cluster and group names
[rgc.names_clu, rgc.names_group] = get_rgc_names();

% Set model weight threshold (restrict number of weights going into model)
weight_threshold = 0.001; % 0.001:weak; 0.2:strong
weight_threshold_strong = 0.2;

% Plot clusters or groups
clu_vs_group = 'group'; % plot clusters

% Save pars
save_figs = true;

fprintf('Plotting results for model type: %s\n', model_type);

%% Get subpopulation of dLGN cells
if ~isempty(get_subpop_)
    [model_units, model_w, model_y_hat, model_corr, model_n_train, model_n_valid, model_rmse, model_y_train] = ...
        get_subpop(get_subpop_, model_units, model_w, model_y_hat, model_corr,model_n_train, model_n_valid, model_rmse, model_y_train);
end

%% Deal with cells w no model
% TODO: find better way - for now done after single cell plotting, assuming
% that modelling script will not exclude cells, but only save idx of
% missing ones ~ current solution is unsatisfying mix of methods.
[idx, ~] = find(isnan(model_corr(:,1)) == 1 | model_corr(:,1) == 0); % idx = 147;
fprintf('Cell(s) with no model: %s\n', mat2str(idx));

model_units(idx,:,:)    = [];
model_w(idx,:,:)        = [];
model_y_hat(idx,:,:)    = [];
model_corr(idx,:)       = [];
model_n_train(idx,:)    = [];
model_n_valid(idx,:)    = [];
model_rmse(idx,:)       = [];
model_y_train(idx,:,:)  = [];
% model_y_valid(idx,:,:)  = [];

%% Normalize weights

[nneurons, ntypes, nrepeats] = size(model_w);
model_w_norm = model_w; % DEBUG TEST
for ineuron = 1 : nneurons
    for irepeat = 1 : nrepeats
%         this_max = max(squeeze(model_w(ineuron, :,irepeat))); % normalize to max weight
        this_sum = sum(abs(squeeze(model_w(ineuron, :, irepeat)))); % normalize to sum of weights
%         norm_w(ineuron,:,irepeat) = model_w(ineuron, :,irepeat)/this_max;
        model_w_norm(ineuron,:,irepeat) = model_w(ineuron, :,irepeat) / this_sum;
    end
end

%% Compute means across cross-validation repeats
model_units_mean = mean(model_units,3);
model_w_norm_mean = mean(model_w_norm,3);
model_w_norm_sd = std(model_w_norm,[],3);
model_units_y_hat_mean = mean(model_y_hat,3);
model_corr_mean = mean(model_corr,2);
model_rmse_mean = mean(model_rmse,2);

% model_corr_mean = mean(model_corr,2); % redundant
% mean_mean_w = mean(mean(model_w_norm,3)); % unused
% mean_sd = mean(std(model_w_norm,[],3)); % unused

%% Plotting

% Summary fig containing other plotting function outputs as subplots
margin = 3; % paper margin
h = 29.7 - margin; % A4 paper height - margin
w = 21 - margin; % A4 paper width - margin
fig_all = figure(figPars, 'Position', [17 5 w h]);
fig_all.Position = [17 5 w*1.5 h*1.5];

%% Plot example unit vs model response time series + model weights
% TODO: re-write plot function and cell selection code

% Set example cells
man_vs_auto = 'man'; % manual or automatic selection
if ~isempty(get_subpop_), man_vs_auto = 'auto', end; % force auto-selection of only subpop selected
if strcmp(man_vs_auto, 'auto') % automatic (by correlation)
    n_cells = 4; % how many cells to plot
    [~,cell_idx] = sort(model_corr,'descend');
    cell_idx = cell_idx(1:n_cells);
elseif strcmp(man_vs_auto, 'man') % manual
    cell_idx = [153, 280, 452, 133];%[154, 281, 453, 133]; % [18, 133, 378, 809, 37, 153, 280, 132]
end

% Plot
[fig_units, lg_units] = plot_units_model(cell_idx, rgc, model_clu_idx, model_units_mean, model_units_y_hat_mean,...
    model_corr_mean, model_w_norm_mean, model_w_norm_sd, weight_threshold,...
    model_type, figPars, axPars, axParsAux, h, color_ind);

% Copy subplot into summary figure
figHandles = findall([fig_units], 'Type', 'axes');
copy = copyobj([figHandles], fig_all);
copy2 = copyobj([figHandles(end), lg_units], fig_all); % workaround to also copy legend w assoc axis

% Some silly workaround to reposition plot
for i = 1:length(copy)
    copy(i).Position = copy(i).Position + [0, 14, 0, 0];
end
copy2(1).Position = copy2(1).Position + [0, 14, 0, 0];

%% Plot correlation histogram
% TODO: implement for varying thresholds

[fig_corr, ax_corr] = plot_corr_hist(model_corr_mean, axPars, axParsAux);

% Copy subplot into summary figure and adjust
copy = copyobj([ax_corr], fig_all);
set(copy(1), axPars, 'Position', [17.5017   36.1692    7.0555    3.6336]); %[2 h-pos(i) 9 2]);
copy.Children(1).Position(2) = copy(1).YLim(2);

%% Plot RMSE histogram

% [fig_rmse, ax_rmse, lg_rmse] = plot_rmse_hist(model_rmse_mean, axPars, axParsAux);
[fig_rmse, ax_rmse] = plot_rmse_hist(model_rmse_mean, axPars, axParsAux);

% Copy subplot into summary figure and adjust
copy = copyobj([ax_rmse], fig_all); % Copy subplot into summary figure
set(copy(1), axPars, 'Position', [17.4978   31.5736    7.0203    3.0339]);
copy.Children(1).Position(2) = copy(1).YLim(2);

%% Plot # of cells for reconstruction = RGC>dLGN convergece
% Do for standard weight threshold and strict weight threshold

% standard weight threshold
[fig_conv, ax_conv] = plot_convergence_hist(model_w_norm, weight_threshold, axPars, axParsAux);

copy = copyobj([ax_conv], fig_all); % Copy subplot into summary figure
set(copy(1), axPars, 'Position', [16.0514   26.1761    3.9511    2.9281], 'XLabel', []);
copy.Children(1).Position(2) = copy(1).YLim(2);
copy.Children(1).Position(1) = copy(1).XLim(1);

% strict weight threshold
[fig_conv2, ax_conv2] = plot_convergence_hist(model_w_norm, weight_threshold_strong, axPars, axParsAux);

copy = copyobj([ax_conv2], fig_all); % Copy subplot into summary figure
set(copy(1), axPars, 'Position', [21.2019   26.1761    3.7394    2.8928], 'XLabel', [], 'YLabel', []);
copy.Children(1).Position(2) = copy(1).YLim(2);
copy.Children(1).Position(1) = copy(1).XLim(1);

%% Plot number of weights as a function of threshold

[fig_thresh, ax_thresh, cbar_thresh] = plot_weights_per_threshold(model_w, model_w_norm, axPars);

copy = copyobj([ax_thresh, cbar_thresh], fig_all); % Copy subplot into summary figure
% Need to reset colormap for some reason
% colormap(copy(1), flipud(gray));
cmap = [1 1 1 ; flipud(bone(cbar_thresh.Limits(2)))]; % add ones to represent NaNs as white
cmap = [cmap(1:2, :); cmap(0.1*cbar_thresh.Limits(2):end,:)]; % clip cmap to show values close to white as more grey
colormap(copy(1), cmap);
set(copy(1), axPars, 'Position', [15.5   18.3   8    6.5]);
copy(1).XLabel.FontSize = 15;

%% Plot model weight distribution
% TODO: clean up functions

% Show if plotting clusters or groups
fprintf('Plotting %s model weight distribution:\n', clu_vs_group);

%%% Plot soft threshold
% Plot mean weights
mean_vs_percent = 'mean'; % plot mean weight
[fig_weight_mean, ax_weight] = plot_weight_distribution(model_w_norm, ...
    model_w_norm_mean, model_corr_mean, weight_threshold, rgc, model_clu_idx, clu_vs_group, mean_vs_percent, axPars, axParsAux);

copy = copyobj([ax_weight], fig_all); % Copy subplot into summary figure
set(copy(1), axPars, 'Position', [1.9050   13.8994   24.0242    2.3636], 'XLabel', [], 'XTickLabels', []);
copy(1).YLabel.FontSize = axParsAux.LabelFontSize;

% Plot percent
mean_vs_percent = 'percent'; % plot mean weight
[fig_weight_prct, ax_weight_prct] = plot_weight_distribution(model_w_norm, ...
    model_w_norm_mean, model_corr_mean, weight_threshold, rgc, model_clu_idx, clu_vs_group, mean_vs_percent, axPars, axParsAux);

copy = copyobj([ax_weight_prct], fig_all); % Copy subplot into summary figure
set(copy(1), axPars, 'Position', [1.8736   11.2631   23.9889    1.7992], 'XLabel', [], 'XTickLabels', 1:1:49, 'Title', []);
copy(1).YLabel.FontSize = axParsAux.LabelFontSize;

%%% Plot strong threshold
% Plot mean weights
mean_vs_percent = 'mean'; % plot mean weight
[fig_weight_mean, ax_weight] = plot_weight_distribution(model_w_norm, ...
    model_w_norm_mean, model_corr_mean, weight_threshold_strong, rgc, model_clu_idx, clu_vs_group, mean_vs_percent, axPars, axParsAux);

copy = copyobj([ax_weight], fig_all); % Copy subplot into summary figure
set(copy(1), axPars, 'Position', [1.8736    7.8058   24.0242    1.7286], 'XLabel', [], 'XTickLabels', []);
copy(1).YLabel.FontSize = axParsAux.LabelFontSize;

% Plot percent
mean_vs_percent = 'percent'; % plot mean weight
[fig_weight_prct, ax_weight_prct] = plot_weight_distribution(model_w_norm, ...
    model_w_norm_mean, model_corr_mean, weight_threshold_strong, rgc, model_clu_idx, clu_vs_group, mean_vs_percent, axPars, axParsAux);

copy = copyobj([ax_weight_prct], fig_all); % Copy subplot into summary figure
set(copy(1), axPars, 'Position', [1.8736    5.0189   23.9183    1.7639], 'XTickLabelRotation', 45, 'Title', []);
copy(1).YLabel.FontSize = axParsAux.LabelFontSize;


%% Save figures

% Save example unit fig
if save_figs
    fig_dir = '../results/example_cells/';
    fname = sprintf('example_cells_%s_%s', model_type, get_subpop_);
    if (~exist(fig_dir, 'dir')); mkdir(fig_dir); end
    saveas(fig_units, fullfile(fig_dir, [fname, '.eps']), 'epsc');
    saveas(fig_units, fullfile(fig_dir, [fname, '.png']));
    fprintf('Saving %s\n', fullfile(fig_dir, [fname, '.eps']))
    fprintf('Saving %s\n', fullfile(fig_dir, [fname, '.png']))
end

% Save correlation fig
if save_figs
    fig_dir = '../results/mean_corr/';
    fname = sprintf('mean_corr_%s_%s', model_type, get_subpop_);
    if (~exist(fig_dir, 'dir')); mkdir(fig_dir); end
    saveas(fig_corr, fullfile(fig_dir, [fname, '.eps']), 'epsc');
    saveas(fig_corr, fullfile(fig_dir, [fname, '.png']));    
end

% Save RMSE fig
if save_figs
    fig_dir = '../results/mean_rmse/';
    fname = sprintf('mean_rmse_%s_%s', model_type, get_subpop_);
    if (~exist(fig_dir, 'dir')); mkdir(fig_dir); end
    saveas(fig_rmse, fullfile(fig_dir, [fname, '.eps']), 'epsc');
    saveas(fig_rmse, fullfile(fig_dir, [fname, '.png']));    
end

% Save convergence fig
if save_figs
    fig_dir = '../results/num_types_histograms/';
    fname = sprintf('num_cells_%s_%s_thresh-%s', model_type, get_subpop_, num2str(weight_threshold));
    if (~exist(fig_dir, 'dir')); mkdir(fig_dir); end     
    saveas(fig_conv, fullfile(fig_dir, [fname, '.eps']), 'epsc');
    saveas(fig_conv, fullfile(fig_dir, [fname, '.png']));
    fprintf('Saving %s\n', fullfile(fig_dir, [fname, '.eps']))
    fprintf('Saving %s\n', fullfile(fig_dir, [fname, '.png']))
end

% Save weight threshold fig
if save_figs
    fig_dir = '../results/weights_per_threshold/';
    fname = sprintf('weights_per_threshold_%s_%s', model_type, get_subpop_);
    if (~exist(fig_dir, 'dir')); mkdir(fig_dir); end     
    saveas(fig_thresh, fullfile(fig_dir, [fname, '.eps']), 'epsc');
    saveas(fig_thresh, fullfile(fig_dir, [fname, '.png']));
    fprintf('Saving %s\n', fullfile(fig_dir, [fname, '.eps']))
    fprintf('Saving %s\n', fullfile(fig_dir, [fname, '.png']))
end
 
% Save summary fig
if save_figs
    fig_dir = '../results/summary_fig/';
    fname = sprintf('summary_fig_%s', model_type);
    
    fname = sprintf('summary_fig_%s_%s', model_type, get_subpop_);

    if (~exist(fig_dir, 'dir')); mkdir(fig_dir); end
    saveas(fig_all, fullfile(fig_dir, [fname, '.eps']), 'epsc');
    saveas(fig_all, fullfile(fig_dir, [fname, '.png']));
    fprintf('Saving %s\n', fullfile(fig_dir, [fname, '.eps']))
    fprintf('Saving %s\n', fullfile(fig_dir, [fname, '.png']))
end


% % Save fig
% if save_figs
%     fig_dir = '../results/weight_type_distribution/';
%     if (~exist(fig_dir, 'dir')); mkdir(fig_dir); end
%     fig_name = sprintf('weights_%s_%s_%s_thresh-%s', model_type, clu_vs_group, mean_vs_percent, num2str(weight_threshold));
%     saveas(fig_weight_mean, fullfile(fig_dir, [fig_name, '.eps']), 'epsc');
%     saveas(fig_weight_mean, fullfile(fig_dir, [fig_name, '.png']));    
% end


return
%% Plot cluster convergence

% Plot convergence matrix
[fig, ax, convMat_norm] = plot_conv_matrix(model_w_norm, model_w_norm_mean, weight_threshold, model_clu_idx, axPars, axParsAux);

% Plot convergence pairs per cluster
[fig, ax] = plot_conv_pairs(convMat_norm, model_clu_idx, maxRGC, axPars, axParsAux)
