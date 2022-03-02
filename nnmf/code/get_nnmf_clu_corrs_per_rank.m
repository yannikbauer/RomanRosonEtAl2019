%% Get NNMF-cluster max corrs
% Performs same NNMF analysis as main_dev.m but with a focus on getting the
% maxCorrs between NNMF and clusters as a function of rank for plotting

clear all

max_corrs = {}; % initialize array
for irank = 2:33
    close all;
    %% Parameters

    % Set data set typ
    data_type = 'crossval'; % OPTIONS: 'crossval', 'sparse'
    % Set rank par
    rank = irank;%4; % number of components = k; []=best k component-model (33 comps)

    % Directories
    addpath(genpath('../../Functions/'));
    data_dir = '../data';
    results_dir = '../results';
    
    % Save
    save_fig = true;

    %% Load data

    if strcmp(data_type, 'sparse')
        % sparse components
        load('data/nnmf_best.mat') % Miro's 25 sparse components
        A = A'; % A' = h | FxT (components x time)
        Y = Y'; % Y' = w | NxF (cells x components)

    elseif strcmp(data_type, 'crossval')
        % Load Miro's model to get dLGN PSTH data (not found elsewhere)
        load('data/nnmf_best.mat');
        clear A D_res k_best X Y;
        X = psth.psth';

        % Load cross-validated components
        if isempty(rank)
            fname = 'nnmf_cv_model.mat';
        elseif ~isempty(rank)
            fname = sprintf('nnmf_cv_model_%s.mat', num2str(rank));
        end
        load(fullfile(data_dir, fname)); % cross-validated components
        A = Vt;
        Y = U;
        k_best = size(U,2);
    end
    %
    %% COMPUTE FEATURES AND WEIGHTS
    disp('Computing features and weights.')

    % normalize the weights to range [0, 1]
    Y_norm = zeros(size(Y));
    for iunit = 1:size(Y,1)
        Y_norm(iunit,:) = (Y(iunit,:)-min(Y(iunit,:)))/(max(Y(iunit,:))-min(Y(iunit,:)));
    end

    %% Get NNMF components dendrogram

    d = pdist(A, 'euclidean');
    tree = linkage(d,'ward');
    leafOrder = optimalleaforder(tree,d, 'CRITERIA', 'adjacent');
    [~,~,OUTPERM] = dendrogram(tree, 0, 'Reorder', leafOrder, 'Orientation', 'left', ... 
                                'ColorThreshold','default');

    %% Correlate NNMF features with dLGN-p RGC type mean responses

    [max_corrs{irank}, fig] = corr_nnmf_rgc(OUTPERM, psth, A);

end

% Remove empty cells 
% max_corrs = max_corrs(find(~cellfun(@isempty, max_corrs)))

% Get mean and SD of max_corrs per rank
max_corr_mean = zeros(length(max_corrs),1);
max_corr_sd = zeros(length(max_corrs),1);
max_corr_se = zeros(length(max_corrs),1);
for i = 1:length(max_corrs)
    max_corr_mean(i) = mean(max_corrs{i});
    max_corr_sd(i) = std(max_corrs{i});    
    max_corr_se(i) = std(max_corrs{i}) / sqrt(length(max_corrs{i}));
end

%% Plot mean max corrs as a function of rank
close all;
figure;
errorbar(max_corr_mean, max_corr_se, 'linewidth', 1)
hold on
plot(max_corr_mean, 'k', 'linewidth', 2)
xlabel('Rank')
ylabel('mean(max(correlations))')
title('Average maxCorrs between NNMF and clusters as a function of rank')
legend('mean + SE')

% Save fig
%     save_fig = false;
if save_fig
    % save figure
    fname = ['nnmf_', data_type, '_mean_max_corrs'];        
    ffname = fullfile(results_dir, 'nnmf_num_components', fname);
    fprintf('Saving figure to %s ...\n', ffname);
    hgexport(gcf, ffname, hgexport('factorystyle'),'Format','eps');
    saveas(gcf, ffname, 'png');
%     close(gcf)
end
