% clear

%% load data
load('../data/weights.mat')

%% set dimensions
[nneurons, ntypes, nrepeats] = size(var_w);

thresholds = 0:0.01:1;
%% explore 1st neuron
norm_w = var_w;
for ineuron = 1 : nneurons
    for irepeat = 1 : nrepeats
%         this_max = max(squeeze(var_w(ineuron, :,irepeat))); % normalize to max weight
        this_sum = sum(squeeze(var_w(ineuron, :,irepeat))); % normalize to sum of weights
%         norm_w(ineuron,:,irepeat) = var_w(ineuron, :,irepeat)/this_max;
        norm_w(ineuron,:,irepeat) = var_w(ineuron, :,irepeat)/this_sum;
    end
end

%% plot figure of normalized weights
figure; imagesc(mean(norm_w,3), [0 1])
colorbar

%% go through the thresholds
edges = 0:ntypes+1;
edges = edges + 0.5;
counts = nan(numel(edges)-1, numel(thresholds));
for ithreshold = 1 : numel(thresholds)
    threshold = thresholds(ithreshold);
    % sum: sum over the RGC types for counting, mean over repeats
    counts(:, ithreshold) = histcounts(mean(sum(norm_w > threshold,2),3), edges);
end

%% plot figure with number of clusters as a function of threshold
figure;
imagesc(1:ntypes, thresholds, counts');
ax = gca();
ylabel('Weight threshold')
xlabel('Number of RGC types')
axis square
box off
colormap(flipud(gray));
[figPar, axPar] = aux_setPlotPars();
set(ax,axPar);
set(ax, 'YDir', 'normal');

cbar = colorbar;
title(cbar,'# cells')
cbar.FontSize = axPar.FontSize;
cbar.TickDirection = 'out';

% Save figures
% print(gcf, '-depsc', 'thresholdedWeights.eps')
% print(gcf, '-dpdf', 'thresholdedWeights.pdf')

%% explore a bit more
figure; plot(mean(sum(norm_w > 0.9,2),3))
figure; plot(mean(sum(norm_w > 0.4,2),3))
figure; plot(mean(sum(norm_w > 0.6,2),3))
figure; plot(mean(sum(norm_w > 0.1,2),3))

%% correlation / consistency of weights across repeats
ineuron = 10;
figure; imagesc(squeeze(var_w(ineuron,:,:)));
