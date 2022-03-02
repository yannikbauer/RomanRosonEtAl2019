


load('/Users/miroslav/Google Drive/BBE_BusseBerensEuler/code/modelling/workspace/weights.mat');
w_mean   = mean(var_w,3);


%% strongest weight / sum of remaining weights

w_max    = zeros(size(w_mean,1),1);
w_sum    = zeros(size(w_mean,1),1);
w_ratio  = zeros(size(w_mean,1),1);

for i = 1:size(w_mean,1)
    tmp = w_mean(i,:);
    w_max(i) = max(tmp);
    w_sum(i) = sum(tmp) - w_max(i);
    w_ratio(i) = w_max(i)/w_sum(i);
end

% remove NaNs
w_ratio(isnan(w_ratio)) = [];

% plot some stuff %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

edges = [0:0.1:2.5 inf];
n_all = histcounts(w_ratio,edges);
n_all = n_all*100/length(w_ratio);
edges(1) = [];
edges(end) = 2.6;

minval = min(w_ratio);
maxval = max(w_ratio);

% plot histo
figure; hold on
b3 = bar(edges,n_all,'histc');
set(b3, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', [1 1 1], 'LineWidth',1.5);
set(gca,'Xtick',minval:2:maxval,'XLim',[minval-1 maxval+1])

% plot mean
x = zeros(100,1);
x(:,1) = median(w_ratio);
plot(x,linspace(0,ceil(max(n_all)),100)','--k','LineWidth',1.5)

% display
disp('Convergence')
fprintf('mean = %.2f\n',mean(w_ratio));
fprintf('meadian = %.2f\n',median(w_ratio));

set(gca,'TickDir','out','FontSize',11,'Xlim',[0 3],'XTick',edges(0.4:5:end))
xlabel(gca,'Weight ratio','FontSize',14)
ylabel(gca,'Fraction of dLGN cells (%)','FontSize',14)
title(gca, 'Strong vs weak weight ratio')


%% Use threshold

thr = 0.1:0.1:1.4;
w_hist = zeros(size(w_mean,1),length(thr));

for i = 1:length(thr)
    for j = 1:size(w_mean,1)
        w_hist(j,i) = nnz(w_mean(j,:) > thr(i));
    end
end

w_max = max(max(w_hist));
edges = 0:1:w_max;
w_plot = zeros(size(w_hist,2),length(edges));
for i = 1:size(w_hist,2)
    w_plot(i,:) = histc(w_hist(:,i),edges)*100/size(w_mean,1);    
end

figure;imagesc(w_plot)
set(gca,'YtickLabel',thr(2:2:end))
ylabel('threshold','FontSize',14);
xlabel('# of weights','FontSize',14);
colorbar
title('Fraction of neurons with weights > threshold')



































