%% rgc2lgn_lin.m
% Models dLGN responses to the chirp stimulus as a combination of RGC mean
% cluster responses. 
% 
% TODO: 
% update help text; 
% - reinstate defaults after debugging" rndNum=1000; do all cells

clear all;close all
%% Load data

load('../data/lgn_s70_f20.mat') % all lgn cells; filtered with sigma = 40 ms and output freq = 20 Hz, used as standard
load('../data/g6rgc.mat') % GCaMP6f RGC traces
load('../data/o1rgc.mat') % OGB-1 RGC traces
load('../data/kerns.mat') % GCaMP6f and OGB-1 convolution kernels


%% SET PARAMETERS

fname      = 'workspace';  % figure folder name
doSave     = false;      % saves the workspace

filename   = 'lin_range_o1g6';
normtype   = 2;      % 1 = mean normalize; 2 = range normalize
minCellNum = 2;      % min # of dLGN-p RGCs assigned to a cluster required for including cluster into model; 0 = all clusters
indicator  = 3;      % 1 = GCaMP6f; 2 = OGB1 3 = OGB1 restricted
nClust     = 49;     % # of RGC clusters

rndNum = 4;       % trial shuffling
ratio  = 0.5;        % training vs. validation set ratio;


%% Compute means & resample data

% Select indicator
if indicator == 1 || indicator == 3 % GCaMP6f
    rgc = g6rgc;
    kernel = g6Kern;
elseif indicator == 2  % OGB1
    rgc = o1rgc;
    kernel = o1Kern;
elseif indicator == 3  % OGB1 restricted
    rgc = g6rgc;
    kernel = o1Kern;
end

% get RGC distribution
rgc_dist = zeros(nClust,1);
for iunit=1:nClust
    rgc_dist(iunit) = nnz(rgc.cluIdx == iunit);
end

cluIdx  = 1:1:nClust;
cluIdx  = cluIdx(rgc_dist >= minCellNum);
nLGN    = size(lgn.chirp,2);
nRGC    = nnz(rgc_dist >= minCellNum);
tRGC    = length(rgc.chirp_time);
nTrials = size(lgn.chirp_trials,3);

if indicator == 3
    rgc = o1rgc;
end

% downsample single dlgn traces
lgn_mean_trials = zeros(nLGN,tRGC,nTrials);
for iunit=1:nLGN
    for j=1:nTrials
        trial = lgn.chirp_trials(iunit,:,j);
        if any(isnan(trial))
            break
        end
        lgn_mean_trials(iunit,:,j) = interp1(lgn.chirp_time',trial,rgc.chirp_time','linear');
    end
end

% compute RGC mean
rgc_mean = zeros(tRGC,nRGC);
for iunit=1:nRGC
    rgc_mean(:,iunit) = mean(rgc.chirp(:,rgc.cluIdx == cluIdx(iunit)),2);
end


%% Convolve LGN cells
% NOTE: The convolved vs. non-convolved responses can be plotted in './plot_convolved_response.m'

lgn_conv_trials = zeros(size(lgn_mean_trials));
for iunit=1:nLGN
    for j=1:nTrials
        
        if ~any(lgn_mean_trials(iunit,:,j)) == 1
            break
        end
        
        % convolve trials
        tmp = conv(lgn_mean_trials(iunit,:,j), kernel);
        
        % cut the trace, so it matches RGC
        lgn_conv_trials(iunit,:,j) = tmp(1:size(lgn_conv_trials,2)); 
        
        % correct for convolution artifact (first three and the last datapoint        
        lgn_conv_trials(iunit,1,j)   = mean(lgn_conv_trials(iunit,5:7,j));
        lgn_conv_trials(iunit,2,j)   = median(lgn_conv_trials(iunit,5:7,j));
        lgn_conv_trials(iunit,3,j)   = mean(lgn_conv_trials(iunit,6:8,j));
        lgn_conv_trials(iunit,4,j)   = median(lgn_conv_trials(iunit,6:8,j));
    end
end

save('../data/lgn_f7_convolved.mat', 'lgn_conv_trials');
save('../data/lgn_f7_mean_trials.mat', 'lgn_mean_trials');

%%
% release some memory
clear lgn_mean_trials rgc_dist

%% Get linear model of LGN activity properties from RGCs activity properties

% initialize vars
var_corr_lin  = zeros(nLGN,rndNum);
var_w         = zeros(nLGN,nRGC,rndNum);
var_rmse_lin  = zeros(nLGN,rndNum);
var_y_train   = zeros(nLGN,size(rgc.chirp,1),rndNum);
var_y_valid   = zeros(nLGN,size(rgc.chirp,1),rndNum);
var_units     = zeros(nLGN,size(rgc.chirp,1),rndNum);
var_yhat      = zeros(nLGN,size(rgc.chirp,1),rndNum);
var_nTrain    = cell(nLGN,1);                  
var_nValid    = cell(nLGN,1);                  

myPool = parpool(4); % parallel for-loop for 4-cored iMac
for iunit = 1:1:2%dLGN    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COMPUTE MODEL %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fprintf('.... processing unit %d out of %d\n', iunit, nLGN)
    
    % reduce dimension
    unit = zeros(size(lgn_conv_trials(iunit,:,:),2),size(lgn_conv_trials(iunit,:,:),3));
    unit(:,:) = lgn_conv_trials(iunit,:,:);
    
    % remove empty trials
    unit(:,sum(unit)==0) = [];
    
    % repeat due to trial randomization
    par_corr_lin = zeros(rndNum,1);
    par_w        = zeros(nRGC,rndNum);
    par_rmse_lin = zeros(rndNum,1);
    par_y_train  = zeros(size(unit,1),rndNum);
    par_y_valid  = zeros(size(unit,1),rndNum);
    par_units    = zeros(size(unit,1),rndNum);
    par_yhat     = zeros(size(unit,1),rndNum);
    par_nTrain   = zeros(length(1:floor(size(unit,2)*ratio)),rndNum);
    par_nValid   = zeros(length(floor(size(unit,2)*ratio)+1:size(unit,2)),rndNum);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    parfor iter = 1:rndNum % shuffle trials        
        
        trials = randperm(size(unit,2));
        par_nTrain(:,iter) = trials(1:floor(size(unit,2)*ratio));
        par_nValid(:,iter) = trials(floor(size(unit,2)*ratio)+1:end);
        
        y_train = mean(unit(:, par_nTrain(:,iter)),2);
        y_valid = mean(unit(:,par_nValid(:,iter)),2);
        unit_mean = mean(unit,2);
        
        if normtype == 1 % mean normalize
            y_train_norm = y_train - mean(y_train(1:8));
            y_train_norm = y_train_norm / max(abs(y_train_norm));
            
            y_valid_norm = y_valid - mean(y_valid(1:8));
            y_valid_norm = y_valid_norm / max(abs(y_valid_norm));
            
            unit_n = unit_mean - mean(unit_mean(1:8));
            unit_n = unit_n / max(abs(unit_n));
            
        elseif normtype == 2 % normalize to 1
            y_train_norm = (y_train-min(y_train))/(max(y_train)-min(y_train));
            y_valid_norm = (y_valid-min(y_valid))/(max(y_valid)-min(y_valid));
            unit_n       = (unit_mean-min(unit_mean))/(max(unit_mean)-min(unit_mean));
        end
        
        % model parameters
        opts1  = optimset('display','off','MaxIter',1000,'Algorithm','interior-point');
        
        % linear constraints
        lb = [-Inf zeros(1, nRGC)];
        
        % construct design matrix
        X = [ones(size(rgc_mean,1),1) rgc_mean];
        
        % predict weights for the training data
        try
            [w, resnorm, res] = lsqlin(X,y_train_norm,[],[],[],[],lb,[],[],opts1);
        catch
            continue
        end
        % Construct linear prediction
        lin = X*w;
        
        % compute correlation between model and data & save vars
        par_corr_lin(iter)  = corr(y_valid_norm,lin); % compute the corr on validation data
        par_w(:,iter)       = w(2:end);
        par_rmse_lin(iter)  = sqrt(mean((y_valid_norm - lin).^2));       
        par_y_train(:,iter) = y_train_norm;
        par_y_valid(:,iter) = y_valid_norm;
        par_units(:,iter)   = unit_n;
        par_yhat(:,iter)    = lin;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    % assign variables
    var_corr_lin(iunit,:) = par_corr_lin;
    var_w(iunit,:,:)      = par_w;
    var_rmse_lin(iunit,:)   = par_rmse_lin;
    var_y_train(iunit,:,:)  = par_y_train;
    var_y_valid(iunit,:,:)  = par_y_valid;    
    var_units(iunit,:,:)    = par_units;    
    var_yhat(iunit,:,:)     = par_yhat;
    var_nTrain{iunit}     = {par_nTrain};
    var_nValid{iunit}     = {par_nValid};
end

delete(myPool)

if doSave
    if(exist(fname, 'dir') ~= 7)
        mkdir(fname)
    end
    save(fullfile(fname,filename))
end




