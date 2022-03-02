%% rgc2lgn_model.m
% Models dLGN responses to the chirp stimulus as a combination of RGC mean
% cluster responses.
%
% Explanation from Roman_Roson_et_al_2018:
% [...], we combined the dLGN-p RGC dataset (Fig. 2) and the dLGN dataset (Fig. 3)
% to study how the dLGN responses are computed from the retinal output.
% We first accounted for the differences in recording methods by convolving 
% the dLGN spiking responses with the OGB-1 Ca2+ indicator kernel (Supp. Fig. 2).
% We then used a linear model – constrained to have non-negative weights – 
% to predict dLGN responses as a sum of weighted RGC inputs (Fig. 5a, Supp. Fig. 1).
% For prediction, we used the RGC-all cluster means 3 that were assigned at 
% least two dLGN-p cells. Each dLGN recording consisted of multiple trials 
% (stimulus repetitions, n = 10-30), which were divided into a training and a
% test set (50 prct / 50 prct). The model was evaluated using repeated random 
% sub-sampling cross-validation with 1,000 repetitions, where the weights were
% fitted on the training set in each repetition and prediction quality was 
% assessed on the test set. The reported weights represent mean values across the repeats.
% 
% NOTE: This file is an updated version of rgc2lgn_lin.m, to include the
% option of different models for comparison
%
% TODO: 
% - re-implement parfor loop after debugging

clear all; 
close all;

%% Get times
t1 = datetime('now', 'Format', 'dd-MMM-yyyy HH:mm:ss');
fprintf('START TIME: %s\n', t1);
%% Load data
% NOTE: Data sources: 
%  lgn_s70_f20.m: see analyses/Functions/makeChirpPSTH.m
%  g6rgc.m: created in Jupyter NBs, see analyses/gcamp6-ogb1/code/notebooks
%  o1rgc.m: created in Jupyter NBs, see analyses/gcamp6-ogb1/code/notebooks
%  kerns.m: created in Jupyter NBs, see analyses/gcamp6-ogb1/code/notebooks
load('../data/lgn_s70_f20.mat') % all lgn cells; filtered with sigma = 40 ms and output freq = 20 Hz, used as standard
load('../data/g6rgc.mat') % GCaMP6f RGC traces
load('../data/o1rgc.mat') % OGB-1 RGC traces
load('../data/kerns.mat') % GCaMP6f and OGB-1 convolution kernels


%% SET PARAMETERS

% Data processing pars
normtype   = 2;     % 1 = mean normalize; 2 = range normalize
minCellNum = 2;     % min # of dLGN-p RGCs assigned to a cluster required for including cluster into model; 0 = all clusters
indicator  = 3;     % RGC data set for modelling
                    % 1 = GCaMP6f; 2 = OGB1 3 = OGB1 restricted
                    % ogb1 restricted = OGB1 traces (cleaner than gcamp6
                    % due to higher cell nums), restricted to clusters
                    % found to be dLGN-projecting
nClust     = 49;    % # of RGC clusters

% Cross-validation pars
nCV = 200;      % repeated random trial sub-sampling cross-val number
cv_ratio  = 0.5;       % training vs. validation set ratio;

% Model pars
% Model type OPTIONS:
% - 'lin': linear
% - 'lin_nonneg': linear w non-neg constraint
% - 'lin_nonneg_exp2': linear with non-neg constraint + two-term exponential
% - 'lasso': linear model w elastic net regularization
% - 'glm': Generalized linear model with exponential non-linearity
% - 'lassoglm': Generalized linear model with exponential non-linearity & elastic net regularization
model_type = 'lasso';

% Initialize parallel cluster
% myPool = parpool(4); % parallel for-loop for 4-cored iMac

% Save pars
doSaveModel = true;  % save model + rgc data
save_dir    = '../data/';   % figure folder name
filename = sprintf('model_%s', model_type);
doSaveConv = false; % saves convolved dLGN traces 

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

% Get RGC distribution (cells per cluster)
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

% Downsample single dlgn traces
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

% Compute RGC mean
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
        lgn_conv_trials(iunit,1,j) = mean(lgn_conv_trials(iunit,5:7,j));
        lgn_conv_trials(iunit,2,j) = median(lgn_conv_trials(iunit,5:7,j));
        lgn_conv_trials(iunit,3,j) = mean(lgn_conv_trials(iunit,6:8,j));
        lgn_conv_trials(iunit,4,j) = median(lgn_conv_trials(iunit,6:8,j));
    end
end

if doSaveConv
    fullfile(save_dir,filename)
    save(fullfile(save_dir, 'lgn_f7_convolved.mat'), 'lgn_conv_trials');
    save(fullfile(save_dir, 'lgn_f7_mean_trials.mat'), 'lgn_mean_trials');
end

%% release some memory
clear lgn_mean_trials rgc_dist

%% Model LGN activity properties from RGCs activity properties

fprintf('Running model type: %s\n', model_type);

% Initialize vars
% NOTE: Unfortunately, the output is >2 GB w/o compression, so saving into
% one structure "model" w subfields is very slow to save and load.
model_type      = model_type; % model type
model_corr      = zeros(nLGN,nCV); % correlations
model_w         = zeros(nLGN,nRGC,nCV); % weights
model_rmse      = zeros(nLGN,nCV); % root mean squared error (RMSE)
model_y_train   = zeros(nLGN,size(rgc.chirp,1),nCV); % training set trials
model_y_valid   = zeros(nLGN,size(rgc.chirp,1),nCV); % test set trials
model_units     = zeros(nLGN,size(rgc.chirp,1),nCV); % all trials
model_y_hat     = zeros(nLGN,size(rgc.chirp,1),nCV); % model prediction
model_n_train   = cell(nLGN,1); % training set size
model_n_valid   = cell(nLGN,1); % test set size
model_clu_idx   = cluIdx; % indices of clusters used for model

% Construct design matrix
X = [ones(size(rgc_mean,1),1) rgc_mean];

% Model pars
switch model_type
    case 'lin'
        opts  = optimset('display','off','MaxIter',1000,'Algorithm','interior-point');          
    case 'lin_nonneg'                 
        opts  = optimset('display','off','MaxIter',1000,'Algorithm','interior-point');
        % Linear constraints
        lb = [-Inf zeros(1, nRGC)]; % lower bound: non-neg constraint (1st=-Inf, bec is offset)                                     
    case 'lin_nonneg_exp2'
        opts  = optimset('display','off','MaxIter',1000,'Algorithm','interior-point');
        % Linear constraints
        lb = [-Inf zeros(1, nRGC)]; % lower bound: non-neg constraint (1st=-Inf, bec is offset) 
        warning('off','curvefit:fit:invalidStartPoint')                    
    case 'lasso'        
    case 'glm'
        warning('off', 'stats:glmfit:IllConditioned');  
    case 'lassoglm'
        
end

% Loop through dLGN cells (n = 815)
for iunit = 1:1:nLGN    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COMPUTE MODEL %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fprintf('... processing unit %d out of %d\n', iunit, nLGN)
    tic
    parfor_progress(nCV); % display progress
    
    % reduce dimension
    unit = zeros(size(lgn_conv_trials(iunit,:,:),2),size(lgn_conv_trials(iunit,:,:),3));
    unit(:,:) = lgn_conv_trials(iunit,:,:);
    
    % remove empty trials
    unit(:,sum(unit)==0) = [];
    
    % Initialize temporary variables for parallel processes
    tmp_corr     = zeros(nCV,1);
    tmp_w        = zeros(nRGC,nCV);
    tmp_rmse     = zeros(nCV,1);
    tmp_y_train  = zeros(size(unit,1),nCV);
    tmp_y_valid  = zeros(size(unit,1),nCV);
    tmp_units    = zeros(size(unit,1),nCV);
    tmp_yhat     = zeros(size(unit,1),nCV);
    tmp_n_train   = zeros(length(1:floor(size(unit,2)*cv_ratio)),nCV);
    tmp_n_valid   = zeros(length(floor(size(unit,2)*cv_ratio)+1:size(unit,2)),nCV);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Repetead random sub-sampling cross-validation
    parfor iter = 1:nCV % avoid parfor for DEBUGGING
        
%         fprintf('.'); % Progress feedback
        parfor_progress; 
        
        % Split trials into randomly shuffle training and test set
        trials = randperm(size(unit,2));
        tmp_n_train(:,iter) = trials(1:floor(size(unit,2)*cv_ratio));
        tmp_n_valid(:,iter) = trials(floor(size(unit,2)*cv_ratio)+1:end);
        
        % Use trial means
        y_train = mean(unit(:, tmp_n_train(:,iter)),2);
        y_valid = mean(unit(:,tmp_n_valid(:,iter)),2);
        unit_mean = mean(unit,2); % all trials for later plotting of unit mean
        
       % Normalize
        if normtype == 1 % mean normalize
            y_train_norm = y_train - mean(y_train(1:8));
            y_train_norm = y_train_norm / max(abs(y_train_norm));
            
            y_valid_norm = y_valid - mean(y_valid(1:8));
            y_valid_norm = y_valid_norm / max(abs(y_valid_norm));
            
            unit_n = unit_mean - mean(unit_mean(1:8));
            unit_n = unit_n / max(abs(unit_n));
            
        elseif normtype == 2 % range normalize to [0, 1]
            y_train_norm = (y_train-min(y_train))/(max(y_train)-min(y_train));
            y_valid_norm = (y_valid-min(y_valid))/(max(y_valid)-min(y_valid));
            unit_n       = (unit_mean-min(unit_mean))/(max(unit_mean)-min(unit_mean));
        end
        
        %% Run model
        
        % Pick model type
        switch model_type
            case 'lin'
                try
                    % Linear model (vanilla type)
                    [w, resnorm, res] = lsqlin(X,y_train_norm,[],[],[],[],[],[],[],opts);

                    % Construct linear prediction
                    y_hat = X*w;

                catch
                    % Fit cannot be computed
                    warning('rgc2lgn_model: Fit cannot be computed for unit: %d; repeat: %d', iunit, iter)
                    continue
                end
            
            case 'lin_nonneg'                 
                try
                    % Linear model w non-neg constraint
                    [w, resnorm, res] = lsqlin(X,y_train_norm,[],[],[],[],lb,[],[],opts);

                    % Construct linear prediction
                    y_hat = X*w;

                catch
                    % Fit cannot be computed
                    warning('rgc2lgn_model: Fit cannot be computed for unit: %d; repeat: %d', iunit, iter)
                    continue
                end
                
            case 'lin_nonneg_exp2'
                try
                    % Linear-exponential model with non-neg constraint                              
                    [w, resnorm, res] = lsqlin(X,y_train_norm,[],[],[],[],lb,[],[],opts);

                    % Construct linear prediction
                    lin = X*w;

                    % Non-linear fitoptions (object needs to be created in
                    % here for parfor to work)
                    F = fitoptions('METHOD','NonlinearLeastSquares','Algorithm',...
                        'Trust-Region','Robust','LAR','Display','off','Lower',[0 0 0 0],...
                        'TolFun', 1e-04);

                    % Fit two-term exponential on top of linear model
                    f_exp2 = fit(lin, y_train_norm, 'exp2', F);

                    % Make prediction
                    y_hat = f_exp2.a * exp(lin * f_exp2.b) + f_exp2.c * exp(lin * f_exp2.d); % two-term exponential
                    
                catch                                    
                    % Fit cannot be computed
                    warning('rgc2lgn_model: Fit cannot be computed for unit: %d; repeat: %d', iunit, iter)
                    continue
                end

            case 'lasso'
                try
                    % Linear model w elastic net regularization
                    % Alpha = ratio of L1/L2-norm, 1=lasso, 0=ridge, inbetween=elastic net
                    [B, FitInfo] = lasso(X, y_train_norm, 'CV', 5, 'Alpha', 1, 'MaxIter', 1000, 'NumLambda', 50, 'RelTol', 1e-3);                    
                    w = B(:,FitInfo.Index1SE); % Use largest lambda value such that the MSE is within one SE of the minimum SE
                    B0 = FitInfo.Intercept(FitInfo.Index1SE);
                    y_hat = B0 + X*w;                    
                    
                catch
                    % Fit cannot be computed
                    warning('rgc2lgn_model: Fit cannot be computed for unit: %d; repeat: %d', iunit, iter)
                    continue
                end
            case 'glm'
                try
                    % Fit GLM with exponential non-linearity
                    % Modified from source: J. Pillow, GLM Tutorial (tutorial1_PoissonGLM.m) 
                    % https://github.com/pillowlab/GLMspiketraintutorial/blob/master/tutorial1_PoissonGLM.m
                    w = glmfit(X, y_train_norm, 'Poisson', 'link', 'log', 'constant', 'on');
                    B0 = w(1);
                    w = w(2:end);                    

                    % Compute predicted spike rate on training data
                    y_hat = exp(B0 + X*w);
                catch
                    % Fit cannot be computed
                    warning('rgc2lgn_model: Fit cannot be computed for unit: %d; repeat: %d', iunit, iter)
                    continue
                end
            case 'lassoglm'
                try
                    % Fit GLM with exponential non-linearity and elastic net regularization
                    % Poisson distribution has default link function 'log'
                    % Alpha = ratio of L1/L2-norm, 1=lasso, 0=ridge, inbetween=elastic net
                    [B, FitInfo] = lassoglm(X, y_train_norm, 'poisson', 'CV', 5, 'Alpha', 1, 'NumLambda', 50, 'RelTol', 1e-3, 'MaxIter', 1e3);
                    w = B(:,FitInfo.Index1SE); % Use largest lambda value such that the MSE is within one SE of the minimum SE
                    B0 = FitInfo.Intercept(FitInfo.Index1SE);
                    
                    y_hat = exp(B0 + X*w);
                catch
                    % Fit cannot be computed
                    warning('rgc2lgn_model: Fit cannot be computed for unit: %d; repeat: %d', iunit, iter)
                    continue                    
                end
%                 try
%                    
%                    % Options
%                    epsilon = 0.001;
%                    opts = struct('alpha', 0, 'lambda', epsilon);
%                    fit = glmnet(X, y_train_norm, 'poisson', opts);
%                    
%                    
%                    y_hat = glmPredict(fit);
%                    
%                    
%                 catch
%                     % Fit cannot be computed
%                     warning('rgc2lgn_model: Fit cannot be computed for unit: %d; repeat: %d', iunit, iter)
%                     continue
%                 end
                
        end

        % compute correlation between model and data & save vars
        tmp_corr(iter)      = corr(y_valid_norm,y_hat); % compute the corr on validation data
        tmp_w(:,iter)       = w(2:end);
        tmp_rmse(iter)      = sqrt(mean((y_valid_norm - y_hat).^2));       
        tmp_y_train(:,iter) = y_train_norm;
        tmp_y_valid(:,iter) = y_valid_norm;
        tmp_units(:,iter)   = unit_n;
        tmp_yhat(:,iter)    = y_hat;
        
%% DEBUGGING
%         [w2, resnorm2, res2] = lsqlin(X,y_train_norm,[],[],[],[],[],[],[],opts);
%         yhat2 = X*w2;
%         plot(unit_n);
%         hold on;
%         plot(yhat);
%         plot(yhat2);
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    % Assign variables
    model_corr(iunit,:)      = tmp_corr;
    model_w(iunit,:,:)       = tmp_w;
    model_rmse(iunit,:)      = tmp_rmse;
    model_y_train(iunit,:,:) = tmp_y_train;
    model_y_valid(iunit,:,:) = tmp_y_valid;    
    model_units(iunit,:,:)   = tmp_units;
    model_y_hat(iunit,:,:)   = tmp_yhat;
    model_n_train{iunit}     = {tmp_n_train}; % TODO: change format
    model_n_valid{iunit}     = {tmp_n_valid}; % TODO: change format
    
    fprintf('Elapsed time for unit: %s \n', datestr(datenum(0,0,0,0,0,toc), 'HH:MM:SS'))
    
    % Intermediate data save
%     save(fullfile(save_dir,sprintf([filename, '_1-73'])))
end

%% Note and save cells w no model
% TODO: think of better way
% [idx, ~] = find(isnan(model_corr(:,1)) == 1); % idx = 147;
% fprintf('Cell with no model: %s/n', mat2str(idx))
% 
% % Save into file
% try
%     load(fullfile(save_dir, 'no_model_cells.m'));   
% catch
%     no_model_cells = struct();    
% end
% switch model_type
%     case 'lin'
%         no_model_cells.lin = idx;
%     case 'lin_nonneg'
%         no_model_cells.lin_nonneg = idx;
%     case 'lin_nonneg_exp2'
%         no_model_cells.lin_nonneg_exp2 = idx;
%     case 'lasso'
%         no_model_cells.lasso = idx;
% end
% 
% save(fullfile(save_dir,'no_model_cells'), 'no_model_cells');

% model_units(idx,:,:)    = [];
% model_w(idx,:,:)        = [];
% model_y_hat(idx,:,:)    = [];
% model_corr(idx,:)       = [];
% model_n_train(idx,:)    = [];
% model_n_valid(idx,:)    = [];
% model_rmse(idx,:)       = [];
% model_y_train(idx,:,:)  = [];
% model_y_valid(idx,:,:)  = [];

%% Clear parpool
% delete(myPool)

%% Save
if doSaveModel
    %%%%%%%%%%%%%%%%%%%%%% DEBUG
%     filename = 'model_lin_nonneg_exp2_1-200';
    %%%%%%%%%%%%%%%%%%%%%%
    fprintf('Saving model to %s\n', fullfile(save_dir,filename));
    if(exist(save_dir, 'dir') ~= 7)
        mkdir(save_dir)
    end
    save(fullfile(save_dir,filename), '-regexp', 'model_', 'rgc')
end

%% Get times
t2 = datetime('now', 'Format', 'dd-MMM-yyyy HH:mm:ss');
dt = between(t1, t2);
fprintf('ELAPSED TOTAL TIME: %s\n', dt);

return;
%% DEBUGGING: concatenate model data sets
% list=['model_units','model_y']
% list=cellstr({'model_units','model_y'})
% save(fullfile(save_dir,'model_test'), list)
% save(fullfile('../data/','model_test'), '-regexp', 'model_', 'rgc')
% clear all
% load(fullfile('../data/','model_test'))
% load(fullfile('../data/', 'model_lin'))

filename = 'model_lin_nonneg_exp2_201-815'
load(fullfile(save_dir,filename), '-regexp', 'model_', 'rgc');

% Concatenate split datasets
model_corr2      = model_corr;
model_w2         = model_w    ;     
model_rmse2      = model_rmse  ;   
model_y_train2   = model_y_train;   
model_y_valid2   = model_y_valid ;  
model_units2     = model_units    ; 
model_y_hat2     = model_y_hat     ;
model_n_train2   = model_n_train   ;
model_n_valid2   = model_n_valid   ;
% model_clu_idx   = cluIdx; % indices of clusters used for model

filename = 'model_lin_nonneg_exp2_1-200'
load(fullfile(save_dir,filename), '-regexp', 'model_', 'rgc');

a=200;
b=201;

model_corr      = [model_corr(1:a,:); model_corr2(b:815,:)];
model_w         = [model_w(1:a,:,:); model_w2(b:815,:,:)];
model_rmse      = [model_rmse(1:a,:); model_rmse2(b:815,:)];
model_y_train   = [model_y_train(1:a,:,:); model_y_train2(b:815,:,:)];
model_y_valid   = [model_y_valid(1:a,:,:); model_y_valid2(b:815,:,:)];
model_units     = [model_units(1:a,:,:); model_units2(b:815,:,:)];
model_y_hat     = [model_y_hat(1:a,:,:); model_y_hat2(b:815,:,:)];
model_n_train   = [model_n_train(1:a,:); model_n_train2(b:815,:)];
model_n_valid   = [model_n_valid(1:a,:); model_n_valid2(b:815,:)];

% model_w=reshape(model_w,[815,33,1000]);
% model_y_train=reshape(model_y_train,[815,249,1000]);
% model_y_valid=reshape(model_y_valid,[815,249,1000]);
% model_units=reshape(model_units,[815,249,1000]);
% model_y_hat=reshape(model_y_hat,[815,249,1000]);

clear model_corr2 model_w2 model_rmse2 model_y_train2 model_y_valid model_units2 model_y_hat2 model_n_train2 model_n_valid2 model_y_valid2 model_type2;

filename = sprintf('model_%s', model_type);
fprintf('Saving model to %s\n', fullfile(save_dir,filename));
save(fullfile(save_dir,filename), '-regexp', 'model_', 'rgc');

