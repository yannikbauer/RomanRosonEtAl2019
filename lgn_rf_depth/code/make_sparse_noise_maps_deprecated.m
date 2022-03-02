% make_sparse_noise_maps()
% Script calling Miro's code for generating RF maps from MUA data.
% NOTES: 
% - modified sparseNoiseExpByChan(key, RFtwin, zcutoff , file_ext, plotFlag, datDir, figDir, matDir, saveRaw)
%   to include datDir as an input argument
% - experiments not clean:
%   > for first series, some error w calculation of n_gratings per condition,
%    due to diff numbers in each experiment
%   > not all expts in list are populated


clear all
%% Project path
cd(fileparts(which(mfilename))); % changes working dir to current file
addpath(genpath('../')); % add entire project path, incl. general code used across analyses
addpath('../../code/m/brainpics')

%% Check DataJoint connection
if ~exist('DJ_DIRS')
    startup_lmu;
end

%% Parameter setup

RFtwin = [0 0.15];
zcutoff = 8;
file_ext = 'envl.dat';
plotFlag = true;
datDir = '../data/';%'/Volumes/lab/users/yannik/projects/ret2dlgn/dlgn_rf_depth/data/';
figDir = '../results/rf_maps/';
matDir = '../results/rf_maps/';
saveRaw = false;

%% Get experiment keys

% % Hard-coded test keys
key(1).mouse_counter = 132;
key(1).series_num = 12;
key(1).exp_num = 1;
key(2).mouse_counter = 132;
key(2).series_num = 12;
key(2).exp_num = 6;
% 
key(1).mouse_counter = 86;
key(1).series_num = 3;
key(1).exp_num = 4;


key(1).mouse_counter = 58;
key(1).series_num = 6;
key(1).exp_num = 1;

% Load key
load('../data/sn_exp.mat');
sn_exp = rmfield(sn_exp, {'mouse_id', 'series_date', 'fid'}); % rm unnec fields


%% Create sparse noise maps by series
% Logic: sparseNoiseExpByChan can take multiple experiments within the same
% series to create a combined RF map, so we generate a key by series

% Get all unique mice an loop through each
[mouse_counter, idx] = unique([sn_exp.mouse_counter].', 'rows');
for i = 1:length(mouse_counter)
        
    % Get mouse with all series as struct
    mouse = sn_exp([sn_exp.mouse_counter] == mouse_counter(i));    
    
    % Now get all unique series within that mouse and loop through those
    [series_num, idx] = unique([mouse.series_num].', 'rows');    
    for j = 1:length(series_num)
        
        % Get series key as input to sparseNoiseExpByChan
        key = mouse([mouse.series_num] == series_num(j));        
        
        fprintf('Doing mouse: %i, series: %i \n', [mouse_counter(i), series_num(j)]);
        
%         sparseNoiseExpByChan(key, RFtwin, zcutoff , file_ext, plotFlag, datDir, figDir, matDir, saveRaw);
        
        % DEBUGGING: some expt have exceptions
        % TODO: implement try-catch for that
        if i > 1
%             sparseNoiseExpByChan(key2, RFtwin, zcutoff , file_ext, plotFlag, datDir, figDir, matDir, saveRaw);
            sparseNoiseExpByChan(key, RFtwin, zcutoff , file_ext, plotFlag, figDir, matDir, saveRaw);
        end
        
    end

end











