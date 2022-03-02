function [PSTH,outTimeStamps,PSTH_ori_samp,unit_mask] = computePSTH(unitList, sigma, dt, timeLims, outFreq, norm)
% computePSTH(unitList, sigma, dt, timeLims, outFreq) computes PSTHs for 
% a part of or the whole chirp stimulus
%
%   unitList    = cotains the list of keys%
%   sigma       = size of the gaussion kernel in ms [default: 50]       
%   dt          = size of the resolution in sec [default: 0.001]
%   timeLims    = PSTHs from timeLims(1) to timeLims(2) in sec 
%                 [default [0 32.05]]
%   outFreq 	= sampling frequency of output smoothed signal 
%                 [default: 20]
%       
% [PSTH]                = computePSTH(...) returns resampled and smoothed 
%                         N-by-M matrix  where N is the number of variables 
%                         and M is the number of observations 
% [~, ~]                = computePSTH(...) returns time stamps
% [~, ~, PSTH_ori_samp] = computePSTH(...) returns smoothed N-by-M matrix 
%                         where N is the number of variables and M is the 
%                         number of observations 
%
if nargin < 5 || isempty(outFreq)
    outFreq = 20;
end
if nargin < 4 || isempty(timeLims) || size(timeLims,2) ~= 2
    timeLims = [0 32.05];
end
if nargin < 3 || isempty(dt)
    dt = 0.001;
end
if nargin < 2 && isempty(sigma)
    sigma = 50;
    dt = 0.001;
    timeLims = [0 32.05];
    outFreq = 20;
end
if nargin < 1
    error('computePSTH:Insuficient inputs!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% apply mask in the case of missing spikes
unit_mask = ones(length(unitList),1);

% Select units of interest
units = 1:length(unitList);

% PSTH
% Compute maximum bin edge
onset = timeLims(1)+dt/2; % dt/2 for taking the center of an edge and not the border of an edge
maxdur = timeLims(2);

edges = onset:dt:maxdur; 

% initialize
PSTH = nan(numel(units), ceil(maxdur/(1/outFreq) - onset/(1/outFreq)));
outTimeStamps = nan(1, ceil(maxdur/(1/outFreq) - onset/(1/outFreq)));
PSTH_ori_samp = nan(numel(units), length(edges));

for unit = 1:numel(units)
    % get all spike times in TrialSpikeExtra as cell array
    spikeTimes = fetchn(data.TrialSpikesExtra(unitList(unit)),'spike_times');
    
    % Convert data format
    spikeTimes = cellfun(@transpose,spikeTimes,'un',0);
    nTrials = size(spikeTimes,1);
    psth = zeros(nTrials,length(edges)); % contains the filtered response
    
    % at least 5 trials should have spikes
    if nnz(~cellfun('isempty',spikeTimes)) < 5 
        unit_mask(unit) = 0;
        continue
    end
    
    % Compute counts for every trial, averages over trials and0.1:dt:maxdur spike rate (Hz)    
    for i = 1:numel(spikeTimes)       
        psth(i,:) = histc(spikeTimes{i},edges); % counts
        psth(i,:) = psth(i,:)*(1/dt); % spike rates (Hz)
    end
    
    meanPSTH = mean(psth); % averages per bin across trials
    
    % Smooth and resample the data
    if ~isempty(sigma)       % 1280
       [PSTH(unit,:), outTimeStamps(1,:)] = SmoothAndResampleC(meanPSTH, edges, [], 'sec', [], sigma, outFreq);  % Smooth and resample
       [PSTH_ori_samp(unit,:), ~] = SmoothAndResampleC(meanPSTH, edges, [], 'sec', [], sigma);        % Smooth
    else
        % not smoothed but resampled
        outTimeStamps = 1/outFreq : 1/outFreq : maxdur;
        PSTH(unit,:) = interp1((1:size(psth,2))*dt,meanPSTH,outTimeStamps);
        
        PSTH_ori_samp(unit,:) = meanPSTH;
    end
    
    % Normalize PSTHMatrix based on the whole length of the stimulus
    if norm == 1
        PSTH(unit,:) = (PSTH(unit,:)-min(PSTH(unit,:)))/(max(PSTH(unit,:))-min(PSTH(unit,:)));
        PSTH_ori_samp(unit,:) = (PSTH_ori_samp(unit,:)-min(PSTH_ori_samp(unit,:)))/(max(PSTH_ori_samp(unit,:))-min(PSTH_ori_samp(unit,:)));
    end
    
end % end of unit loop

% apply mask
PSTH = PSTH(unit_mask == 1,:);
PSTH_ori_samp = PSTH_ori_samp(unit_mask == 1,:);

PSTH = PSTH';
PSTH_ori_samp = PSTH_ori_samp';
end % end of function


















