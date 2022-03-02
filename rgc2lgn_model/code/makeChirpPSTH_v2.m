%% Make chirp PSTHs
% v2 is based on original function makeChirpPSTH.m located in
% analyses/clustering/code, and includes convolution with calcium kernels

% PSTH - resampled and smoothed
% PSTH_ori_samp - smoothed
% TS - PSTH timestamps
%
function [PSTH, outTimeStamps, PSTH_ori_samp] = makeChirpPSTH(unitList, sigma, dt, timeLims, outFreq)
% makeChirpPSTH computes PSTHs for a part of or the whole chirp stimulus
%
%   [PSTHMatrix, PSTH_ori_samp] = makeChirpPSTH(unitList, sigma, dt, onset, offset, outFreq)
%
%   unitList - cotains the list of keys
%
%   sigma - size of the gaussion kernel in ms
%       [default: 50]
%
%   dt - size of the resolution in sec
%       [default: 0.001]
%
%   timeLims - PSTHs from timeLims(1) to timeLims(2) in sec
%       [default [0 32.05]]
%
%   outFreq - sampling frequency of output smoothed signal
%       [default: 20]
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
    error('Insuficient inputs!')
end


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
    
    % Compute counts for every trial
    for i = 1:numel(spikeTimes)
        psth(i,:) = histc(spikeTimes{i},edges); % counts
        psth(i,:) = psth(i,:)*(1/dt); % spike rates (Hz)
    end
    
    
    %%%%%%%%%%%%
    % Convolve LGN cells
    % NOTE: The convolved vs. non-convolved responses can be plotted in './plot_convolved_response.m'

    psth_conv = zeros(size(psth));
    for itrial=1:nTrials
%         for j=1:nTrials
% 
%         if ~any(psth(itrial,:,j)) == 1
%             break
%         end

        % convolve trials
        tmp = conv(psth(itrial,:), kernel);

        % cut the trace, so it matches RGC
        psth_conv(itrial,:,j) = tmp(1:size(psth_conv,2)); 

        % correct for convolution artifact (first three and the last datapoint        
        psth_conv(itrial,1,j) = mean(psth_conv(itrial,5:7,j));
        psth_conv(itrial,2,j) = median(psth_conv(itrial,5:7,j));
        psth_conv(itrial,3,j) = mean(psth_conv(itrial,6:8,j));
        psth_conv(itrial,4,j) = median(psth_conv(itrial,6:8,j));
%         end
    end
    
    
    
    %%%%%%%%%%%%%%%%%
    
    % Average over trials and0.1:dt:maxdur spike rate (Hz)
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
    
    % Normalize PSTHMatrix
    PSTH(unit,:) = (PSTH(unit,:)-min(PSTH(unit,:)))/(max(PSTH(unit,:))-min(PSTH(unit,:)));
    PSTH_ori_samp(unit,:) = (PSTH_ori_samp(unit,:)-min(PSTH_ori_samp(unit,:)))/(max(PSTH_ori_samp(unit,:))-min(PSTH_ori_samp(unit,:)));
    
end % end of unit loop

PSTH = PSTH';
PSTH_ori_samp = PSTH_ori_samp';





end % end of function


















