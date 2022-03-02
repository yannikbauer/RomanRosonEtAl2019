%% Creates convolved dLGN traces

clear all;
close all;

%% Get dLGN unit list
% GET SAME KEYS USED FOR DLGN UNITS IN OTHER ANALYSES OF ROMAN_ROSON_ET_AL_2018
corr_p      = 0.001;         % 0.001; corr_p represents the wilcoxon ransum correlation for between/withing segments
qi          = 0.05;          % 0.05; Philipps quality index
cluster_qi  = 3;             % 3; cluster quality index of likely single- vs multiunit
r           = 0.65;          % 0.7; regression index
type        = 3;             % 1 = chirp units; 2 = SFTFORI units; 3 = both

unitList = getAllUnits(corr_p, qi, cluster_qi, r, type);

%% Make PSTH

% sigma = 5;
% outFreq = 160;
dt = 0.001;
timeLims = [0 32.05];
outFreq = 1000;

t = 1/outFreq : 1/outFreq : timeLims(2);

disp('Computing PSTH')
[PSTH, ts, ~] = makeChirpPSTH(unitList, [], dt, timeLims, outFreq);
PSTH = PSTH';
% function [PSTH, outTimeStamps, PSTH_ori_samp] = makeChirpPSTH(unitList, sigma, dt, timeLims, outFreq)

%% Convolve PSTH and full-resolution PSTH

% Load Calcium kernels
load('../data/kerns.mat') % GCaMP6f and OGB-1 convolution kernels
kernel.trace = o1Kern; % Use OGB-1 kernel
kernel.fps = 7.8251462420517992; % Kernel sampling rate
kernel.dt = 1/kernel.fps; % kernel time steps
kernel.dur = kernel.dt*(size(kernel.trace,2)-1); % kernel duration
kernel.t = 0:kernel.dt:kernel.dur; % Kernel time points (ms)

% Resample kernel to sampling rate of PSTH
t_interp = 0 : 1/outFreq : kernel.dur;
kernel.trace = interp1(kernel.t, kernel.trace, t_interp);

% Convolve traces with Calcium kernel
PSTH_conv = zeros(size(PSTH));
for iunit = 1:size(PSTH_conv,1)
    %PSTH_conv(iunit,:) = conv(PSTH(iunit,:), kernel.trace, 'full');
    tmp = conv(PSTH(iunit,:), kernel.trace, 'full');
    PSTH_conv(iunit,:) = tmp(1:size(PSTH_conv,2));
end

% Normalize traces
for iunit = 1:size(PSTH_conv,1)
    PSTH_conv(iunit,:) = PSTH_conv(iunit,:) - min(PSTH_conv(iunit,:));
    PSTH_conv(iunit,:) = PSTH_conv(iunit,:) / max(PSTH_conv(iunit,:));
end


%% Smooth the full-resolution PSTH

sigma = 50;

PSTH_smooth = nan(size(PSTH,1), length(t));
for iunit = 1:size(PSTH,1)    
    [PSTH_smooth(iunit,:), t_smooth] = SmoothAndResampleC(PSTH(iunit,:), t, [], 'sec', 'gaussian', sigma); % smooth
%     else
%         % not smoothed but resampled
%         outTimeStamps = 1/outFreq : 1/outFreq : maxdur;
%         PSTH(unit,:) = interp1((1:size(psth,2))*dt,meanPSTH,outTimeStamps);
%         
%         PSTH_ori_samp(unit,:) = meanPSTH;
%     end
end

% Normalize traces
for iunit = 1:size(PSTH_smooth,1)
    PSTH_smooth(iunit,:) = PSTH_smooth(iunit,:) - min(PSTH_smooth(iunit,:));
    PSTH_smooth(iunit,:) = PSTH_smooth(iunit,:) / max(PSTH_smooth(iunit,:));
end

%% Resample convolved PSTH
t_interp = 0 : kernel.dt : timeLims(2);
PSTH_conv_res = nan(size(PSTH_conv,1), length(t_interp));
for iunit = 1:size(PSTH_conv,1)
    PSTH_conv_res(iunit,:) = interp1(t, PSTH_conv(iunit,:), t_interp);
end

% Normalize traces
for iunit = 1:size(PSTH_conv_res,1)
    PSTH_conv_res(iunit,:) = PSTH_conv_res(iunit,:) - min(PSTH_conv_res(iunit,:));
    PSTH_conv_res(iunit,:) = PSTH_conv_res(iunit,:) / max(PSTH_conv_res(iunit,:));
end


%% Save data
% Will be loaded in plot_convolved_response.m

save_dir    = '../data/';
save(fullfile(save_dir, 'lgn_psth.mat'), 'PSTH', 't');
save(fullfile(save_dir, 'lgn_psth_s50_f1000.mat'), 'PSTH_smooth', 't');
save(fullfile(save_dir, 'lgn_psth_convolved.mat'), 'PSTH_conv', 't');
save(fullfile(save_dir, 'lgn_psth_convolved_resampled.mat'), 'PSTH_conv_res', 't_interp');




