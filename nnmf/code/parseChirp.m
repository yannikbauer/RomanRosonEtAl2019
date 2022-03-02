function psth = parseChirp(unitList, dt, timeLims)
% psth = parseChirp(unitList,dt,timeLims) computes psth the chirp 
% stimulus. 
%
%   unitList    = list with unit keys
%   dt          = data sampling
%   timeList    = are the time limits of the stimulus
%
%   [psth] = parseChirp(...) structure containing the psth data for steps,
%                  freq and the con component of the chirp stimulus.
%                  Moreover, it contains pshts for the whole stimulus and
%                  a psth suitable for plotting
%
if nargin < 3
    timeLims = [0 32.05]; 
end
if nargin < 2
    dt = 0.001;
end
if nargin < 1
    unitList = getAllUnits(0.001,0.05,3,3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nFUNCTION: COMPUTE PSTHS\n')

% LOW SAMPLING (STEPS) - NORMALIZED TO 1
sigma = 40;  
outFreq = 20; %7.8 
norm_chirp = 0;

[PSTH, ts] = computePSTH(unitList, sigma, dt, timeLims, outFreq, norm_chirp);
PSTH = PSTH(ts >= timeLims(1) & ts <= timeLims(2),:);
ts = ts(ts >= timeLims(1) & ts <= timeLims(2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2 FULLSTIM - STEPS - NORMALIZED TO 1 
psth.psth     = PSTH';
psth.tl       = timeLims;
psth.ts       = ts;
psth.sigma    = sigma;
psth.outfreq  = outFreq;

disp('...done')


















