%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display total number of cells %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PARAMETERS
corr_p      = 0.001;         % 0.001; corr_p represents the wilcoxon ransum correlation for between/withing segments
qi          = 0.05;          % 0.05; Philipps quality index
cluster_qi  = 3;             % 3; cluster quality index of likely single- vs multiunit
r           = 0.7;           % 0.7; regression index

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% GET ALL UNITS
% get the chirp keys but exclude all experiments before '2014-03-14' and all ntsr-mice
chirp_keys = (data.Series('series_date > "2014-03-14"') & ...
    data.Experiments('exp_name LIKE "%chirp%"')); %- data.Mice('mouse_id LIKE "%ntsr%"');

% Exclude all experiments before '2014-03-14'
date_keys = fetch(data.Series(chirp_keys),'series_date');
mask_date = datenum(fetchn(data.Series(chirp_keys),'series_date')) > datenum('2014-03-14', 'yyyy-mm-dd');

% Contains all valid chirp keys
chirp_keys = date_keys(mask_date);
chirp_keys = rmfield(chirp_keys,'series_date');

% get all units
exp_keys = fetch(data.Experiments(chirp_keys) & 'exp_name LIKE "%chirp%"');
all_units = fetch(data.Units(exp_keys),'unit_id');


%% GET CHIRP UNITS

% Create a table of units that respond satisfactorily to chirp stimulus
chirp_units = fetch(miro.ChirpQuality(chirp_keys) ...
    & data.ClusterInfo(sprintf('quality <= %d',cluster_qi)) ...
    & miro.ChirpQuality(sprintf('corr_p <= %f',corr_p)) ...
    & miro.ChirpQuality(sprintf('berens_qi >= %f',qi)));


%% Get SFTFOri UNITS
sftfori_units = fetch(food.SFTFOriUnitStability(sprintf('r>=%f',r)) & data.ClusterInfo(sprintf('quality <= %f',cluster_qi)));

%% GET CHIRP & SFTFORI UNITS
unitList = getAllUnits(corr_p,qi,cluster_qi,r,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  DISPLAY STATISTICS
fprintf('\nUNIT STATS:\n');
fprintf('----------\n');
fprintf('all units:  %d\n',numel(all_units));
fprintf('chirp units:  %d\n',numel(chirp_units));
fprintf('SFTFOri units:  %d\n',numel(sftfori_units));
fprintf('chirp & SFTFOri units:  %d\n',numel(unitList));














































