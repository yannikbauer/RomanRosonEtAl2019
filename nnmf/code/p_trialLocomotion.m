function [speed,loco_ts,trial_mean_speed,trial_max_speed,trial_percent_moving,trial_dir] = p_trialLocomotion(myKey)

samplingrate = 30000;
trial_movement_threshold = 1;

% Fetch on and offset of trial
[onset, offset] = fetchn(data.ExtraTrials(myKey),'trial_onset','trial_offset');

% Fetch all spikes of trial and set timing relative to onset
loco_key = myKey;
loco_key = rmfield(loco_key, 'trial_num');

try
    [t, speed, direction] = fetch1(data.Locomotion(loco_key), 'timestamps', 'speed', 'direction');
    
    % convert to seconds
    onset = onset/samplingrate;
    offset = offset/samplingrate;
    t = t/samplingrate;
    
    index_of_trials = (t >= onset) & (t <= offset);
    
    if(all(~index_of_trials))
        error('TrialMovement:loco', 'Trial smaller 11 ms or corrupted TreadmillData found, aborting...');
    end
    
    speed = speed(index_of_trials);
    
    if isnan(direction)
        trial_dir = nan;
    else
        trial_dir = direction(index_of_trials);
    end
    
    loco_ts = t(index_of_trials) - onset;
    
    % Get some means and stds
    trial_mean_speed = mean(speed);
    trial_max_speed = max(speed);
    
    % Mark movement and calculate percentage
    trial_percent_moving = sum(speed > trial_movement_threshold) / length(speed);
    
catch
    speed = [];
    loco_ts = [];
    trial_mean_speed = [];
    trial_max_speed = [];
    trial_percent_moving = [];
    trial_dir = [];
end
end
