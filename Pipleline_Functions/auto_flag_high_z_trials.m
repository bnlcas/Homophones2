function is_high_z_trial = auto_flag_high_z_trials(ERPs, varargin)
% This function takes the standard homophone ERP structure and generates a
% num_trials boolean array that is true if a trial has events with high
% z-scores
%
% Inputs:
% ERPs - standard homophone ERP structure
%
% Variable Inputs:
% 1 - high_z_thresh: the ECoG z-score threshold of what is considered
% a high z-score data point - (default = 10)
%
% 2 - num_chan_thresh: (the threshold number of channels with high z-score
% time points required to consider a trial bad - (default = 1)
%
% 3 - num_timept_thresh: (the threshold number of timepoints with high z-score
% time points required to consider a channel to be a bad channel - (default = 1)
% 
% Outputs:
%
% high_z_trials - num_trials x 1 Boolean Array (ON if trial matches
% conditions set to be bad)
%
% Created by Ben Lucas 10-20-16

%% Note on the assebly of ECoG for a trial:
% ecog_prime and ecog_target overlap - to prevent redundancy (messes with
% the thresholds), we will clip the prime to exclude data after
% presentation of the stimulus, and the target to exclude data from before
% the 500 ms pre stimulus
%
% Additionally, it is assumed that trials that have high z-events in either
% Primes or Targets are bad and these will be rejected (same treatment is
% given to manually rejected data):



%% Function Parameters:
high_z_thresh = 10;
num_chan_thresh = 1;
num_timept_thresh = 1;
if length(varargin) > 0
    high_z_thresh = varargin{1};
end
if length(varargin)>1
    num_chan_thresh = varargin{2};
end
if length(varargin)>1
    num_timept_thresh = varargin{3};
end

time_axis = ERPs.time_axis;
is_good_channel = true(size(ERPs.ecog_primes,1),1);
is_good_channel(ERPs.BadChans) = false;

is_high_z_trial = false(size(ERPs.is_good_trial)); % Copy form of the manual list of good trials

%% loop through each trial to detmine if it is clear of high-z:
% (Cannot be vectorized because duration of each trial is variable...
ch_high = false(sum(is_good_channel),length(is_high_z_trial));
for k = 1:length(is_high_z_trial)
    pre_prime_plocked = (time_axis < (ERPs.prime_duration(k)+0.1));    % time points of prime locked that occur before end of prime
    post_prime_tlocked = (time_axis >= -0.5);                % time points of target locked that occur after end of prime

    ecog_prime = squeeze(ERPs.ecog_targets(is_good_channel,pre_prime_plocked,k));
    ecog_target = squeeze(ERPs.ecog_targets(is_good_channel,post_prime_tlocked,k));

    ecog_trial = cat(2,ecog_prime, ecog_target); % combined ecog for the whole trial

    is_high_z = abs(ecog_trial) > high_z_thresh;

    high_times = sum(is_high_z,2); % total timepts above thresh in channels
    ch_is_high_z = (high_times >= num_timept_thresh);
    
    high_chans = sum(ch_is_high_z); %number of channels with >= threshold bad time poitns
    is_high_z_trial(k) = (high_chans >= num_chan_thresh);
    
    ch_high(find(ch_is_high_z),k) = true;
end

% figure; imagesc(ch_high); ylabel('Good Channel'); xlabel('Trial');
% title('Distribution of High Z-Channels for Each Trial')
a = 1;
end

