function [spectral_erps, time_axis, freq_axis] = create_spectrogram_erps(ERPs, data_dir, subj, varargin);
% returns a chans x 40 band x timepoints maxtrix of spectral
% ERPs for the 40 bands of hilbert transformed data averaged across trials
% 
% Inputs:
% ERPs - standard ERP data structure for a subject  in questions
%
% data_dir - directory for the folder containing the ecog data for the
% various blocks of Homophone Experiments for a subject (eg:
% data_dir = 'userdata/data_dump/EC125/ecog_data';
%
% subj - (string) subject name [data_dir \filesep subj '_' block], contains
% block data
%
% Variable Inputs:
% 1 - select_trials (boolean for list of trials to choose from)
%
% 2 - erp_twin (2 point array for start and stop time relative to event for
% erps (ex [-2 4])
%
% 3 - use_target_locked (default = true), false - use prime locked events
% 
% Outputs:
% spectral_erps - n_chans x n_bands x n_timepts array of the average
% intensity across select trials for a given frequency band, channel and time point
%
% time_axis - n_timept array listing the time relative to event onset
%
% freq_axis - n_band array listing the center frequency for each band
%
%

%% Set Variable Inputs:
include_trial = ERPs.is_good_trial;
erp_twin = [-1 2]; % time window of erps
erp_zwin = [-0.5 0]; % time window to z-score erps (pre prime)
use_target_lock = true;
if length(varargin) > 0
    if ~isempty(varargin{1})
        include_trial = varargin{1};
    end
end
if length(varargin) > 1
    if ~isempty(varargin{2})
        erp_twin = varargin{2};
    end
end
if length(varargin) > 2
    if ~isempty(varargin{3})
        use_target_lock = varargin{3};
    end
end
if use_target_lock
    event_start_times = ERPs.target_starttimes;
else
    event_start_times = ERPs.prime_starttimes;
end

%% Initialize Outputs:
fs = 400;
fs_ds = 100; % downsample frequency
time_axis = linspace(erp_twin(1), erp_twin(2), round(diff(erp_twin)*fs_ds));
load([data_dir filesep subj '_' ERPs.data_block{1} filesep 'cfs_4to200_40band.mat']);
freq_axis = cfs;

spectral_erps = zeros(size(ERPs.ecog_targets,1), length(time_axis), length(freq_axis));


%% Loop through trials and add to spectral data:
for k = 1:length(include_trial)
    if include_trial(k)
        spectral_erps_trial = zeros(size(ERPs.ecog_targets,1), length(time_axis), length(freq_axis));

        hilb_dir = [data_dir filesep subj '_' ERPs.data_block{k} filesep 'HilbAA_4to200_40band'];
        dat_twin = (event_start_times(k) + erp_twin)*1000;
        zscore_twin = (ERPs.prime_starttimes(k) + erp_zwin)*1000; % (in ms)
        for i = 1:4
            for j = 1:64
                [dat_event, fs] = readhtk([hilb_dir filesep 'Wav' num2str(i) num2str(j) '.htk'], dat_twin);
                dat_zscore = readhtk([hilb_dir filesep 'Wav' num2str(i) num2str(j) '.htk'], zscore_twin);
                dat_event_z = gdivide(gsubtract(dat_event,mean(dat_zscore,2)), std(dat_zscore,[],2));
                for h = 1:size(dat_event_z,1)
                    dat_ds = resample(dat_event_z(h,:), 1, 4); % magic constants - happy times (just dont change anything (EVER!))
                    if length(dat_ds) ~= length(time_axis)
                        if length(dat_ds) > length(time_axis)
                            dat_ds((length(time_axis)+1):end) = [];
                        else
                            dat_ds(length(time_axis)) = 0;
                        end
                    end
                    spectral_erps_trial((j+64*(i-1)),:,h) = dat_ds; %dat_event_z(1:length(time_axis));
                end
            end
        end
        spectral_erps = spectral_erps+spectral_erps_trial;
    end
end
   
%% Average ERPs:
spectral_erps = spectral_erps/(sum(include_trial));

%% Zero Bad Channels:
spectral_erps(ERPs.BadChans,:,:) = 0;

end