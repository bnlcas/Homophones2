function [] = create_channel_spectrogram_erps(ERPs, data_dir, subj, varargin);
% [spectral_erps, time_axis, freq_axis] = create_channel_spectrogram_erps(ERPs, data_dir, subj, varargin);
% saves a band x timepoints maxtrix x trials matrix of spectral data
% for the 40 bands of hilbert transformed data averaged across trials
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
% Defaults, could be variable inputs:
% twin = [-1 1];
% 
% saveout = true;
%
% select_trials = all of them
% use_target_locked (default = true), false - use prime locked events
% 
% Outputs:
% Saved to [data_dir filesep SpectralERPs]
%
% spectral_erps_ch_(number) (256x - one for each channel)
% each is 40x(timepts)x(n_trials) and each band is znormalized to [-0.5 1]
% prestim
%
% time_axis - n_timept array listing the time relative to event onset
%
% freq_axis - n_band array listing the center frequency for each band
%
%

%% Set Variable Inputs:
%include_trial = ERPs.is_good_trial;
erp_twin = [-1 1]; % time window of erps
erp_zwin = [-0.5 0]; % time window to z-score erps (pre prime)

event_start_times = ERPs.target_starttimes;
trial_start_times = ERPs.prime_starttimes;

%% Initialize Outputs:
fs_in = 400;
fs_ds = 100; % downsample frequency
time_axis = linspace(erp_twin(1), erp_twin(2), round(diff(erp_twin)*fs_ds));

load([data_dir filesep subj '_' ERPs.data_block{1} filesep 'cfs_4to200_40band.mat']);
freq_axis = cfs;

out_dir = [data_dir filesep 'spectrogram_data'];
mkdir(out_dir)
save([out_dir filesep 'freq_axis.mat'], 'freq_axis');
save([out_dir filesep 'time_axis.mat'], 'time_axis');

n_chans = 256;
n_trials = length(ERPs.is_good_trial);
%% Make spectral ERPs for each channel
for n = 1:n_chans
    spectral_erps_ch = zeros(length(time_axis),length(freq_axis), n_trials);
    %% loop through trials
    if (sum(ERPs.BadChans == n) == 0)
        for k = 1:n_trials
            % trial parameters:
            hilb_dir = [data_dir filesep subj '_' ERPs.data_block{k} filesep 'HilbAA_4to200_40band'];
            dat_twin = (event_start_times(k) + erp_twin)*1000;
            zscore_twin = (trial_start_times(k) + erp_zwin)*1000; % (in ms)
            i = ceil(n/64);
            j = mod(n,64);
            if j == 0
                j = 64;
            end
            % Load data and zscore
            [dat_event, fs_in] = readhtk([hilb_dir filesep 'Wav' num2str(i) num2str(j) '.htk'], dat_twin);
            dat_zscore = readhtk([hilb_dir filesep 'Wav' num2str(i) num2str(j) '.htk'], zscore_twin);
            dat_event_z = gdivide(gsubtract(dat_event,mean(dat_zscore,2)), std(dat_zscore,[],2));       

            % Downsample each band to 100 hz
            for h = 1:size(dat_event_z,1)
                dat_ds = resample(dat_event_z(h,:), 1, 4); % magic constants - happy times (just dont change anything (EVER!))
                spectral_erps_ch(:,h,k) = dat_ds(1:length(time_axis));
            end
        end
    end
    save([out_dir filesep 'spectral_erps_ch_' num2str(n) '.mat'], 'spectral_erps_ch');

end

end