function reaction_times = get_response_times(evnt)
%% Attempt to find the reaction times for an evnt structure:
block = evnt(1).block;

%dat_dir = '/Users/changlab/Documents/data/EC123/data';
%audio_file = [dat_dir '/EC123_' block '/Analog/ANIN1.wav'];
audio_dir = [evnt(1).dpath filesep 'Analog'];
if exist([audio_dir filesep 'ANIN1.wav']) == 2
    [audio, fs_audio] = audioread([audio_dir filesep 'ANIN1.wav']);
else
    [audio, fs_audio] = readhtk([audio_dir filesep 'ANIN1.htk']);
end
%% Load event, load audio

%[audio, fs_audio] = audioread(audio_file);


%% order start times
start_times = zeros(1,length(evnt));
for j = 1:length(evnt) 
    start_times(j) = evnt(j).StartTime;
end
[~,order] = sort(start_times);
evnt = evnt(order);

max_resp_time = 5;
onset_thresh = 0.04;
num_events = length(evnt)/2;
reaction_times = zeros(num_events, 1);
for i = 1:num_events
    prompt_time = evnt(2*i).StopTime + 0.5;
    audio_dat = audio(round(prompt_time*fs_audio):round((prompt_time+max_resp_time)*fs_audio));
    audio_dat_env = smooth(abs(hilbert(audio_dat)),fs_audio/10);
    audio_dat_env(1:round(fs_audio/10)) = 0;
    onset_ind = find(audio_dat_env > onset_thresh,1);
    if isempty(onset_ind)
        reaction_times(i) = max_resp_time;
    else
        reaction_times(i) = onset_ind/fs_audio;
    end
end


end
    

