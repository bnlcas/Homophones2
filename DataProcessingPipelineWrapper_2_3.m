%% Wrapper to generate ERP structure for data.
%% MARK 2:
% Generates standard Homophone ERPs structures for a given subject
% Similar to the functionality of DataProcessingPipeline (_1)
% and Builds upon the existing framework
% 
% New Features:
% Generates ERPs for Several Frequency Bands of data (Alpha, Theta,...)
%
% Reads Data Directly from
% Automatically labels trials with High Z data points
%
% Automatically Detects Folders and Fills in Gaps if necessary
% 
% Reads data from 40 band 400 hz pre-processed data folders to create
% erps for several bands of data
% 
% Designed to be more seemless in its function and also provide status
% updates/prompt user when necessary.
%
%
% Version 2.1:
% Automatically add a .grd field to the Main ERP structure
% that contains information about the tiling orientation of this grid
% also plot using an updated standard plot function
% 
% Version 2.3
% unifies the plot functions of 2.2 and grd functions of 2.1
%
%% Wrapper Mark 1 Example Data:
% Must list the director for the subject's data
% This folder should contain the following folders:
%
% Event_Audio - audio files for homophones
%
% behavior - contains the experiment data files from the task
% This should also contain the file 'homophone_list_full.mat'
% (This is important for determining dominant/subordinance)
%
% ecog_data - contains subfolders for each block containing the the Hilbert
% Transfored data, bad trials, bad channels and audio
% % Each block folder should be formatted
% ECXXX_BYY (Ex: EC125_B100)
% 
% % root_dir = '/Users/changlab/Documents/data/EC125';  % Subject directory described above
% % subj = 'EC125';
% % blocks = {'B46', 'B10', 'B35', 'B38', 'B40', 'B29', 'B31'};



%function [] = DataProcessingPipelineWrapper_2(rdoot_dir)

root_dir = '/Users/changlab/Documents/data/EC133';  % Subject directory described above



%% Establish Directories:
tmp = strsplit(root_dir, filesep);
if strcmp(tmp{end},'')
    tmp(end) =[]; % kill tail
    root_dir = strjoin(tmp, filesep);
end
subj = tmp{end};
fprintf('Subject ID: %s, ', subj); fprintf('\n')


%% Get List of Blocks:
if exist([root_dir filesep 'ecog_data']) == 7
    dat_files = dir([root_dir filesep 'ecog_data']);
    blocks = {};
    for k = 1:length(dat_files)
        tmp = strsplit(dat_files(k).name,'_');
        if length(tmp) > 1
            if strcmpi(tmp{1},subj)
                if strncmp(tmp{2},'B',1);
                    blocks = [blocks, tmp(2)];
                end
            end
        end
    end
    fprintf(' %s, ', blocks{:}); fprintf('\n')
else
    error('ECoG Data not found - Create folder "ecog_data" with block data')
end

%blocks = {'B19', 'B21', 'B26', 'B31', 'B35', 'B42', 'B49'};

%% Behavior Files:
if exist([root_dir filesep 'behavior']) ~= 7
    error('"behavior" folder not found - Create folder with behavioral data')
end
if exist([root_dir filesep 'Event_Audio']) ~= 7
    error('Stimulus Audio Data Not Found - Create folder Homophone Stimuli Audio files')
end





%% Evnt Files:
%% Combine behavior data and generate homophone evnt files
create_evnt_structs = false;
if exist([root_dir filesep 'events']) == 7
    exists_evnt_block = false(1,length(blocks));
    for k = 1:length(blocks)
        exists_evnt_block(k) = (2 == exist([root_dir filesep 'events' filesep blocks{k} '_evnt.mat']));
    end
    missing_blocks = blocks(~exists_evnt_block);
    if sum(~exists_evnt_block) > 0
        fprintf('event structure not detected for %s \n', missing_blocks{:});
    end
    if sum(exists_evnt_block) ~= length(blocks)
        reply = input('Generate Missing event structures? (\''Y''/''N'') \n (The quote markes must be entered for whatever reason...)\n');
        if strcmpi(reply,'y')
            create_evnt_structs = true;
        elseif strcmpi(reply,'n')
            create_evnt_structs = false;
        else
            create_evnt_structs = true;
            fprintf('Wrong Answer, Generating Evnt Structs anyway')
        end
    end
else
    create_evnt_structs = true;
    missing_blocks = blocks;
    fprintf('event structures not detected\n\n');
end
    
if create_evnt_structs
    fprintf('generating evnt structures... \n Will Take awhile...\n')
    anin_chan = 2; % varies, typically 2, on EC131, this was set to 3
    homophone_process_data(subj, root_dir, missing_blocks,anin_chan); % Generates evnt Strucutres
end

% Manually Clean up evnt_files
evnt_dir = [root_dir filesep 'events'];
manually_modify_evnts = false;
for k = 1:length(blocks)
    if 2 == exist([evnt_dir filesep blocks{k} '_evnt.mat'])
        load([evnt_dir filesep blocks{k} '_evnt.mat'])
        if mod(length(evnt),48) ~= 0
            fprintf('data missing in block %s - possible concern\n', blocks{k});
            if mod(length(evnt),2) ~= 0
                fprintf('Trial missing in block %s - Must Manually delete data\n',blocks{k});
                manually_modify_evnts = true;
            end
        end
    end
end
if manually_modify_evnts
   d = dialog('Position',[300 300 250 150],'Name','Event Size Error');
    txt = uicontrol('Parent',d,...
               'Style','text',...
               'Position',[20 80 210 40],...
               'String','Trial Data Missing from Event Files Manually Delete Missing Pair Note the Change and Rerun');
    btn = uicontrol('Parent',d,...
           'Position',[89 20 70 25],...
           'String','Continue',...
           'Callback','delete(gcf)');
    uiwait(d);
end

localize_evnt_files(root_dir, evnt_dir); % necessary for preprocessing

    
%% Manually Identify Mis-Identifications
% % % Fix bug on EC123_B22 (trail 53 dropped)
% % load([root_dir filesep 'events' filesep 'B31_evnt.mat']);
% % for i = 1:length(evnt) starts(i) = evnt(i).StartTime; end
% % [~,order] = sort(starts);
% % evnt = evnt(order);
% % evnt(53) = [];
% %
% % % Fix bug on EC125_B31 (trial 121 dropped)
% % load([root_dir filesep 'events' filesep 'B31_evnt.mat']);
% % for i = 1:length(evnt) starts(i) = evnt(i).StartTime; end
% % [~,order] = sort(starts);
% % evnt = evnt(order);
% % evnt(121) = [];

% % % Fix bug on EC131_B21 (trial 135 dropped)
% % load([root_dir filesep 'events' filesep 'B21_evnt.mat']);
% % for i = 1:length(evnt) starts(i) = evnt(i).StartTime; end
% % [~,order] = sort(starts);
% % evnt = evnt(order);
% % evnt(135) = [];

%% Save modified evnt files with mis-ID's removed
% % save([root_dir filesep 'events' filesep 'B31_evnt.mat']);
% % save([root_dir filesep 'events' filesep 'B21_evnt.mat']);








%% Generate 40 Band ECoG Data:
% detect is there are already 40 band data folders for each block,
% otherwise make them:
for k = 1:length(blocks)
    block_dir = [root_dir filesep 'ecog_data' filesep subj '_' blocks{k}];
    generate_40_band = (exist([block_dir filesep 'HilbAA_4to200_40band']) ~= 7);
    if generate_40_band
        fprintf('Block % s Hilbert Data (40 band) Not Found.\nGenerating 40 Band Hilbert Data... \nThis will take a while...\n',blocks{k});
        [~,cfs] = transformData(block_dir,'outFlag',2,'CARflag',0);
        fprintf('\n');
        save([block_dir filesep 'cfs_4to200_40band.mat'], 'cfs');
    end
end







%% Create ERP structure from the individual blocks and frequencies:
band_freqs = [4 7; 8 15; 18 30; 33 55; 70 150]; % list of of freuqnecies bounds for the various bands 
band_names = {'theta', 'alpha', 'beta', 'lg', 'hg'}; % names of the bands in the preceding matrix (band_names{k} is found between band_freqs(k,1) and band_freqs(k,2)...

save_erps = true;           % saves the ERP structures once created
plot_test_erps = true;          % plot the ERPs for Subject - useful for debugging
clear_when_saved = true;    % clears erps after generation - saves space
for k = 1:length(band_names);
    EPR_name = ['ERPs_' band_names{k}];
    fprintf('Generating %s band ERP structure\n', band_names{k});
    %% Make unified ERP Structure for blocks of Subject in band_frequencies (k)
    % load center frequencies:
    load([block_dir filesep 'cfs_4to200_40band.mat'])
    band_inds = find((cfs >= band_freqs(k,1)) & (cfs <= band_freqs(k,2))); % indecies of center frequencies in the range
    ERPs = Make_unified_ERP_struct(evnt_dir, band_inds);
    ERPs.BadChans(ERPs.BadChans > size(ERPs.ecog_primes,1)) = []; % kill bad channels with no data.
    
    %% Auto Flag high z-ecents
    high_z_thresh = 10; % z-score >10 considered high
    num_chan_thresh = 1; % 1 or more high channels suffienct to be called high trial
    num_timept_thresh = 4; % 1 or more time points in a trial sufficient to be called high channel
    is_high_z_trial = auto_flag_high_z_trials(ERPs, high_z_thresh, num_chan_thresh, num_timept_thresh); % flag trials with 1 ch and 1 timept with z>10)
    ERPs.is_high_z_trial = is_high_z_trial;

    
    %% Add Lexical Stats
    current_func = 'DataProcessingPipelineWrapper_2_3.m';
    function_dir = which(current_func); 
    lex_stats_dir = [function_dir(1:(length(function_dir) - length(current_func))) 'Pipleline_Functions' filesep 'Linguisitic_Data'];
    ERPs = append_linguistic_stats(ERPs, lex_stats_dir);
    
    
    %% Add Grid Tiling Data structure for Plotting Purposes:
    intial_corner = 'br'; tiling_dir = 'y';
    warning('Grid Plot Tiling varies by Patient - Please Confirm that this is correct') 
    ERPs.grd = grid_subplot_layout('br','y', 'plot_test_grid', true);
    % EC123:
    % ERPs.grd = grid_subplot_layout('br','y', 'plot_test_grid', true);
    % EC125:
    % ERPs.grd = grid_subplot_layout('br','y', 'plot_test_grid', true);
    % EC131:
    % ERPs.grd = grid_subplot_layout('bl','x', 'plot_test_grid', true, 'split_grid_flip', 'vertical');
    
    %% Plot Test Data
    if plot_test_erps
        is_good = (ERPs.is_good_trial == 1);    % Trials with no bad time segments
        is_sub = (ERPs.is_related_subordinant == 1) & is_good;  % subordinate sense primed
        is_dom = (ERPs.is_related_dominant == 1) & is_good; % dominante sense primed
        is_rel = is_sub | is_dom;                       % Related prime
        is_unr = is_good & ~is_sub &~is_dom;            % Rando prime
        is_homophone = is_good & strcmpi(ERPs.target_names,'mouse'); % target homophone is mouse

        ecog1 = ERPs.ecog_targets(:,:,is_rel & is_homophone);
        ecog2 = ERPs.ecog_targets(:,:,is_unr & is_homophone);
        is_sig = ttest2(ecog1, ecog2, 'Alpha', 0.01, 'Dim',3); % t-significant electrodes for pval < 0.001
        PlotECogGrid_std_v3(ERPs, [-1 1], ecog2, ecog1, 'SigPts',is_sig, 'Normalize', false);


        %% Get list of significant channels for some comparisons
        is_sig_chan_rel = find_sig_chans_homophones(ERPs, is_rel, is_unr);
        is_sig_chan_rel_mouse = find_sig_chans_homophones(ERPs, is_rel & is_homophone, is_unr & is_homophone);
    end


    %% Save?
    if save_erps
        structure_dir = [root_dir filesep 'ecog_data' filesep 'ERP_structure'];
        if (exist(structure_dir) ~= 7)
            mkdir(structure_dir);
        end
        save([structure_dir filesep subj '_' EPR_name '.mat'], 'ERPs','-v7.3'); % Is freaking huge...
    end

    %% Clear used data:
    if clear_when_saved
        clear ERPs
    end
    
end




