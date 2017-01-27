function [] = localize_evnt_files(root_dir, evnt_dir)
%% this function loads a directory of event files that were generated
% on a different machine and changes the file directories to include
% the localize directory locations specified in root_dir


%% get list of evnt files:
event_files = dir(evnt_dir); % list of all possible event structures
is_file = false(size(event_files));
for k = 1:length(event_files)
    file_name_parts = strsplit(event_files(k).name,'_');
    if length(file_name_parts) ==2
        if strcmpi(file_name_parts{2}, 'evnt.mat') & strncmpi(file_name_parts{1}(1),'B',1)
            is_file(k) = true;
        end
    end
end
event_files(~is_file) = [];


%% Loop through and change:
for k = 1:length(event_files)
    load([evnt_dir filesep event_files(k).name]);
    evnt = dpaths2local(evnt, root_dir); % Localize dpathes and wnames
    save([evnt_dir filesep event_files(k).name], 'evnt');
end

end