function evnt_out = dpaths2local(evnt, root_dir)
%% This function takes evnt structures generated on a different machine
% and ammends the event structures to the local directory
evnt_out = evnt;
for i = 1:length(evnt)
    % Split dpath (check for mac or window format:
    dpath_mac = strsplit(evnt(i).dpath,'/');
    wname_mac = strsplit(evnt(i).wname,'/');
    
    dpath_win = strsplit(evnt(i).dpath,'\');
    wname_win = strsplit(evnt(i).wname,'\');
%     if length(wname_mac) < 2
%         a =1 ;
%     end
    if length(dpath_win) > length(dpath_mac)
        dpath_stem_dir = dpath_win;
        wname_stem_dir = wname_win;
    else
        dpath_stem_dir = dpath_mac;
        wname_stem_dir = wname_mac;
    end
    if isempty(dpath_stem_dir{end}) % kill haning empty cell (is extra slash added)
        dpath_stem_dir(end) = [];
    end
    if isempty(wname_stem_dir{end})
        wname_stem_dir(end) = [];
    end
    
    %dpath_stem_dir = strjoin(dpath_stem_dir, filesep); % combine end of directory
    %dpath_out = [root_dir filesep dpath_stem_dir];
    %dpath_stem_dir = strjoin(dpath_stem_dir((end-1):end), filesep);
    dpath_out = [root_dir filesep 'ecog_data' filesep dpath_stem_dir{end}];

    
    %wname_stem_dir = strjoin(wname_stem_dir, filesep);
    wname_out = [root_dir filesep 'Event_Audio' filesep wname_stem_dir{end}];
    
    evnt_out(i).dpath = dpath_out;
    evnt_out(i).wname = wname_out;
end



end
