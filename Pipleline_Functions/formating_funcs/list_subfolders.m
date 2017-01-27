function sub_folder_list = list_subfolders(root_dir)
%% This function takes a directory location as its input and returns a cell array
% listing the names of the subfolders found in that array.

files_dir = dir(root_dir);
files_dir(1:2) = [];
file_names = cell(size(files_dir));
is_folder = true(size(files_dir));
for k = 1:length(files_dir)
    is_folder(k) = (exist([root_dir filesep files_dir(k).name]) == 7);
    file_names{k} = files_dir(k).name;
end
file_names(~is_folder) = [];

sub_folder_list = file_names;

end