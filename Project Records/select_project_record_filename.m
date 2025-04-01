function project_record_filename = select_project_record_filename()
%% Basic function to allow user to select project record filename.
% This file needs to reside in same location as ProjectRecord*.json files.

% Author - Whitney Griggs
% October 20, 2022

%% Find all Project Records
startPath = pwd;
current_file_loc = mfilename('fullpath');
[filepath, ~, ~] = fileparts(current_file_loc);

cd(filepath);
rawfiles = dir('ProjectRecord*.json');
if ~isempty(rawfiles)
    project_record_filenames = {rawfiles.name};
    
    cd(startPath);
    
    %% Select Project Record to use
    
    % get Project for which this is a dataset
    [indx,~] = listdlg('PromptString','Which project record do you want to use?',...
        'ListString',project_record_filenames,...
        'SelectionMode','single',...
        'Name','Project Record Selection', ...
        'ListSize', [600 300]);
    project_record_filename = project_record_filenames{indx};
else
    error('No appropriate project record filename was chosen.');
end

