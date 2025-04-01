%% Set up code repository and data directory
% This script will let the user set the path to the downloaded data. The
% data can be found on CaltechDATA at https://doi.org/10.22002/f3y3k-em558


%% 1. Add folder and subfolder to matlab search path
% Add the parent folder above the current folder and all
% subfolders of that parent folder to the MATLAB search path. Then save the
% search path.

folder = fileparts(mfilename('fullpath')); 
addpath(genpath(folder));
rmpath(genpath(fullfile(folder, '.git')));
savepath;

%% 2. Set the path to the downloaded and extracted data
% If data is not extracted, then extract before doing this step. Keep same
% data organization and file structure as in .zip file.

path = specify_data_path;
fprintf('\nSpecified "%s" as the path to the downloaded data.\n', path);



