function path = get_user_data_path(varargin)
% get_user_data_path Retrieve path to the user's data
%
% INPUTS:
%   varargin:
%       path_type:      String; Available options are: 
%                       {'root', 'doppler', 'rois', 
%                       'anatomicalpolygons','sulcusmap', 'output'}
% 
% Outputs:
%   path:               string; Absolute path to where data is stored

%% Variable argument input parser
p = inputParser;
p.addParameter('pathType','root', @ischar)
p.parse(varargin{:});
variable_inputs = p.Results;

% we don't care about case-sensitivity
variable_inputs.pathType = lower(variable_inputs.pathType); 

%% Check if root has already been set
folder = fileparts(mfilename('fullpath')); 
folder_root = fileparts(folder);
data_path_file =  fullfile(folder_root, 'data_path.mat');
if isfile(data_path_file)
    % If path has already been specified, then load that.
    load(data_path_file, 'base_path');
else
    % If path has not been set, then ask user to specify a location
    base_path = specify_data_path;
end

%% create end of path
switch variable_inputs.pathType
    
    % DOPPLER
    case 'doppler'
        append_path = 'doppler';
    % ROIs and anatomical polygons and sulcus maps
    case 'rois'
        append_path = 'ROIs';
    case 'anatomicalpolygons'
        append_path = fullfile('ROIs', 'Anatomical_Polygons');
    case 'sulcusmap'
        append_path = fullfile('ROIs', 'SulcusMaps');
    case 'output'
        append_path = 'output';
    case 'glm'
        append_path = fullfile('output', 'glm');
    case 'decoding'
        append_path = fullfile('output', 'decoding');
    case 'across session analyses'
        append_path = fullfile('output', 'across session analyses');
    case 'root'
        append_path = '';
    otherwise
        error('path type problem: path not set!');
end

%% create full path 
path = fullfile(base_path,append_path);

end