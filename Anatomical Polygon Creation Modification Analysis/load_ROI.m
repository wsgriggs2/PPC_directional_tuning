function [combined_ROI_mask] = load_ROI(Session, Run, mask_size, varargin)
%% Load ROIs for specified Session and Run
%
% Inputs:
%   Session - scalar;
%   Run     - scalar;
%   mask_size - [yPix xPix]; Desired size of the resulting ROI mask
%   varargin
%       'verbose' - Do you want to plot things?
%       'angiogram' - If no ROIs already exist, this is required to
%                     generate new ROIs.
%       'automatic_loading' - Allows user to use batch processing where all
%                             ROIs in most recently created matching file
%                             are used.
% Author: Whitney Griggs


% Parse input arguments
p = inputParser;
p.CaseSensitive = false;
p.addParameter('verbose', false)
p.addParameter('angiogram', []);
p.addParameter('automatic_loading', false);
p.parse(varargin{:});
inputs = p.Results;

pixelsize = 0.1;
xPix = mask_size(2);
yPix = mask_size(1);

%% Ask if user wants to load a predefined polygon ROI, define a new one, or use entire image
if ~inputs.automatic_loading
    user_choice = questdlg('Choose one of the following options',...
        'ROI choice', 'Load existing ROI(s)','Define new ROI(s)','Use entire image','Use entire image');
    
    
    user_choice = string(user_choice);
else
    user_choice = 'Load existing ROI(s)';
end

%% If load predefined polygon ROI
% If multiple polygon ROIs, then create a mask that is the union of all
% polygon ROIs
if strcmp(user_choice, 'Load existing ROI(s)')
    if inputs.automatic_loading
        % Find most recent matching file
        wildcard_string = fullfile(get_user_data_path('pathType', 'AnatomicalPolygons'), sprintf('*BCI_ROIs_S%dR%d_*.mat', Session, Run));
        allFileInfo = dir(wildcard_string);
        file_datenums = cellfun(@datenum, {allFileInfo.date});
        [~, ind] = max(file_datenums);
        load_fullfile = fullfile(allFileInfo(ind).folder, allFileInfo(ind).name);
        load(load_fullfile, 'ROI_table');
        
    else
        [loadfile, loadpath] = uigetfile(fullfile(get_user_data_path('pathType', 'AnatomicalPolygons'), sprintf('*BCI_ROIs_S%dR%d_*.mat', Session, Run)));
        if all(~loadfile) || all(~loadpath)
            % Create new ROIs
            warning('No file selected, so letting user create new ROI');
            ROI_table = save_new_ROIs(inputs.angiogram, Session, Run);
        else
            load_fullfile = fullfile(loadpath, loadfile);
            try
                
                load(load_fullfile, 'ROI_table');
            catch
                error('Problem loading the ROI file - %s', load_fullfile);
            end
        end
    end
    
    
    % Find union of the specified ROI masks
    combined_ROI_mask = select_ROIs(ROI_table, yPix, xPix, ...
        'verbose', inputs.verbose, ...
        'angiogram', inputs.angiogram, ...
        'automatic_loading', inputs.automatic_loading);
    
    
    
    %% If define new one
    % This should support multiple polygon ROIs, not just one
elseif strcmp(user_choice, 'Define new ROI(s)')
    if isempty(inputs.angiogram)
        error('Cannot create new ROIs without angiogram. Check inputs to this function');
    end
    
    % Create new ROIs
    ROI_table = save_new_ROIs(inputs.angiogram, Session, Run);
    
    % Find union of specified ROI masks
    combined_ROI_mask = select_ROIs(ROI_table, yPix, xPix);
    
    %% If use entire_image, then return ROI_mask as ones
elseif isempty(user_choice) || strcmp(user_choice, 'Use entire image')
    combined_ROI_mask = true(mask_size);
end

if inputs.verbose
    X_img_mm = pixelsize/2 + (0:xPix-1)*pixelsize; % So that the center of the image is at true halfway point, not shifted.
    Z_img_mm = pixelsize/2 + (0:yPix-1)*pixelsize;
    
    
    figure;
    imagesc(X_img_mm, Z_img_mm, combined_ROI_mask); colorbar;
    title('Logical mask for union of selected ROIs');
    xlabel('mm');
    ylabel('Relative depth (mm)');
    
end