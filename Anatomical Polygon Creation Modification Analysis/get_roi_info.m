function ROI_table = get_roi_info(Session, Run, angiogram, UF, varargin)

%% This function will return the ROI table for a given session/run combo.
% If no roi table exists for a given session/run combo, it will allow you
% to create and save the ROI table.

% Added ability to have custom filanem on June 20, 2024 to allow loading of
% more restrictive ROI.

%% Input parser
p = inputParser;
addOptional(p, 'region_names',{'LIP', 'PRR'});
addOptional(p, 'pixelsize', 0.1);
addOptional(p, 'custom_filename', '');
parse(p, varargin{:});
result = p.Results;
region_names = result.region_names;
pixelsize = result.pixelsize;

%%
X_img_mm = pixelsize/2 + (0:size(angiogram,2)-1)*pixelsize; % So that the center of the image is at true halfway point, not shifted.
Z_img_mm = pixelsize/2 + (0:size(angiogram,1)-1)*pixelsize + UF.Depth(1);

% Get size of image
[yPix, xPix] = size(angiogram);

%% Check if ROIs already exist for this session/run
% Define filename and location
if ~isempty(result.custom_filename)
    ROI_filename = result.custom_filename;
    data_path = get_user_data_path('pathType','rois');
    ROI_fullfile = fullfile(data_path, ROI_filename);
elseif ~isempty(Session) && ~isempty(Run)
    ROI_filename = strcat('Anatomical_ROI_S', num2str(Session), '_R', num2str(Run), '.mat');
    data_path = get_user_data_path('pathType','rois');
    ROI_fullfile = fullfile(data_path, ROI_filename);
else
    ROI_fullfile = '';
end

% Check for file existence
if isfile(ROI_fullfile)
    create_new_ROI_response = questdlg('Would you like to define new ROIs for this session/run?','create_new_ROIs','yes','no','no');
    create_new_ROIs = strcmp(create_new_ROI_response, 'yes');
else
    % create the ROIs if they do not exist already
    create_new_ROIs = true;
end

if create_new_ROIs
    
    % Display anatomical angiogram
    figure;
    plotDuplexImage(X_img_mm, Z_img_mm, NaN(size(angiogram)), angiogram,...
        'nonlinear_bg',2, ...
        'showColorbar', false);
    
    AnatomicalPolygon_info = cell(length(region_names), 2);
    for region = 1:length(region_names)
        title(['Define region - ' region_names{region}]);
        [~, x, y] = roipoly; %Select polygon anatomical polygon
        
        % Define conversion from world coordinates (mm) into image
        % coordinates(px)
        x_limits = xlim;
        x_spacing = (x_limits(2)-x_limits(1))/size(angiogram,2);
        y_limits = ylim;
        y_spacing = (y_limits(2)-y_limits(1))/size(angiogram,1);
        
        % Use coordinate conversion to convert vertices into image coordinates
        % y_limits/y_spacing used for initial offset (e.g. starting at 3
        % mm instead of 0 mm).
        x_px = x/x_spacing;
        y_px = y/y_spacing - y_limits(1)/y_spacing;
        
        %If polygon goes out of figure bounds, define edge to be figure bound.
        x(x>x_limits(2))=x_limits(2); x(x<x_limits(1))=x_limits(1);
        y(y>y_limits(2))=y_limits(2); y(y<y_limits(1))=y_limits(1);
        
        x_px(x_px>xPix)=xPix; x_px(x_px<0)=0;
        y_px(y_px>yPix)=yPix; y_px(y_px<0)=0;
        
        %Draw anatomical polygon onto subplot
        hold on;
        hline = plot(x, y,...
            ['g' '.-'],...
            'MarkerSize', 15);
        hold off;
        
        % Store the anatomical polygons
        AnatomicalPolygon_info{region, 1} = [x_px, y_px];
        
        % Smooth ROI
        smoothed_points = fnplt(cscvn([AnatomicalPolygon_info{region, 1}(:, 1)' AnatomicalPolygon_info{region, 1}(1, 1);
            AnatomicalPolygon_info{region, 1}(:, 2)'  AnatomicalPolygon_info{region, 1}(1, 2)]));
        
        % If smoothed ROI goes out of figure bounds, define edge to be
        % figure bound.
        x_smoothed = smoothed_points(1, :);
        y_smoothed = smoothed_points(2, :);
        x_smoothed(x_smoothed>xPix)=xPix; x_smoothed(x_smoothed<0)=0;
        y_smoothed(y_smoothed>yPix)=yPix; y_smoothed(y_smoothed<0)=0;
        
        AnatomicalPolygon_info{region, 1} = [x_smoothed' y_smoothed'];
    end
    
    % Create the rotation matrices for each region. This can then
    % be used by downstream functions to automatically split into dorsal
    % and ventral components.
    for region = 1:length(region_names)
        % Define shape of LIP
        roi_poly = polyshape(AnatomicalPolygon_info{region, 1}(:, 1), AnatomicalPolygon_info{region, 1}(:, 2));

        % Rotate polygon so that the y-axis of the ROI is as vertical as
        % possible.
        app = interactive_rotate_polygon(roi_poly);

        pauseCount=0;
        %Mandatory wait period because other the code
        %keeps running forward without waiting for user input
        %from the polygon rotation GUI.
        while isempty(app.rotation_tform)
            pause(1);
            pauseCount=pauseCount+1;
            if pauseCount == 2000
                error('Rotation_tform was not created within 2000 seconds. Stopping code');
                break;
            end
        end
        AnatomicalPolygon_info{region, 2} = app.rotation_tform;
        close(app.UIFigure);
        clear app;
    end
    
    % Create a table of the ROIs, thus linking the boundaries with their
    % respective ROI name.
    ROI_table = table(AnatomicalPolygon_info(:, 1), cell2mat(AnatomicalPolygon_info(:, 2)), 'RowNames', region_names,'VariableNames', {'Boundary', 'Rotation'});
    
    % Save ROIs if desired
    if ~isempty(Session) && ~isempty(Run)
        save_ROI_response = questdlg('Would you like to save these new ROIs for this session/run?','save_new_ROIs','yes','no','no');
        save_ROIs = strcmp(save_ROI_response, 'yes');
    else
        save_ROIs = false;
    end

    if save_ROIs
        save(ROI_fullfile, 'ROI_table',...
            'Session', 'Run', ...
            'angiogram', 'UF');
    end
else
    
    % If ROI is already been defined previously, then load that.
    load(ROI_fullfile, 'ROI_table');
end
