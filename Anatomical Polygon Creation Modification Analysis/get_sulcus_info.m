function sulcus_vertices = get_sulcus_info(Session, Run, angiogram, UF, varargin)

%% This function will return the sulcus information for a given session/run combo.
% If no sulcus info exists for a given session/run combo, it will allow you
% to create and save the sulcus info.

%% Input parser
p = inputParser;
addOptional(p, 'pixelsize', 0.1);
addOptional(p, 'verbose', false);
addOptional(p, 'use_existing_sulcus_map', false);
parse(p, varargin{:});
result = p.Results;
pixelsize = result.pixelsize;

%%
X_img_mm = pixelsize/2 + (0:size(angiogram,2)-1)*pixelsize; % So that the center of the image is at true halfway point, not shifted.
Z_img_mm = pixelsize/2 + (0:size(angiogram,1)-1)*pixelsize + UF.Depth(1);

% Get size of image
[yPix, xPix] = size(angiogram);

%% Check if sulcus ROIs already exist for this session/run
% Define filename and location
if ~isempty(Session) && ~isempty(Run)
    ROI_filename = strcat('Sulcus_ROI_S', num2str(Session), '_R', num2str(Run), '.mat');
    data_path = get_user_data_path('pathType','sulcusmap');
    sulcus_fullfile = fullfile(data_path, ROI_filename);
else
    sulcus_fullfile = '';
end

% Check for file existence
if isfile(sulcus_fullfile)
    create_new_sulcus_response = questdlg('Would you like to define new sulcus for this session/run?','create_new_sulcus','yes','no','no');
    create_new_sulcus = strcmp(create_new_sulcus_response, 'yes');
else
    % create the ROIs if they do not exist already
    create_new_sulcus = true;
end

if create_new_sulcus
    % Display anatomical angiogram
    figure;
    plotDuplexImage(X_img_mm, Z_img_mm, NaN(size(angiogram)), angiogram,...
        'nonlinear_bg',2, ...
        'showColorbar', false);
    title('Draw out borders of the sulcus and brain surface');
    %Select open polygon using drawpolyline.
    polyline_handle = drawpolyline;
    VerticePositions = polyline_handle.Position;
    
    %Draw anatomical polygon onto subplot
    hold on;
    hline = plot(VerticePositions(:, 1), VerticePositions(:, 2),...
        ['g' '.-'],...
        'MarkerSize', 15);
    hold off;
    
    % Smooth polyline
    cs = cscvn([VerticePositions(:, 1)';
        VerticePositions(:, 2)']);
    smoothed_points = fnplt(cs);
    
    % Define conversion from world coordinates (mm) into image coordinates(px)
    x = smoothed_points(1, :);
    y = smoothed_points(2,:);
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
    
    smoothed_points_px = [x_px' y_px'];
    
    % Add the remaining vertices to finish out the sulcus polygon
    [~, min_ind] = min(smoothed_points_px(:, 1));
    if min_ind > size(smoothed_points_px, 1)/2
        smoothed_points_px = smoothed_points_px(end:-1:1, :);
    end
    smoothed_points_px = [0 smoothed_points_px(1, 2);
        smoothed_points_px;
        xPix smoothed_points_px(end, 2);
        xPix 0;
        0 0];
    
    % Define polygon - Useful for visualization
    sulcus_poly = polyshape(smoothed_points_px(:, 1), smoothed_points_px(:, 2));
    figure;
    plot(sulcus_poly);
    xlim([0 xPix]);
    ylim([0 yPix]);
    
    % Create mask from the sulcus polygon
    sulcus_vertices = sulcus_poly.Vertices;
    sulcus_vertices(isnan(sulcus_vertices)) = 0;
    
    sulcus_mask = poly2mask(sulcus_vertices(:, 1), sulcus_vertices(:, 2), yPix, xPix);
    if result.verbose
        figure;
        imagesc(sulcus_mask);
    end
    
    % Save ROIs if desired
    if ~isempty(Session) && ~isempty(Run)
        save_sulcus_response = questdlg('Would you like to save the sulcus drawing for this session/run?','save_new_sulcus','yes','no','no');
        save_sulcus = strcmp(save_sulcus_response, 'yes');
    else
        save_sulcus = false;
    end
    
    if save_sulcus
        save(sulcus_fullfile, 'sulcus_vertices', 'sulcus_mask',...
            'Session', 'Run', ...
            'angiogram', 'UF');
    end
else
    
    % If ROI is already been defined previously, then load that.
    load(sulcus_fullfile, 'sulcus_vertices', 'angiogram');
    [yPix, xPix] = size(angiogram);
    sulcus_mask = poly2mask(sulcus_vertices(:, 1), sulcus_vertices(:, 2), yPix, xPix);
    if result.verbose
        figure;
        imagesc(sulcus_mask);
    end
end
