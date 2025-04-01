function ROI_table = save_new_ROIs(angiogram, Session, Run)
%% Create and save new ROIs
% Create new ROIs using roipoly, perform any desired rotations on these
% polygons, then save the result in a table with user-specified ROI title.
% Inputs:
%   angiogram - anatomical image that you wish to draw ROIs onto
%   Session - scalar; for creating filename for saving
%   Run     - scalar; for creating filename for saving
% Outputs:
%   ROI_table - Copy of the ROI table that is saved to disk. Has 3 columns.
%                   ROI name
%                   ROI vertices (in image space)
%                   Rotation amount - Useful for calculating dorsal/ventral
%                                     automatically
%
% Written by Whitney 2022/02/08

    [yPix, xPix] = size(angiogram);
    
    pixelsize = 0.1;
    X_img_mm = pixelsize/2 + (0:size(angiogram,2)-1)*pixelsize; % So that the center of the image is at true halfway point, not shifted.
    Z_img_mm = pixelsize/2 + (0:size(angiogram,1)-1)*pixelsize;


    figure;
    plotDuplexImage(X_img_mm, Z_img_mm, NaN(size(angiogram)), angiogram,...
        'nonlinear_bg',2, ...
        'showColorbar', false);
    
    createAnotherROI = true;
    ROI_num = 0;
    while createAnotherROI
        ROI_num = ROI_num + 1;
        title(['Define region - ' num2str(ROI_num)]);
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
        anatomical_polygon_vertices = [x_px, y_px];
        
        % Smooth ROI
        smoothed_points = fnplt(cscvn([anatomical_polygon_vertices(:, 1)' anatomical_polygon_vertices(1, 1);
            anatomical_polygon_vertices(:, 2)'  anatomical_polygon_vertices(1, 2)]));
        
        temp_cell = inputdlg('Type in name of this ROI');
        ROI_names{ROI_num} = temp_cell{1};
        AnatomicalPolygon_info{ROI_num, 1} = smoothed_points';
        
        % Ask user if they want to add another ROI
        new_roi_answer = questdlg('Do you want to create another ROI?',...
        'Add another ROI?', 'Yes','No','No');
        if isempty(new_roi_answer)
            createAnotherROI = false;
        else
            createAnotherROI = strcmp(new_roi_answer, 'Yes');
        end
    end
    
    % Create the rotation matrices for each region. This can then
    % be used by downstream functions to automatically split into dorsal
    % and ventral components.
    for ROI = 1:length(ROI_names)
        % Define shape of LIP
        roi_poly = polyshape(AnatomicalPolygon_info{ROI, 1}(:, 1), AnatomicalPolygon_info{ROI, 1}(:, 2));

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
        AnatomicalPolygon_info{ROI, 2} = app.rotation_tform;
        close(app.UIFigure);
        clear app;
    end
    
    % Create a table of the ROIs, thus linking the boundaries with their
    % respective ROI name.
    ROI_table = table(AnatomicalPolygon_info(:, 1), cell2mat(AnatomicalPolygon_info(:, 2)), 'RowNames', ROI_names,'VariableNames', {'Boundary', 'Rotation'});
    
    % Save ROIs if desired
    if ~isempty(Session) && ~isempty(Run)
        save_ROI_response = questdlg('Would you like to save these new ROIs for this session/run?','save_new_ROIs','yes','no','no');
        save_ROIs = strcmp(save_ROI_response, 'yes');
    else
        save_ROIs = false;
    end

    if save_ROIs       
        savepath = fullfile(get_user_data_path('pathType', 'AnatomicalPolygons'), ...
        sprintf('BCI_ROIs_S%dR%d_generated_%s.mat', ...
        Session, Run, datetime('today', format='yyyyMMdd')));
        save(savepath, 'ROI_table',...
            'Session', 'Run', ...
            'angiogram');
    end