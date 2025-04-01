%% Script to plot results of `AnalyzeAcrossPlanes_8tgt`

%% Load workspace
clear; close all; clc
[loadfile, loadpath] = uigetfile(fullfile(get_user_data_path('pathType', 'across session analyses'), '*8dir_tuningCurve_data*.mat'));
load_fullfile = fullfile(loadpath, loadfile);

load(load_fullfile);
pixelsize = 0.1;

PositiveNegativeColormap = load('PercentChangeColormap.mat');
PositiveNegativeColormap = PositiveNegativeColormap.PercentChangeColormap;

%% User parameters
disk_filter_size = 1;


%% Create a movie
pvalue_threshold = 0.01; % For 1-way anova

use_transparency = true;

min_windows = min(nWindows(nWindows > 0));
% Scale effectsize for transparency
effectsize_1D = [];


fix_start_end = mean(task_timing(:, 1:2));
memory_end = mean(task_timing(:, 3));
movement_end = mean(task_timing(:, 4));

for timenum = 1:min_windows
    for plane = 1:length(ap_planes)
        effectsize_anova = effectsize_anova_time_cell{plane};
        effectsize_anova_timepoint = squeeze(effectsize_anova(:, :, timenum));
        effectsize_1D = [effectsize_1D; reshape(effectsize_anova_timepoint, [], 1)];
    end
end
transparency_axis = max(abs(effectsize_1D), [], 'all') * 0.8;

% Generate figure showing saturation combined with circular colormap
number_saturation_elements = 1000;
number_colors = size(circle_colormap, 1);

circle_colormap_repmat = repmat(flipud(circle_colormap), 1, 1, number_saturation_elements);
circle_colormap_repmat = permute(circle_colormap_repmat, [1 3 2]);
alpha_map = repmat(linspace(0, 1, number_saturation_elements), number_colors, 1, 1);

ValuesToDisplay = [-pi pi];
if use_transparency
    filename = sprintf('%s_serial_3d_tuningcurve_rho_anovaThreshold_%0.2f_withTransparency_DateGenerated_%s.mp4', Monkey, pvalue_threshold, datetime('today', format='yyyyMMdd'));
else
    filename = sprintf('%s_serial_3d_tuningcurve_rho_anovaThreshold_%0.2f_DateGenerated_%s.mp4', Monkey, pvalue_threshold, datetime('today', format='yyyyMMdd'));
end
data_save_root = get_user_data_path('pathType', 'across session analyses');
video_fullfile = fullfile(data_save_root, 'movies', filename);

vidfile = VideoWriter(video_fullfile,'MPEG-4');
vidfile.FrameRate = 1;
open(vidfile);

fig1 = figure('units','normalized','outerposition',[0 0 0.5 1]);
t = tiledlayout('flow');
for timenum = 1:min_windows
    for plane = 1:length(ap_planes)
        nexttile(plane)
        ap_plane = ap_planes(plane);
        iDop_tuningcurve_theta = iDop_tuningcurve_theta_cell{plane};
        iDop_tuningcurve_rho = iDop_tuningcurve_rho_cell{plane};
        effectsize_anova = effectsize_anova_time_cell{plane};
        effectsize_anova_timepoint = effectsize_anova(:, :, timenum);
        anova_across_time = qvalueFDR_anova_time_cell{plane};
        anova_timepoint = anova_across_time(:, :, timenum);
        
        iDop_plot = iDop_tuningcurve_theta(:, :, timenum);
        iDop_plot2 = iDop_tuningcurve_rho(:, :, timenum);
        
        % Threshold on 1-way ANOVA
        anova_qvalue_mask = anova_timepoint < pvalue_threshold;
        iDop_plot(~anova_qvalue_mask) = NaN;
        iDop_plot2(~anova_qvalue_mask) = NaN;
        
        
        pixelsize = 0.1;
        
        % Transparency map
        transparency_map = nan(size(iDop_plot));
        if use_transparency
            % Set only the voxels that are shown to be different
            scaled_effectsize = effectsize_anova_timepoint/transparency_axis;
            scaled_effectsize(scaled_effectsize>1) = 1;
            if any(anova_qvalue_mask, 'all')
                transparency_map(anova_qvalue_mask) = scaled_effectsize(anova_qvalue_mask);
            end
        end
        
        
        X_img_mm = pixelsize/2 + (0:size(iDop_tuningcurve_theta,2)-1)*pixelsize; % So that the center of the image is at true halfway point, not shifted.
        Z_img_mm = pixelsize/2 + (0:size(iDop_tuningcurve_theta,1)-1)*pixelsize + UF{plane}.Depth(1);
        
        plotDuplexImage(X_img_mm, Z_img_mm, iDop_plot, angiogram{plane},...
            'colormap2use', circle_colormap, 'nonlinear_bg',2, ...
            'showColorbar', false, ...
            'AutoColorBarLimits',ValuesToDisplay,...
            'transparencyMap', transparency_map, ...
            'showColorbar', false,...
            'ColorBarTitle', 'direction', ...
            'FontSize', 9,...
            'darkestShade', 0.25);
        if mod(ap_plane, 1) == 0
            title(sprintf('AP Plane: %d', ap_plane));
        else
            title(sprintf('AP Plane: %0.2f', ap_plane));
        end
    end
    
    %Specify the location of the colorbar
    nexttile(length(ap_planes)+1)
    
    % Specify the legend
    if use_transparency
        im = image( linspace(0, transparency_axis, number_saturation_elements), linspace(-pi, pi, number_colors), circle_colormap_repmat);
        im.AlphaData = alpha_map;
        xlabel(sprintf('Effect size (%s)', anova_effectsize_measure));
        ylabel('Direction (rad)');
        axis square
        yticks([-pi, -pi/2, 0, pi/2, pi])
        yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
        set(gca,'YAxisLocation','right')
        
        
    else
        % Turn off the axes
        axis square;
        axis off;
        % Draw the colorbar
        c = colorbar('FontSize',9,...
            'Location','north',...
            'Ticks', [-pi, -pi/2, 0, pi/2, pi], ...
            'TickLabels', {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
        % Specify the colormap
        colormap(gca, circle_colormap)
        % Define the colormap limits.
        clim([-pi pi]);
        % Label the colorbar
        c.Label.String = 'direction (rad)';
    end
    
    time_aligned_to_cue = round(timenum - 1 + fix_start_end(1));
    if time_aligned_to_cue < 0
        title(t, sprintf('Preferred direction maps:\nThreshold 1-way ANOVA q < %0.2f\n Time since cue: %d sec \n Fixation period', pvalue_threshold, time_aligned_to_cue));
    elseif time_aligned_to_cue == 0
        title(t, sprintf('Preferred direction maps:\nThreshold 1-way ANOVA q < %0.2f\n Time since cue: %d sec \n Cue period', pvalue_threshold, time_aligned_to_cue));
    elseif time_aligned_to_cue > 0 && time_aligned_to_cue < memory_end
        title(t, sprintf('Preferred direction maps:\nThreshold 1-way ANOVA q < %0.2f\n Time since cue: %d sec \n Memory period', pvalue_threshold, time_aligned_to_cue));
    elseif time_aligned_to_cue > memory_end
        title(t, sprintf('Preferred direction maps:\nThreshold 1-way ANOVA q < %0.2f\n Time since cue: %d sec \n Movement period', pvalue_threshold, time_aligned_to_cue));
    end
    
    F = getframe(fig1);
    writeVideo(vidfile, F);
    
    
end
close(vidfile)


%% Align to composite functional angiogram to set angiogram for AP Plane
% This prevents the weird overlapping lines situtation from using the
% composite functional angiogram.

realign_angiogram = true;

if realign_angiogram
    aligned_angiogram = align_data_to_session_angiogram(angiogram, session_run_list);
    
    save_fullfile = fullfile(data_save_root, sprintf('%s_aligned_angiogram_dateGenerated_%s.mat', Monkey, datetime('today', format='yyyyMMdd')));
    save(save_fullfile, 'aligned_angiogram', 'Monkey');
else
    [loadfile, loadpath] = uigetfile(fullfile(get_user_data_path('pathType', 'across session analyses'), '*aligned_angiogram*.mat'));
    load_fullfile = fullfile(loadpath, loadfile);
    load(load_fullfile, 'aligned_angiogram');
    
    if size(aligned_angiogram, 1) > size(ap_planes)
        aligned_angiogram = aligned_angiogram(~remove_ind);
    end
    
end


%% % Define LIP and MIP ROIs and mask area outside of LIP/MIP ROIs.
make_new_rois = false;

if realign_angiogram || make_new_rois
    figure;
    pixelsize = 0.1;
    region_names = {'LIP', 'MIP'};
    
    % Preallocate space
    AnatomicalPolygon_info = cell(length(ap_planes), length(region_names));
    
    for plane = 1:length(ap_planes)
        % Load in slot data
        ap_plane = ap_planes(plane);
        anatomical_angiogram = aligned_angiogram{plane};
        
        xPixels = size(anatomical_angiogram,2);
        yPixels = size(anatomical_angiogram,1);
        
        % Cut all data to the x and z size of the anatomical angiogram. Since
        % some planes anatomical data are at angles. Use any column/row which
        % has no NaNs. Can later revise to keep any row with at least one
        % real number if that is cleaner (image-wise).
        anatomical_angiogram(anatomical_angiogram==0) = NaN;
        row_ind = all(isnan(anatomical_angiogram), 2);
        col_ind = all(isnan(anatomical_angiogram), 1);
        
        anatomical_angiogram(row_ind, :) = [];
        anatomical_angiogram(:, col_ind) = [];
        
        X_img_mm = pixelsize/2 + (0:xPixels-1)*pixelsize; % So that the center of the image is at true halfway point, not shifted.
        Z_img_mm = pixelsize/2 + (0:yPixels-1)*pixelsize + UF{plane}.Depth(1) + nnz(row_ind)*pixelsize;
        
        
        % Display anatomical angiogram
        plotDuplexImage(X_img_mm, Z_img_mm, NaN(size(anatomical_angiogram)), anatomical_angiogram,...
            'nonlinear_bg',2, ...
            'showColorbar', false, ...
            'darkestShade', 0.25);
        
        for region = 1:length(region_names)
            f = msgbox(sprintf('Define %s ROI', region_names{region}));
            [~, x, y] = roipoly; %Select polygon anatomical polygon
            
            % Define conversion from world coordinates (mm) into image
            % coordinates(px)
            x_limits = xlim;
            x_spacing = (x_limits(2)-x_limits(1))/size(anatomical_angiogram,2);
            y_limits = ylim;
            y_spacing = (y_limits(2)-y_limits(1))/size(anatomical_angiogram,1);
            
            % Use coordinate conversion to convert vertices into image coordinates
            % y_limits/y_spacing used for initial offset (e.g. starting at 3
            % mm instead of 0 mm).
            x_px = x/x_spacing;
            y_px = y/y_spacing - y_limits(1)/y_spacing;
            
            %If polygon goes out of figure bounds, define edge to be figure bound.
            x(x>x_limits(2))=x_limits(2); x(x<x_limits(1))=x_limits(1);
            y(y>y_limits(2))=y_limits(2); y(y<y_limits(1))=y_limits(1);
            
            x_px(x_px>xPixels)=xPixels; x_px(x_px<0)=0;
            y_px(y_px>yPixels)=yPixels; y_px(y_px<0)=0;
            
            %Draw anatomical polygon onto subplot
            hold on;
            hline = plot(x, y,...
                ['g' '.-'],...
                'MarkerSize', 15);
            hold off;
            
            % Store the anatomical polygons
            AnatomicalPolygon_info{plane, region} = [x_px, y_px];
            
            % Smooth ROI
            smoothed_points = fnplt(cscvn([AnatomicalPolygon_info{plane, region}(:, 1)' AnatomicalPolygon_info{plane, region}(1, 1);
                AnatomicalPolygon_info{plane, region}(:, 2)'  AnatomicalPolygon_info{plane, region}(1, 2)]));
            
            AnatomicalPolygon_info{plane, region} = smoothed_points';
        end
        
    end
    
    save_fullfile = fullfile(data_save_root, 'across session analyses', sprintf('%s_LIP_ROIs_dateGenerated_%s.mat', Monkey, datetime('today', format='yyyyMMdd')));
    save(save_fullfile, 'AnatomicalPolygon_info', 'region_names', 'Monkey');
else
    [loadfile, loadpath] = uigetfile(fullfile(get_user_data_path('pathType', 'across session analyses'), '*LIP_ROIs*.mat'));
    load_fullfile = fullfile(loadpath, loadfile);
    load(load_fullfile, 'AnatomicalPolygon_info', 'region_names');
    
    if size(AnatomicalPolygon_info, 1) > size(ap_planes)
        AnatomicalPolygon_info = AnatomicalPolygon_info(~remove_ind, :);
    end
end


%% Create LIP, MIP, and d/v splitted ROIs for each plane
% Based upon Liu et al. 2010, LIPd/v division is at 53% of sulcal depth.
% Will split into two groups. 0-53, 53 - 100 sulcal depth. This prevents
% overlap of the two regions. We are estimating bottom of LIP as bottom of
% sulcus, but this may be inaccurate. This is largely based upon 2015
% Calabrese atlas. Lewis and Van Essen 2000 shows the VIP as starting just
% above the bottom of sulcus.

% Using 1 - 0.53 because the polygon is flipped from the anatomical
% orientation. Therefore 1 - 0.53 represents the distance to the LIPd/v
% divide from the "top" of the polygon as defined in pixel space.
LIPdv_sulcal_percent = .47;

% For dorsal/ventral MIP, using the same division.
% 1) Define/Load in ROI
% 2) Create shape that is ROI size/shape.
% 3) Rotate shape until it is roughly vertical. Use the same rotation for
%    LIP and MIP.
% 4) Define top 53% as LIPd and define bottom 47% as LIPv
% 4.5) Using same divisions for MIP
% 5) Map back to original orientation
% 6) Apply mask

[LIPd_mask, LIPv_mask, MIPd_mask, MIPv_mask, LIP_mask, MIP_mask] = deal(cell(length(ap_planes), 1));

if realign_angiogram || make_new_rois
    
    ROI_rotations = cell(length(ap_planes), 1);
    for plane = 1:length(ap_planes)
        % Define shape of LIP
        LIP_poly = polyshape(AnatomicalPolygon_info{plane, 1}(:, 1), AnatomicalPolygon_info{plane, 1}(:, 2));
        
        % Rotate polygon so that the y-axis of the ROI is as vertical as
        % possible.
        app = interactive_rotate_polygon(LIP_poly);
        
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
        ROI_rotations{plane} = app.rotation_tform;
        close(app.UIFigure);
        clear app;
    end
    
    save_fullfile = fullfile(data_save_root, sprintf('%s_ROI_rotation_matrices_dateGenerated_%s.mat', Monkey, datetime('today', format='yyyyMMdd')));
    save(save_fullfile, 'ROI_rotations', 'Monkey');
else
    [loadfile, loadpath] = uigetfile(fullfile(get_user_data_path('pathType', 'across session analyses'), '*ROI_rotation_matrices*.mat'));
    load_fullfile = fullfile(loadpath, loadfile);
    load(load_fullfile, 'ROI_rotations');
    
    if size(ROI_rotations, 1) > size(ap_planes)
        ROI_rotations = ROI_rotations(~remove_ind);
    end
end

[LIPd_poly_cell, LIPv_poly_cell, MIPd_poly_cell, MIPv_poly_cell] = deal(cell(length(ap_planes), 1));
rows_cols_to_cut = cell(length(ap_planes), 2);
imaging_plane_sizes = zeros(length(ap_planes), 2);


for plane = 1:length(ap_planes)
    % Load in slot data
    ap_plane = ap_planes(plane);
    anatomical_angiogram = aligned_angiogram{plane};
    
    % Cut all data to the x and z size of the anatomical angiogram. Since
    % some planes anatomical data are at angles. Use any column/row which
    % has no NaNs. Can later revise to eliminate any row with at least one
    % NaN if that is cleaner.
    anatomical_angiogram(anatomical_angiogram==0) = NaN;
    row_ind = all(isnan(anatomical_angiogram), 2);
    col_ind = all(isnan(anatomical_angiogram), 1);
    
    anatomical_angiogram(row_ind, :) = [];
    anatomical_angiogram(:, col_ind) = [];
    
    % Instead of repeatedly loading the anatomical angiogram just to get
    % these measure, just save them into a variable.
    rows_cols_to_cut{plane, 1} = row_ind;
    rows_cols_to_cut{plane, 2} = col_ind;
    
    imaging_plane_sizes(plane, 1) = size(anatomical_angiogram,2); %xPixels
    imaging_plane_sizes(plane, 2) = size(anatomical_angiogram,1); %yPixels
    
    % Define shape of LIP
    LIP_poly = polyshape(AnatomicalPolygon_info{plane, 1}(:, 1), AnatomicalPolygon_info{plane, 1}(:, 2));
    
    
    % Define shape of MIP
    MIP_poly = polyshape(AnatomicalPolygon_info{plane, 2}(:, 1), AnatomicalPolygon_info{plane, 2}(:, 2));
    
    % Rotate polygon
    LIP_rotated = rotate(LIP_poly, ROI_rotations{plane});
    MIP_rotated = rotate(MIP_poly, ROI_rotations{plane});
    
    % Define polygons in mm space, not pixel space
    X_img_mm = pixelsize/2 + (0:size(anatomical_angiogram,2)-1)*pixelsize; % So that the center of the image is at true halfway point, not shifted.
    Z_img_mm = pixelsize/2 + (0:size(anatomical_angiogram,1)-1)*pixelsize + UF{plane}.Depth(1) + nnz(row_ind)*pixelsize;
    
        
    % Create masks for each area
    LIP_mask{plane} = poly2mask(LIP_poly.Vertices(:, 1), LIP_poly.Vertices(:, 2), ...
        imaging_plane_sizes(plane, 2), imaging_plane_sizes(plane, 1));
    MIP_mask{plane} = poly2mask(MIP_poly.Vertices(:, 1), MIP_poly.Vertices(:, 2), ...
        imaging_plane_sizes(plane, 2), imaging_plane_sizes(plane, 1));    
end


%% Specify which regions we want to analyze
% use a GUI to select which ones you want to load
inclusive_region_names = [region_names(:)', 'all'];
[roi_indx,tf] = listdlg('PromptString', 'Select region(s) to analyze:',...
    'SelectionMode', 'multiple',...
    'ListString', inclusive_region_names);

%% Visualize strongest tuning during given period -
% Standardized across monkeys to have same scale for figure purposes
strength_clim = [0 0.75];
circ_std_clim = [1 3.25];
fscore_clim = [0.4 2.55];
kurtosis_clim = [-1 1];
skewness_clim = [-0.3 0.3];
    

% find absolute strongest tuning in specified window
window_of_interest = memory_end;
pvalue_threshold = 0.01;

ROI_mask = 'None'; % Other options are 'LIP', 'MIP', 'None'

% Define smoothing filter
if disk_filter_size > 0
    disk_filter = fspecial('disk', disk_filter_size);
end


baseline = fix_start_end;

ValuesToDisplay = [-pi pi];

timeValuesToDisplay = window_of_interest;

index_of_interest = round(window_of_interest - baseline(1));
% Handle single timepoint index
if isscalar(index_of_interest)
    index_of_interest = [index_of_interest index_of_interest];
end

fig_tuning = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
t_tuning = tiledlayout('flow');
fig_strength = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
t_strength = tiledlayout('flow');
fig_Fscore = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
t_Fscore = tiledlayout('flow');
fig_std = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
t_std = tiledlayout('flow');
fig_skewness = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
t_skewness = tiledlayout('flow');
fig_kurtosis = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
t_kurtosis = tiledlayout('flow');

fig_laterality = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
t_laterality = tiledlayout('flow');


for plane = 1:length(ap_planes)
    
    % Load in slot data
    ap_plane = ap_planes(plane);
    iDop_tuningcurve_theta = iDop_tuningcurve_theta_cell{plane};
    iDop_tuningcurve_rho = iDop_tuningcurve_rho_cell{plane};
    F_value = F_cell{plane};

    iDop_circ_stats = iDop_tuningcurve_circ_stats{plane};
    iDop_pvalue_8dir = qvalueFDR_F_cell{plane};
    anatomical_angiogram = aligned_angiogram{plane};
    
    % Extract circular statistics of interest
    timepoints_to_extract = length(iDop_circ_stats);
    [circ_std, skewness, kurtosis] = deal(NaN(size(iDop_circ_stats{1}.circ_std, 1),...
        size(iDop_circ_stats{1}.circ_std, 2), ...
        timepoints_to_extract));
    for i = 1:timepoints_to_extract
        circ_std(:, :, i) = iDop_circ_stats{1}.circ_std_alt; % This one is the one used by Wikipedia. Not bounded.
        skewness(:, :, i) = iDop_circ_stats{1}.skewness_alt; %Fisher standardized version
        kurtosis(:, :, i) = iDop_circ_stats{1}.kurtosis_alt; %Fisher standardized version 
    end


    
    % Laterality index
    % Ignore up/down
    % Formula
    % From "Sustained Activity in Topographic Areas of Human Posterior 
    % Parietal Cortex during Memory-Guided Saccades
    % Denis Schluppeck, Clayton E. Curtis, Paul W. Glimcher, and David J. Heeger
    % Same formula is also used by Kagan et al. 2010
    contra_ind = [1 2 8];
    ipsi_ind = [4 5 6];

    response = iDop_cohenD_8dir_cell{plane};

    contra_response = mean(response(:, :, :, contra_ind), 4);
    ipsi_response = mean(response(:, :, :, ipsi_ind), 4);
    laterality_response = (contra_response - ipsi_response) ./ (abs(contra_response) + abs(ipsi_response));
    
    % Cut all data to the x and z size of the anatomical angiogram.
    row_ind = rows_cols_to_cut{plane, 1};
    col_ind = rows_cols_to_cut{plane, 2};
    iDop_tuningcurve_theta(row_ind, :, :) = [];
    iDop_tuningcurve_theta(:, col_ind, :) = [];
    iDop_tuningcurve_rho(row_ind, :,:) = [];
    iDop_tuningcurve_rho(:, col_ind, :) = [];
    F_value(row_ind,:) = [];
    F_value(:, col_ind) = [];
    anatomical_angiogram(row_ind, :) = [];
    anatomical_angiogram(:, col_ind) = [];
    iDop_pvalue_8dir(row_ind, :) = [];
    iDop_pvalue_8dir(:, col_ind) = [];
    
    circ_std(row_ind, :,:) = [];
    circ_std(:, col_ind, :) = [];
    skewness(row_ind, :, :) = [];
    skewness(:, col_ind, :) = [];
    kurtosis(row_ind, :, :) = [];
    kurtosis(:, col_ind, :) = [];


    laterality_response(row_ind, :, :) = [];
    laterality_response(:, col_ind, :) = [];
    
    % cut to portion of interest
    iDop_tuningcurve_theta = iDop_tuningcurve_theta(:, :, index_of_interest(1):min(index_of_interest(2), size(iDop_tuningcurve_theta, 3)));
    iDop_tuningcurve_rho = iDop_tuningcurve_rho(:, :, index_of_interest(1):min(index_of_interest(2), size(iDop_tuningcurve_rho, 3)));
    circ_std = circ_std(:, :, index_of_interest(1):min(index_of_interest(2), size(circ_std, 3)));
    skewness = skewness(:, :, index_of_interest(1):min(index_of_interest(2), size(skewness, 3)));
    kurtosis = kurtosis(:, :, index_of_interest(1):min(index_of_interest(2), size(kurtosis, 3)));
    laterality_response = laterality_response(:, :, index_of_interest(1):min(index_of_interest(2), size(laterality_response, 3)));

    
    [max_strength, iDop_tuningcurve_theta_maxInd] = max(iDop_tuningcurve_rho, [], 3);
    [iDop_tuningcurve_theta_max, circ_std_max, ...
        skewness_max, kurtosis_max, laterality_max] = deal(zeros(size(iDop_tuningcurve_theta, 1), size(iDop_tuningcurve_theta, 2)));

    
    for yrange = 1:size(iDop_tuningcurve_theta_maxInd, 1)
        for xrange = 1:size(iDop_tuningcurve_theta_maxInd, 2)
            iDop_tuningcurve_theta_max(yrange, xrange) = iDop_tuningcurve_theta(yrange, xrange, iDop_tuningcurve_theta_maxInd(yrange, xrange));
            circ_std_max(yrange, xrange) = circ_std(yrange, xrange, iDop_tuningcurve_theta_maxInd(yrange, xrange));
            skewness_max(yrange, xrange) = skewness(yrange, xrange, iDop_tuningcurve_theta_maxInd(yrange, xrange));
            kurtosis_max(yrange, xrange) = kurtosis(yrange, xrange, iDop_tuningcurve_theta_maxInd(yrange, xrange));
            laterality_max(yrange, xrange) = laterality_response(yrange, xrange, iDop_tuningcurve_theta_maxInd(yrange, xrange));
        end
    end
    
    max_time = iDop_tuningcurve_theta_maxInd;
    
    % Apply mask where we ignore voxels with pvalue for any of the
    % directions that is less than pvalue_threshold.
    
    anova_pvalue_ind = iDop_pvalue_8dir < pvalue_threshold;
    
    iDop_tuningcurve_theta_max(~anova_pvalue_ind) = NaN;
    max_strength(~anova_pvalue_ind) = NaN;
    F_value(~anova_pvalue_ind) = NaN;
    max_time(~anova_pvalue_ind) = NaN;
    
    circ_std_max(~anova_pvalue_ind) = NaN;
    skewness_max(~anova_pvalue_ind) = NaN;
    kurtosis_max(~anova_pvalue_ind) = NaN;
    
    laterality_max(~anova_pvalue_ind) = NaN;

    
    if disk_filter_size > 0
        iDop_tuningcurve_theta_max_smoothed = imfilter(iDop_tuningcurve_theta_max, disk_filter);
        max_time_smoothed = imfilter(max_time, disk_filter);
        max_strength_smoothed = imfilter(max_strength, disk_filter);
        F_value_smoothed = imfilter(F_value, disk_filter);
        circ_std_max_smoothed = imfilter(circ_std_max, disk_filter);
        skewness_max_smoothed = imfilter(skewness_max, disk_filter);
        kurtosis_max_smoothed = imfilter(kurtosis_max, disk_filter);
        
        laterality_response_smoothed = imfilter(laterality_max, disk_filter);
     
    else
        iDop_tuningcurve_theta_max_smoothed = iDop_tuningcurve_theta_max;
        max_time_smoothed = max_time;
        max_strength_smoothed = max_strength;
        F_value_smoothed = F_value;
        circ_std_max_smoothed = circ_std_max;
        skewness_max_smoothed = skewness_max;
        kurtosis_max_smoothed = kurtosis_max;
        
        laterality_response_smoothed = laterality_max;
    end
    


    
    xPixels = size(anatomical_angiogram,2);
    yPixels = size(anatomical_angiogram,1);
    
    switch ROI_mask
        case 'LIP'
            % Generate mask of voxels for each specific LIP AnatomicalPolygon
            roi_mask = poly2mask(AnatomicalPolygon_info{plane, 1}(:,1), AnatomicalPolygon_info{plane, 1}(:,2), yPixels, xPixels);
        case 'MIP'
            % Generate mask of voxels for each specific MIP AnatomicalPolygon
            roi_mask = poly2mask(AnatomicalPolygon_info{plane, 2}(:,1), AnatomicalPolygon_info{plane, 2}(:,2), yPixels, xPixels);
        case 'LIP_and_MIP'
            % Generate mask of voxels for each specific LIP AnatomicalPolygon
            LIP_mask = poly2mask(AnatomicalPolygon_info{plane, 1}(:,1), AnatomicalPolygon_info{plane, 1}(:,2), yPixels, xPixels);
            
            % Generate mask of voxels for each specific MIP AnatomicalPolygon
            MIP_mask = poly2mask(AnatomicalPolygon_info{plane, 2}(:,1), AnatomicalPolygon_info{plane, 2}(:,2), yPixels, xPixels);
            
            % Combine the two masks
            roi_mask = LIP_mask | MIP_mask;
            clear LIP_mask MIP_mask
        case 'None'
            roi_mask = true(yPixels, xPixels);
    end
    
    % Apply ROI mask
    iDop_tuningcurve_theta_max_smoothed(~roi_mask) = NaN;
    max_time_smoothed(~roi_mask) = NaN;
    max_strength_smoothed(~roi_mask) = NaN;
    F_value_smoothed(~roi_mask) = NaN;
    circ_std_max_smoothed(~roi_mask) = NaN;
    skewness_max_smoothed(~roi_mask) = NaN;
    kurtosis_max_smoothed(~roi_mask) = NaN;
    
    laterality_response_smoothed(~roi_mask) = NaN;
    
    X_img_mm = pixelsize/2 + (0:size(anatomical_angiogram,2)-1)*pixelsize; % So that the center of the image is at true halfway point, not shifted.
    Z_img_mm = pixelsize/2 + (0:size(anatomical_angiogram,1)-1)*pixelsize + UF{plane}.Depth(1) + nnz(row_ind)*pixelsize;
    
    % Tuning
    
    % Add transparency layer
    
    transparency_map = isnan(iDop_tuningcurve_theta_max_smoothed);
    
    nexttile(t_tuning, plane);
    plotDuplexImage(X_img_mm, Z_img_mm, iDop_tuningcurve_theta_max_smoothed, anatomical_angiogram,...
        'colormap2use', circle_colormap, 'nonlinear_bg',2, ...
        'showColorbar', false, ...
        'AutoColorBarLimits',ValuesToDisplay,...
        'darkestShade', 0.25, ...
        'transparencyMap', ~transparency_map);
    if mod(ap_plane, 1) == 0
        title(sprintf('AP Plane: %d', ap_plane));
    else
        title(sprintf('AP Plane: %0.2f', ap_plane));
    end
   
    
    % Strength
    nexttile(t_strength, plane);
    plotDuplexImage(X_img_mm, Z_img_mm, max_strength_smoothed, anatomical_angiogram,...
        'nonlinear_bg',2, ...
        'showColorbar', false, ...
        'colormap2use', magma,...
        'AutoColorBarLimits', strength_clim,...
        'darkestShade', 0.25);
    if mod(ap_plane, 1) == 0
        title(sprintf('AP Plane: %d', ap_plane));
    else
        title(sprintf('AP Plane: %0.2f', ap_plane));
    end

    nexttile(t_Fscore, plane);
    plotDuplexImage(X_img_mm, Z_img_mm, log10(F_value_smoothed), anatomical_angiogram,...
        'nonlinear_bg',2, ...
        'showColorbar', false,...
        'AutoColorBarLimits', fscore_clim,...
        'colormap2use', magma,...
        'darkestShade', 0.25);
    if mod(ap_plane, 1) == 0
        title(sprintf('AP Plane: %d', ap_plane));
    else
        title(sprintf('AP Plane: %0.2f', ap_plane));
    end
    
    % Circular std
    nexttile(t_std, plane);
    plotDuplexImage(X_img_mm, Z_img_mm, circ_std_max_smoothed, anatomical_angiogram,...
        'nonlinear_bg',2, ...
        'showColorbar', false,...
        'AutoColorBarLimits', circ_std_clim,...
        'colormap2use', magma,...
        'darkestShade', 0.25);
    if mod(ap_plane, 1) == 0
        title(sprintf('AP Plane: %d', ap_plane));
    else
        title(sprintf('AP Plane: %0.2f', ap_plane));
    end
    
    % Skewness
    abs_max_kurtosis_value = max(abs(skewness_max_smoothed), [], 'all', 'omitnan');
    
    nexttile(t_skewness, plane);
    plotDuplexImage(X_img_mm, Z_img_mm, skewness_max_smoothed, anatomical_angiogram,...
        'nonlinear_bg',2, ...
        'showColorbar', true, ...
        'colormap2use', PositiveNegativeColormap,...
        'showColorbar', false, ...
        'AutoColorBarLimits', skewness_clim,...
        'darkestShade', 0.25);
    if mod(ap_plane, 1) == 0
        title(sprintf('AP Plane: %d', ap_plane));
    else
        title(sprintf('AP Plane: %0.2f', ap_plane));
    end
    
    % Kurtosis
    abs_max_kurtosis_value = max(abs(kurtosis_max_smoothed), [], 'all', 'omitnan');
    nexttile(t_kurtosis, plane);
    plotDuplexImage(X_img_mm, Z_img_mm, kurtosis_max_smoothed, anatomical_angiogram,...
        'nonlinear_bg',2, ...
        'colormap2use', PositiveNegativeColormap,...
        'showColorbar', false, ...
        'AutoColorBarLimits', kurtosis_clim,...
        'darkestShade', 0.25);
    if mod(ap_plane, 1) == 0
        title(sprintf('AP Plane: %d', ap_plane));
    else
        title(sprintf('AP Plane: %0.2f', ap_plane));
    end
    
    % Laterality
    % By definition, goes from -1 to 1
    %abs_max_laterality_value = max(abs(laterality_response_smoothed), [], 'all', 'omitnan');
    
    nexttile(t_laterality, plane);
    plotDuplexImage(X_img_mm, Z_img_mm, laterality_response_smoothed, anatomical_angiogram,...
        'nonlinear_bg',2, ...
        'colormap2use', PositiveNegativeColormap,...
        'showColorbar', false, ...
        'AutoColorBarLimits', [-1 1],...
        'darkestShade', 0.25);
    if mod(ap_plane, 1) == 0
        title(sprintf('AP Plane: %d', ap_plane));
    else
        title(sprintf('AP Plane: %0.2f', ap_plane));
    end
end

nexttile(t_tuning)
% Turn off the axes
axis off;
% Draw the colorbar
c = colorbar('FontSize',18,'Location','layout');
% Specify the colormap
colormap(gca,circle_colormap)
% Define the colormap limits.
clim(ValuesToDisplay);
% Label the colorbar
c.Label.String = 'Preferred direction (rad)';
c.Ticks = -pi:pi/2:pi;
TargetDirectionsRadiansNegPiToPosPi = {'-\pi','-\pi/2','0','\pi/2','\pi'};
c.TickLabels = TargetDirectionsRadiansNegPiToPosPi;

nexttile(t_strength)
% Turn off the axes
axis off;
% Draw the colorbar
c = colorbar('FontSize',18,'Location','layout');
% Specify the colormap
colormap(gca,magma)
% Define the colormap limits.
clim(strength_clim);
c.Ticks = strength_clim(1):diff(strength_clim)/3:strength_clim(2);
% Label the colorbar
c.Label.String = sprintf('Tuning strength (%s)', 'Cohens d');

nexttile(t_Fscore)
% Turn off the axes
axis off;
% Draw the colorbar
c = colorbar('FontSize',18,'Location','layout');
% Specify the colormap
colormap(gca,magma)
% Define the colormap limits.
clim(fscore_clim);
c.Ticks = [fscore_clim(1) 0.5:0.5:fscore_clim(2) fscore_clim(2)];
% Label the colorbar
c.Label.String = 'log(F score)';

nexttile(t_std)
% Turn off the axes
axis off;
% Draw the colorbar
c = colorbar('FontSize',18,'Location','layout');
% Specify the colormap
colormap(gca,magma)
% Define the colormap limits.
clim(circ_std_clim);
c.Ticks = [1:circ_std_clim(2) circ_std_clim(2)];
% Label the colorbar
c.Label.String = 'Circular std';

nexttile(t_kurtosis)
% Turn off the axes
axis off;
% Draw the colorbar
c = colorbar('FontSize',18,'Location','layout');
% Specify the colormap
colormap(gca,PositiveNegativeColormap)
% Define the colormap limits.
clim(kurtosis_clim);
c.Ticks = kurtosis_clim(1):diff(kurtosis_clim)/4:kurtosis_clim(2);
% Label the colorbar
c.Label.String = 'Kurtosis';

nexttile(t_skewness)
% Turn off the axes
axis off;
% Draw the colorbar
c = colorbar('FontSize',18,'Location','layout');
% Specify the colormap
colormap(gca,PositiveNegativeColormap)
% Define the colormap limits.
clim(skewness_clim);
c.Ticks = skewness_clim(1):diff(skewness_clim)/4:skewness_clim(2);
% Label the colorbar
c.Label.String = 'Skewness';


nexttile(t_laterality)
% Turn off the axes
axis off;
% Draw the colorbar
c = colorbar('FontSize',18,'Location','layout');
% Specify the colormap
colormap(gca,PositiveNegativeColormap)
% Define the colormap limits.
clim([-1 1]);
% Label the colorbar
c.Label.String = 'Laterality index';
c.Ticks = -1:0.5:1;


title(t_tuning, sprintf('Strongest preferred direction (rad) - q<%0.2f', pvalue_threshold));
title(t_strength, sprintf('Strength of tuning - q<%0.2f', pvalue_threshold));
title(t_Fscore, sprintf('F-score - q<%0.2f', pvalue_threshold));
title(t_std, sprintf('Circular standard deviation - q<%0.2f', pvalue_threshold));
title(t_skewness, sprintf('Skewness - q<%0.2f', pvalue_threshold));
title(t_kurtosis, sprintf('Kurtosis - q<%0.2f', pvalue_threshold));
title(t_laterality, sprintf('Laterality index - q<%0.2f', pvalue_threshold));

provide_save_option = false;
if provide_save_option
    save_resp = questdlg('Would you like to save these figures just generated?',...
            'save','yes','no','yes');
    save_figures = strcmp(save_resp,'yes');
else
    save_figures = false;
end
if save_figures
    save_directional_analysis_figure(fig_tuning, 'tuned_direction', Monkey);
    save_directional_analysis_figure(fig_strength, 'max_strength', Monkey); %Used to be 'peak_strength` but that was confusing given other findpeak options with similar name
    save_directional_analysis_figure(fig_Fscore, 'fscore', Monkey);
    save_directional_analysis_figure(fig_std, 'std', Monkey);
    save_directional_analysis_figure(fig_skewness, 'skewness', Monkey);
    save_directional_analysis_figure(fig_kurtosis, 'kurtosis', Monkey);
    save_directional_analysis_figure(fig_laterality, 'laterality', Monkey);
end

%% Extract data and analyze it from LIP and MIP
% Define smoothing filter
if disk_filter_size > 0
    disk_filter = fspecial('disk', disk_filter_size);
end

% Define time region of interest.
window_of_interest = memory_end;
baseline = fix_start_end;

index_of_interest = round(window_of_interest - baseline(1));
% Handle single timepoint index
if isscalar(index_of_interest)
    index_of_interest = [index_of_interest index_of_interest];
end

% Initialize variables
[LIPd_preferred_direction, LIPd_peak_tuning_time, ...
    LIPv_preferred_direction, LIPv_peak_tuning_time, ...
    MIPd_preferred_direction, MIPd_peak_tuning_time, ...
    MIPv_preferred_direction, MIPv_peak_tuning_time] = deal(NaN(length(ap_planes), 2));

switch Monkey
    case 'L'
        plane_positions = [-1.667:6.667/4:5];
        plane_positions = ap_planes;
    case 'P'
        plane_positions = linspace(-3, 5, 8);
        plane_positions = ap_planes;
end

% Create record for voxels at level of dorsal/ventral subregions ('dv')
[concatenated_voxels_theta_dv, concatenated_voxels_time_dv, region_record_dv, plane_record_dv] = deal([]);

% Create record for voxels at level of entire region ('overall')
[concatenated_voxels_theta_overall, concatenated_voxels_time_overall, ...
    concatenated_voxels_percent_depth_overall, region_record_overall, ...
    plane_record_overall, concatenated_voxels_pixel_depth_overall] = deal([]);

for plane = 1:length(ap_planes)
    % Load in slot data
    ap_plane = ap_planes(plane);
    iDop_tuningcurve_theta = iDop_tuningcurve_theta_cell{plane};
    iDop_tuningcurve_rho = iDop_tuningcurve_rho_cell{plane};
    iDop_pvalue_F = qvalueFDR_F_cell{plane};
    anatomical_angiogram = aligned_angiogram{plane};
    
    % Cut all data to the x and z size of the anatomical angiogram. Since
    % some planes anatomical data are at angles. Use any column/row which
    % has no NaNs. Can later revise to eliminate any row with at least one
    % real number if that is cleaner.
    anatomical_angiogram(anatomical_angiogram==0) = NaN;
    row_ind = all(isnan(anatomical_angiogram), 2);
    col_ind = all(isnan(anatomical_angiogram), 1);
    
    iDop_tuningcurve_theta(row_ind, :, :) = [];
    iDop_tuningcurve_theta(:, col_ind, :) = [];
    
    iDop_tuningcurve_rho(row_ind, :, :) = [];
    iDop_tuningcurve_rho(:, col_ind, :) = [];
    
    
    if any(contains(inclusive_region_names(roi_indx), 'all'))
        % Already did this under FDR correction for image subset. Only do if
        % entire image is being analyzed.
        
        iDop_pvalue_F(row_ind, :) = [];
        iDop_pvalue_F(:, col_ind) = [];
    end
    
    
    anatomical_angiogram(row_ind, :) = [];
    anatomical_angiogram(:, col_ind) = [];
    
    % cut to portion of interest
    iDop_tuningcurve_theta = iDop_tuningcurve_theta(:, :, index_of_interest(1):index_of_interest(2));
    iDop_tuningcurve_rho = iDop_tuningcurve_rho(:, :, index_of_interest(1):index_of_interest(2));
    
    [max_strength, iDop_tuningcurve_theta_maxInd] = max(iDop_tuningcurve_rho, [], 3);
    
    iDop_tuningcurve_theta_max = zeros(size(iDop_tuningcurve_theta, 1), size(iDop_tuningcurve_theta, 2));
    
    for yrange = 1:size(iDop_tuningcurve_theta_maxInd, 1)
        for xrange = 1:size(iDop_tuningcurve_theta_maxInd, 2)
            iDop_tuningcurve_theta_max(yrange, xrange) = iDop_tuningcurve_theta(yrange, xrange, iDop_tuningcurve_theta_maxInd(yrange, xrange));
        end
    end
    
    max_time = iDop_tuningcurve_theta_maxInd;
    
    % Apply mask where we ignore voxels with pvalue for any of the
    % directions that is less than pvalue_threshold.
    F_pvalue_ind = iDop_pvalue_F < pvalue_threshold;
    
    
    statistical_threshold = F_pvalue_ind;
    
    iDop_tuningcurve_theta_max(~statistical_threshold) = NaN;
    max_strength(~statistical_threshold) = NaN;
    max_time(~statistical_threshold) = NaN;
    
    if disk_filter_size > 0
        iDop_tuningcurve_theta_max_smoothed = imfilter(iDop_tuningcurve_theta_max, disk_filter);
        max_time_smoothed = imfilter(max_time, disk_filter);
        max_strength_smoothed = imfilter(max_strength, disk_filter);
    else
        iDop_tuningcurve_theta_max_smoothed = iDop_tuningcurve_theta_max;
        max_time_smoothed = max_time;
        max_strength_smoothed = max_strength;
    end
    
    
    xPixels = size(max_time_smoothed, 2);
    yPixels = size(max_time_smoothed,1);
    
    
    % Generate mask of voxels for each mask
    LIP_ROI_mask = LIP_mask{plane};
    MIP_ROI_mask = MIP_mask{plane};
    
    LIP_combined_mask = LIP_ROI_mask & statistical_threshold;
    MIP_combined_mask = MIP_ROI_mask & statistical_threshold;
    
    % Find depth of every (statistically significant) voxel within LIP and
    % MIP
    % Two approaches.
    % 1) Absolute depth from top of LIP (in mm).
    % 2) Relative depth within LIP (percent).
    
    % Create mask without NaN values included
    nan_mask = isnan(iDop_tuningcurve_theta_max_smoothed);
    LIP_combined_mask = LIP_combined_mask & ~nan_mask;
    MIP_combined_mask = MIP_combined_mask & ~nan_mask;
    
    if nan_mask ~= isnan(max_time_smoothed)
        error('There is a mismatch in the number and location of NaNs. Check code.');
    end
    
    % Get relative depth of every voxel
    [LIP_percent_depth, LIP_pixel_depth] = get_percent_depth(LIP_combined_mask);
    [MIP_percent_depth, MIP_pixel_depth] = get_percent_depth(MIP_combined_mask);
    
    
    % Label each voxel's region and plane for later analysis
    % LIP
    significant_voxel_num_LIP = nnz(LIP_combined_mask);
    concatenated_voxels_theta_overall = [concatenated_voxels_theta_overall; iDop_tuningcurve_theta_max_smoothed(LIP_combined_mask)];
    concatenated_voxels_time_overall = [concatenated_voxels_time_overall; max_time_smoothed(LIP_combined_mask)];
    concatenated_voxels_percent_depth_overall = [concatenated_voxels_percent_depth_overall; LIP_percent_depth];
    concatenated_voxels_pixel_depth_overall = [concatenated_voxels_pixel_depth_overall; LIP_pixel_depth];
    region_record_overall = [region_record_overall; repmat('LIP', significant_voxel_num_LIP, 1)];
    plane_record_overall = [plane_record_overall; repmat(plane, significant_voxel_num_LIP, 1)];
    
    % MIP
    significant_voxel_num_MIP = nnz(MIP_combined_mask);
    concatenated_voxels_theta_overall = [concatenated_voxels_theta_overall; iDop_tuningcurve_theta_max_smoothed(MIP_combined_mask)];
    concatenated_voxels_time_overall = [concatenated_voxels_time_overall; max_time_smoothed(MIP_combined_mask)];
    concatenated_voxels_percent_depth_overall = [concatenated_voxels_percent_depth_overall; MIP_percent_depth];
    concatenated_voxels_pixel_depth_overall = [concatenated_voxels_pixel_depth_overall; MIP_pixel_depth];
    region_record_overall = [region_record_overall; repmat('MIP', significant_voxel_num_MIP, 1)];
    plane_record_overall = [plane_record_overall; repmat(plane, significant_voxel_num_MIP, 1)];  
end

%% Create inds for different regions
% For the 'overall' record
LIP_ind = all(region_record_overall == 'LIP', 2);
MIP_ind = all(region_record_overall == 'MIP', 2);


%% Bee swarm plot for relative depth within LIP and MIP
switch Monkey 
    case 'L'
        EBZ_mm = ap_planes * -5/3 + 9*5/3;
        x_ticks = -5:5/3:5/3;
    case 'P'
        EBZ_mm = ap_planes * -5/3 + 12;
        x_ticks = -8:5/3:2;
end

% use a GUI to select which planes you want
[indx,tf] = listdlg('PromptString','Select Planes to use:',...
    'SelectionMode','multiple',...
    'ListString',string(ap_planes), ...
    'ListSize', [600 300]);

plane_ind = ismember(plane_record_overall, indx);

marker_size = 20;

figure;
tld = tiledlayout(2, 2);

% For Direction
% LIP
nexttile(1);


swarm_LIP_theta = swarmchart(EBZ_mm(plane_record_overall(LIP_ind & plane_ind)), 180/pi*concatenated_voxels_theta_overall(LIP_ind & plane_ind), ...
    marker_size, ...
    concatenated_voxels_pixel_depth_overall(LIP_ind & plane_ind)*pixelsize,...
    'filled',...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
% Add colorbar
clrbar_LIP_theta = colorbar;
clrbar_LIP_theta.Label.String = 'Depth (mm)';
clim([0 11]);
xlimits = xlim;

% Add line of best fit
x = EBZ_mm(plane_record_overall(LIP_ind & plane_ind));
y = 180/pi*concatenated_voxels_theta_overall(LIP_ind & plane_ind);
mdl = fitlm(x,y);
mdl
anova(mdl,'summary')

coefficients = polyfit(x, y, 1);
xFit = linspace(min(x)-1.66, max(x)+1.66, 1000);
yFit = polyval(coefficients, xFit);
hold on;
plot(xFit, yFit, 'r-');
hold off;
fprintf('Line of best fit (p=%0.4f for slope)\n y (deg) = %0.2fx (mm) + %0.2f (deg)\n', mdl.Coefficients.pValue(2), coefficients)


% Add title and axe labels
title('Directional preference across planes - LIP');
ylabel('Directional preference (rad)');
xlabel('Position relative to EBZ (mm)');

yticks(180/pi * [-pi:pi/4:pi]);
xlim(xlimits);
xticks(x_ticks)


% scale to standardized size
ylim(180/pi*[-pi pi]);

% MIP
nexttile();
swarm_MIP_theta = swarmchart(EBZ_mm(plane_record_overall(MIP_ind & plane_ind)), 180/pi*concatenated_voxels_theta_overall(MIP_ind & plane_ind), ...
    marker_size, ...
    concatenated_voxels_percent_depth_overall(MIP_ind & plane_ind),...
    'filled',...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);


% Add title and axe labels
title('Directional preference across planes - MIP');
ylabel('Directional preference (rad)');
xlabel('Position relative to EBZ (mm)');

% Add colorbar
clrbar_MIP_theta = colorbar;
clrbar_MIP_theta.Label.String = 'Percent depth';

yticks(180/pi * [-pi:pi/4:pi]);
xticks(x_ticks);


% scale to standardized size
ylim(180/pi*[-pi pi]);


%% Histogram across all planes for distribution of RFs in LIP
figure;

bin_edges = [-180:5:180];

nexttile();
histogram(180/pi*concatenated_voxels_theta_overall(LIP_ind), bin_edges, 'Normalization', 'pdf');
xlim([-180 180]);
title("Distribution of LIP voxel's peak tuning");
ylabel('Probability');
xlabel('Degree');


% Fit distribution
pd = fitdist(180/pi*concatenated_voxels_theta_overall(LIP_ind), 'Normal');

pdf_values = pdf(pd, bin_edges);
hold on;
plot(bin_edges, pdf_values);
hold off;

% Add patch to show expected preferred range between -90 and +90
ylimits = ylim;
hContralateral = patch([-90 -90 ...
    90 90], ...
    [ylimits(1) ylimits(2) ylimits(2) ylimits(1)], 'k', ...
    'facecolor', 'k', 'edgecolor', 'none');
hContralateral.FaceAlpha = 0.2;

contralateral_preference_ind = (-pi/2 < concatenated_voxels_theta_overall(LIP_ind)) & (concatenated_voxels_theta_overall(LIP_ind) < pi/2);
fprintf('Percent of tuned voxels with peak tuning within [-90 90] degrees: %0.1f%%\n', 100*nnz(contralateral_preference_ind)/nnz(LIP_ind));


%% Histogram across all planes for depth of activated voxels
figure;
bin_edges = [0:0.04:1];

histogram(concatenated_voxels_percent_depth_overall(LIP_ind), bin_edges, 'Normalization', 'pdf');
xlabel('Percent depth in sulcus (<0.47 is dorsal LIP)');

% Fit distribution
pd = fitdist(concatenated_voxels_percent_depth_overall(LIP_ind), 'Normal');

pdf_values = pdf(pd, bin_edges);
hold on;
plot(bin_edges, pdf_values);

% Draw line to represent LIP d/v boundary
ylimits = ylim;
empirical_boundary = 0.47;
hold on;
line([empirical_boundary empirical_boundary], ...
    ylimits);
hold off;
xlim([0 1]);



%% Beeswarm with color and y-axis swapped
% Using absolute depth

marker_size = 10;

figure;
tld = tiledlayout(1, 2);

% For Direction
% LIP
nexttile();

jitter_width = 2 * min(diff(unique(ap_planes(plane_record_overall(LIP_ind)))));

swarm_LIP_theta = swarmchart(ap_planes(plane_record_overall(LIP_ind)), concatenated_voxels_pixel_depth_overall(LIP_ind)*pixelsize, ...
    marker_size, ...
    concatenated_voxels_theta_overall(LIP_ind),...
    'filled', ...
    'MarkerFaceAlpha', 0.5, ...
    'MarkerEdgeAlpha', 0.5, ...
    'XJitterWidth', jitter_width);
colormap(circle_colormap);
% Add colorbar
clrbar_LIP_theta = colorbar;
clrbar_LIP_theta.Label.String = 'Directional preference (rad)';


% Add title and axe labels
title('Directional preference across planes - LIP');
ylabel('Depth (mm)');
xlabel('Plane');
set(gca, 'YDir','reverse')

% For Direction
% LIP - As percent depth
nexttile();

jitter_width = 2 * min(diff(unique(ap_planes(plane_record_overall(LIP_ind)))));

swarm_LIP_theta = swarmchart(ap_planes(plane_record_overall(LIP_ind)), concatenated_voxels_percent_depth_overall(LIP_ind), ...
    marker_size, ...
    concatenated_voxels_theta_overall(LIP_ind),...
    'filled', ...
    'MarkerFaceAlpha', 0.5, ...
    'MarkerEdgeAlpha', 0.5, ...
    'XJitterWidth', jitter_width);
colormap(circle_colormap);
% Add colorbar
clrbar_LIP_theta = colorbar;
clrbar_LIP_theta.Label.String = 'Directional preference (rad)';

% Add line dividing dorsal/ventral LIP
hold on;
xlimits = xlim;
line(xlimits, [0.47 0.47], 'LineStyle', '--');
hold off;
set(gca, 'YDir','reverse')


% Add title and axe labels
title('Directional preference across planes - LIP');
ylabel('Depth (%)');
xlabel('Plane');

%% 1D lines for each plane that represent tuning and depth
% X axis - Plane
% Y axis - Relative Depth (percent)
% color - Tuning

% How many bins do we want to split depth into?
num_bins = 50;

% Define smoothing filter
disk_filter = fspecial('disk', 1);

% Define time region of interest.
window_of_interest = memory_end;
baseline = fix_start_end;

index_of_interest = round(window_of_interest - baseline(1));
% Handle single timepoint index
if numel(index_of_interest) == 1
    index_of_interest = [index_of_interest index_of_interest];
end

[LIP_tuning_1D_in_percent, MIP_tuning_1D_in_percent] = deal(NaN(num_bins, length(ap_planes)));
bins_in_percent = linspace(0, 1, num_bins+1);

for plane = 1:length(ap_planes)
    % Load in slot data
    ap_plane = ap_planes(plane);
    iDop_tuningcurve_theta = iDop_tuningcurve_theta_cell{plane};
    iDop_tuningcurve_rho = iDop_tuningcurve_rho_cell{plane};
    iDop_pvalue_8dir = qvalueFDR_F_cell{plane};
    anatomical_angiogram = aligned_angiogram{plane};
    
    % Cut all data to the x and z size of the anatomical angiogram.
    row_ind = rows_cols_to_cut{plane, 1};
    col_ind = rows_cols_to_cut{plane, 2};
    iDop_tuningcurve_theta(row_ind, :, :) = [];
    iDop_tuningcurve_theta(:, col_ind, :) = [];
    iDop_tuningcurve_rho(row_ind, :,:) = [];
    iDop_tuningcurve_rho(:, col_ind, :) = [];
    anatomical_angiogram(row_ind, :) = [];
    anatomical_angiogram(:, col_ind) = [];
    iDop_pvalue_8dir(row_ind, :) = [];
    iDop_pvalue_8dir(:, col_ind) = [];
    
    % cut to portion of interest
    iDop_tuningcurve_theta = iDop_tuningcurve_theta(:, :, index_of_interest(1):index_of_interest(2));
    iDop_tuningcurve_rho = iDop_tuningcurve_rho(:, :, index_of_interest(1):index_of_interest(2));
    
    [max_strength, iDop_tuningcurve_theta_maxInd] = max(iDop_tuningcurve_rho, [], 3);
    
    iDop_tuningcurve_theta_max = zeros(size(iDop_tuningcurve_theta, 1), size(iDop_tuningcurve_theta, 2));
    
    for yrange = 1:size(iDop_tuningcurve_theta_maxInd, 1)
        for xrange = 1:size(iDop_tuningcurve_theta_maxInd, 2)
            iDop_tuningcurve_theta_max(yrange, xrange) = iDop_tuningcurve_theta(yrange, xrange, iDop_tuningcurve_theta_maxInd(yrange, xrange));
        end
    end
    
    max_time = iDop_tuningcurve_theta_maxInd;
    
    % Apply mask where we ignore voxels with pvalue for any of the
    % directions that is less than pvalue_threshold.
    pvalue_mask = iDop_pvalue_8dir < pvalue_threshold;
    pvalue_ind = squeeze(any(pvalue_mask, 4));
    
    F_pvalue_ind = iDop_pvalue_F < pvalue_threshold;
    
    statistical_threshold = pvalue_ind;
    
    iDop_tuningcurve_theta_max(~statistical_threshold) = NaN;
    max_strength(~statistical_threshold) = NaN;
    max_time(~statistical_threshold) = NaN;
    
    iDop_tuningcurve_theta_max_smoothed = imfilter(iDop_tuningcurve_theta_max, disk_filter);
    max_time_smoothed = imfilter(max_time, disk_filter);
    max_strength_smoothed = imfilter(max_strength, disk_filter);
    
    
    xPixels = size(max_time_smoothed, 2);
    yPixels = size(max_time_smoothed,1);
    
    
    % Generate mask of voxels for each mask
    LIP_ROI_mask = LIP_mask{plane};
    MIP_ROI_mask = MIP_mask{plane};
    
    % Combine masks with pvalue thresholding
    LIP_combined_mask = LIP_ROI_mask & statistical_threshold;
    MIP_combined_mask = MIP_ROI_mask & statistical_threshold;
    
    
    % Find depth of every (statistically significant) voxel within LIP and
    % MIP
    % Two possible approaches.
    % 1) Absolute depth from top of LIP (in mm).
    % 2) Relative depth within LIP (percent).
    % Using method 2 because it normalizes across different animals and
    % planes. Can easily modify this to use Method 1 though.
    
    % Create mask without NaN values included
    nan_mask = isnan(iDop_tuningcurve_theta_max_smoothed);
    LIP_combined_mask = LIP_combined_mask & ~nan_mask;
    MIP_combined_mask = MIP_combined_mask & ~nan_mask;
    
    if nan_mask ~= isnan(max_time_smoothed)
        error('There is a mismatch in the number and location of NaNs. Check code.');
    end
    
    % Get relative depth of every voxel
    LIP_percent_depth = get_percent_depth(LIP_combined_mask);
    MIP_percent_depth = get_percent_depth(MIP_combined_mask);
    
    LIP_tuning = iDop_tuningcurve_theta_max_smoothed(LIP_combined_mask);
    MIP_tuning = iDop_tuningcurve_theta_max_smoothed(MIP_combined_mask);
    
    LIP_strength = max_strength_smoothed(LIP_combined_mask);
    MIP_strength = max_strength_smoothed(MIP_combined_mask);
    
    for a = 1:num_bins
        voxel_ind = LIP_percent_depth > bins_in_percent(a) & LIP_percent_depth < bins_in_percent(a+1);
        LIP_tuning_1D_in_percent(a, plane) = mean(LIP_tuning(voxel_ind), 'all');
        LIP_strength_1D_in_percent(a, plane) = mean(LIP_strength(voxel_ind), 'all');
        
        voxel_ind = MIP_percent_depth > bins_in_percent(a) & MIP_percent_depth < bins_in_percent(a+1);
        MIP_tuning_1D_in_percent(a, plane) = mean(MIP_tuning(voxel_ind), 'all');
        MIP_strength_1D_in_percent(a, plane) = mean(MIP_strength(voxel_ind), 'all');
        
    end
end

% Normalize to the strongest
LIP_strength_1D_norm_in_percent = LIP_strength_1D_in_percent./max(LIP_strength_1D_in_percent, [], 'all');
MIP_strength_1D_norm_in_percent = MIP_strength_1D_in_percent./max(MIP_strength_1D_in_percent, [], 'all');


figure;
circle_colormap_to_use = circle_colormap;
tld = tiledlayout(1, 2);
nexttile()
imAlpha=LIP_strength_1D_norm_in_percent;
imAlpha(isnan(LIP_tuning_1D_in_percent))=0;
ax1 = imagesc(ap_planes, bins_in_percent, LIP_tuning_1D_in_percent, 'AlphaData',imAlpha);
clim([-pi pi]);
colormap(circle_colormap_to_use);
%c = colorbar;
%c.Label.String = ('Preferred direction (rad)');
ylabel('Percent Depth (%)');
xlabel('Plane number');

nexttile();
transparency_axis = max(LIP_strength_1D_in_percent, [], 'all');

% Generate figure showing saturation combined with circular colormap
number_saturation_elements = 1000;
number_colors = size(circle_colormap_to_use, 1);

circle_colormap_repmat = repmat(flipud(circle_colormap_to_use), 1, 1, number_saturation_elements);
circle_colormap_repmat = permute(circle_colormap_repmat, [1 3 2]);
alpha_map = repmat(linspace(0, 1, number_saturation_elements), number_colors, 1, 1);

im = image( linspace(0, transparency_axis, number_saturation_elements), linspace(-pi, pi, number_colors), circle_colormap_repmat);
im.AlphaData = alpha_map;
xlabel(sprintf('Effect size (%s)', anova_effectsize_measure));
ylabel('Direction (rad)');
axis square
yticks([-pi, -pi/2, 0, pi/2, pi])
yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'}); % Check this is correct. There was something weird
set(gca,'YAxisLocation','right')



%% 1D lines for each plane that represent tuning and depth
% X axis - Plane
% Y axis - Absolute depth (mm)
% color - Tuning


% How many bins do we want to split depth into?
num_bins = 50;


% Define smoothing filter
if disk_filter_size > 0
    disk_filter = fspecial('disk', 1);
end

% Define time region of interest.
window_of_interest = memory_end;
baseline = fix_start_end;

index_of_interest = round(window_of_interest - baseline(1));
% Handle single timepoint index
if isscalar(index_of_interest)
    index_of_interest = [index_of_interest index_of_interest];
end

[LIP_tuning_1D_in_mm, MIP_tuning_1D_in_mm] = deal(NaN(num_bins, length(ap_planes)));
bins_in_mm = linspace(0, 12, num_bins+1);

for plane = 1:length(ap_planes)
    % Load in slot data
    ap_plane = ap_planes(plane);
    iDop_tuningcurve_theta = iDop_tuningcurve_theta_cell{plane};
    iDop_tuningcurve_rho = iDop_tuningcurve_rho_cell{plane};
    iDop_pvalue_8dir = qvalueFDR_F_cell{plane};
    anatomical_angiogram = aligned_angiogram{plane};
    
    
    % Cut all data to the x and z size of the anatomical angiogram.
    row_ind = rows_cols_to_cut{plane, 1};
    col_ind = rows_cols_to_cut{plane, 2};
    iDop_tuningcurve_theta(row_ind, :, :) = [];
    iDop_tuningcurve_theta(:, col_ind, :) = [];
    iDop_tuningcurve_rho(row_ind, :,:) = [];
    iDop_tuningcurve_rho(:, col_ind, :) = [];
    anatomical_angiogram(row_ind, :) = [];
    anatomical_angiogram(:, col_ind) = [];
    iDop_pvalue_8dir(row_ind, :) = [];
    iDop_pvalue_8dir(:, col_ind) = [];
    
    % cut to portion of interest
    iDop_tuningcurve_theta = iDop_tuningcurve_theta(:, :, index_of_interest(1):index_of_interest(2));
    iDop_tuningcurve_rho = iDop_tuningcurve_rho(:, :, index_of_interest(1):index_of_interest(2));
    
    [max_strength, iDop_tuningcurve_theta_maxInd] = max(iDop_tuningcurve_rho, [], 3);
    
    iDop_tuningcurve_theta_max = zeros(size(iDop_tuningcurve_theta, 1), size(iDop_tuningcurve_theta, 2));
    
    for yrange = 1:size(iDop_tuningcurve_theta_maxInd, 1)
        for xrange = 1:size(iDop_tuningcurve_theta_maxInd, 2)
            iDop_tuningcurve_theta_max(yrange, xrange) = iDop_tuningcurve_theta(yrange, xrange, iDop_tuningcurve_theta_maxInd(yrange, xrange));
        end
    end
    
    max_time = iDop_tuningcurve_theta_maxInd;
    
    % Apply mask where we ignore voxels with pvalue for any of the
    % directions that is less than pvalue_threshold.
    pvalue_mask = iDop_pvalue_8dir < pvalue_threshold;
    pvalue_ind = squeeze(any(pvalue_mask, 4));
    
    F_pvalue_ind = iDop_pvalue_F < pvalue_threshold;
    
    statistical_threshold = pvalue_ind;
    
    iDop_tuningcurve_theta_max(~statistical_threshold) = NaN;
    max_strength(~statistical_threshold) = NaN;
    max_time(~statistical_threshold) = NaN;
    
    if disk_filter_size > 0
        iDop_tuningcurve_theta_max_smoothed = imfilter(iDop_tuningcurve_theta_max, disk_filter);
        max_time_smoothed = imfilter(max_time, disk_filter);
        max_strength_smoothed = imfilter(max_strength, disk_filter);
    else
        iDop_tuningcurve_theta_max_smoothed = iDop_tuningcurve_theta_max;
        max_time_smoothed = max_time;
        max_strength_smoothed = max_strength;
    end
    
    
    xPixels = size(max_time_smoothed, 2);
    yPixels = size(max_time_smoothed,1);
    
    
    % Generate mask of voxels for each mask
    LIP_ROI_mask = LIP_mask{plane};
    MIP_ROI_mask = MIP_mask{plane};
    
    % Combine masks with pvalue thresholding
    LIP_combined_mask = LIP_ROI_mask & statistical_threshold;
    MIP_combined_mask = MIP_ROI_mask & statistical_threshold;
    
    
    % Find depth of every (statistically significant) voxel within LIP and
    % MIP
    % Two approaches.
    % 1) Absolute depth from top of LIP (in mm).
    % 2) Relative depth within LIP (percent).
    
    % Create mask without NaN values included
    nan_mask = isnan(iDop_tuningcurve_theta_max_smoothed);
    LIP_combined_mask = LIP_combined_mask & ~nan_mask;
    MIP_combined_mask = MIP_combined_mask & ~nan_mask;
    
    if nan_mask ~= isnan(max_time_smoothed)
        error('There is a mismatch in the number and location of NaNs. Check code.');
    end
    
    % Get relative depth of every voxel
    [LIP_percent_depth, LIP_pixel_depth] = get_percent_depth(LIP_combined_mask);
    [MIP_percent_depth, MIP_pixel_depth] = get_percent_depth(MIP_combined_mask);
    
    LIP_mm_depth = LIP_pixel_depth*pixelsize;
    MIP_mm_depth = MIP_pixel_depth*pixelsize;
    
    LIP_tuning = iDop_tuningcurve_theta_max_smoothed(LIP_combined_mask);
    MIP_tuning = iDop_tuningcurve_theta_max_smoothed(MIP_combined_mask);
    
    LIP_strength = max_strength_smoothed(LIP_combined_mask);
    MIP_strength = max_strength_smoothed(MIP_combined_mask);
    
    
    for a = 1:num_bins
        voxel_ind = LIP_mm_depth > bins_in_mm(a) & LIP_mm_depth < bins_in_mm(a+1);
        LIP_tuning_1D_in_mm(a, plane) = mean(LIP_tuning(voxel_ind), 'all');
        LIP_strength_1D_in_mm(a, plane) = mean(LIP_strength(voxel_ind), 'all');
        
        voxel_ind = MIP_mm_depth > bins_in_mm(a) & MIP_mm_depth < bins_in_mm(a+1);
        MIP_tuning_1D_in_mm(a, plane) = mean(MIP_tuning(voxel_ind), 'all');
        MIP_strength_1D_in_mm(a, plane) = mean(MIP_strength(voxel_ind), 'all');
        
    end
end

% Normalize to the strongest
LIP_strength_1D_norm_in_mm = LIP_strength_1D_in_mm./max(LIP_strength_1D_in_mm, [], 'all');
MIP_strength_1D_norm_in_mm = MIP_strength_1D_in_mm./max(MIP_strength_1D_in_mm, [], 'all');


figure;
circle_colormap_to_use = circle_colormap;
tld = tiledlayout(1, 2);
nexttile()
imAlpha=LIP_strength_1D_norm_in_mm;
imAlpha(isnan(LIP_tuning_1D_in_mm))=0;
ax1 = imagesc(ap_planes, bins_in_mm, LIP_tuning_1D_in_mm, 'AlphaData',imAlpha);
clim([-pi pi]);
colormap(circle_colormap_to_use);
%c = colorbar;
%c.Label.String = ('Preferred direction (rad)');
ylabel('Depth (mm)');
xlabel('Plane number');

nexttile();
transparency_axis = max(LIP_strength_1D_in_percent, [], 'all');

% Generate figure showing saturation combined with circular colormap
number_saturation_elements = 1000;
number_colors = size(circle_colormap_to_use, 1);

circle_colormap_repmat = repmat(flipud(circle_colormap_to_use), 1, 1, number_saturation_elements);
circle_colormap_repmat = permute(circle_colormap_repmat, [1 3 2]);
alpha_map = repmat(linspace(0, 1, number_saturation_elements), number_colors, 1, 1);

im = image( linspace(0, transparency_axis, number_saturation_elements), linspace(-pi, pi, number_colors), circle_colormap_repmat);
im.AlphaData = alpha_map;
xlabel(sprintf('Effect size (%s)', anova_effectsize_measure));
ylabel('Direction (rad)');
axis square
yticks([-pi, -pi/2, 0, pi/2, pi])
yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'}); % Check this is correct. There was something weird
set(gca,'YAxisLocation','right')



%% Use approximate slice thickness

% use a GUI to select which planes you want
[indx,tf] = listdlg('PromptString','Select Planes to use:',...
    'SelectionMode','multiple',...
    'ListString',string(ap_planes), ...
    'ListSize', [600 300]);

planes_to_use = ap_planes(indx);

relative_or_absolute_depth = 'absolute';

switch relative_or_absolute_depth
    case 'absolute'
        tuning_matrix_to_use = LIP_tuning_1D_in_mm;
        strength_matrix_to_use = LIP_strength_1D_norm_in_mm;
        bins_to_use = bins_in_mm;
        depth_label = 'Depth (mm)';
    case 'relative'
        tuning_matrix_to_use = LIP_tuning_1D_in_percent;
        strength_matrix_to_use = LIP_strength_1D_norm_in_percent;
        bins_to_use = bins_in_percent;
        depth_label = 'Depth (%)';
end


switch Monkey
    case 'P'
        ap_start = 3.66;
    case 'L'
        ap_start = 1.66;
end


ap_planes_mm = ap_planes*5/3; % in mm now
ap_planes_mm = ap_start - (ap_planes_mm - min(ap_planes_mm));

ap_range = range(ap_planes_mm) + 0.4; %Plane of view is 200 um thicker on both ends


plane_edges = NaN(length(ap_planes), 2);
for plane = 1:length(planes_to_use)
    current_plane = planes_to_use(plane);
    plane_indx = ap_planes == current_plane;
    plane_edges(plane_indx, :) = [ap_planes_mm(plane_indx)-0.2 ap_planes_mm(plane_indx)+0.2];
end

ap_resolution = 500;
ap_subsample = linspace(min(plane_edges, [], 'all'), max(plane_edges, [], 'all'), ap_resolution);

[LIP_strength_adjusted, LIP_tuning_adjusted] = deal(NaN(num_bins, ap_resolution));

for i = 1:length(ap_subsample)
    current_ap_position = ap_subsample(i);
    plane_ind = plane_edges(:, 1) <= current_ap_position & plane_edges(:, 2) >= current_ap_position;
    if nnz(plane_ind)
        LIP_tuning_adjusted(:, i) = mean(tuning_matrix_to_use(:, plane_ind), 2, 'omitnan');
        LIP_strength_adjusted(:, i) = mean(strength_matrix_to_use(:, plane_ind), 2, 'omitnan');
    end
end

% Plot the adjusted image
figure;
imagesc(ap_subsample, bins_in_mm, LIP_tuning_adjusted, 'AlphaData', LIP_strength_adjusted);
xlabel('Plane position relative to EBZ');
ylabel(depth_label);
clim([-pi pi]);
colormap(circle_colormap_to_use);
clb = colorbar;
clb.Ticks = -pi:pi/4:pi;