%% Analyze_8tgt_withGLM
% Calculate pixelwise metrics, such as tuning strength, laterality,
% circular standard deviation, kurtosis, skewness
%
% Uses previously generated GLM results for statistical masking

%% User settings
clear; close all; clc;

% Project record name
project_record_filename = 'ProjectRecord_paper.json';

% Specify the pvalue to threshold maps at
pvalue_threshold = 1e-5;
pvalue_glm_threshold = pvalue_threshold;

% Motion correction method
mc_method = 'normcorre';

if mc_method
    mc_method_string = mc_method;
else
    mc_method_string = 'false';
end

% Size of patch to use for certain calculations . Smaller takes longer,
% but has a more fine grained resolution.
e = 1;

% Specify the names of the 8 different directions
TrialsToUse = {'0 deg','45 deg','90 deg','135 deg','180 deg','- 135 deg','- 90 deg','- 45 deg'};

% Effect size measure
effect_size_measure_pairwise = 'Cohen d'; % Another option is 'percent change'

% Do we want to save the output
savePlots = false;
savePath_base = fullfile(get_user_data_path('pathType', 'glm'));

%% Load specific colormaps
PositiveNegativeColormap = load('PercentChangeColormap.mat');
PositiveNegativeColormap = PositiveNegativeColormap.PercentChangeColormap;

%% Load data
% Load trial-aligned iDop data
loadDopplerData(whos,'multiple', 'false', ...
    'mc_method', mc_method, ...
    'project_record_filename', project_record_filename);

%% Load in GLM masks for each coronal plane
start_dir = pwd;

% Logic for filenaming
% Complete filename up until the date generated
partial_filename = sprintf('glm_data_%sgen*.mat', sprintf('S%dR%d_', sessionRunList'));

cd(fullfile(get_user_data_path('pathType', 'glm')));

[filename, filepath] = uigetfile(partial_filename, 'Select the GLM data file that matches this session.');

if filepath
    % Load data
    load(fullfile(filepath, filename), 'F', 'pvalues_corrected', ...
        'angiogram', 'sessionRunList', 'analysis_parameters');

else
    error('No GLM data file was chosen. Please try again or create the necessary file using a GLM analysis script.');
end
% Rename the angiogram variable to avoid later conflict and improve clarity
glm_angiogram = angiogram;

cd(start_dir);

%% Preprocess data

% Define parameters for the Gaussian smoothing kernel
if isfield(analysis_parameters, 'FWHM')
    FWHM = analysis_parameters.FWHM; % In voxels;
else
    FWHM = 1;
    warning('default FWHM is %d - make sure this is what you wanted', FWHM);
end

% Define the gaussian filter
if FWHM > 0
    filter_type = 'gaussian';
    sigma = FWHM/(sqrt(8*log(2)));
    filter_size = 2*ceil(2*sigma)+1; % This is the default filter size used in `imgaussfilt`.
else
    filter_type = [];
    sigma = [];
    filter_size = [];
end


% Apply spatial filter to entire session.
iDopP = preProcess(iDop, 'spatialFilter', {filter_type, filter_size, sigma}, 'zscore', false);

if size(sessionRunList, 1) > 1
    session_names = '';
    for i = 1:size(sessionRunList, 1)
        session_names = [session_names sprintf('%d_', sessionRunList(i, 1))];
    end
else
    session_names = sprintf('%d', sessionRunList(1, 1));
end

% Rename the angiogram variable to avoid later conflict and improve clarity
dop_angiogram = angiogram;

% pixels-->mm
pixelsize = 0.1;
X_img_mm = pixelsize/2 + (0:size(iDopP,2)-1)*pixelsize; % So that the center of the image is at true halfway point, not shifted.
Z_img_mm = pixelsize/2 + (0:size(iDopP,1)-1)*pixelsize + UF.Depth(1);

% Get size of image
[yPix, xPix, nWindows, nTrials] = size(iDopP);


%% Ideally the GLM mask and the loaded data already match.
% Let's make sure they match
% Align the GLM mask to match the loaded data.

% Use manual alignment to align glm with the dop data
alignment_tform = align_angiogram(dop_angiogram, glm_angiogram);

%Convert transform into matlab format
% Alignment_tform is passed back from ManualImageAlignment
tform = affine2d(alignment_tform);

% Align image to match previously loaded data
glm_angiogram_registered = imwarp(...
    glm_angiogram, ...
    tform, ...
    'OutputView', imref2d(size(dop_angiogram)));
F_registered = imwarp(...
    F, ...
    tform, ...
    'OutputView', imref2d(size(dop_angiogram)), ...
    'FillValues', NaN);
p_cor_registered = imwarp(...
    pvalues_corrected, ...
    tform, ...
    'OutputView', imref2d(size(dop_angiogram)), ...
    'FillValues', NaN);

%% Parsing behavioral data for trial-averaged data
% Limit trial-averaged analysis to successful trials only
success = cell2mat({behavior.success});
behavior(~success) = [];

% Extract timing information about fixation, memory/cue, and movement
% periods.
[fix, mem, mov] = getEpochs(behavior);

% Specify time periods of interest
baselinePeriod = [-1 1];
window_halfsize = 0;

% Extract information about targets
target = cell2mat({behavior.target});

% Extract information about target position
targetPos = vertcat(behavior.targetPos);

if iscell(targetPos)
    targetPos = cell2mat(targetPos);
end

% Convert from cartesian coordinates to polar coordinates
[angle,~] = cart2pol(targetPos(:,1),targetPos(:,2)); % convert to polar coordinates
angle(angle<0) = angle(angle<0)+2*pi; %Convert to have only positive angles.

% Find unique directions seen.
UniqueAngles = unique(angle);
UniqueAnglesInDegrees = UniqueAngles * 360/(2*pi);

% create an index showing which trials are for which directions.
TargetPosInd = zeros(length(target),1);
for angleInd = 1:length(UniqueAngles)
    TargetPosInd(angle==UniqueAngles(angleInd)) = angleInd;
end

%% Computing statistical map of pixels responsive to task
% With the trial-aligned data

% Pre-allocate space for arrays
[pvalue_dir, tscore_dir, effectsize_dir] = deal(NaN(...
    yPix, ...
    xPix, ...
    length(UniqueAngles)...
    ));

% Specify multiple comparison method to use
multcompmethod = 'FDR'; %'bonferroni' or 'FDR' or 'none'

% Calculate indices for baseline and memory period
% Attempting to find the fUS indices that are closest to the bseline
% period. If using a single point as reference, find the closest. If using
% multiple points, find all points within desired range.

if baselinePeriod(1) == baselinePeriod(2)
    baseline_indices = round(baselinePeriod*coreParams.framerate - fix(1)*coreParams.framerate);
else
    baseline_indices = [(ceil(baselinePeriod(1)*coreParams.framerate - fix(1)*coreParams.framerate)):(floor(baselinePeriod(2)*coreParams.framerate - fix(1)*coreParams.framerate))];
end

% Base memStat off length of memory period
% Round to the nearest integer.
memory_end = mem(2);
memStat = [memory_end-window_halfsize memory_end+window_halfsize];

if memStat(1) == memStat(2)
    memory_indices = round(memStat*coreParams.framerate - fix(1)*coreParams.framerate);
else
    memory_indices = [(ceil(memStat(1)*coreParams.framerate - fix(1)*coreParams.framerate)):(floor(memStat(2)*coreParams.framerate - fix(1)*coreParams.framerate))];
end

% Convert into percent change
iDopP_percentchange = convert_to_percent_change(iDopP, [baseline_indices(1) baseline_indices(end)]);

% Clean up some outliers around image border if motion correction is used.
if coreParams.motionCorrection
    % This is empirically determined for now. It will vary based upon the
    % session. It should match with the approximate maximum shift/rotation
    % of the image.
    borderwidth = 2;

    % Remove outliers
    iDopP_percentchange = remove_outliers(iDopP_percentchange, 2);
end

% Calculate percentChange during memory period
PercentChange_memPeriod = squeeze(mean(iDopP_percentchange(:,:,memory_indices,:),3, 'omitnan'));
PercentChange_memPeriod_angles = NaN(yPix, xPix, length(UniqueAngles));

% Calculate percentChange for each direction during memory period
for i = 1:length(UniqueAngles)
    PercentChange_memPeriod_angles(:,:,i) =  squeeze(mean(PercentChange_memPeriod(:,:,TargetPosInd == i),3));
end

%Obtaining statistics for each unique direction
for position = 1:length(UniqueAngles)
    trial_indices = TargetPosInd==position;

    effectsize_dir(:, :, position) = calculate_voxel_patch_effect_size(...
        iDopP_percentchange(:, :, memory_indices, trial_indices),...
        iDopP_percentchange(:, :, baseline_indices, trial_indices), ...
        e, 'effect_size_measure', effect_size_measure_pairwise);
end

%% Thresholding using the GLM results
pvalue_glm_mask = p_cor_registered > pvalue_glm_threshold | isnan(p_cor_registered);

% Threshold the direction results using the GLM result
pvalue_glm_8Dmask = repmat(pvalue_glm_mask, 1, 1, length(UniqueAngles));
effectsize_thresholded_by_glm = effectsize_dir;
effectsize_thresholded_by_glm(pvalue_glm_8Dmask) = NaN;

% Threshold F-score plot by p-value
F_score_thresholded = F_registered;
F_score_thresholded(pvalue_glm_mask) = NaN;

%% Plot F-score
figure;
cmap = plotDuplexImage(X_img_mm, Z_img_mm, F_score_thresholded, angiogram,...
    'colormap2use',flipud(magma), ...
    'nonlinear_bg',2, ...
    'showColorbar', true, ...
    'ColorBarTitle', 'F-score');
title(sprintf('F-score - S%dR%d\n q < %.0e', Session, Run, pvalue_threshold));

%% Plotting
% runAlready is used in later code blocks to make changes to the plots
% generated in this code block. For a new figure, it should always be
% cleared first.
clear runAlready

% Specify which plots to show
plotESmap = true;

% For plotting the effect size maps
if plotESmap
    % The range of values to display. Anything beyond this range is set
    % to max/min colorvalues
    minValue = min(effectsize_thresholded_by_glm, [],'all', 'omitnan');
    maxValue = max(effectsize_thresholded_by_glm, [],'all', 'omitnan');

    maxAbsoluteValue = max([abs(minValue), abs(maxValue)], [], 'omitnan') * 0.8;
    esValuesToDisplay = [-maxAbsoluteValue maxAbsoluteValue];

    % Open figure for the t-map
    fig_directionalresponse=figure('Name',['S' num2str(Session) ' R' num2str(Run)  ' effect size-map' ...
        ': B: ' num2str(baselinePeriod(1)) ':' num2str(baselinePeriod(end)) ...
        ', M: ' num2str(memStat(1)) ':' num2str(memStat(end))],...
        'NumberTitle','On', ...
        'Units', 'Normalized', ...
        'OuterPosition', [0 0 1 1]);

    % Specify the plot configuration using tight_subplot
    % In future, convert to use tiledlayout instead.
    tld = tiledlayout(3, 4);
    tld.TileSpacing = 'tight';

    % Specify the order of things to plot. This will generate the
    % 8-directions all in a circle in their appropriate positions. 0-degree
    % target is at 0-degrees. 90 is at 90, -90 is at -90. etc.
    figureMapping = [7, 3, 2, 1, 5, 9, 10, 11];

    for i=1:length(UniqueAngles)

        %Specify which axes to plot to
        nexttile(figureMapping(i));

        effectsize_plot = effectsize_thresholded_by_glm(:, :, i);
        pvalue_threshold_to_save =  pvalue_glm_threshold;

        % Generate the colormap of the t-scored activation to each direction
        cmap = plotDuplexImage(X_img_mm, Z_img_mm, effectsize_plot, angiogram,...
            'colormap2use',PositiveNegativeColormap, 'nonlinear_bg',2, ...
            'AutoColorBarLimits',esValuesToDisplay);

        % Remove all axes labels
        axis off
        xlabel('');
        ylabel('');

        obj(i) = scalebar; % Put scalebar onto each subplot

        % Put scalebar text on the up left subplot. Remove scalebar text
        % from other subplots.
        if figureMapping(i) == 1
            obj(i).XUnit = 'mm';
        else
            delete(obj(i).hTextX);
        end
        axisXLim = get(gca,'XLim');
        axisYLim = get(gca,'YLim');
        axisXWidth = diff(axisXLim);
        axisYWidth = diff(axisYLim);

        % Display trial type in each subplot
        useText = false;
        if useText
            text((axisXLim(1)+axisXLim(2))/2,axisYLim(2) - 0.1*axisYWidth,TrialsToUse{i},...
                'Color','black','FontSize',18,'HorizontalAlignment','Center')
        end

        %Define colorbar to be shown
        if i==1
            % Specify the center axes for the colorbar
            nexttile(6);

            % Turn off the axes
            axis off;
            % Draw the colorbar
            c = colorbar('FontSize',18,'Location','east');
            % Specify the colormap
            colormap(gca,PositiveNegativeColormap)
            % Define the colormap limits.
            caxis(esValuesToDisplay);
            % Label the colorbar
            c.Label.String = sprintf('%s', effect_size_measure_pairwise);
        end %End colorbar
    end %End for loop over direction

    if savePlots
        filename_SPM = fullfile(savePath_base, ...
            sprintf('S%sR%d_8dirSPM_MC_%s_thresh_%0.4f_dg%s.svg', ...
            session_names, ...
            Run, ...
            mc_method_string, ...
            pvalue_threshold_to_save, ...
            datetime('today', format='yyyyMMdd')));
        saveas(fig_directionalresponse, filename_SPM, 'svg');

        % Also save as MATLAB figure
        saveas(fig_directionalresponse, strrep(filename_SPM, '.svg', '.fig'), 'fig');
    end
end %End plotESmap


%% %%%%%%%%%%%%%%%%%%%%%%%%%% Tuning Curve Analysis %%%%%%%%%%%%%%%%%%%%%%%%
sizeOfVoxel = 4;

% Use hard-coded positions or dynamic positioning
useDynamicPositioning = false;

if useDynamicPositioning
    figure(fig_directionalresponse);

    % Create a button to stop the drawpoint feature
    tb = uicontrol(gcf, 'Style', 'togglebutton',...
        'String', 'Stop updating graph');

    % Dynamically decide on point of interest
    pointer = drawpoint(nexttile(figureMapping(2)));
    MillimeterOfInterest = pointer.Position;
else
    % Specify patches of interest.
    switch Session
        case 52
            %S52R1
            MillimeterOfInterest = ...
                [6.95 10.35
                6.95 13.65
                9.15 11.75];

            sizeOfVoxel = repmat(sizeOfVoxel, size(MillimeterOfInterest));
        case 76
            %S76R1
            MillimeterOfInterest = ...
                [5.45 8.65
                6.55 11.85
                8.85 11.85];

            sizeOfVoxel = repmat(sizeOfVoxel, size(MillimeterOfInterest));
        case 102
            %S102R1
            MillimeterOfInterest = ...
                [6.35 4.85;
                6.85 8.55;
                8.85 8.65];

            sizeOfVoxel = repmat(sizeOfVoxel, size(MillimeterOfInterest));
        otherwise
            error('Not defined for this session yet.');
    end

end

% Replicate voxel size if only one size specified. This is a hack to handle
% case where just a single voxel size is specified. In this case, assume we
% want squares of equal size for each ROI.
if isscalar(sizeOfVoxel)
    sizeOfVoxel = repmat(sizeOfVoxel, size(MillimeterOfInterest, 1), 2);
end

% Specify several unique linestyles
linestyles = {'-.',':','-','--','-.', ':', '-', '--'};

% Pre-allocate space
PercentChange_ROI_timeseries = NaN(size(MillimeterOfInterest,1), nWindows, nTrials);
stePercentChange_ROI = NaN(size(MillimeterOfInterest, 1), 1);

meanPercentChange_ROI_direction = NaN(length(UniqueAngles), size(MillimeterOfInterest, 1));

% Set up the tuning curve and time series plot
if ~useDynamicPositioning
    fig_direction_ERA = figure(); cla reset;
    tld = tiledlayout('flow');

    fig_tuningcurve = figure(); cla reset;
end


while true
    % Continously update the point of interest
    if useDynamicPositioning
        MillimeterOfInterest = pointer.Position;
    end
    % For each patch of interest, plot box onto statistical maps, and generate
    % a timeseries lineplot.
    ax_roi = cell(size(MillimeterOfInterest, 1), 1);
    for j = 1:size(MillimeterOfInterest,1)
        % Convert mm to pixel coordinates
        [~, Xidx] = min(abs(X_img_mm-MillimeterOfInterest(j,1)));
        [~, Zidx] = min(abs(Z_img_mm-MillimeterOfInterest(j,2)));

        % Define center of patch in pixel coordinates
        pixelOfInterest = [Zidx Xidx];


        xrange = pixelOfInterest(2)-round(sizeOfVoxel(j, 1)/2):pixelOfInterest(2)+round(sizeOfVoxel(j, 1)/2);
        zrange = pixelOfInterest(1)-round(sizeOfVoxel(j, 2)/2):pixelOfInterest(1)+round(sizeOfVoxel(j, 2)/2);


        % If box boundaries extend beyond image size, restrict to image size
        xrange = xrange(xrange>0 & xrange < xPix);
        zrange = zrange(zrange>0 & zrange < yPix);

        % Convert from pixel back to mm coordinates
        leftside = min(xrange)*0.1;
        rightside = max(xrange) * 0.1;
        bottomside = (yPix - min(zrange))*0.1;
        topside = (yPix - max(zrange))*0.1;

        % Define color of box edges
        boxcolor = 'w';

        % Choose appropriate figure
        figure(fig_directionalresponse);

        % For each direction plot, put rectangle onto t-map, p-value, and/or
        % percent change plots.
        for i = 1:length(UniqueAngles)
            % Add rectangle to t-map
            if plotESmap
                % Specify appropriate axes in appropriate figure
                ax = nexttile(figureMapping(i));
                xlimits = xlim;
                ylimits = ylim;

                % Making sure the box only ends to right or top(bottom) edge,
                % but not beyond. Causes weird dynamic plotting issues
                % otherwise.
                xwidth = min([xlimits(2) - leftside, sizeOfVoxel(j, 1)*0.1, rightside - xlimits(1)]);
                ywidth = min([bottomside, sizeOfVoxel(j, 2)*0.1, ylimits(2) - ylimits(1) - topside]);


                %Add rectangle or update rectangle
                if exist('runAlready')
                    rect1{j,i}.Position = [leftside, ylimits(2) - bottomside, xwidth, ywidth];
                    rect1{j,i}.EdgeColor = boxcolor;
                else
                    rect1{j,i} = rectangle('Position',[leftside, ylimits(2) - bottomside,...
                        xwidth, ywidth],'LineWidth',2,...
                        'EdgeColor',boxcolor);
                end
            end

        end

        %Calculate STE across trials instead of across space
        PercentChange_ROI = squeeze(mean(mean(PercentChange_memPeriod(zrange,xrange,:),1, 'omitnan'),2));

        % For each direction, calculate the mean and ste time-series across
        % trials.

        % To remove extra points for dynamic positioning
        if useDynamicPositioning && exist('meanPercentChange_ROI_direction','var') && size(meanPercentChange_ROI_direction,2)>1
            meanPercentChange_ROI_direction = meanPercentChange_ROI_direction(:,1);
            stePercentChange_ROI_direction = meanPercentChange_ROI_direction(:,1);
        end


        for i = 1:length(UniqueAngles)
            meanPercentChange_ROI_direction(i,j) = mean(PercentChange_ROI(TargetPosInd==i));
            stePercentChange_ROI_direction(i,j) = std(PercentChange_ROI(TargetPosInd==i))/sqrt(length(PercentChange_ROI(TargetPosInd==i)));
        end

        % Average across the entire patch
        PercentChange_ROI_timeseries(j,:,:) = squeeze(mean(mean(iDopP_percentchange(zrange,xrange,:,:),1, 'omitnan'),2));

        % STE across all conditions
        if ~useDynamicPositioning
            stePercentChange_ROI(j) = std(PercentChange_ROI)/sqrt(length(PercentChange_ROI));
            fprintf('STE across all conditions for ROI %d is %0.3f \n', j, stePercentChange_ROI(j));
        end
    end


    runAlready = true;
    % Three color-blind safe colors
    if size(MillimeterOfInterest,1) <= 3
        c = [51 34 136;
            17 119 51;
            221 204 119;] ./ 255;
    elseif size(MillimeterOfInterest,1) <= 5
        % Five color-blind safe colors
        c = [100 143 255;
            120 94 240;
            220 38 127;
            254 97 0;
            255 176 0;] ./ 255;
    else
        % Need more colorblind safe colors, but for now will just use
        % distinguisable colors module

        c = distinguishable_colors(size(MillimeterOfInterest, 1));
    end


    % colorblind optimized colors
    color2 = [51 34 136;
        17 119 51;
        68 170 153;
        136 204 238;
        221 204 119;
        204 102 119;
        170 68 153;
        136 34 85;] ./ 255;


    % Equally spaced around RGB colorwheel
    color2 = [116 255 0;
        0 253 75;
        0 244 249;
        0 63 254;
        117 0 255;
        254 0 205;
        255 0 0;
        254 194 0;] ./ 255;
    color2 = color2 .* 0.85; %Adjust brightness/vibrance/etc

    % 8 direction timecourse plot for each patch of interest
    for j = 1:size(MillimeterOfInterest,1)
        if ~useDynamicPositioning
            figure(fig_direction_ERA);
            ax_roi{j} = nexttile(j);
        else
            ax_roi{4} = nexttile(4); cla reset;
        end


        % Calculate mean timeseries across trials for each direction
        for position = 1:length(UniqueAngles)
            PositionInd = TargetPosInd == position;
            TimeSeries_ROI_mean(j,:,position) = squeeze(mean(PercentChange_ROI_timeseries(j,:, PositionInd),3));
            TimeSeries_ROI_ste(j,:,position) = squeeze(std(PercentChange_ROI_timeseries(j,:, PositionInd),0,3)./ sqrt(nnz(PositionInd)));

        end


        % Get timestamps aligned to memory/cue start
        timepoints = (1:size(TimeSeries_ROI_mean,2))./coreParams.framerate + fix(1);

        % Plot the 8 lines with error bars onto the graph
        hold on;
        for position = 1:length(UniqueAngles)
            ts(10*(j-1)+position) = shadedErrorBar(timepoints, ...
                TimeSeries_ROI_mean(j,:,position),...
                abs(TimeSeries_ROI_ste(j,:,position)),...
                'lineprops',{'color',color2(position,:),'LineStyle','-','LineWidth',1});
        end
        hold off;

        % Label the plot
        xlabel('Trial Time (s)');
        ylabel('\DeltaPower Doppler (%)');

        % add scale bar
        hold on; plot([-2 -1], [-2.5 -2.5], 'k', 'LineWidth',1);
        text(-2, -2, '1 sec', 'HorizontalAlignment', 'left', 'FontSize', 15);

        % Save the ylimit before changing the figure
        ylimits = ylim;

        % Use actual time windows or indices
        type_of_boundary = 'window';
        switch type_of_boundary
            case 'window'
                baseline_time_period = baselinePeriod;
                analysis_time_period = memStat;

            case 'indices'
                baseline_time_period = [timepoints(baseline_indices(1)) timepoints(baseline_indices(end))];
                analysis_time_period = [timepoints(memory_indices(1)) timepoints(memory_indices(end))];
        end

        % If a single timepoint, use a line to ensure visiblity on
        % plot.
        if baseline_time_period(1) == baseline_time_period(2)
            hBase = line([baseline_time_period(1) baseline_time_period(1)], ...
                [-100 100]);
            hBase.Color = [0, 1, 0, 0.2];
        else
            hBase = patch([baseline_time_period(1) baseline_time_period(1) ...
                baseline_time_period(2) baseline_time_period(2)], ...
                [-100 100 100 -100], 'k', 'facecolor', 'g', 'edgecolor', 'none');
            hBase.FaceAlpha = 0.2;
        end

        if analysis_time_period(1) == analysis_time_period(2)
            hAnalysis = line([analysis_time_period(1) analysis_time_period(1)], ...
                [-100 100]);
            hAnalysis.Color = [0, 0, 1, 0.2];
        else
            hAnalysis = patch([analysis_time_period(1) analysis_time_period(1) ...
                analysis_time_period(2) analysis_time_period(2)], ...
                [-100 100 100 -100], 'k', 'facecolor', 'r', 'edgecolor', 'none');
            hAnalysis.FaceAlpha = 0.2;
        end


        % Set y-limit back
        ylim(ylimits)

        %Add legend
        legend_array = [];
        legend_name = {};
        for position = 1:length(UniqueAngles)
            legend_array = [legend_array, ts(10*(j-1)+position).mainLine];

            legend_name{position} = sprintf('%0.1f deg', UniqueAnglesInDegrees(position));
        end

        if j == 1
            legend(legend_array, legend_name)
        end

        % set the type
        xticks([fix(1), 0, mem(end)])
        xticklabels({'fixation','cue', 'move'})
        xlim([timepoints(1) min(timepoints(end),20)])
    end %End end over patch of interest

    % Fit curve to percent change data for the patch
    activation = meanPercentChange_ROI_direction;
    ste = stePercentChange_ROI_direction;

    if ~useDynamicPositioning
        figure(fig_tuningcurve);
    else
        ax_roi{8} = nexttile(8); cla reset;
    end
    for roi = 1:size(activation,2)
        % create sine fit
        sineEqn = 'a*sin(x*pi/4 + b) + c';
        f = fit((1:length(UniqueAngles))',activation(:,roi),'cubicspline');

        % plot error bars (actual data)
        plotErrBar(1:length(UniqueAngles), activation(:,roi), ste(:,roi), c(roi,:))
        % plot fitted sine
        pl(roi) = plot(1:.01:9,f(1:.01:9),'LineWidth',4,'Color',c(roi,:),'LineStyle',linestyles{roi});
    end


    % type setting now...
    line([0 8], [0 0],'LineStyle','--','Color','black','LineWidth',4)
    if ~useDynamicPositioning

        legendLabels = {'a','b','c','d','e', 'f', 'g', 'h', 'i'};
        legend(pl,legendLabels(1:size(MillimeterOfInterest,1)));
    end
    TargetDirectionsRadians = {'0','\pi/4','\pi/2','3\pi/4','\pi','5\pi/4','3\pi/2','7\pi/4'};
    xlim([1 8])
    xticks(1:length(TargetDirectionsRadians))
    xticklabels(TargetDirectionsRadians)
    xtickangle(-45)
    xlabel('Target Position (radians)');
    ylabel('\DeltaPower Doppler (%)')
    title('Voxel patch tuning curve');

    % update figure and process any pending callbacks
    if useDynamicPositioning
        drawnow; %Important for the dynamic updating to work
    end

    % If non-dynamic updating, exit loop on first run
    if ~useDynamicPositioning
        break;
    elseif (get(tb,'Value')==1)    % If stop button is pressed, stop updating the graphs
        delete(pointer)
        break;
    end
end %End while loop

% Resize all axes to same ylimit scale
if size(MillimeterOfInterest, 1) > 1
    linkaxes([ax_roi{:}]);
end


if savePlots && ~useDynamicPositioning
    filename_ERA = fullfile(savePath_base, ...
        sprintf('S%sR%d_EventRelatedAverage_MC_%s_thresh_%0.4f_dg%s.svg', ...
        session_names, ...
        Run, ...
        mc_method_string,...
        pvalue_threshold_to_save, ...
        datetime('today', format='yyyyMMdd')));
    saveas(fig_direction_ERA, filename_ERA, 'svg');
    % Also save as MATLAB figure
    saveas(fig_direction_ERA, strrep(filename_ERA, '.svg', '.fig'), 'fig');

    filename_tuningcurve = fullfile(savePath_base, ...
        sprintf('S%sR%d_TuningCurves_MC_%s_thresh_%0.4f_dg%s.svg', ...
        session_names, ...
        Run, ...
        mc_method_string, ...
        pvalue_threshold_to_save, ...
        datetime('today', format='yyyyMMdd')));
    saveas(fig_tuningcurve, filename_tuningcurve, 'svg');
    % Also save as MATLAB figure
    saveas(fig_tuningcurve, strrep(filename_tuningcurve, '.svg', '.fig'), 'fig');

    filename_ROI_overlay = fullfile(savePath_base, ...
        sprintf('S%sR%d_8dirSPM_withROIoverlay_MC_%s_%s_thresh_%0.4f_dg%s.svg', ...
        session_names, ...
        Run, ...
        mc_method_string, ...
        pvalue_threshold_to_save, ...
        datetime('today', format='yyyyMMdd')));
    saveas(fig_directionalresponse, filename_ROI_overlay, 'svg');
    % Also save as MATLAB figure
    saveas(fig_directionalresponse, strrep(filename_ROI_overlay, '.svg', '.fig'), 'fig');
end



%% How many voxels within each anatomically-defined region are directionally tuned?
% Quantification of directional activity within non-LIP regions.
show_subfigures = true;
switch Session
    case 102
region_names = {'Area 7op', 'Area 7a/b', 'LIP', 'VIP', 'MIP', 'Area 5', 'MP'}; %S102
    case 76
region_names = {'Area 7op', 'Area 7a/b', 'LIP', 'VIP', 'MIP', 'Area 5', 'Area 1/2', 'Area 23c'}; %76
    case 52
region_names = {'Area 7op', 'Area 7a/b', 'LIP', 'VIP', 'MIP', 'Area 5'}; %S52
    otherwise
        error('Region names not defined for %d', Session);
end
num_regions = length(region_names);
filename_suffix = '_manyRegions';


if size(sessionRunList, 1) > 1
    warning('Limited support for saving and loading of ROIs for aligned sessions yet.');
    if all(sessionRunList(:, 1) == sessionRunList(1, 1))
        standard_filename = strcat('Anatomical_ROI_S', num2str(Session), '_R', num2str(sessionRunList(1, 2)), '.mat');
        custom_filename = strrep(standard_filename, '.mat', strcat(filename_suffix, '.mat'));


        ROI_table = get_roi_info(Session, sessionRunList(1, 2), angiogram, UF, ...
            'region_names', region_names, ...
            'pixelsize', pixelsize,...
            'custom_filename', custom_filename);
    else
        ROI_table = get_roi_info([], [], angiogram, UF, ...
            'region_names', region_names, ...
            'pixelsize', pixelsize);
    end
else
    standard_filename = strcat('Anatomical_ROI_S', num2str(Session), '_R', num2str(Run), '.mat');
    custom_filename = strrep(standard_filename, '.mat', strcat(filename_suffix, '.mat'));

    ROI_table = get_roi_info(Session, Run, angiogram, UF, ...
        'region_names', region_names, ...
        'pixelsize', pixelsize, ...
        'custom_filename', custom_filename);
end

%master figure with all regions superimposed with activation map
significant_voxels_in_image = double(~pvalue_glm_mask);
significant_voxels_in_image(significant_voxels_in_image==0) = NaN;
significant_voxels_in_image(1,1) = 0;

figure;
title('Significant activation with anatomical region approximate boundaries');

cmap = plotDuplexImage(X_img_mm, Z_img_mm, significant_voxels_in_image, angiogram,...
    'colormap2use', winter, ...
    'nonlinear_bg',2);
[significant_voxel_count, num_region_voxels, percent_significant_voxels_region] = deal(NaN(num_regions, 1));
warning('off', 'MATLAB:polyshape:repairedBySimplify');
for region = 1:num_regions
    region_indx = strcmp(region_names{region}, ROI_table.Properties.RowNames);

    % Define mask
    region_boundary = cell2mat(ROI_table{region_indx, 'Boundary'});
    region_boundary = region_boundary(1:end-1,:);
    region_mask = poly2mask(region_boundary(:, 1), region_boundary(:, 2), yPix, xPix);

    % Create mask for activated voxels
    significant_voxels_in_region_mask = ~pvalue_glm_mask & region_mask;

    % How many voxels activated?
    significant_voxel_count(region) = sum(significant_voxels_in_region_mask, 'all');
    num_region_voxels(region) = sum(region_mask, 'all');

    percent_significant_voxels_region(region) = significant_voxel_count(region)/num_region_voxels(region)*100;
    fprintf('%0.2f %% of voxels within %s are activated\n', percent_significant_voxels_region(region), region_names{region}); 

    % Define region boundary in mm space
    region_boundary_X = region_boundary(:, 1)*pixelsize;
    region_boundary_Z = region_boundary(:, 2)*pixelsize+min(Z_img_mm);
    region_boundary_mm = polyshape(region_boundary_X, region_boundary_Z);

    % Plot region boundary
    hold on;
    plot(region_boundary_mm, ...
        'EdgeColor', 'w', ...
        'LineWidth', 1, ...
        'FaceAlpha', 0, ...
        'FaceColor', 'w');
    hold off;
end

% Reset to max image size
xlim([min(X_img_mm), max(X_img_mm)]);
ylim([min(Z_img_mm), max(Z_img_mm)]);


%% %%%%%%%%%%%%%%%%%%%%%%%%%% Tuning Curve Analysis %%%%%%%%%%%%%%%%%%%%%%%%
% Center of mass approach

% Make the theta matrix for each voxel in the iDop matrix
theta = repmat(UniqueAngles', yPix * xPix, 1);

% Convert the effect size into the same dimensions as the theta
rho = reshape(effectsize_dir,[],8); % entire matrix
rho_scaled = rho./repmat(max(abs(rho), [], 2), 1, 8);
%rho_scaled = rho;

[x,y] = pol2cart(theta,rho_scaled); %Converts each pixel into x,y representation of polar tuning curve

% Temporarily turning off warnings to reduce spurious output to the command
% window
id1 = 'MATLAB:polyshape:repairedBySimplify';
warning('off',id1);
id2 = 'MATLAB:polyshape:boundary3Points';
warning('off',id2);
h = waitbar(0,'Please wait ... computing center of mass of each pixel');
%Find centroid for the tuning curve at each pixel
[centroidTheta, centroidRho] = deal(NaN(size(x, 1), 1));
for pixel = 1:size(x,1)
    polyIn = polyshape([x(pixel,:)',y(pixel,:)']);
    [centroidX, centroidY] = centroid(polyIn);
    [centroidTheta(pixel), centroidRho(pixel)] = cart2pol(centroidX,centroidY);
    waitbar(pixel/(length(x)));
end
close(h);
% Turn back on the warnings
warning('on',id1);
warning('on',id2);

%Reshape the centroids back to the x,z dimensions of iDop
centroidTheta = reshape(centroidTheta, yPix, xPix);
centroidRho = reshape(centroidRho, yPix, xPix);


%% Plot tuning curves for entire image
pvalue_mask_to_use = pvalue_glm_mask;

% Convert all centroidTheta values to be between [- pi pi]
greaterThanPIind = centroidTheta > pi;
lessThanNegPIind = centroidTheta < - pi;

centroidTheta_plot = centroidTheta;
centroidTheta_plot(greaterThanPIind) = centroidTheta_plot(greaterThanPIind) - 2*pi;
centroidTheta_plot(lessThanNegPIind) = centroidTheta_plot(lessThanNegPIind) + 2*pi;

% Define transparency map for the center of mass approach
centroidRho_plot = centroidRho;
transparency_map = 2 * centroidRho_plot;
transparency_map(transparency_map>1) = 1;

transparency_map(~pvalue_mask_to_use) = 1;

transparency_map(~pvalue_mask_to_use) = 1;


%Threshold centroid plots using p-value threshold ind from above
centroidTheta_plot(pvalue_mask_to_use) = NaN;
centroidRho_plot(pvalue_mask_to_use) = NaN;


% Load colormap to be used
ColorMapToUse = load('circular_colormap.mat');
ColorMapToUse = ColorMapToUse.colorwheel;

% Define labels depending on number of directions
TargetDirections = {'0 deg','45 deg','90 deg','135 deg','180 deg','- 135 deg','- 90 deg','- 45 deg'};
TargetDirectionsRadiansNegPiToPosPi = {'-\pi','-3\pi/4','-\pi/2','-\pi/4','0','\pi/4','\pi/2','3\pi/4','\pi'};

fig_preferred_direction = figure;

tuning_cmap{1} = plotDuplexImage(X_img_mm,Z_img_mm,centroidTheta_plot,angiogram,...
    'colormap2use',ColorMapToUse, ...
    'nonlinear_bg',2,...
    'AutoColorBarLimits',[ - pi pi]);
title('Center of Mass');

cb = colorbar('Ticks',0:1/8:1,'TickLabels',TargetDirectionsRadiansNegPiToPosPi);

set(cb,'YLim',[0 1])
ylabel(cb, 'Target Position');

currentRow = find(ProjectRecord.Session==sessionRunList(1,1) & ProjectRecord.Run==sessionRunList(1,2));


if savePlots
    filename = fullfile(savePath_base, ...
        sprintf('S%sR%d_PreferredDirectionTuning_MC_%s_thresh_%0.4f_dg%s.svg', ...
        session_names, ...
        Run, ...
        mc_method_string, ...
        pvalue_threshold_to_save, ...
        datetime('today', format='yyyyMMdd')));
    saveas(fig_preferred_direction, filename, 'svg');

    % Also save as MATLAB figure
    saveas(fig_preferred_direction, strrep(filename, '.svg', '.fig'), 'fig');
end

%% Visualize tuning strength (Cohen d) from center of mass

% Put centroidTheta in [-pi pi]
greaterThanPIind = centroidTheta > pi;
lessThanNegPIind = centroidTheta < - pi;

tuned_direction = centroidTheta;
tuned_direction(greaterThanPIind) = tuned_direction(greaterThanPIind) - 2*pi;
tuned_direction(lessThanNegPIind) = tuned_direction(lessThanNegPIind) + 2*pi;

% Add circular statistics
theta = repmat(UniqueAngles', yPixels * xPixels, 1);
rho = reshape(effectsize_dir,[],8); % entire matrix
rho_scaled = rho./repmat(max(abs(rho), [], 2), 1, 8);

% Calculate the new circular statistics
circ_stats = get_circular_statistics(theta, rho);
circ_std = circ_stats.circ_std_alt;
skewness = circ_stats.skewness_alt;
kurtosis = circ_stats.kurtosis_alt;

% Laterality index
% Formula
% From "Sustained Activity in Topographic Areas of Human Posterior
% Parietal Cortex during Memory-Guided Saccades
% Denis Schluppeck, Clayton E. Curtis, Paul W. Glimcher, and David J. Heeger
% Same formula is also used by Kagan et al. 2010
contra_ind = [1 2 8];
ipsi_ind = [4 5 6];
response = rho;
contra_response = mean(response(:, contra_ind), 2);
ipsi_response = mean(response(:, ipsi_ind), 2);
laterality_index = (contra_response - ipsi_response) ./ (abs(contra_response) + abs(ipsi_response));
laterality_index = reshape(laterality_index, yPix, xPix);

% Smooth
disk_filter_size = 1;
% Define smoothing filter
if disk_filter_size > 0
    disk_filter = fspecial('disk', disk_filter_size);
end
tuned_direction_smooth = imfilter(tuned_direction, disk_filter);
strength_smooth = imfilter(centroidRho, disk_filter);
fscore_smooth = imfilter(F_registered, disk_filter);
laterality_smooth = imfilter(laterality_index, disk_filter);
circ_std_smooth = imfilter(circ_std, disk_filter);
skewness_smooth = imfilter(skewness, disk_filter);
kurtosis_smooth = imfilter(kurtosis, disk_filter);

tuned_direction_plot = threshold_matrix(tuned_direction_smooth, pvalue_mask_to_use);
strength_plot = threshold_matrix(strength_smooth, pvalue_mask_to_use);
fscore_plot = threshold_matrix(fscore_smooth, pvalue_mask_to_use);
laterality_plot = threshold_matrix(laterality_smooth, pvalue_mask_to_use);
circ_std_plot = threshold_matrix(circ_std_smooth, pvalue_mask_to_use);
skewness_plot = threshold_matrix(skewness_smooth, pvalue_mask_to_use);
kurtosis_plot = threshold_matrix(kurtosis_smooth, pvalue_mask_to_use);

fig_statsOverlay = figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);
tld = tiledlayout('flow');
% Tuned direction
nexttile(1);
plotDuplexImage(X_img_mm,Z_img_mm,tuned_direction_plot,angiogram,...
    'colormap2use',ColorMapToUse, ...
    'nonlinear_bg',2,...
    'AutoColorBarLimits',[ - pi pi],...
    'showColorbar', true);
axis off;
title('Center of Mass');

% Tuned strength (cohen d)
max_value = max(centroidRho_plot, [], 'all');
nexttile(2);
plotDuplexImage(X_img_mm,Z_img_mm,strength_plot,angiogram,...
    'colormap2use',magma, ...
    'nonlinear_bg',2,...
    'showColorbar', true, ...
    'AutoColorBarLimits',[0 max_value]);
axis off;
title('Tuning strength');

% F-score
nexttile(3);
plotDuplexImage(X_img_mm,Z_img_mm,fscore_plot, angiogram,...
    'colormap2use',magma, ...
    'nonlinear_bg',2, ...
    'showColorbar', true, ...
    'ColorBarTitle', 'F-score');
axis off;
title('F-score');

% circular std
nexttile(4);
plotDuplexImage(X_img_mm,Z_img_mm,circ_std_plot,angiogram,...
    'colormap2use',magma, ...
    'nonlinear_bg',2,...
    'showColorbar', true);
axis off;
title('Circular STD');

% skewness
abs_max = max(abs(skewness_plot), [], 'all');
nexttile(5);
plotDuplexImage(X_img_mm,Z_img_mm,skewness_plot,angiogram,...
    'colormap2use',PositiveNegativeColormap, ...
    'nonlinear_bg',2,...
    'showColorbar', true, ...
    'AutoColorBarLimits', [-abs_max abs_max]);
axis off;
title('Skewness');

% kurtosis
abs_max = max(abs(kurtosis_plot), [], 'all');
nexttile(6);
plotDuplexImage(X_img_mm,Z_img_mm,kurtosis_plot,angiogram,...
    'colormap2use',PositiveNegativeColormap, ...
    'nonlinear_bg',2,...
    'showColorbar', true, ...
    'AutoColorBarLimits', [-abs_max abs_max]);
axis off;
title('Kurtosis');

% laterality
abs_max = max(abs(laterality_plot), [], 'all');
nexttile(7);
plotDuplexImage(X_img_mm, Z_img_mm, laterality_plot, angiogram,...
    'colormap2use',PositiveNegativeColormap, ...
    'nonlinear_bg',2,...
    'showColorbar', true, ...
    'AutoColorBarLimits', [-abs_max abs_max]);
axis off;
title('Laterality');

title(tld, sprintf('S%dR%d summary statistics overlaid\n q < %.0e', Session, Run, pvalue_threshold));

if savePlots
    filename = fullfile(savePath_base, ...
        sprintf('S%sR%d_SummaryStatsOverlay_MC_%s_thresh_%0.4f_dg%s.svg', ...
        session_names, ...
        Run, ...
        mc_method_string, ...
        pvalue_threshold_to_save, ...
        datetime('today', format='yyyyMMdd')));
    saveas(fig_statsOverlay, filename, 'svg');

    % Also save as MATLAB figure
    saveas(fig_statsOverlay, strrep(filename, '.svg', '.fig'), 'fig');
end

%% Report how many LIP voxels and non-LIP voxels contralateral preferring
num_LIP_voxels = sum(~isnan(laterality_plot) & LIP_mask, 'all');
num_nonLIP_voxels = sum(~isnan(laterality_plot) & ~LIP_mask, 'all');
percent_contralateral_preferring_LIP = sum(laterality_plot(LIP_mask)>0, 'all')/num_LIP_voxels;
percent_contralateral_preferring_nonLIP = sum(laterality_plot(~LIP_mask)>0, 'all')/num_nonLIP_voxels;
percent_contralateral_preferring = sum(laterality_plot(:)>0)/sum(~isnan(laterality_plot), 'all');

fprintf('%0.2f %% of voxels are contralateral preferring\n', percent_contralateral_preferring*100);
fprintf('%0.2f %% of voxels OUTside of LIP are contralateral preferring\n', percent_contralateral_preferring_nonLIP*100);
fprintf('%0.2f %% of voxels INside LIP are contralaterally preferring\n', percent_contralateral_preferring_LIP*100);


% Show how many voxels within each region
show_figure = true;
if show_figure
    figure;
    tld = tiledlayout('flow');
    nexttile()
    histogram(laterality_plot(LIP_mask))
    title('LIP - laterality');

    nexttile()
    histogram(laterality_plot(~LIP_mask))
    title('nonLIP - laterality');
end

%% Save results to a file
filename = fullfile(savePath_base, ...
    sprintf('S%sR%d_8tgt_analysis_mc_%s_dg%s.mat', ...
    session_names, ...
    Run, ...
    mc_method_string, ...
    datetime('today', format='yyyyMMdd')));

save(filename, 'pvalue_threshold', 'pvalue_glm_threshold', ...
    'centroidTheta', 'centroidRho',...
    'p_cor_registered', 'pvalue_dir', ...
    'F_registered',...
    'effectsize_dir', ...
    'angiogram', 'Session', 'Run', ...
    'mc_method', 'memStat', 'baselinePeriod');

%% extra functions
function plotErrBar(x,m,s,c)
% uses the matlab built in bar plot & adds an error bar
% x is the x position to plot
% m is the mean (bar y pos), s is the error (bar length)
% c is the color
for i = 1:length(x)
    hold on
    plot([x(i) x(i)],[m(i)-s(i) m(i)+s(i)],'LineWidth',1,'Color',c)
end
end

% Threshold function
function thresholded_mat = threshold_matrix(mat, threshold_mask)
% Helper function to quickly threshold and mask matrices, such as for
% thresholding of pvalue mask
thresholded_mat = mat;
thresholded_mat(threshold_mask) = NaN;
end

