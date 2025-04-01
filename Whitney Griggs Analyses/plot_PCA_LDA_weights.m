%% Visualize PCA+LDA weights
% Feburary 2, 2025 WG
% Calculate and plot PCA+LDA weights
close all

%% hard coded settings

% preprocessing
timeGain = false;
diskFilter = 2;
zScore = false;
detrend = true;

saveMe = true;

% Motion correction method (if any)
motion_correction_method = 'normcorre';
% Defining this string for saving of figures/data later
if isempty(motion_correction_method)
    motion_correction_string = 'None';
else
    motion_correction_string = motion_correction_method;
end


%% Select which sessions we want to iterate over
project_record_filename = 'ProjectRecord_paper.json';
ProjectRecord = load_json_as_table(project_record_filename);
session_run_list = specify_sessions_of_interest(...
    'project_record_filename', project_record_filename);


%% get data, preprocess, get epochs & dimension sizes
loadDopplerData(whos, 'multiple', false, ...
    'mc_method', motion_correction_method, ...
    'project_record_filename', project_record_filename, ...
    'SessionRunList', session_run_list);
if detrend
    iDopP = detrend_sliding_window(iDop, 50);
else
    iDopP = iDop;
end
iDopP = preProcess(iDopP, 'timeGain', timeGain, 'spatialFilter', {'disk', diskFilter}, 'zScore', zScore);


[fix, memory_end, mov_end] = getEpochs(behavior);
fixation_timepoints = ceil(abs(fix(1))*coreParams.framerate);
memory_timepoints = ceil(abs(memory_end(2))*coreParams.framerate);

[yPix, xPix, nWindows, nTrials] = size(iDopP);

pixelsize = 0.1;
X_img_mm = pixelsize/2 + (0:size(angiogram,2)-1)*pixelsize; % So that the center of the image is at true halfway point, not shifted.
Z_img_mm = pixelsize/2 + (0:size(angiogram,1)-1)*pixelsize + UF.Depth(1);
trialTime = fix(1):1/coreParams.framerate:fix(1)+nWindows/coreParams.framerate-1/coreParams.framerate;

% Clear this variable to reduce memory usage
clear iDop

%% Assigning labels to multi-targets
target = cell2mat({behavior.target});
targetPos = vertcat(behavior.targetPos);
if iscell(targetPos)
    targetPos = cell2mat(targetPos);
end
[angle,distance] = cart2pol(targetPos(:,1),targetPos(:,2)); % convert to polar coordinates
angle(angle<0) = angle(angle<0)+2*pi; %Convert to have only positive angles.

UniqueAngles = unique(angle);
TargetPosInd = zeros(length(target),1);
nClasses = length(UniqueAngles);
for position = 1:nClasses
    TargetPosInd(angle==UniqueAngles(position)) = position;
end

% create the training data labels (2 x n_trials) - x, y columns
% For horizontal and vertical axes, use the following labels
% 1 = Negative (down or left)
% 2 = At center (along vertical or horizontal axis)
% 3 (or 2 if not using middle points) = Positive (up or right)

label_strings = {'Negative','Center','Positive'};
train_labels = NaN(size(targetPos,1),2);

for dimension = 1:2
    train_labels(targetPos(:, dimension) < 0, dimension) = 1;

    train_labels(targetPos(:, dimension) == 0, dimension) = 2;
    train_labels(targetPos(:, dimension) > 0, dimension) = 3;
end

%% Compute PCA+LDA results for 8 classes using only single timepoints
% Single timepoint only

lda_weighting_reshaped_singleTimepoints = zeros(yPix, xPix, 7*7, size(iDopP, 3));
for timepoint = 1:size(iDopP, 3)
    % Will create an upper diagonal matrix showing the boundaries for the
    % between classes.
    trainData = flattenDoppler2D(iDopP, timepoint);  % nImages x (nPixels*yPixels)
    train_labels_8class = TargetPosInd;

    % Z-score data
    trainData = zscore(trainData);

    % PCA
    [pcaCoefficients, pcaScores, ~, ~, explained, pcaCenters] = pca(trainData);
    explainedVarianceToKeepAsFraction = 95/100;
    numComponentsToKeep = find(cumsum(explained)/sum(explained) >= explainedVarianceToKeepAsFraction, 1);
    pcaCoefficients = pcaCoefficients(:,1:numComponentsToKeep);

    trainPredictors = pcaScores(:,1:numComponentsToKeep);   % (PCA transformed trainData)

    % LDA
    MdlLinear = fitcdiscr(trainPredictors,train_labels_8class); % No cross-validation here

    % Visualize LDA in PCA space
    counter = 1;
    for class1 = 1:7
        for class2 = 2:8
            if class2 > class1
                coefficients_constant = MdlLinear.Coeffs(class1, class2).Const; %Don't know how to handle this coefficient currently
                coefficients_linear = MdlLinear.Coeffs(class1, class2).Linear;

                pca_lda_weightings = pcaCoefficients * coefficients_linear;

                lda_weighting_reshaped_singleTimepoints(:, :, counter, timepoint) = reshape(pca_lda_weightings, yPix, xPix);
            end
            counter = counter+1;
        end
    end

end

%%  Compute PCA+LDA weights for 8 classes using all timepoints
% Do this for all timepoints
PCA_timeOfInterest = fixation_timepoints:size(iDopP, 3);

% Will create an upper diagonal matrix showing the boundaries for the
% between classes.
trainData = flattenDoppler2D(iDopP, PCA_timeOfInterest);  % nImages x (nPixels*yPixels)
train_labels_8class = TargetPosInd;

% Remap so that it is consistent with the order of the confusion matrix.
class_remapping = [4 3 2 1 8 7 6 5];
remapped_train_labels = class_remapping(train_labels_8class);

% Z-score data
trainData = zscore(trainData);

% PCA
[pcaCoefficients, pcaScores, ~, ~, explained, pcaCenters] = pca(trainData);
explainedVarianceToKeepAsFraction = 95/100;
numComponentsToKeep = find(cumsum(explained)/sum(explained) >= explainedVarianceToKeepAsFraction, 1);
pcaCoefficients = pcaCoefficients(:,1:numComponentsToKeep);

trainPredictors = pcaScores(:,1:numComponentsToKeep);   % (PCA transformed trainData)

% 8-class LDA
MdlLinear = fitcdiscr(trainPredictors,remapped_train_labels); % No cross-validation here

% Visualize 8-class LDA in PCA space
lda_weighting_reshaped = zeros(yPix, xPix*length(PCA_timeOfInterest), 7*7);
counter = 1;
for class1 = 1:7
    for class2 = 2:8
        if class2 > class1
            coefficients_constant = MdlLinear.Coeffs(class1, class2).Const; %Don't know how to handle this coefficient currently
            coefficients_linear = MdlLinear.Coeffs(class1, class2).Linear;

            pca_lda_weightings = pcaCoefficients * coefficients_linear ;%+ coefficients_constant;

            lda_weighting_reshaped(:, :, counter) = reshape(pca_lda_weightings, yPix, xPix*length(PCA_timeOfInterest));
        end
        counter = counter+1;
    end
end

%% Visualize 2 x 3-class LDA in PCA space
% 2 x 3-way LDA
MdlLinear_horz = fitcdiscr(trainPredictors,train_labels(:, 1)); % No cross-validation here
MdlLinear_vert = fitcdiscr(trainPredictors,train_labels(:, 2)); % No cross-validation here

% Order for consistency is UpperLeft (UL), UMiddle, URight, Right, BR, BM, BL, L
horz_class_lookup_table = [1 2 3 3 3 2 1 1];
vert_class_lookup_table = [3 3 3 2 1 1 1 2];
lda_weighting_multicoder_reshaped = zeros(yPix, xPix*length(PCA_timeOfInterest), 7*7);
counter = 1;
for class1 = 1:7
    for class2 = 2:8
        if class2 > class1

            % Horz
            horz_class1 = horz_class_lookup_table(class1);
            horz_class2 = horz_class_lookup_table(class2);
            if horz_class1 ~= horz_class2
                coefficients_linear = MdlLinear_horz.Coeffs(horz_class1, horz_class2).Linear;
                pca_lda_weightings_horz = pcaCoefficients * coefficients_linear ;%+ coefficients_constant;
                use_horz = true;
            else
                pca_lda_weightings_horz = ones(size(pcaCoefficients, 1), 1);
                use_horz = false;
            end

            % Vert
            vert_class1 = vert_class_lookup_table(class1);
            vert_class2 = vert_class_lookup_table(class2);
            if vert_class1 ~= vert_class2
                coefficients_linear = MdlLinear_vert.Coeffs(vert_class1, vert_class2).Linear;
                pca_lda_weightings_vert = pcaCoefficients * coefficients_linear ;%+ coefficients_constant;
                use_vert = true;
            else
                pca_lda_weightings_vert = ones(size(pcaCoefficients, 1), 1);
                use_vert = false;
            end

            % Combine horz and vert
            pca_lda_weightings = pca_lda_weightings_horz .* pca_lda_weightings_vert;

            % To get similar scale, need to square root if both use_vert
            % and use_horz are true. Otherwise, those combinations will be
            % approximately a square smaller than combinations with only
            % one horz or vert class.
            if use_vert && use_horz
                pca_lda_weightings = sqrt(abs(pca_lda_weightings)) .* sign(pca_lda_weightings);
            end

            lda_weighting_multicoder_reshaped(:, :, counter) = reshape(pca_lda_weightings, yPix, xPix*length(PCA_timeOfInterest));
        end
        counter = counter+1;
    end
end


%% Visualize PCA+LDA results for 8 classes
% Single timepoint only


absolute_max = max(abs(lda_weighting_reshaped_singleTimepoints), [], 'all');

for timepoint = 1:size(lda_weighting_reshaped_singleTimepoints, 4)
    % Plotting the visualization
    fig_pcaProjection = figure('units','normalized','outerposition',[0 0 1 1]);
    tld = tiledlayout(7, 7);

    counter = 1;
    custom_colormap = [[ones(100, 1) linspace(0, 1, 100)' linspace(0, 1, 100)']; [linspace(1, 0, 100)' linspace(1, 0, 100)' ones(100, 1)]];

    display_percentage = 0.05;

    fontsize = 12;

    for class1 = 1:7
        for class2 = 2:8
            if class2 > class1
                nexttile(counter);
                plotDuplexImage(X_img_mm, Z_img_mm, lda_weighting_reshaped_singleTimepoints(:, :, counter, timepoint), angiogram, ...
                    'displayStrongest', display_percentage,...
                    'colormap2use',custom_colormap, ...
                    'nonlinear_bg',2, ...
                    'AutoColorBarLimits',[-absolute_max absolute_max], ...
                    'FontSize', fontsize, ...
                    'showColorbar', false);
                %axis off;
            end

            if class1 == 3 && class2 == 2
                nexttile(counter);
                % Turn off the axes
                axis off;
                % Draw the colorbar
                c = colorbar('FontSize',fontsize,'Location','east');
                % Specify the colormap
                colormap(gca,custom_colormap)
                % Define the colormap limits.
                caxis([-absolute_max absolute_max]);
                % Label the colorbar
                c.Label.String = 'PCA+LDA weightings';
                c.Ticks = [-absolute_max 0 absolute_max];
            end
            counter = counter+1;
        end
    end

    tld.TileSpacing = 'none';
    tld.Padding = 'tight';
    title(tld, sprintf('%0.2f%% strongest (pos and neg) PCA+LDA weightings - \nS%dR%d - Index %d', display_percentage*100, Session, Run, timepoint));

    save_to_file = false;
    if save_to_file && saveMe

        suggested_fname = sprintf('pca_projection_maps_S%dR%d_gen%s.svg', Session, Run, datetime('today', format='yyyyMMdd'));

        [filename, filepath] = uiputfile(fullfile(get_user_data_path('pathType', 'decoding'), suggested_fname));
        if any([filename filepath])
            saveas(fig_pcaProjection, fullfile(filepath, filename));
        end
    end
end

%%  Visualize PCA+LDA results for 8 classes
% Since we do not analyze timepoints before fixation, cut the trialTime to
% reflect this
trialTime_cut = trialTime(fixation_timepoints:end);

absolute_max = max(abs(lda_weighting_reshaped), [], 'all');

% Plotting the visualization
save_to_file = true;
if save_to_file && saveMe

    data_save_root = fullfile(get_user_data_path('pathType', 'decoding'), 'movies');
    suggested_fname = sprintf('pca_projection_video_S%dR%d_gen%s.mp4', Session, Run, datetime('today', format='yyyyMMdd'));

    [filename, filepath] = uiputfile(fullfile(data_save_root, suggested_fname));
    if any([filename filepath])
        video_fullfile = fullfile(filepath, filename);

        vidfile = VideoWriter(video_fullfile,'MPEG-4');
        vidfile.FrameRate = 1;
        open(vidfile);
    end
end


fig_pcaProjection = figure('units','normalized','outerposition',[0 0 1 1]);
% Assume that we have a PCA frame for every timepoint in the trial. Used to
% only use post fixation time, but this is cleaner. Not making backwards
% compatible since we did not document in those data files which timepoints
% were used.
for timepoint = 1:length(PCA_timeOfInterest)
    clf;
    tld = tiledlayout(7, 7);

    counter = 1;
    custom_colormap = [[ones(100, 1) linspace(0, 1, 100)' linspace(0, 1, 100)']; [linspace(1, 0, 100)' linspace(1, 0, 100)' ones(100, 1)]];

    display_percentage = 0.05;

    fontsize = 12;

    for class1 = 1:7
        for class2 = 2:8
            if class2 > class1
                nexttile(counter);
                plotDuplexImage(X_img_mm, Z_img_mm, lda_weighting_reshaped(:, (1+((timepoint-1)*xPix)):xPix*timepoint, counter), angiogram, ...
                    'displayStrongest', display_percentage, ...
                    'colormap2use',custom_colormap, ...
                    'nonlinear_bg',2, ...
                    'AutoColorBarLimits',[-absolute_max absolute_max], ...
                    'FontSize', fontsize, ...
                    'showColorbar', false);
                axis off;
            end

            if class1 == 1 && class2 == class1 + 1
                axis on;
            end

            if class1 == 3 && class2 == 2
                nexttile(counter);
                % Turn off the axes
                axis off;
                % Draw the colorbar
                c = colorbar('FontSize',fontsize,'Location','east');
                % Specify the colormap
                colormap(gca,custom_colormap)
                % Define the colormap limits.
                caxis([-absolute_max absolute_max]);
                % Label the colorbar
                c.Label.String = 'PCA+LDA weightings';
                c.Ticks = [-absolute_max 0 absolute_max];
            end
            counter = counter+1;
        end
    end

    tld.TileSpacing = 'none';
    tld.Padding = 'tight';

    if trialTime_cut(timepoint) < 0
        title(tld, sprintf('%0.2f%% strongest (pos and neg) PCA+LDA weightings - \nS%dR%d - \nTime since cue %0.1f sec\n Fixation period', display_percentage*100, Session, Run, trialTime_cut(timepoint)));
    elseif trialTime_cut(timepoint) >= 0 && trialTime_cut(timepoint) <= memory_end(2)
        title(tld, sprintf('%0.2f%% strongest (pos and neg) PCA+LDA weightings - \nS%dR%d - \nTime since cue %0.1f sec\n Memory period', display_percentage*100, Session, Run, trialTime_cut(timepoint)));
    elseif trialTime_cut(timepoint) > memory_end(2) && trialTime_cut(timepoint) <= mov_end(2)
        title(tld, sprintf('%0.2f%% strongest (pos and neg) PCA+LDA weightings - \nS%dR%d - \nTime since cue %0.1f sec\n Movement period', display_percentage*100, Session, Run, trialTime_cut(timepoint)));
    else
        title(tld, sprintf('%0.2f%% strongest (pos and neg) PCA+LDA weightings - \nS%dR%d - \nTime since cue %0.1f sec\n ITI', display_percentage*100, Session, Run, trialTime_cut(timepoint)));
    end
    if exist('vidfile', 'var')
        F = getframe(fig_pcaProjection);
        writeVideo(vidfile, F);
    end
    drawnow;
end

if exist('vidfile', 'var')
    close(vidfile);
    clear vidfile
end


%%  Visualize PCA+LDA results for 2 x 3-class LDA multicoder
absolute_max = max(abs(lda_weighting_multicoder_reshaped), [], 'all');

% Plotting the visualization
save_to_file = true;
if save_to_file && saveMe

    data_save_root = fullfile(get_user_data_path('pathType', 'decoding'), 'movies');
    suggested_fname = sprintf('pca_projection_multicoder_video_S%dR%d_gen%s.mp4', Session, Run, datetime('today', format='yyyyMMdd'));

    [filename, filepath] = uiputfile(fullfile(data_save_root, suggested_fname));
    if any([filename filepath])
        video_fullfile = fullfile(filepath, filename);

        vidfile = VideoWriter(video_fullfile,'MPEG-4');
        vidfile.FrameRate = 1;
        open(vidfile);
    end
end


fig_pcaProjection = figure('units','normalized','outerposition',[-0.4188 0.3493 0.4133 0.902]);
% Assume that we have a PCA frame for every timepoint in the trial. Used to
% only use post fixation time, but this is cleaner. Not making backwards
% compatible since we did not document in those data files which timepoints
% were used.
for timepoint = 1:length(PCA_timeOfInterest)
    clf;
    tld = tiledlayout(7, 7);

    counter = 1;
    custom_colormap = [[ones(100, 1) linspace(0, 1, 100)' linspace(0, 1, 100)']; [linspace(1, 0, 100)' linspace(1, 0, 100)' ones(100, 1)]];

    display_percentage = 0.05;

    fontsize = 12;


    for class1 = 1:7
        for class2 = 2:8
            if class2 > class1
                nexttile(counter);

                normalized_values = lda_weighting_multicoder_reshaped(:, (1+((timepoint-1)*xPix)):xPix*timepoint, counter);
                normalized_values = normalized_values./max(abs(normalized_values), [], 'all');
                plotDuplexImage(X_img_mm, Z_img_mm, lda_weighting_multicoder_reshaped(:, (1+((timepoint-1)*xPix)):xPix*timepoint, counter), angiogram, ...
                    'displayStrongest', display_percentage, ...
                    'colormap2use',custom_colormap, ...
                    'nonlinear_bg',2, ...
                    'AutoColorBarLimits',[-absolute_max absolute_max], ...
                    'FontSize', fontsize, ...
                    'showColorbar', false);
                axis off;
            end

            if class1 == 1 && class2 == class1 + 1
                axis on;
            end

            if class1 == 3 && class2 == 2
                nexttile(counter);
                % Turn off the axes
                axis off;
                % Draw the colorbar
                c = colorbar('FontSize',fontsize,'Location','east');
                % Specify the colormap
                colormap(gca,custom_colormap)
                % Define the colormap limits.
                caxis([-absolute_max absolute_max]);
                % Label the colorbar
                c.Label.String = 'PCA+LDA weightings';
                c.Ticks = [-absolute_max 0 absolute_max];
            end
            counter = counter+1;
        end
    end

    tld.TileSpacing = 'none';
    tld.Padding = 'tight';

    if trialTime_cut(timepoint) < 0
        title(tld, sprintf('%0.2f%% strongest (pos and neg) PCA+LDA weightings - \nS%dR%d - \nTime since cue %0.1f sec\n Fixation period', display_percentage*100, Session, Run, trialTime_cut(timepoint)));
    elseif trialTime_cut(timepoint) >= 0 && trialTime_cut(timepoint) <= memory_end(2)
        title(tld, sprintf('%0.2f%% strongest (pos and neg) PCA+LDA weightings - \nS%dR%d - \nTime since cue %0.1f sec\n Memory period', display_percentage*100, Session, Run, trialTime_cut(timepoint)));
    elseif trialTime_cut(timepoint) > memory_end(2) && trialTime_cut(timepoint) <= mov_end(2)
        title(tld, sprintf('%0.2f%% strongest (pos and neg) PCA+LDA weightings - \nS%dR%d - \nTime since cue %0.1f sec\n Movement period', display_percentage*100, Session, Run, trialTime_cut(timepoint)));
    else
        title(tld, sprintf('%0.2f%% strongest (pos and neg) PCA+LDA weightings - \nS%dR%d - \nTime since cue %0.1f sec\n ITI', display_percentage*100, Session, Run, trialTime_cut(timepoint)));
    end
    if exist('vidfile', 'var')
        F = getframe(fig_pcaProjection);
        writeVideo(vidfile, F);
    end
    drawnow;
end

if exist('vidfile', 'var')
    close(vidfile);
    clear vidfile
end


%% Average each image into single number for each timepoint
n_timepoints = length(PCA_timeOfInterest);
activity_1D_time = NaN(7, 7, n_timepoints);

for timepoint = 1:n_timepoints
    counter = 1;
    for class1 = 1:7
        for class2 = 2:8
            if class2 > class1
                activity_1D_time(class1, class2, timepoint) = mean(abs(lda_weighting_reshaped(:, (1+((timepoint-1)*xPix)):xPix*timepoint, counter)), 'all');
            end
            counter = counter+1;
        end
    end
end

absolute_max = max(activity_1D_time, [], 'all');
absolute_min = min(activity_1D_time, [], 'all');

fig_pcaProjection = figure('units','normalized','outerposition',[ 0.1 0.1 0.8 0.8]);
tld = tiledlayout(7, 7);

counter = 1;

mov_start_time = mov_end(1)*coreParams.framerate;

for class1 = 1:7
    for class2 = 2:8
        if class2 > class1
            nexttile(counter);
            plot(0:n_timepoints-1, squeeze(activity_1D_time(class1, class2, :)));

            ylim([absolute_min absolute_max]);
            xlim([0 n_timepoints]);

            % shade fixation & memory period
            hMem = patch([0 0 mov_start_time mov_start_time], [absolute_min absolute_max absolute_max absolute_min], ...
                'k', 'facecolor','g','edgecolor','none');
            hMem.FaceAlpha = 0.2;
            xticks([0 mov_start_time])
            xticklabels({'Cue','Go'})

        end

        if class1 == 1 && class2 == class1 + 1

        else
            set(gca, 'XTickLabel', []);
            set(gca, 'XTick', []);
            set(gca, 'YTickLabel', []);
            set(gca, 'YTick', []);
        end

        counter = counter+1;
    end
end
%
tld.TileSpacing = 'none';
tld.Padding = 'tight';


%% Average each class-class timeseries into single image
weighting_to_use = '2x3class'; % options are '2x3class' or '8class'

switch weighting_to_use
    case '2x3class'
        lda_weighting_reshaped_to_use = lda_weighting_multicoder_reshaped;
    case '8class'
        lda_weighting_reshaped_to_use = lda_weighting_reshaped;
end


n_timepoints = size(lda_weighting_reshaped_to_use, 2)/xPix;
extracted_images_time = NaN(yPix, xPix, n_timepoints, 7, 8);

for timepoint = 1:n_timepoints

    counter = 1;

    for class1 = 1:7
        for class2 = 2:8
            if class2 > class1
                extracted_images_time(:, :, timepoint, class1, class2) = lda_weighting_reshaped_to_use(:, (1+((timepoint-1)*xPix)):xPix*timepoint, counter);
            end

            counter = counter+1;
        end
    end
end

average_across_time = squeeze(mean(extracted_images_time, 3));

fig_pcaProjection = figure('units','normalized','outerposition',[  0 0 1 1]);

tld = tiledlayout(8, 8);
counter = 1;
display_percentage = 0.05;
fontsize = 12;
absolute_max = max(abs(average_across_time), [], 'all');

% Hard-coding for publication purposes
for class1 = 1:8
    for class2 = 1:8
        if class2 > class1
            nexttile(counter);
            plotDuplexImage(X_img_mm, Z_img_mm, average_across_time(:, :, class1, class2), angiogram, ...
                'displayStrongest', display_percentage, ...
                'colormap2use',custom_colormap, ...
                'nonlinear_bg',2, ...
                'AutoColorBarLimits',[-absolute_max absolute_max], ...
                'FontSize', fontsize, ...
                'showColorbar', false);
            axis off;
        elseif class2 < class1
            nexttile(counter);
            plotDuplexImage(X_img_mm, Z_img_mm, average_across_time(:, :, class2, class1), angiogram, ...
                'displayStrongest', display_percentage, ...
                'colormap2use',custom_colormap, ...
                'nonlinear_bg',2, ...
                'AutoColorBarLimits',[-absolute_max absolute_max], ...
                'FontSize', fontsize, ...
                'showColorbar', false);
            axis off;
        end

        if class1 == 1 && class2 == class1 + 1
            axis on;
        end

        if class1 == 2 && class2 == 2
            nexttile(counter);
            % Turn off the axes
            axis off;
            % Draw the colorbar
            c = colorbar('FontSize',fontsize,'Location','east');
            % Specify the colormap
            colormap(gca,custom_colormap)
            % Define the colormap limits.
            caxis([-absolute_max absolute_max]);
            % Label the colorbar
            c.Label.String = 'PCA+LDA weightings';
            c.Ticks = [-absolute_max 0 absolute_max];
        end
        counter = counter+1;
    end
end

tld.TileSpacing = 'none';
tld.Padding = 'tight';

save_to_file = false;
if save_to_file && saveMe

    suggested_fname = sprintf('pca_projection_maps_multicoder_S%dR%d_gen%s.svg', Session, Run, datetime('now', 'Format', 'yyyyMMdd'));

    [filename, filepath] = uiputfile(fullfile(get_user_data_path('pathType', 'decoding'), suggested_fname));
    if any([filename filepath])
        saveas(fig_pcaProjection, fullfile(filepath, filename));
    end
end

%% Look at conjunction
PCA_LDA_ind = zeros(yPix, xPix);

for class1 = 1:7
    for class2 = 2:8
        if class2 > class1
            current_boundary = average_across_time(:, :, class1, class2);

            PCA_LDA_ind = PCA_LDA_ind + abs(current_boundary);
        end
    end
end

figure;
plotDuplexImage(X_img_mm, Z_img_mm, PCA_LDA_ind, angiogram, ...
    'nonlinear_bg',2, ...
    'colormap2use',magma, ...
    'showColorbar', true, ...
    'displayTop', 0.1);
title('The color represents # of plots where that pixel was within most significant X% of image. ');

% Mask the PCA_LDA index

% For threshold where we just add up absolute value of each boundary
threshold = quantile(abs(PCA_LDA_ind(:)), 0.9);
PCA_LDA_mask = PCA_LDA_ind > threshold;

figure;
tiledlayout('flow');
nexttile
histogram(PCA_LDA_ind);
nexttile
imagesc(PCA_LDA_mask);
axis image;

%% Compute Dice-Sorenson with searchlight and GLM analyses
figure_conjunction = figure('Units', 'normalized', 'OuterPosition', [ 0 0 1 1]);
tld = tiledlayout('flow');
title(tld, sprintf('S%dR%d', session_run_list'));
% Colorblind palette
color1 = [100 143 255]./255;
color2 = [255 176 0]./255;
color3 = [220 38 127]./255;

% Specify the pvalue to threshold maps at
pvalue_threshold = 1e-5;

% LDA + GLM
start_dir = pwd;
% Logic for filenaming
% Complete filename up until the date generated
partial_filename = sprintf('glm_data_%sgen*.mat', sprintf('S%dR%d_', session_run_list'));

cd(fullfile(get_user_data_path('pathType', 'glm')));

% Find most recent matching file
file_listing = dir(partial_filename);


% Load this most recent file
if isempty(file_listing)
    error('No appropriate GLM file exists.');
end
    
[~,most_recent_file_index] = max([file_listing.datenum]);
most_recent_file = file_listing(most_recent_file_index).name;
load(most_recent_file, 'F', 'pvalues_corrected', ...
    'angiogram');

% To handle older files where I did not separately store F and other useful
% variables
if ~exist('F', 'var') || ~exist('pvalues_corrected', 'var')
    load(most_recent_file, 'glm_data');
    F = glm_data.F;
    pvalues_corrected = glm_data.p_corrected;

    % Clear from memory since it is large.
    clear glm_data;
end

% Rename for convenience
glm_angiogram = angiogram;
glm_pvalues = pvalues_corrected;
glm_F = F;

% Define GLM mask
glm_mask = ~(glm_pvalues>pvalue_threshold);

% Return to original directory
cd(start_dir);

% Define axes labels for eventual plots
pixelsize = 0.1;
X_img_mm = pixelsize/2 + (0:size(F,2)-1)*pixelsize; % So that the center of the image is at true halfway point, not shifted.
Z_img_mm = pixelsize/2 + (0:size(F,1)-1)*pixelsize;

% compute sorenson dice coefficient
dice_similarity_glm = dice(glm_mask, PCA_LDA_mask);

% Create conjunction mask where each overlay is a different color
conjunction_mask = nan(size(glm_mask));
conjunction_mask(glm_mask) = 2; % For consistency with other dice-sorenson plots in directional tuning paper, keep this as 2.
conjunction_mask(PCA_LDA_mask) = 1;
conjunction_mask(glm_mask & PCA_LDA_mask) = 3;


nexttile(tld);
plotDuplexImage(X_img_mm, Z_img_mm, conjunction_mask, angiogram,...
    'colormap2use', [color1; color2; color3;], 'nonlinear_bg', 2, ...
    'showColorbar', false, ...
    'colorbarTitle', 'conjunction',...
    'darkestShade', 0.25, ...
    'AutoColorBarLimits', [1 3]);
title(sprintf('GLM vs PCA+LDA \nSorenson Dice Coefficient: %0.10f', ...
    dice_similarity_glm));


% Create specified legend
hold on;
c1 = scatter(nan, nan, 'filled', 'MarkerFaceColor', color1);
c2 = scatter(nan, nan, 'filled', 'MarkerFaceColor', color2);
c3 = scatter(nan, nan, 'filled', 'MarkerFaceColor', color3);
hold off
legend([c1, c2, c3], {'PCA+LDA', 'GLM', 'Both'}, 'Location', 'southwest');


% LDA + searchlight Dice-Sorenson

display_type = 'pvalue_threshold'; % 'bottom_voxels' or 'top_voxels' Or 'pvalue_threshold'
display_accuracy_or_angularError = 'angular_error'; % 'angular_error' or 'accuracy'
radii_to_test = 2;

% Motion correction method
mc_method = 'normcorre';

if mc_method
    mc_method_string = mc_method;
else
    mc_method_string = 'false';
end

partial_filename = sprintf('data_*S%dR%d*.mat', session_run_list');

cd(fullfile(get_user_data_path('pathType', 'decoding')));

% Find most recent matching file
file_listing = dir(partial_filename);
% Load this most recent file
if isempty(file_listing)
    error('No appropriate searchlight file exists.');
end
[~,most_recent_file_index] = max([file_listing.datenum]);
most_recent_file = file_listing(most_recent_file_index).name;

% Load the data
decoding_results = load(most_recent_file);
cd(start_dir);


% Extract `empiric_mean_angular_error` from `result_combined`
[sl_pvalue_angular_error, sl_empiric_mean_angular_error, ...
    sl_accuracy, sl_pvalue] = deal(NaN(size(decoding_results.sl_result_combined)));
for i = 1:size(decoding_results.sl_result_combined, 1)
    for j = 1:size(decoding_results.sl_result_combined, 2)
        for k = 1:size(decoding_results.sl_result_combined, 3)
            try
                sl_empiric_mean_angular_error(i, j, k) = decoding_results.sl_result_combined{i, j, k}.angular_error;
                sl_accuracy(i, j, k) = decoding_results.sl_result_combined{i, j, k}.percentCorrect;
                sl_pvalue(i, j, k) = decoding_results.sl_result_combined{i, j, k}.p;
                sl_pvalue_angular_error(i, j, k) = decoding_results.sl_result_combined{i, j, k}.angular_error_pvalue;
            catch
                sl_empiric_mean_angular_error(i, j, k) = NaN;
                sl_accuracy(i, j, k) = NaN;
                sl_pvalue(i, j, k) = NaN;
                sl_pvalue_angular_error(i, j, k) = NaN;
            end
        end
    end
end


radii_index = find(radii_to_test ==2);



% Apply FDR correction for the p-values
pvalues_1D = reshape(sl_pvalue,[],1);
noNaN_ind = isnan(pvalues_1D);
[Q] = mafdr(pvalues_1D(~noNaN_ind), 'BHFDR', true);
pvalues_corrected = NaN(size(pvalues_1D));
pvalues_corrected(~noNaN_ind) = Q;
matrix_size = size(sl_pvalue);
sl_pvalue_corrected = reshape(pvalues_corrected, matrix_size);

% Apply FDR correction for the angular error p-values
pvalues_1D = reshape(sl_pvalue_angular_error,[],1);
noNaN_ind = isnan(pvalues_1D);
[Q] = mafdr(pvalues_1D(~noNaN_ind), 'BHFDR', true);
pvalues_corrected = NaN(size(pvalues_1D));
pvalues_corrected(~noNaN_ind) = Q;
matrix_size = size(sl_pvalue_angular_error);
sl_pvalue_angular_error_corrected = reshape(pvalues_corrected, matrix_size);

colormap_to_use = inferno;
% Pick appropriate SL dataset

switch display_accuracy_or_angularError
    case 'angular_error'
        pvalue_to_use = sl_pvalue_angular_error_corrected;
        performance_metric_to_use = sl_empiric_mean_angular_error * 180/pi;
        colorbar_title = 'Angular error (deg)';
        unit_measure = ' deg';
        colorscale = flipud(colormap_to_use);

        colorbar_max = 90;
        colorbar_min = 30;
    case 'accuracy'
        pvalue_to_use = sl_pvalue_corrected;
        performance_metric_to_use = sl_accuracy;
        colorbar_title = 'Performance (% correct)';
        unit_measure = '%';
        colorscale = colormap_to_use;

        colorbar_max = 40;
        colorbar_min = 0;
end


% Overlay performance on anatomy
switch display_type
    case 'pvalue_threshold'
        pvalue_mask = pvalue_to_use(:, :, radii_index) < pvalue_threshold;
        performance_masked = squeeze(performance_metric_to_use(:, :, radii_index));
        performance_masked(~pvalue_mask) = NaN;
        %fprintf('thresholded at p<=%d\n', pvalue_threshold);
    case 'top_voxels'
        % Find top X% of voxels
        performance = performance_metric_to_use(:, :, radii_index);
        pvalue = pvalue_to_use(:,:, radii_index);
        top_quantile_threshold = quantile(performance(:), 1-percent_of_voxels_to_keep);
        performance_mask = performance >= top_quantile_threshold;
        performance_masked = performance;
        performance_masked(~performance_mask) = NaN;
        fprintf('keeping top %d%% of voxels\n', percent_of_voxels_to_keep*100);

        % Report the p-value associated with the respective quantile thresholds
        pvalue_of_quantile_threshold = max(pvalue(performance == top_quantile_threshold),[], 'all');
        fprintf('pvalue_threshold - %d\n', pvalue_of_quantile_threshold);
    case 'bottom_voxels'
        % Find top X% of voxels
        performance = performance_metric_to_use(:, :, radii_index);
        pvalue = pvalue_to_use(:,:, radii_index);
        bottom_quantile_threshold = quantile(performance(:), percent_of_voxels_to_keep);
        performance_mask = performance <= bottom_quantile_threshold;
        performance_masked = performance;
        performance_masked(~performance_mask) = NaN;
        title_string = sprintf('keeping top %d%% of voxels\n', percent_of_voxels_to_keep*100);

        % Report the p-value associated with the respective quantile thresholds
        pvalue_of_quantile_threshold = max(pvalue(performance == bottom_quantile_threshold),[], 'all');
        fprintf('pvalue_threshold - %d\n', pvalue_of_quantile_threshold);
end

% Update colorbar limits
colorbar_max = ceil(max([colorbar_max; performance_masked(:)], [], 'all'));
colorbar_min = floor(min([colorbar_min; performance_masked(:)], [], 'all'));

% Create mask from searchlight
searchlight_map = performance_masked;
searchlight_mask = ~isnan(searchlight_map);

% compute sorenson dice coefficient
dice_similarity_searchlight = dice(searchlight_mask, PCA_LDA_mask);

% Create conjunction mask where each overlay is a different color
conjunction_mask = nan(size(searchlight_mask));
conjunction_mask(searchlight_mask) = 2;
conjunction_mask(PCA_LDA_mask) = 1;
conjunction_mask(searchlight_mask & PCA_LDA_mask) = 3;

nexttile(tld);
plotDuplexImage(X_img_mm, Z_img_mm, conjunction_mask, angiogram,...
    'colormap2use', [color1; color2; color3], 'nonlinear_bg', 2, ...
    'showColorbar', false, ...
    'colorbarTitle', 'conjunction',...
    'darkestShade', 0.25, ...
    'AutoColorBarLimits', [1 3]);
title(sprintf('Searchlight vs PCA+LDA \nSorenson Dice Coefficient: %0.10f', ...
    dice_similarity_searchlight));

% Create specified legend
hold on;
c1 = scatter(nan, nan, 'filled', 'MarkerFaceColor', color1);
c2 = scatter(nan, nan, 'filled', 'MarkerFaceColor', color2);
c3 = scatter(nan, nan, 'filled', 'MarkerFaceColor', color3);
hold off
legend([c1, c2, c3], {'PCA+LDA', 'Searchlight', 'Both'}, 'Location', 'southwest');

