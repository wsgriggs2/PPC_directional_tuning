%% Script to run offline decoding analysis for single combination of concatenated sessions/runs
% Author: Whitney Griggs
% Date: 2021


%% Initialize workspace
clear; clc; close all;

startDir = pwd;

%% hard coded settings
% choose K for K-fold validation
K = 10;

% preprocessing
timeGain = false;
zScore = false;
detrend = true;
diskFilter = 2;
motion_correction_method = 'normcorre';

% Do you want to save out results and figures?
saveMe = true;

%% get classifier string
if ~exist('classifier_string','var') || ~exist('validation_string','var')
    [classifierString, validationString] = uigetclassifier;
end


%% get data, preprocess, get epochs & dimension sizes
loadDopplerData(whos, 'multiple', false, ...
    'mc_method', motion_correction_method, ...
    'project_record_filename', 'ProjectRecord_paper.json');
if detrend
    iDopP = detrend_sliding_window(iDop, 50);
else
    iDopP = iDop;
end
iDopP = preProcess(iDopP, 'timeGain', timeGain, 'spatialFilter', {'disk', diskFilter}, 'zScore', zScore);

% Clear this variable to reduce memory usage
clear iDop

% Get timing information from behavior
[fix, memory, mov] = getEpochs(behavior);


[yPix, xPix, nWindows, nTrials] = size(iDopP);
pixelsize = 0.1;
X_img_mm = pixelsize/2 + (0:size(angiogram,2)-1)*pixelsize; % So that the center of the image is at true halfway point, not shifted.
Z_img_mm = pixelsize/2 + (0:size(angiogram,1)-1)*pixelsize + UF.Depth(1);

%% ROI Selection
% If wanted, load user-defined ROIs
ROI_mask = load_ROI(Session, Run, size(angiogram), ...
    'verbose', false, ...
    'angiogram', angiogram);

%% Flatten Doppler for easier masking
% nPixels x nWindows x nTrials
% This allows us to eliminate specific voxels that we don't want. Supports
% non-rectangular ROIs.
iDop_3D = reshape(iDopP, [yPix*xPix, nWindows, nTrials]);

%% implement ROIs, masks, cropping
% first, flatten the ROI mask
ROI_mask_flat = reshape(ROI_mask, yPix*xPix, 1);

% Apply ROI mask
iDop_crop = iDop_3D(ROI_mask_flat, :, :);

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

% Decide if we want to use the center positions.
useMiddlePosition = true;

% create the training data labels (2 x n_trials) - x, y columns
% For horizontal and vertical axes, use the following labels
% 1 = Negative (down or left)
% 2 = At center (along vertical or horizontal axis)
% 3  = Positive (up or right)

label_strings = {'Negative','Center','Positive'};

train_labels = NaN(size(targetPos,1),2);

for dimension = 1:2
    train_labels(targetPos(:, dimension) < 0, dimension) = 1;
    train_labels(targetPos(:, dimension) == 0, dimension) = 2;
    train_labels(targetPos(:, dimension) > 0, dimension) = 3;
end

%% run cross validation for each sample in the memory epoch
num_repetitions = 100000;
try
    load(sprintf('angular_error_null_distribution_%dtgts_%dreps.mat', nClasses, num_repetitions));
    null_mean_angular_error = null_mean_angular_error(:, 1:length(train_labels));
catch
    choice = questdlg(sprintf('This combo of %d Classes and %d num_repetitions did not exist. Do you want to create now?', nClasses, num_repetitions),...
        'Create new null distribution?','yes','no','yes');
    if strcmp(choice,'yes')
        [~, null_mean_angular_error] = generate_null_distribution_angular_error(nClasses, num_repetitions);
    else
        error('The distribution does not exist currently. Aborting.\n');
    end
    null_mean_angular_error = null_mean_angular_error(:, 1:length(train_labels));
end

% Determine dimensionality of cPCA subspace (if used)
% For now, using (number of classes - 1), i.e m = 2 for 3 classes
m = length(unique(train_labels)) - 1;

% allocate memory
[result_horz, confusion_horz, cp_horz] = deal(cell(1,nWindows));   % one cell for each data point in memory epoch
[p_horz, percentCorrect_horz, t_horz] = deal(NaN(1,nWindows));     % one array point for each point in mem epoch

[result_vert, confusion_vert, cp_vert] = deal(cell(1,nWindows));   % one cell for each data point in memory epoch
[p_vert, percentCorrect_vert] = deal(NaN(1,nWindows));     %

[result_combined, confusion_combined, cp_combined] = deal(cell(1,nWindows));   % one cell for each data point in memory epoch
[p_combined, percentCorrect_combined, empiric_mean_angular_error] = deal(NaN(1,nWindows));     %


fprintf('\n ||| ||| ||| --- TESTING SESSION %i RUN %i --- ||| ||| ||| \n\n',Session,Run)

% set up timing for train/test sets
fixTime = ceil(abs(fix(1))*coreParams.framerate);
memTime = ceil(abs(memory(end))*coreParams.framerate);

disp('Time with for loop');
tic
variance_to_keep = 95;
for i = 1:nWindows
    % this is the meat of this whole loop where we continue to pile on data
    % (concatenate it, really) as we progress through time. This is on the
    % underlying assumption that the classifier is going to use data from
    % whatever epoch of the task we're currently in to classify data from
    % that same epoch.

    % use all of the data during the fixation period...
    if i <= fixTime
        testTimeI = 1:i;
    % or use all the data during the memory delay and movement
    else
        testTimeI = fixTime:i;
    end

    % prepare the training data for this time window
    trainData = flattenDoppler2D(iDop_crop, testTimeI);  % nImages x (nPixels*yPixels)

    % cross validated performance is calculated here
    [cp_horz{i}, p_horz(i), cp_vert{i}, p_vert(i), cp_combined{i}, p_combined(i), custom_counting_matrix] = crossValidate_multicoder(trainData, train_labels, ...
        'classificationMethod', classifierString,...
        'validationMethod', validationString, ...
        'K', K, ...
        'm', m,...
        'variance_to_keep', variance_to_keep);

    % grabbing useful values
    percentCorrect_horz(i) = cp_horz{i}.CorrectRate*100;
    confusion_horz{i} = cp_horz{i}.CountingMatrix;

    percentCorrect_vert(i) = cp_vert{i}.CorrectRate*100;
    confusion_vert{i} = cp_vert{i}.CountingMatrix;

    percentCorrect_combined(i) = cp_combined{i}.CorrectRate*100;
    confusion_combined{i} = custom_counting_matrix;

    % Reconstruct a record for final performance. This
    % is not the actual order of predictions since that
    % is not easily accessible. This is correct number
    % of each prediction though.
    [actual_class, predicted_class] = deal([]);
    for true_class_num = 1:9
        for predicted_class_num = 1:9
            actual_class = [actual_class; true_class_num*ones(custom_counting_matrix(predicted_class_num,true_class_num), 1)];
            predicted_class = [predicted_class; predicted_class_num*ones(custom_counting_matrix(predicted_class_num, true_class_num), 1)];
        end
    end
    % Calculate mean_angular_error for each timepoint
    [empiric_mean_angular_error(i), ~, pvalue_empiric(i), ste_angular_error(i)] = calculate_angular_error(predicted_class, actual_class, ...
        'ignore_middle_points', false, ...
        'num_repetitions', num_repetitions, ...
        'null_mean_angular_error', null_mean_angular_error);

    % save the data to results
    result_horz{i}.percentCorrect = percentCorrect_horz(i);
    result_horz{i}.confusion = confusion_horz{i};
    result_horz{i}.p = p_horz(i);
    result_horz{i}.t = i-1;

    result_vert{i}.percentCorrect = percentCorrect_vert(i);
    result_vert{i}.confusion = confusion_vert{i};
    result_vert{i}.p = p_vert(i);
    result_vert{i}.t = i-1;

    result_combined{i}.percentCorrect = percentCorrect_combined(i);
    result_combined{i}.confusion = confusion_combined{i};
    result_combined{i}.p = p_combined(i);
    result_combined{i}.t = i-1;
    result_combined{i}.angular_error = empiric_mean_angular_error(i);
    result_combined{i}.ste_angular_error = ste_angular_error(i);
    result_combined{i}.angular_error_pvalue = pvalue_empiric(i);


    t_horz(i) = i-1;
    i
end
toc

%% save out the best results
[bestPercentCorrect_horz, bestIndex_horz] = max(percentCorrect_horz);
bestConfusion_horz = confusion_horz{bestIndex_horz};
bestP_horz = p_horz(bestIndex_horz);

[bestPercentCorrect_vert, bestIndex_vert] = max(percentCorrect_vert);
bestConfusion_vert = confusion_vert{bestIndex_vert};
bestP_vert = p_vert(bestIndex_vert);

[bestPercentCorrect_combined, bestIndex_combined] = max(percentCorrect_combined);
bestConfusion_combined = confusion_combined{bestIndex_combined};
bestP_combined = p_combined(bestIndex_combined);

%% Plot confusion matrices for horizontal, vertical, and combined
f_confusion = figure('units','normalized','OuterPosition',[0 0 1 1]);
t = tiledlayout(1,3);
title(t, sprintf('S%dR%d - Confusion charts for multi-target multicoder using %s', Session, Run, classifierString), 'FontSize', 16, 'FontWeight', 'bold')
nexttile
imagesc_confusion_matrix(bestConfusion_horz(1:end-1,:)', ...
    'label_strings', label_strings, ...
    'title', 'Horizontal prediction', ...
    'FontColor', 'k', ...
    'show_colorbar', false, ...
    'superimpose_text', 'percentage', ...
    'row_normalize', true);
set(gcf,'color','w');
axis image;

nexttile;
imagesc_confusion_matrix(bestConfusion_vert(1:end-1,:)', ...
    'label_strings', label_strings, ...
    'title', 'Vertical prediction', ...
    'FontColor', 'k', ...
    'show_colorbar', false, ...
    'superimpose_text', 'percentage', ...
    'row_normalize', true);
set(gcf,'color','w');
axis image


nexttile;

combined_label_strings = {'Ipsilateral Down', 'Center Down', 'Contralateral Down', ...
    'Ipsilateral Center',  'Middle', 'Contralateral Center', ...
    'Ipsilateral Up', 'Center Up', 'Contralateral Up'};

reordered_label_strings = combined_label_strings([7 8 9 6 3 2 1 4 5]);
picture_label_strings = {char(8598), char(8593), char(8599), char(8594), ...
    char(8600), char(8595), char(8601), char(8592), char(8857)};


reordered_confusion_mat = bestConfusion_combined([7 8 9 6 3 2 1 4 5], [7 8 9 6 3 2 1 4 5]);


imagesc_confusion_matrix(reordered_confusion_mat(:, 1:8)', ...
    'label_strings', picture_label_strings, ...
    'title', 'Combined prediction', ...
    'FontColor', 'k', ...
    'show_colorbar', true, ...
    'superimpose_text', 'percentage', ...
    'row_normalize', true);
set(gcf,'color','w');
axis image;

% Save to file if requested
save_to_file = false;
if save_to_file && saveMe

    suggested_fname = sprintf('confusion_charts_S%dR%d_gen%s.svg', Session, Run, datetime('today', format='yyyyMMdd'));

    [filename, filepath] = uiputfile(fullfile(get_user_data_path('pathType', 'decoding'), suggested_fname));

    if any([filename filepath])
        saveas(f_confusion, fullfile(filepath, filename));
    end
end

%% ------------------------------ DISPLAYING / PLOTTING / SAVING ------------------------------ %%
% plot percent correct
trialTime = fix(1):1/coreParams.framerate:fix(1)+nWindows/coreParams.framerate-1/coreParams.framerate;


f_accuracy_across_time = figure();
set(f_accuracy_across_time,'Position',[800 500 700 500])

% Plot horizontal performance across time
subplot(1, 3, 1)
chance_level = 100/size(unique(train_labels(:,1)),1);
titleString = ['S' num2str(Session) 'R' num2str(Run) ', ' classifierString ', horizontal-only ' validationString];

plot_decoder_performance_across_time(trialTime, percentCorrect_horz, ...
    p_horz, chance_level, ...
    mov(1)*coreParams.framerate, ...
    'p_threshold', 0.05, ...
    'title', titleString, ...
    'ylim', [0 100]);

% For vertical plots
subplot(1, 3, 2)
chance_level = 100/size(unique(train_labels(:,2)),1);
titleString = ['S' num2str(Session) 'R' num2str(Run) ', ' classifierString ', vertical-only ' validationString];

plot_decoder_performance_across_time(trialTime, percentCorrect_vert, ...
    p_vert, chance_level, ...
    mov(1)*coreParams.framerate, ...
    'p_threshold', 0.05, ...
    'title', titleString, ...
    'ylim', [0 100]);

% For combined plot
subplot(1, 3, 3);

chance_level = 100/(size(unique(train_labels(:,2)),1)*size(unique(train_labels(:,1)),1));

% For directional tuning paper, using modified conservative estimate of
% chance level due to class imbalance, i.e., ideal chance is 11.11% but due
% to no middle category, chance is actually closer to 12.5%
chance_level = 100/length(unique(TargetPosInd));

titleString = ['S' num2str(Session) 'R' num2str(Run) ', ' classifierString ', horz+vert combined ' validationString];

plot_decoder_performance_across_time(trialTime, percentCorrect_combined, ...
    p_combined, chance_level, ...
    mov(1)*coreParams.framerate, ...
    'p_threshold', 0.05, ...
    'title', titleString, ...
    'ylim', [0 75]);

% Print output to console
% display selected final results for horizontal
fprintf('\n ||| FINAL RESULTS (selected) ||| \n')
fprintf('\n\tbest accuracy: %2.1f%%\n', bestPercentCorrect_horz)
fprintf('\tbest p: %i\n', bestP_horz)
fprintf('\tbest results are from time index %i (trial time %2.1f s)\n\n', bestIndex_horz, bestIndex_horz+fix(1))

% display selected final results for vertical
fprintf('\n ||| FINAL RESULTS (selected) ||| \n')
fprintf('\n\tbest accuracy: %2.1f%%\n', bestPercentCorrect_vert)
fprintf('\tbest p: %i\n', bestP_vert)
fprintf('\tbest results are from time index %i (trial time %2.1f s)\n\n', bestIndex_vert, bestIndex_vert+fix(1))

% display selected final results for combined
fprintf('\n ||| FINAL RESULTS (selected) ||| \n')
fprintf('\n\tbest accuracy: %2.1f%%\n', bestPercentCorrect_combined)
fprintf('\tbest p: %i\n', bestP_combined)
fprintf('\tbest results are from time index %i (trial time %2.1f s)\n\n', bestIndex_combined, bestIndex_combined+fix(1))

% Save to file if requested
save_to_file = false;
if save_to_file && saveMe

    suggested_fname = sprintf('accuracy_S%dR%d_gen%s.svg', Session, Run, datetime('today', format='yyyyMMdd'));

    [filename, filepath] = uiputfile(fullfile(get_user_data_path('pathType', 'decoding'), suggested_fname));
    if any([filename filepath])
        saveas(f_accuracy_across_time, fullfile(filepath, filename));
    end
end

%% Plot angular error
f_angular_error = figure();
set(f_angular_error,'Position',[800 500 700 500])
titleString = ['S' num2str(Session) 'R' num2str(Run) ', ' classifierString ', Angular Error ' validationString];

% plot percent correct
plot_angular_error_across_time(trialTime, empiric_mean_angular_error, ...
    pvalue_empiric, pi/2, ...
    mov(1)*coreParams.framerate, ...
    'p_threshold', 0.05, ...
    'title', titleString, ...
    'ylim', [0 100*pi/180]);

% Save to file if requested
save_to_file = false;
if save_to_file && saveMe

    suggested_fname = sprintf('angular_error_S%dR%d_gen%s.svg', Session, Run, datetime('today', format='yyyyMMdd'));

    [filename, filepath] = uiputfile(fullfile(get_user_data_path('pathType', 'decoding'), suggested_fname));
    if any([filename filepath])
        saveas(summary_fig, fullfile(filepath, filename));
    end
end

%% Plot angular error variance
titleString = ['S' num2str(Session) 'R' num2str(Run) ', ' classifierString ', Angular Error STE ' validationString];

figure;
% plot percent correct
plot_angular_error_across_time(trialTime, ste_angular_error, ...
    pvalue_empiric, pi/2, ...
    mov(1)*coreParams.framerate, ...
    'p_threshold', 0.05, ...
    'title', titleString, ...
    'ylim', [0 10*pi/180]);
ylimits = ylim;
yticks(0:1:ylimits(2));


%% Show confusion matrix and time plots in a single figure;
% Imagesc, subplots (or tiled layouts), and these other plots do not play
% well together. 

summary_fig = figure('units','normalized','outerposition',[0 0 1 1]);
tld = tiledlayout(1, 3);
chance_level = 100/(size(unique(train_labels(:,2)),1)*size(unique(train_labels(:,1)),1));
% For directional tuning paper, using modified conservative estimate of
% chance level due to class imbalance, i.e., ideal chance is 11.11% but due
% to no middle category, chance is actually closer to 12.5%
chance_level = 100/length(unique(TargetPosInd));


titleString = ['S' num2str(Session) 'R' num2str(Run)];

plot_decoder_performance_across_time(trialTime, percentCorrect_combined, ...
    p_combined, chance_level, ...
    mov(1)*coreParams.framerate, ...
    'p_threshold', 0.05, ...
    'title', titleString, ...
    'ylim', [0 75], ...
    'tiledlayout', tld);

plot_angular_error_across_time(trialTime, empiric_mean_angular_error, ...
    pvalue_empiric, pi/2, ...
    mov(1)*coreParams.framerate, ...
    'p_threshold', 0.05, ...
    'title', titleString, ...
    'ylim', [0 100*pi/180], ...
    'tiledlayout', tld);

imagesc_confusion_matrix(reordered_confusion_mat(:, 1:8)', ...
    'label_strings', picture_label_strings, ...
    'title', '', ...
    'FontColor', 'k', ...
    'tiledlayout', tld);
set(gcf,'color','w');

% Save to file if requested
save_to_file = false;
if save_to_file && saveMe

    suggested_fname = sprintf('summary_decoding_results_S%dR%d_gen%s.svg', Session, Run, datetime('today', format='yyyyMMdd'));

    [filename, filepath] = uiputfile(fullfile(get_user_data_path('pathType', 'decoding'), suggested_fname));
    if any([filename filepath])
        saveas(summary_fig, fullfile(filepath, filename));
    end
end

%% Searchlight analysis
% Use the last time point in the trial and run prediction at small group level. Then
% display a map showing prediction ability of individual small groups.

num_repetitions = 100000; % For the angular error
sl_K = K;
sl_validationString = 'kFold';

split_on_sulcus = true;
remove_sulcus_pixels = true;

if split_on_sulcus || remove_sulcus_pixels
    sulcus_vertices = get_sulcus_info(Session, Run, angiogram, UF, ...
        'verbose', true);

    sulcus_pixels = poly2mask(sulcus_vertices(:, 1), sulcus_vertices(:, 2), yPix, xPix);
end

radii_to_test = [1 2];

[columnsInImage, rowsInImage] = meshgrid(1:xPix, 1:yPix);

% For speed purposes, load and pass the null distribution.
% To make smaller, we can shorten to length of trainLabels. This
% will still be longer than necessary, but sets an upper bound on
% possible # of trials.
n_trials = length(train_labels);
null_mean_angular_error = load_null_mean_angular_error(nClasses, num_repetitions, n_trials);

[sl_percentCorrect_horz, sl_percentCorrect_vert, ...
    sl_percentCorrect_combined, sl_pvalue_combined, ...
    sl_angularError_combined, sl_angularError_pvalue_combined] = deal(NaN(size(angiogram, 1), size(angiogram, 2), length(radii_to_test)));

[sl_confusion_horz, sl_confusion_vert, sl_confusion_combined, sl_result_horz, sl_result_vert, sl_result_combined] = deal(cell(size(angiogram, 1), size(angiogram, 2), length(radii_to_test)));

timeOfInterest = fixTime:nWindows;

trainData = flattenDoppler2D(iDopP, timeOfInterest);  % nImages x (nPixels*yPixels)
if size(trainData, 1) == 1 && size(trainData, 2) > 1
    trainData = trainData';
end


p = gcp; % Initialize the parallel pool.
fprintf('Searchlight analysis started at: %s\n', datestr(now()));
for k = 1:length(radii_to_test)
    radius = radii_to_test(k);
    fprintf('Searchlight analysis for radius %0.2f ... ', radius);

    ppm = ProgressBar(yPix * xPix, taskname='Searchlight analysis', no_log=true);


    parfor i = 1:xPix
        for j = 1:yPix
            % Check if (i, j) is within sulcus. If so, then skip.
            if (~split_on_sulcus && ~remove_sulcus_pixels) || ~sulcus_pixels(j, i)
                % prepare the training data for this time window
                ROI_pixels = (rowsInImage - j).^2 ...
                    + (columnsInImage - i).^2 <= radius.^2;

                % If sulcus splits circle pixels, then only keep half with the
                % center pixel.
                if split_on_sulcus
                    ROI_sulcus_overlap = ROI_pixels & sulcus_pixels;
                    if nnz(ROI_sulcus_overlap)
                        ROI_sulcus_exclude = ROI_pixels & ~sulcus_pixels;
                        binary_comps = bwlabel(ROI_sulcus_exclude);
                        center_pixel_comp = binary_comps(j, i);
                        ROI_pixels = binary_comps == center_pixel_comp;
                    end
                elseif remove_sulcus_pixels
                    ROI_pixels = ROI_pixels & ~sulcus_pixels;
                end

                % Define mask
                mask_ind = false(yPix, xPix);
                mask_ind(ROI_pixels) = true;

                % Convert mask into 1D
                mask_ind_1D = reshape(mask_ind, 1, []);
                mask_ind_1D = repmat(mask_ind_1D, 1, length(timeOfInterest));

                % Apply mask to training data
                trainData_sl = trainData(:,mask_ind_1D);

                % cross validated performance is calculated here
                [cp_horz, p_horz, cp_vert, p_vert, cp_combined, p_combined, custom_counting_matrix] = crossValidate_multicoder(trainData_sl, train_labels, ...
                    'classificationMethod', classifierString, ...
                    'validationMethod', 'kFold', ...
                    'K', sl_K, ...
                    'm', m,...
                    'variance_to_keep', 95);

                % Angular error calculation
                % Reconstruct a record for final performance. This
                % is not the actual order of predictions since that
                % is not easily accessible. This is correct number
                % of each prediction though.
                [actual_class, predicted_class] = deal([]);
                for true_class_num = 1:9
                    for predicted_class_num = 1:9
                        actual_class = [actual_class; true_class_num*ones(custom_counting_matrix(predicted_class_num,true_class_num), 1)];
                        predicted_class = [predicted_class; predicted_class_num*ones(custom_counting_matrix(predicted_class_num, true_class_num), 1)];
                    end
                end
                [empiric_mean_angular_error, ~, pvalue_empiric] = calculate_angular_error(predicted_class, actual_class, ...
                    'ignore_middle_points', false, ...
                    'num_repetitions', num_repetitions, ...
                    'null_mean_angular_error', null_mean_angular_error);

                % grabbing useful values
                sl_percentCorrect_horz(j,i, k) = cp_horz.CorrectRate*100;
                sl_confusion_horz{j,i, k} = cp_horz.CountingMatrix;

                sl_percentCorrect_vert(j,i, k) = cp_vert.CorrectRate*100;
                sl_confusion_vert{j,i, k} = cp_vert.CountingMatrix;

                sl_percentCorrect_combined(j,i, k) = cp_combined.CorrectRate*100;
                sl_pvalue_combined(j, i, k) = p_combined;
                sl_confusion_combined{j,i, k} = custom_counting_matrix;
                sl_angularError_combined(j, i, k) = empiric_mean_angular_error;
                sl_angularError_pvalue_combined(j, i, k) = pvalue_empiric;

                % save the data to results
                sl_result_horz{j,i, k}.percentCorrect =sl_percentCorrect_horz(j,i,k);
                sl_result_horz{j,i, k}.confusion = sl_confusion_horz{j,i,k};
                sl_result_horz{j,i, k}.p = p_horz;

                sl_result_vert{j,i, k}.percentCorrect = sl_percentCorrect_vert(j,i,k);
                sl_result_vert{j,i, k}.confusion = sl_confusion_vert{j,i,k};
                sl_result_vert{j,i, k}.p = p_vert;

                sl_result_combined{j,i, k}.percentCorrect = sl_percentCorrect_combined(j,i,k);
                sl_result_combined{j,i, k}.confusion = sl_confusion_combined{j,i,k};
                sl_result_combined{j,i, k}.p = p_combined;
                sl_result_combined{j, i, k}.angular_error = sl_angularError_combined(j, i, k);
                sl_result_combined{j, i, k}.angular_error_pvalue = sl_angularError_pvalue_combined(j, i, k);
                ppm.count;
            end
        end
    end
    delete(ppm);

    fprintf('Completed %s\n', datestr(now()));
end

%% Visualize searchlight analysis
colorbar_limits = [0 ceil(max(sl_percentCorrect_combined,[], 'all'))];

fig_searchlight = figure('units','normalized','outerposition',[0 0 1 1]);
tld = tiledlayout('flow');
circle_center = [ 2 max(Z_img_mm)-2];
for k = 1:length(radii_to_test)
    radius = radii_to_test(k);

    % Apply FDR correction for the p-values
    pvalues_1D = reshape(sl_pvalue_combined(:, :, k),[],1);
    nan_ind = isnan(pvalues_1D);
    [Q] = mafdr(pvalues_1D(~nan_ind), 'BHFDR', true);
    pvalues_corrected = NaN(size(pvalues_1D));
    pvalues_corrected(~nan_ind) = Q;
    matrix_size = [size(sl_pvalue_combined, 1) size(sl_pvalue_combined, 2)];
    sl_pvalue_combined_corrected = reshape(pvalues_corrected, matrix_size);

    % Overlay performance on anatomy
    pvalue_mask = sl_pvalue_combined_corrected < 0.001;

    performance_masked = squeeze(sl_percentCorrect_combined(:, :, k));
    performance_masked(~pvalue_mask) = NaN;

    nexttile(k);

    cmap = plotDuplexImage(X_img_mm, Z_img_mm, performance_masked, angiogram,...
        'colormap2use',magma, 'nonlinear_bg',2, ...
        'AutoColorBarLimits',colorbar_limits, ...
        'showColorbar', true, ...
        'ColorbarTitle', 'Performance (% correct)');
    title(sprintf('Searchlight analysis (radius = %0.2f) \nS%dR%d', radius, Session, Run));

    hold on;
    circle_handle = viscircles(circle_center,radius*pixelsize, 'EnhanceVisibility', false, ...
        'Color', 'w');
end

save_to_file = false;
if save_to_file && saveMe

    suggested_fname = sprintf('searchlight_results_S%dR%d_gen%s.svg', Session, Run, datetime('today', format='yyyyMMdd'));

    [filename, filepath] = uiputfile(fullfile(get_user_data_path('pathType', 'decoding'), suggested_fname));
    if any([filename filepath])
        saveas(fig_searchlight, fullfile(filepath, filename));
    end
end


%% save to file if desired
saveMe = questdlg('Would you like to save these results to file?',...
    'save results','yes','no','no');
saveMe = strcmp(string(saveMe),'yes');


if saveMe
    % create the filename string
    suggested_fname = sprintf('data_S%dR%d_mc%s_gen%s.mat', Session, Run, motion_correction_method, datetime('today', format='yyyyMMdd'));
    save_fullfile = fullfile(get_user_data_path('pathType', 'decoding'), suggested_fname);

    % move to the output folder
    cd(fullfile(get_user_data_path('pathType','Output'),'decoding'));

    % save the results file
    fprintf('Saving %s to %s...',suggested_fname,pwd)
    save(save_fullfile, 'classifierString', 'validationString', 'cp_horz','result_horz','percentCorrect_horz','confusion_horz','p_horz','t_horz','ROI_mask',...
        'cp_vert','result_vert','percentCorrect_vert','confusion_vert','p_vert', 'detrend',...
        'cp_combined', 'result_combined', 'percentCorrect_combined', 'confusion_combined', 'p_combined', ...
        'diskFilter', 'sessionRunList', ...
        'sl_result_combined',...
        'radii_to_test','sl_K', 'sl_validationString');
    fprintf('done.\n')
else
    disp('Session data not saved!')
end
cd(startDir)
