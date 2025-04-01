%% Examine all different timewindows for searchlight analysis
% January 14, 2025 WG

%% Initialize workspace
clear; clc; close all;

startDir = pwd;

% choose K for K-fold validation
sl_K = 10;
% preprocessing
timeGain = false;
zScore = false;
detrend = true;

diskFilter = 2;

% set parameters for computing the mask (if requested)
saveMe = false;

motion_correction_method = 'normcorre';

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


    [fix, memory_end, mov_end] = getEpochs(behavior);



[yPix, xPix, nWindows, nTrials] = size(iDopP);

%% ROI Selection
ROI_mask = load_ROI(Session, Run, size(angiogram), ...
    'verbose', true', ...
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
% 3 = Positive (up or right)
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


fprintf('\n ||| ||| ||| --- TESTING SESSION %i RUN %i --- ||| ||| ||| \n\n',Session,Run)

% set up timing for train/test sets
fixTime = ceil(abs(fix(1))*coreParams.framerate);
memTime = ceil(abs(memory_end(end))*coreParams.framerate);

variance_to_keep = 95;

%% Searchlight analysis
% Ran prediction at single voxel level. Then
% display a map showing prediction ability of individual pixels.

timepoints_to_test = 1:nWindows;
n_timepoints = length(timepoints_to_test);

split_on_sulcus = true;
remove_sulcus_pixels = true;

if split_on_sulcus || remove_sulcus_pixels
    sulcus_vertices = get_sulcus_info(Session, Run, angiogram, UF, ...
        'verbose', true, ...
        'use_existing_sulcus_map', true);

    sulcus_pixels = poly2mask(sulcus_vertices(:, 1), sulcus_vertices(:, 2), yPix, xPix);
end

radii_to_test = [2];

[columnsInImage, rowsInImage] = meshgrid(1:xPix, 1:yPix);

% For speed purposes, load and pass the null distribution.
% To make smaller, we can shorten to length of trainLabels. This
% will still be longer than necessary, but sets an upper bound on
% possible # of trials.
n_trials = length(train_labels);
null_mean_angular_error = load_null_mean_angular_error(nClasses, num_repetitions, n_trials);

[sl_percentCorrect_horz, sl_percentCorrect_vert, ...
    sl_percentCorrect_combined, sl_pvalue_combined, ...
    sl_angularError_combined, sl_angularError_pvalue_combined] = deal(NaN(size(angiogram, 1), size(angiogram, 2), n_timepoints, length(radii_to_test)));

[sl_confusion_horz, sl_confusion_vert, sl_confusion_combined, sl_result_horz, sl_result_vert, sl_result_combined] = deal(cell(size(angiogram, 1), size(angiogram, 2), n_timepoints, length(radii_to_test)));


for timepoint = 1:n_timepoints
    timepoint
    trial_timepoint = timepoints_to_test(timepoint);
    % use all of the data during the fixation period...
    if trial_timepoint <= fixTime
        timeOfInterest = 1:trial_timepoint;
        % or use all the data during the memory delay and movement
    else
        timeOfInterest = fixTime:trial_timepoint;
    end


    trainData = flattenDoppler2D(iDopP, timeOfInterest);  % nImages x (nPixels*yPixels)
    if size(trainData, 1) == 1 && size(trainData, 2) > 1
        trainData = trainData';
    end


    p = gcp; % Initialize the parallel pool.
    fprintf('Searchlight analysis started at: %s\n', datestr(now()));
    for k = 1:length(radii_to_test)
        radius = radii_to_test(k);
        fprintf('Searchlight analysis for radius %0.2f ... \n', radius);

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

                    mask_ind = false(yPix, xPix);
                    mask_ind(ROI_pixels) = true;

                    mask_ind_1D = reshape(mask_ind, 1, []);
                    mask_ind_1D = repmat(mask_ind_1D, 1, length(timeOfInterest));

                    trainData_sl = trainData(:,mask_ind_1D);


                    % cross validated performance is calculated here
                    [cp_horz, p_horz, cp_vert, p_vert, cp_combined, p_combined, custom_counting_matrix] = crossValidate_multicoder(trainData_sl, train_labels, ...
                        'classificationMethod', classifierString,...
                        'validationMethod', 'kFold', 'K', sl_K, ...
                        'm', m,...
                        'variance_to_keep', variance_to_keep);

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
                    [empiric_mean_angular_error, ~, pvalue_empiric, ste_empiric] = calculate_angular_error(predicted_class, actual_class, ...
                        'ignore_middle_points', false, ...
                        'num_repetitions', num_repetitions, ...
                        'null_mean_angular_error', null_mean_angular_error);

                    % grabbing useful values
                    sl_percentCorrect_horz(j,i, timepoint, k) = cp_horz.CorrectRate*100;
                    sl_confusion_horz{j,i, timepoint, k} = cp_horz.CountingMatrix;

                    sl_percentCorrect_vert(j,i, timepoint, k) = cp_vert.CorrectRate*100;
                    sl_confusion_vert{j,i, timepoint, k} = cp_vert.CountingMatrix;

                    sl_percentCorrect_combined(j,i, timepoint, k) = cp_combined.CorrectRate*100;
                    sl_pvalue_combined(j, i, timepoint, k) = p_combined;
                    sl_confusion_combined{j,i, timepoint, k} = custom_counting_matrix;
                    sl_angularError_combined(j, i, timepoint, k) = empiric_mean_angular_error;
                    sl_angularError_pvalue_combined(j, i, timepoint, k) = pvalue_empiric;

                    % save the data to results
                    sl_result_horz{j,i, timepoint, k}.percentCorrect =sl_percentCorrect_horz(j,i, timepoint, k);
                    sl_result_horz{j,i, timepoint, k}.confusion = sl_confusion_horz{j,i,timepoint, k};
                    sl_result_horz{j,i, timepoint, k}.p = p_horz;

                    sl_result_vert{j,i, timepoint, k}.percentCorrect = sl_percentCorrect_vert(j,i,timepoint, k);
                    sl_result_vert{j,i, timepoint, k}.confusion = sl_confusion_vert{j,i,timepoint, k};
                    sl_result_vert{j,i, timepoint, k}.p = p_vert;

                    sl_result_combined{j,i, timepoint, k}.percentCorrect = sl_percentCorrect_combined(j,i,timepoint, k);
                    sl_result_combined{j,i, timepoint, k}.confusion = sl_confusion_combined{j,i,timepoint, k};
                    sl_result_combined{j,i, timepoint, k}.p = p_combined;
                    sl_result_combined{j, i, timepoint, k}.angular_error = sl_angularError_combined(j, i, timepoint, k);
                    sl_result_combined{j, i, timepoint, k}.angular_error_pvalue = sl_angularError_pvalue_combined(j, i, timepoint, k);
                    
                    ppm.count;
                end
            end
        end
        delete(ppm);
        fprintf('Completed %s\n', datestr(now()));
    end
end

%% Extract `empiric_mean_angular_error` from `result_combined`
% This is unnecessary if running entire script, but added so that you can
% load previously saved data instead.
[sl_pvalue_angular_error, sl_empiric_mean_angular_error, ...
    sl_accuracy, sl_pvalue] = deal(NaN(size(sl_result_combined)));
for i = 1:size(sl_result_combined, 1)
    for j = 1:size(sl_result_combined, 2)
        for k = 1:size(sl_result_combined, 3)
            for m = 1:size(sl_result_combined, 4)
                try
                    sl_empiric_mean_angular_error(i, j, k, m) = sl_result_combined{i, j, k, m}.angular_error;
                    sl_accuracy(i, j, k, m) = sl_result_combined{i, j, k, m}.percentCorrect;
                    sl_pvalue(i, j, k, m) = sl_result_combined{i, j, k, m}.p;
                    sl_pvalue_angular_error(i, j, k, m) = sl_result_combined{i, j, k, m}.angular_error_pvalue;
                catch
                    sl_empiric_mean_angular_error(i, j, k, m) = NaN;
                    sl_accuracy(i, j, k, m) = NaN;
                    sl_pvalue(i, j, k, m) = NaN;
                    sl_pvalue_angular_error(i, j, k, m) = NaN;
                end
            end
        end
    end
end

%% Visualize searchlight analysis
pvalue_threshold = 1e-3;

n_timepoints = size(sl_result_combined, 3);

trialTime = fix(1):1/coreParams.framerate:fix(1)+nWindows/coreParams.framerate-1/coreParams.framerate;


pixelsize = 0.1;
X_img_mm = pixelsize/2 + (0:size(angiogram,2)-1)*pixelsize; % So that the center of the image is at true halfway point, not shifted.
Z_img_mm = pixelsize/2 + (0:size(angiogram,1)-1)*pixelsize + UF.Depth(1);

colorbar_limits = [0 ceil(max(sl_accuracy,[], 'all'))];

circle_center = [ 2 max(Z_img_mm)-2];
for k = 1:length(radii_to_test)

    pvalues_1D = reshape(sl_pvalue(:, :, :, k),[],1);
    nan_ind = isnan(pvalues_1D);
    [Q] = mafdr(pvalues_1D(~nan_ind), 'BHFDR', true);
    pvalues_corrected = NaN(size(pvalues_1D));
    pvalues_corrected(~nan_ind) = Q;
    matrix_size = [size(sl_pvalue, 1) size(sl_pvalue, 2) size(sl_pvalue, 3)];
    sl_pvalue_combined_corrected_all = reshape(pvalues_corrected, matrix_size);

    radius = radii_to_test(k);
    % Plot searchlight analysis results
    fig_searchlight = figure('units','normalized','outerposition',[0 0 1 1]);
    tld = tiledlayout('flow');
        

    for timepoint = 1:n_timepoints
        trial_timepoint = timepoints_to_test(timepoint);

        when_in_trial = trialTime(trial_timepoint);
        if when_in_trial < 0
            task_state = 'fixation';
        elseif when_in_trial == 0
            task_state = 'cue';
        elseif when_in_trial > 0 && when_in_trial < memory_end(2)
            task_state = 'memory';
        elseif when_in_trial > memory_end(2) && when_in_trial < mov_end(2)
            task_state = 'movement';
        elseif when_in_trial > mov_end(2)
            task_state = 'ITI';
        end


        sl_pvalue_combined_corrected = sl_pvalue_combined_corrected_all(:, :, timepoint);

    
        % Overlay performance on anatomy
        pvalue_mask = sl_pvalue_combined_corrected < pvalue_threshold;

        
        performance_masked = squeeze(sl_accuracy(:, :, timepoint, k));
        performance_masked(~pvalue_mask) = NaN;
    
        nexttile(tld, timepoint);
        if timepoint ==n_timepoints
            show_colorbar = true;
        else
            show_colorbar = false;
        end
        cmap = plotDuplexImage(X_img_mm, Z_img_mm, performance_masked, angiogram,...
            'colormap2use',magma, 'nonlinear_bg',2, ...
            'AutoColorBarLimits',colorbar_limits, ...
            'showColorbar', show_colorbar, ...
            'ColorbarTitle', 'Performance (% correct)');
        axis off;
        title(sprintf('%0.2f sec - %s', when_in_trial, task_state));
    
        hold on;
        circle_handle = viscircles(circle_center,radius*pixelsize, 'EnhanceVisibility', false, ...
            'Color', 'w');
        hold off;
    end
    title(tld, sprintf('Searchlight analysis S%dR%d (radius = %0.2f)', Session, Run, radius))
end

save_to_file = true;
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
    suggested_fname = sprintf('searchlight_data_S%dR%d_mc%s_gen%s.mat', Session, Run, motion_correction_method, datetime('today', format='yyyyMMdd'));
    save_fullfile = fullfile(get_user_data_path('pathType', 'decoding'), suggested_fname);

    % save the results file
    fprintf('Saving %s to %s...',suggested_fname,pwd)
    save(save_fullfile, 'classifierString', 'validationString', 'ROI_mask',...
        'detrend', 'num_repetitions',...
         'diskFilter', 'sessionRunList', ...
        'sl_result_combined', 'timepoints_to_test',... 
        'radii_to_test', 'sl_K', 'validationString');
    fprintf('done.\n')
else
    disp('Session data not saved!')
end

