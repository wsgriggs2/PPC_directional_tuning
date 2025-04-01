%% Script to run offline decoding analysis for all possible combinations of concatenated sessions/runs
% First does cross-validation on each single session, then performs
% training on one session and testing on all other sessions. It iterates
% over each combination of having one session as training and others as
% testing.
%
% Author: Whitney Griggs
% Date: 2021


close all; clc; clear;

startDir = pwd;
%% hard coded settings
% choose K for K-fold validation
K = 10;
% preprocessing
timeGain = false;
diskFilter = 2;
zScore = false;
detrend = true;
% set parameters for computing the mask (if requested)
project_record_filename = 'ProjectRecord_paper.json';

%% get classifier string
if ~exist('classifier_string','var') || ~exist('validation_string','var')
    [classifierString, validationString] = uigetclassifier;
end

%% Specify the sessions we want to analyze
ProjectRecord = load_json_as_table(project_record_filename);
session_run_list = specify_sessions_of_interest(...
    'project_record_filename', project_record_filename);

% Group by session
session_id = NaN(size(session_run_list, 1), 1);

unique_sessions = unique(session_run_list(:, 1));
for session_num = 1:length(unique_sessions)
    session_id(session_run_list(:, 1) == unique_sessions(session_num)) = session_num;
end

%% Load all data
num_sessions = length(unique_sessions);
iDop_cell = cell(num_sessions, 1);

for session = 1:num_sessions
    current_session_ind = session_id == session;

    session_num = session_run_list(current_session_ind,1);
    run_nums = session_run_list(current_session_ind,2);
    indx = ismember(ProjectRecord.Session, session_num) & ismember(ProjectRecord.Run, run_nums);

    Monkey = ProjectRecord{indx, 'Monkey'};
    % If using RT BMI data, then "mistakes" can be just issues with the
    % decoder, so it is worth using all trials that make it to memory
    % period.
    if any(ismember(ProjectRecord{indx, 'Project'}, {'RT BMI'}))
        load_all_attempted_trials = true;
    else
        load_all_attempted_trials = false;
    end

    data = loadDopplerData([], 'multiple', true, ...
        'project_record_filename', project_record_filename, ...
        'sessionRunList', session_run_list(current_session_ind, :), ...
        'load_all_attempted_trials', load_all_attempted_trials, ...
        'mc_method', 'normcorre');

    iDop_cell{session} = data;
end

clear iDop behavior UF Session Run yPixels xPixels nWindows nTrials
clear angiogram coreParams sessionRunList

%% Align all data
define_new_alignment = true;

align_all_to_one_session = false;
if define_new_alignment
    % Do we want to align all to a single session or create the pairwise
    % alignments so that the training set is always its native orientation and
    % each test set is aligned to the new train set?

    % Create pairwise transformation matrix
    transform_matrices = cell(num_sessions);
    if align_all_to_one_session
        % This has to be defined by the user. Typically is can be first, but may
        % vary depending on image sizes and such.
        align_to_session_num = select_session_from_selection(session_run_list);

        % Extract base data to be aligned to
        session_ind = unique_sessions == align_to_session_num;
        base_data = iDop_cell{session_ind};
        base_angiogram = makeAngiogram(base_data.iDop);


        for session = setdiff(unique_sessions, align_to_session_num)

            current_session_ind = unique_sessions == session;
            current_session_data = iDop_cell{current_session_ind};

            % Align test to specified session
            [current_session_data.iDop, transform_matrices{session_ind, current_session_ind}] = align_iDop(base_data.iDop, current_session_data.iDop);

            % Overwriting to save computer memory
            iDop_cell{current_session_ind} = current_session_data;
        end
    else
        counter = 1;
        for session1 = 1:num_sessions
            for session2 = setdiff(1:num_sessions, session1)

                session1_data = iDop_cell{session1};
                session2_data = iDop_cell{session2};

                % Find the current alignment between the sessions but don't
                % apply until later.
                [~, transform_matrices{session1, session2}] = align_iDop(session1_data.iDop, session2_data.iDop);

                fprintf('Finished %d/%d\n', counter, num_sessions^2-num_sessions);
                counter = counter+1;
            end
        end
    end

else
    base_filename = '*pairwise_session_decoding_performance_*.mat';
    [filename, pathname] = uigetfile(fullfile(get_user_data_path('pathType','decoding'),base_filename), 'Load across session performance');

    % save the results file
    fprintf('Loading %s...', fullfile(pathname, filename));
    load(fullfile(pathname, filename), 'transform_matrices');
end

%% Preprocess each session
for session = 1:num_sessions
    current_session_data = iDop_cell{session};

    if detrend
        dop_detrend = detrend_sliding_window(current_session_data.iDop, 50);
    else
        dop_detrend = current_session_data.iDop;
    end

    % Overwriting to save computer memory
    current_session_data.iDopP = preProcess(dop_detrend, 'timeGain', timeGain, 'spatialFilter', {'disk', diskFilter}, 'zScore', zScore);
    iDop_cell{session} = current_session_data;

    figure;
    tld = tiledlayout(1, 2);
    nexttile()
    imagesc(makeAngiogram(current_session_data.iDop));
    axis image;
    nexttile()
    imagesc(makeAngiogram(current_session_data.iDopP));
    axis image;
    title(tld, sprintf('S%d', unique_sessions(session)));

end

%% Load null mean angular error distribution
num_repetitions = 100000;
nClasses = 8;
try
    load(sprintf('angular_error_null_distribution_%dtgts_%dreps.mat', nClasses, num_repetitions), 'null_mean_angular_error');
    null_mean_angular_error_complete = null_mean_angular_error;
catch
    choice = questdlg(sprintf('This combo of %d Classes and %d num_repetitions did not exist. Do you want to create now?', nClasses, num_repetitions),...
        'Create new null distribution?','yes','no','yes');
    if strcmp(choice,'yes')
        [~, null_mean_angular_error] = generate_null_distribution_angular_error(nClasses, num_repetitions);
    else
        error('The distribution does not exist currently. Aborting.\n');
    end
    null_mean_angular_error_complete = null_mean_angular_error;
end

%% Load GLM mask (if desired)
% For directional tuning paper, add ability to load GLM mask from a single
% session and use that across all sessions. This is a type of
% data leak for the self session, but not the other sessions.

user_choice = questdlg('Choose one of the following options', 'use GLM mask?',...
    'GLM mask', 'Use entire image', 'Use entire image');

user_choice = string(user_choice);
if strcmp(user_choice, 'GLM mask')
    [glm_session, glm_run] = select_session_from_selection(session_run_list);
    ROI_mask = load_glm_mask([glm_session glm_run], 1e-5);
    glm_session_ind = find(unique_sessions == glm_session);
else
    ROI_mask = NaN(1, 1);
end

%% Perform decoding on each session by itself
% We can technically skip this code block and just load previously analyzed
% data, but I want to run it here to ensure that the same analysis steps
% are being done to the data. I.e. don't want a difference in the analysis
% steps to explain why we get worse performance on decoding other sessions.

[across_session_performance] = cell(num_sessions);
[best_accuracy, best_angular_error, pvalues_accuracy, pvalues_angular_error] = deal(NaN(num_sessions));

for session = 1:num_sessions
    current_session_data = iDop_cell{session};

    trainSession = current_session_data.sessionRunList(1, 1);
    trainRun = current_session_data.sessionRunList(1, 2);

    iDop = current_session_data.iDopP;
    behavior = current_session_data.behavior;
    [yPix, xPix, nWindows, nTrials] = size(iDop);

    iDop_3D = reshape(iDop, [yPix*xPix, nWindows, nTrials]);

    % implement ROIs, masks, cropping
    % Warp the ROI_mask to the current session
    if strcmp(user_choice, 'GLM mask')
        iDop_train_angiogram = makeAngiogram(iDop);
        train_session_ind = session;
        if train_session_ind ~= glm_session_ind
            ROI_mask_warp =  imwarp(...
                ROI_mask, ...
                transform_matrices{train_session_ind, glm_session_ind}, ...
                'OutputView', imref2d(size(iDop_train_angiogram)));
        else
            ROI_mask_warp = ROI_mask;
        end

        % first, flatten the ROI mask
        ROI_mask_flat = reshape(ROI_mask_warp, yPix*xPix, 1);

        % Apply ROI mask
        iDop_3D = iDop_3D(ROI_mask_flat, :, :);
    end

    % Assigning labels to multi-targets
    [fix, mem, mov] = getEpochs(behavior);
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

    label_strings = {'Negative','Center','Positive'};

    train_labels = NaN(size(targetPos,1),2);
    for dimension = 1:2
        train_labels(targetPos(:, dimension) < 0, dimension) = 1;
        train_labels(targetPos(:, dimension) == 0, dimension) = 2;
        train_labels(targetPos(:, dimension) > 0, dimension) = 3;
    end

    % Determine dimensionality of cPCA subspace (if used)
    % For now, using (number of classes - 1), i.e m = 2 for 3 classes
    m = length(unique(train_labels)) - 1;

    % allocate memory
    [result_horz, confusion_horz, cp_horz] = deal(cell(1,nWindows));   % one cell for each data point in memory epoch
    [p_horz, percentCorrect_horz, t_horz] = deal(NaN(1,nWindows));     % one array point for each point in mem epoch

    [result_vert, confusion_vert, cp_vert] = deal(cell(1,nWindows));   % one cell for each data point in memory epoch
    [p_vert, percentCorrect_vert, t_vert] = deal(NaN(1,nWindows));     %

    [result_combined, confusion_combined, cp_combined] = deal(cell(1,nWindows));   % one cell for each data point in memory epoch
    [p_combined, percentCorrect_combined, empiric_mean_angular_error, pvalue_empiric] = deal(NaN(1,nWindows));

    fprintf('\n ||| ||| ||| --- TESTING SESSION %i RUN %i --- ||| ||| ||| \n\n',trainSession,trainRun)

    % set up timing for train/test sets
    coreParams = current_session_data.coreParams;
    fixTime = ceil(abs(fix(1))*coreParams.framerate);
    memTime = ceil(abs(mem(end))*coreParams.framerate);

    h = waitbar(0,'Please wait ... testing ability to decode');

    variance_to_keep = 95;

    % To speed up the analysis, ignore times before the cue onset.
    start_timepoint = fixTime;

    for timepoint = start_timepoint:nWindows
        % this is the meat of this whole loop where we continue to pile on data
        % (concatenate it, really) as we progress through time. This is on the
        % underlying assumption that the classifier is going to use data from
        % whatever epoch of the task we're currently in to classify data from
        % that same epoch.
        % use all of the data during the fixation period...
        if timepoint <= fixTime
            testTimeI = 1:timepoint;
            % or use all the data during the memory delay and movement
        else
            testTimeI = fixTime:timepoint;
        end

        % prepare the training data for this time window
        trainData = flattenDoppler2D(iDop_3D, testTimeI);  % nImages x (nPixels*yPixels)

        % cross validated performance is calculated here
        [cp_horz{timepoint}, p_horz(timepoint), cp_vert{timepoint}, p_vert(timepoint), cp_combined{timepoint}, p_combined(timepoint), custom_counting_matrix] = crossValidate_multicoder(trainData, train_labels, ...
            'classificationMethod', classifierString,...
            'validationMethod', validationString,...
            'K', K, ...
            'm', m,...
            'variance_to_keep', variance_to_keep);

        % grabbing useful values
        percentCorrect_horz(timepoint) = cp_horz{timepoint}.CorrectRate*100;
        confusion_horz{timepoint} = cp_horz{timepoint}.CountingMatrix;

        percentCorrect_vert(timepoint) = cp_vert{timepoint}.CorrectRate*100;
        confusion_vert{timepoint} = cp_vert{timepoint}.CountingMatrix;

        percentCorrect_combined(timepoint) = cp_combined{timepoint}.CorrectRate*100;
        confusion_combined{timepoint} = custom_counting_matrix;

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
        [empiric_mean_angular_error(timepoint), ~, pvalue_empiric(timepoint)] = calculate_angular_error(predicted_class, actual_class, ...
            'ignore_middle_points', false, ...
            'num_repetitions', num_repetitions, ...
            'null_mean_angular_error', null_mean_angular_error_complete(:, 1:length(train_labels)));
        % Calculate mean_angular_error for each timepoint


        % save the data to results
        result_combined{timepoint}.percentCorrect = percentCorrect_combined(timepoint);
        result_combined{timepoint}.confusion = confusion_combined{timepoint};
        result_combined{timepoint}.p = p_combined(timepoint);
        result_combined{timepoint}.t = timepoint-1;
        result_combined{timepoint}.angular_error = empiric_mean_angular_error;
        result_combined{timepoint}.angular_error_pvalue = pvalue_empiric;

        t_horz(timepoint) = timepoint-1;
        waitbar(timepoint / nWindows)
    end
    close(h);

    % Extract performance
    [bestPercentCorrect, bestIndex] = max(percentCorrect_combined);
    [bestAngularError, bestIndex_angularError] = min(empiric_mean_angular_error);
    bestConfusion_combined = confusion_combined{bestIndex};
    bestP_combined = p_combined(bestIndex);
    bestP_angular_error = pvalue_empiric(bestIndex_angularError);

    across_session_performance{session, session} = result_combined;
    best_accuracy(session, session) = bestPercentCorrect;
    best_angular_error(session, session) = bestAngularError;


    pvalues_accuracy(session, session) = bestP_combined;
    pvalues_angular_error(session, session) = bestP_angular_error;

end

%% Define some things for the confusion matrices
combined_label_strings = {'Ipsilateral Down', 'Center Down', 'Contralateral Down', ...
    'Ipsilateral Center',  'Middle', 'Contralateral Center', ...
    'Ipsilateral Up', 'Center Up', 'Contralateral Up'};
reordered_label_strings = combined_label_strings([7 8 9 6 3 2 1 4 5]);


%% Now do pairwise comparison
debug = false;

for train_session = 1:num_sessions
    fprintf('Training session - %d\n', unique_sessions(train_session));
    train_session_data = iDop_cell{train_session};

    % Extract some things
    iDopP_train_orig = train_session_data.iDopP;
    behavior_train_orig = train_session_data.behavior;
    coreParams_train_orig = train_session_data.coreParams;


    iDop_train_angiogram = makeAngiogram(iDopP_train_orig);

    % Compile some things
    [fix_train, mem_train, mov_train] = getEpochs(behavior_train_orig);

    for test_session = setdiff(1:num_sessions, train_session)
        fprintf(' -- Testing session - %d\n', unique_sessions(test_session));

        % Redefine these in case they got changed.
        iDopP_train = iDopP_train_orig;
        coreParams_train = coreParams_train_orig;
        behavior_train = behavior_train_orig;

        test_session_data = iDop_cell{test_session};
        iDopP_test = test_session_data.iDopP;
        behavior_test = test_session_data.behavior;
        coreParams_test = test_session_data.coreParams;
        test_session_run_list = test_session_data.sessionRunList;

        if ~align_all_to_one_session
            % Perform alignment here since it has not been done yet
            iDopP_test =  imwarp(...
                iDopP_test, ...
                transform_matrices{train_session, test_session}, ...
                'OutputView', imref2d(size(iDop_train_angiogram)));
        end

        [fix_test, mem_test, mov_test] = getEpochs(behavior_test);

        % cut to same timepoints
        if round(fix_train(1)) ~= round(fix_test(1))
            min_time_ind = min(abs(round(fix_train(1))), abs(round(fix_test(1))));
            iDopP_train = iDopP_train(:, :, -round(fix_train(1)*coreParams_train.framerate)-(min_time_ind*coreParams_train.framerate)+1:end, :);
            iDopP_test = iDopP_test(:, :, -round(fix_test(1)*coreParams_test.framerate)-(min_time_ind*coreParams_test.framerate)+1:end, :);
        end


        % If different framerates
        if coreParams_train.framerate ~= coreParams_test.framerate
            % Convert to 1 Hz.
            if coreParams_train.framerate > 1
                iDopP_train = iDopP_train(:, :, 1:2:end, :);
                coreParams_train.framerate = 1;
            elseif coreParams_test.framerate > 1
                iDopP_test = iDopP_test(:, :, 1:2:end, :);
                coreParams_test.framerate = 1;
            end
        end

        % If still different lengths, then cut to shortest
        if size(iDopP_train, 3) ~= size(iDopP_test, 3)
            min_time_length = min(size(iDopP_train, 3), size(iDopP_test, 3));
            iDopP_train = iDopP_train(:, :, 1:min_time_length, :);
            iDopP_test = iDopP_test(:, :, 1:min_time_length, :);
        end

        [yPix_test, xPix_test, nWindows_test, nTrials_test] = size(iDopP_test);
        [yPix_train, xPix_train, nWindows_train, nTrials_train] = size(iDopP_train);


        %% Flatten Doppler for easier masking
        % nPixels x nWindows x nTrials
        % This allows us to eliminate specific voxels that we don't want. Supports
        % non-rectangular ROIs.
        iDop_3D_train = reshape(iDopP_train, [yPix_train*xPix_train, nWindows_train, nTrials_train]);
        iDop_3D_test = reshape(iDopP_test, [yPix_test*xPix_test, nWindows_test, nTrials_test]);

        %% Apply ROI mask (if desired)
        % Warp the ROI_mask to the current session
        if strcmp(user_choice, 'GLM mask')

            if train_session ~= glm_session_ind
                ROI_mask_warp =  imwarp(...
                    ROI_mask, ...
                    transform_matrices{train_session, glm_session_ind}, ...
                    'OutputView', imref2d(size(iDop_train_angiogram)));
            else
                ROI_mask_warp = ROI_mask;
            end

            % first, flatten the ROI mask
            ROI_mask_flat = reshape(ROI_mask_warp, yPix_train*xPix_train, 1);

            % Apply ROI mask
            iDop_3D_train = iDop_3D_train(ROI_mask_flat, :, :);
            iDop_3D_test = iDop_3D_test(ROI_mask_flat, :, :);
        end

        %% Assigning labels to multi-targets
        [targetPosInd_train, train_labels, label_strings_train] = extract_direction_labels(behavior_train);
        [targetPosInd_test, test_labels, label_strings_test] = extract_direction_labels(behavior_test);

        %% run cross validation for each sample in the memory epoch
        % Determine dimensionality of cPCA subspace (if used)
        % For now, using (number of classes - 1), i.e m = 2 for 3 classes
        m = length(unique(train_labels)) - 1;

        % allocate memory
        [result_horz, confusion_horz, cp_horz] = deal(cell(1,nWindows_train));   % one cell for each data point in memory epoch
        [p_horz, percentCorrect_horz, t_horz] = deal(NaN(1,nWindows_train));     % one array point for each point in mem epoch

        [result_vert, confusion_vert, cp_vert] = deal(cell(1,nWindows_train));   % one cell for each data point in memory epoch
        [p_vert, percentCorrect_vert, t_vert] = deal(NaN(1,nWindows_train));     %

        [result_combined, confusion_combined, cp_combined] = deal(cell(1,nWindows_train));   % one cell for each data point in memory epoch
        [p_combined, percentCorrect_combined, t_combined] = deal(NaN(1,nWindows_train));     %

        % set up timing for train/test sets
        fixTime = ceil(abs(fix_train(1))*coreParams_train.framerate);

        test_labels_combined = test_labels(:,1) + length(unique(test_labels(:,1))) * (test_labels(:,2)-1);

        num_combined_labels = length(unique(test_labels(:, 1))) * length(unique(test_labels(:, 2)));

        % To speed up the analysis, ignore times before the cue onset.
        start_timepoint = fixTime;

        h = waitbar(0,'Please wait ... testing ability to decode');

        for timepoint = start_timepoint:nWindows_train
            % this is the meat of this whole loop where we continue to pile on data
            % (concatenate it, really) as we progress through time. This is on the
            % underlying assumption that the classifier is going to use data from
            % whatever epoch of the task we're currently in to classify data from
            % that same epoch.

            % use all of the data during the fixation period...
            if timepoint <= fixTime
                testTimeI = 1:timepoint;
                % or use all the data during the memory delay and movement
            else
                testTimeI = fixTime:timepoint;
            end

            % Initialize empty variable
            countingmatrix_custom = zeros(num_combined_labels);

            % prepare the training data for this time window
            trainData = flattenDoppler2D(iDop_3D_train, testTimeI);  % nImages x (nPixels*yPixels)
            testData = flattenDoppler2D(iDop_3D_test, testTimeI);

            trainData_norm = zscore(trainData);
            testData_norm = zscore(testData);

            class_horz = classifyDoppler(trainData_norm, ...
                train_labels(:,1), ...
                testData_norm, ...
                'method', classifierString, ...
                'm', m,...
                'variance_to_keep', 95);
            class_vert = classifyDoppler(trainData_norm, ...
                train_labels(:,2) , ...
                testData_norm, ...
                'method', classifierString, ...
                'm', m,...
                'variance_to_keep', 95);

            if all(~isnan(test_labels_combined))
                class_combined = class_horz + length(unique(test_labels(:,1))) * (class_vert-1);
            end
            cp_horz{timepoint} = classperf(test_labels(:, 1), class_horz);
            cp_vert{timepoint} = classperf(test_labels(:, 2), class_vert);
            if all(~isnan(test_labels_combined))
                cp_combined{timepoint} = classperf(test_labels_combined, class_combined);
                countingmatrix_custom = update_counting_matrix(countingmatrix_custom, class_combined, test_labels_combined);
            end

            % Angular error calculation
            % Reconstruct a record for final performance. This
            % is not the actual order of predictions since that
            % is not easily accessible. This is correct number
            % of each prediction though.
            [actual_class, predicted_class] = deal([]);
            for true_class_num = 1:9
                for predicted_class_num = 1:9
                    actual_class = [actual_class; true_class_num*ones(countingmatrix_custom(predicted_class_num,true_class_num), 1)];
                    predicted_class = [predicted_class; predicted_class_num*ones(countingmatrix_custom(predicted_class_num, true_class_num), 1)];
                end
            end
            [empiric_mean_angular_error, ~, pvalue_empiric] = calculate_angular_error(predicted_class, actual_class, ...
                'ignore_middle_points', false, ...
                'num_repetitions', num_repetitions, ...
                'null_mean_angular_error', null_mean_angular_error_complete(:, 1:length(test_labels)));


            % grabbing useful values
            percentCorrect_horz(timepoint) = cp_horz{timepoint}.CorrectRate*100;
            nCorrect_horz = sum(diag(cp_horz{timepoint}.CountingMatrix));
            nCounted_horz = sum(cp_horz{timepoint}.CountingMatrix(:));
            chance_horz = 1/length(unique(class_horz));
            p_horz(timepoint) = binomialTest(nCorrect_horz, nCounted_horz, chance_horz, 'one');
            confusion_horz{timepoint} = cp_horz{timepoint}.CountingMatrix;

            percentCorrect_vert(timepoint) = cp_vert{timepoint}.CorrectRate*100;
            nCorrect_vert = sum(diag(cp_vert{timepoint}.CountingMatrix));
            nCounted_vert = sum(cp_vert{timepoint}.CountingMatrix(:));
            chance_vert = 1/length(unique(class_vert));
            p_vert(timepoint) = binomialTest(nCorrect_vert, nCounted_vert, chance_vert, 'one');
            confusion_vert{timepoint} = cp_vert{timepoint}.CountingMatrix;

            percentCorrect_combined(timepoint) = cp_combined{timepoint}.CorrectRate*100;
            nCorrect_combined = sum(diag(cp_combined{timepoint}.CountingMatrix));
            nCounted_combined = sum(cp_combined{timepoint}.CountingMatrix(:));
            chance_combined = 1/length(unique(class_combined));
            p_combined(timepoint) = binomialTest(nCorrect_combined, nCounted_combined, chance_combined, 'one');
            confusion_combined{timepoint} = countingmatrix_custom;

            if debug
                figure;

                reordered_confusion_mat = countingmatrix_custom([7 8 9 6 3 2 1 4 5], [7 8 9 6 3 2 1 4 5]);

                imagesc_confusion_matrix(reordered_confusion_mat(:, 1:8)', ...
                    'label_strings', reordered_label_strings, ...
                    'title', '', ...
                    'FontColor', 'k');
                title(sprintf('train S%d\ntest S%d\nTimepoint %d', unique_sessions(train_session), unique_sessions(test_session), timepoint));
            end
            % save the data to results
            result_combined{timepoint}.percentCorrect = percentCorrect_combined(timepoint);
            result_combined{timepoint}.confusion = confusion_combined{timepoint};
            result_combined{timepoint}.p = p_combined(timepoint);
            result_combined{timepoint}.t = timepoint-1;
            result_combined{timepoint}.angular_error = empiric_mean_angular_error;
            result_combined{timepoint}.angular_error_pvalue = pvalue_empiric;

            t_horz(timepoint) = timepoint-1;

            waitbar(timepoint / nWindows_train);
        end

        try
            close(h);
        catch
            warning('Waitbar was closed inadvertently');
        end

        %% Concatenate performance information for each session

        [bestPercentCorrect, bestIndex] = max(percentCorrect_combined);
        bestConfusion_combined = confusion_combined{bestIndex};
        [bestAngularError, bestIndex_angular_error] = min(empiric_mean_angular_error);
        bestP_combined = p_combined(bestIndex);
        bestP_angular_error = pvalue_empiric(bestIndex_angular_error);

        across_session_performance{train_session, test_session} = result_combined;
        best_accuracy(train_session, test_session) = bestPercentCorrect;
        best_angular_error(train_session, test_session) = bestAngularError;


        pvalues_accuracy(train_session,test_session) = bestP_combined;
        pvalues_angular_error(train_session, test_session) = bestP_angular_error;
    end

end

%% Save across session performance
detrend_string = '';
if detrend
    detrend_string = '_withDetrend';
else
    detrend_string = '';
end


if align_all_to_one_session
    alignment_string = sprintf('_align_all_to_S%d', align_to_session_num);
else
    alignment_string = '_pairwise_alignment';
end

analysis_parameters.K = K;
analysis_parameters.timeGain = timeGain;
analysis_parameters.diskFilter = diskFilter;
analysis_parameters.zScore = zScore;
analysis_parameters.detrend = detrend;

defaultName = sprintf('%s_pairwise_session_decoding_performance_%s%s%s_dategen_%s.mat', ...
    Monkey, classifierString, detrend_string, alignment_string,...
    datetime('today', format='yyyyMMdd'));
[filename, pathname] = uiputfile('*.mat', 'Save pairwise performance as',fullfile(get_user_data_path('pathType','across session analyses'),defaultName));


% save the results file
fprintf('Saving %s to %s...',filename,pwd)
if ~exist('align_to_session_num', 'var')
    align_to_session_num = NaN;
end

if ~exist('glm_session', 'var')
    glm_session = NaN;
end

if ~exist('ROI_mask', 'var')
    ROI_mask = NaN;
end
save(fullfile(pathname, filename), 'session_run_list', ...
    'across_session_performance', ...
    'best_accuracy', 'best_angular_error', 'analysis_parameters', ...
    'pvalues_accuracy', 'pvalues_angular_error', ...
    'transform_matrices', 'align_to_session_num', ...
    'align_all_to_one_session', 'Monkey', 'user_choice', 'glm_session', ...
    'ROI_mask');
fprintf('done.\n')


%% Helper function
function [session, run] = select_session_from_selection(session_run_list)
% Helper function to choose one session/run from a previously specified
% session/run list

run_strings = cell(size(session_run_list,1),1);
for i = 1:size(run_strings,1)
    Session = session_run_list(i, 1);
    Run = session_run_list(i, 2);

    run_strings{i} = sprintf('S%dR%d', Session, Run);
end

% Create the dialog box
[indx,tf] = listdlg('PromptString','Select Session/Run(s) to load:',...
    'SelectionMode','single',...
    'ListString',run_strings, ...
    'ListSize', [600 300]);

session = session_run_list(indx, 1);
run = session_run_list(indx, 2);
end
