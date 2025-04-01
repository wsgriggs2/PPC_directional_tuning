%% Script to run decoding analysis across all specified sessions independently

clear
close all

startDir = pwd;
project_record_filename = 'ProjectRecord_paper.json';
ProjectRecord = load_json_as_table(project_record_filename);

sessions_to_use = [ProjectRecord{1:end, 'Session'} ProjectRecord{1:end, 'Run'}];

bad_sessions = [];

% Retrieve the number of unique sessions
[unique_sessions, ia] = unique(sessions_to_use(:, 1));
num_of_sessions = length(unique_sessions);

%% hard coded settings
% choose K for K-fold validation
K = 10;
% preprocessing
timeGain = false;
diskFilter = 2;
zScore = false;
detrend = true;

use_motionCorrection = true;

%% get classifier string
classifierString = 'PCA+LDA';
validationString = 'kFold';

%% Loop over sessions

% preallocation
across_session_results = cell(1, num_of_sessions);
task_timing = nan(num_of_sessions, 4);
[nWindows_session, nTrials_session] = deal(nan(num_of_sessions, 1));

% For bootstrapping
num_repetitions = 100000;
nClasses = 8;
% Get null distribution for angular error
try
    load(sprintf('angular_error_null_distribution_%dtgts_%dreps.mat', nClasses, num_repetitions));
    null_mean_angular_error_full = null_mean_angular_error;
catch
    choice = questdlg(sprintf('This combo of %d Classes and %d num_repetitions did not exist. Do you want to create now?', nClasses, num_repetitions),...
        'Create new null distribution?','yes','no','yes');
    if strcmp(choice,'yes')
        [~, null_mean_angular_error_full] = generate_null_distribution_angular_error(nClasses, num_repetitions);
    else
        error('The distribution does not exist currently. Aborting.\n');
    end
end


% Iterate over each session
for session_num = 1:num_of_sessions
    % Extract current session number
    Session = unique_sessions(session_num);
    Run_1 = sessions_to_use(ia(session_num), 2);
    
    % Find all runs of that session
    session_ind = sessions_to_use(:, 1) == Session;
    session_run_list = sessions_to_use(session_ind, :);
    
    % Find what monkey this session is from
    ProjectRecord_indx = intersect(find(ProjectRecord.Session==Session),...
        find(ProjectRecord.Run==Run_1));
    Monkey = ProjectRecord{ProjectRecord_indx, 'Monkey'};
    
    % Load data from a single session
    [iDop, behavior, coreParams, angiogram, UF] = ...
        concatenateRuns(session_run_list, ...
        use_motionCorrection, ...
        'manual_alignment', true, ...
        'project_record_filename', project_record_filename, ...
        'load_all_attempted_trials', false);
    
    % Apply any preprocessing necessary
    if detrend
        iDopP = detrend_sliding_window(iDop, 50);
    else
        iDopP = iDop;
    end
    iDopP = preProcess(iDopP, 'timeGain', timeGain, 'spatialFilter', {'disk', diskFilter}, 'zScore', zScore);
    
    % Get task timing for each session
    [fix, memm, mov] = getEpochs(behavior);
    task_timing(session_num,:) = [fix(1) fix(2) memm(2) mov(2)];
    
    
    [yPix, xPix, nWindows_session(session_num), nTrials_session(session_num)] = size(iDopP);
    
    %% create a clear mask if we need to (necessary for parfor)
    if ~exist('flatMask','var')  || length(flatMask)~=(xPix*yPix)
        flatMask = ones(1,yPix*xPix);
    end
    
    
    %% Assigning labels to multi-targets
    show_figure = false;

    target = cell2mat({behavior.target});
    targetPos = vertcat(behavior.targetPos);
    if iscell(targetPos)
        targetPos = cell2mat(targetPos);
    end
    [angle,distance] = cart2pol(targetPos(:,1),targetPos(:,2)); % convert to polar coordinates
    angle(angle<0) = angle(angle<0)+2*pi; %Convert to have only positive angles.
    
    UniqueAngles = unique(angle);
    TargetPosInd = zeros(length(target),1);
    for position = 1:length(UniqueAngles)
        TargetPosInd(angle==UniqueAngles(position)) = position;
    end
   
    
    % create the training data labels (2 x n_trials) - x, y columns
    % For horizontal and vertical axes, use the following labels
    % 1 = Negative (down or left)
    % 2 = At center (along vertical or horizontal axis)
    % 3 = Positive (up or right)

        label_strings = {'Negative','Center','Positive'};

    if show_figure
    figure;
    tiledlayout(1,2);
    end
    
    
    train_labels = NaN(size(targetPos,1),2);
    
    for dimension = 1:2
        train_labels(targetPos(:, dimension) < 0, dimension) = 1;
            train_labels(targetPos(:, dimension) == 0, dimension) = 2;
            train_labels(targetPos(:, dimension) > 0, dimension) = 3;

        if show_figure
        nexttile;
        histogram(train_labels(:,dimension));
        title(sprintf('Number of trials for each training label: dimension %d',dimension));
        end
    end
    
    
    % Get null distribution for angular error
    if nClasses ~= 8
        % Throw an error if not 8 classes. We loaded a specific null
        % distribution. If this is not always true, then revise code to
        % load the bootstrap distribution for each new session.
        error('This script assumes there are 8 classes. Find the bug and fix.');
    end
    null_mean_angular_error = null_mean_angular_error_full(:, 1:length(train_labels));
    
    
    try
        %% run cross validation for each sample in the memory epoch
        
        % Determine dimensionality of cPCA subspace (if used)
        % For now, using (number of classes - 1), i.e m = 2 for 3 classes
        m = length(unique(train_labels)) - 1;
        
        % allocate memory
        [result_horz, confusion_horz, cp_horz] = deal(cell(1,nWindows_session(session_num)));   % one cell for each data point in memory epoch
        [p_horz, percentCorrect_horz, t_horz] = deal(NaN(1,nWindows_session(session_num)));     % one array point for each point in mem epoch
        
        [result_vert, confusion_vert, cp_vert] = deal(cell(1,nWindows_session(session_num)));   % one cell for each data point in memory epoch
        [p_vert, percentCorrect_vert, t_vert] = deal(NaN(1,nWindows_session(session_num)));     %
        
        [result_combined, confusion_combined, cp_combined] = deal(cell(1,nWindows_session(session_num)));   % one cell for each data point in memory epoch
        [p_combined, percentCorrect_combined, t_combined] = deal(NaN(1,nWindows_session(session_num)));     %
        
        % Create string for number of runs being used
        run_string = [];
        for i = 1:size(session_run_list, 1)
            if i == size(session_run_list, 1)
                run_string = [run_string sprintf('%d', session_run_list(i, 2))];
            else
                run_string = [run_string sprintf('%d ', session_run_list(i, 2))];
            end
        end
        fprintf('\n ||| ||| ||| --- TESTING SESSION %i RUN [%s] --- ||| ||| ||| \n\n',Session, run_string)
        
        % set up timing for train/test sets
        fixTime = ceil(abs(fix(1))*coreParams.framerate);
        
        thePool = gcp;
        H = ProgressBar(nWindows_session(session_num), taskname='Running BCI (cross-validated): ');
        for i = 1:nWindows_session(session_num)
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
            trainData = flattenDoppler2D(iDopP, testTimeI);  % nImages x (nPixels*yPixels)
            
            % cross validated performance is calculated here
            [cp_horz{i}, p_horz, cp_vert{i}, p_vert, cp_combined{i}, p_combined, custom_counting_matrix] = crossValidate_multicoder(trainData, train_labels, ...
                'classificationMethod', classifierString, ...
                'validationMethod', validationString, 'K', K, ...
                'm', m,...
                'variance_to_keep', 95);
            
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
            [empiric_mean_angular_error(i), ~, pvalue_empiric(i)] = calculate_angular_error(predicted_class, actual_class, ...
                'ignore_middle_points', false, ...
                'num_repetitions', num_repetitions, ...
                'null_mean_angular_error', null_mean_angular_error);
            
            % save the data to results
            result_horz{i}.percentCorrect = cp_horz{i}.CorrectRate*100;
            result_horz{i}.confusion =  cp_horz{i}.CountingMatrix;
            result_horz{i}.p = p_horz;
            result_horz{i}.t = i-1;
            
            result_vert{i}.percentCorrect = cp_vert{i}.CorrectRate*100;
            result_vert{i}.confusion = cp_vert{i}.CountingMatrix;
            result_vert{i}.p = p_vert;
            result_vert{i}.t = i-1;
            
            result_combined{i}.percentCorrect = cp_combined{i}.CorrectRate*100;
            result_combined{i}.confusion = custom_counting_matrix;
            result_combined{i}.p = p_combined;
            result_combined{i}.t = i-1;
            result_combined{i}.angular_error = empiric_mean_angular_error(i);
            result_combined{i}.angular_error_pvalue = pvalue_empiric(i);
            
            % Increment the parallel pool counter
            H.count(); 
            
        end
    catch
        fprintf('Session - %d, Run(s) - [%s] \n', Session, run_string)
        warning('This session broke for some reason \n');
        bad_sessions = [bad_sessions; [repmat(session_num, size(session_run_list, 1), 1) session_run_list];];
    end
    across_session_results{session_num} = {result_horz, result_vert, result_combined, session_run_list, Monkey, cp_horz, cp_vert, cp_combined};
end

%% Saving
root_savepath = get_user_data_path('pathType', 'decoding');
filename = sprintf('AcrossSessionResults_8dir_allMonkeys_dateGenerated_%s.mat', datetime('today', format='yyyyMMdd'));
save(fullfile(root_savepath, filename), 'across_session_results', 'bad_sessions', 'task_timing', 'nWindows_session', 'nTrials_session', 'detrend', 'use_motionCorrection');

