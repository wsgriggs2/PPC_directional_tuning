%% Script to create a file for `plot_across_session_offline_decoding_results`
% Previously had to run `BCI_CV_direction_multicoder_acrossSessions` but
% this was inefficient when just a single file was changed. This script now
% extracts the relevant info to allow just single sessions/runs to be
% updated as needed.
% Whitney Griggs January 25, 2025

%% Specify which files to load
project_record_filename = 'ProjectRecord_paper.json';
ProjectRecord = load_json_as_table(project_record_filename);
session_run_list = specify_sessions_of_interest(...
    'project_record_filename', project_record_filename);

n_files = size(session_run_list, 1);

bad_sessions = [];
detrend = true;

use_motionCorrection = true;
%% Now iterate and extract relevant content from each file
across_session_results = cell(1, n_files);
task_timing = nan(n_files, 4);
[nWindows_session, nTrials_session] = deal(nan(n_files, 1));

home_dir = pwd;

for filenum = 1:n_files
    try
    reduced_session_run_list = session_run_list(filenum, :);

    current_session = session_run_list(filenum, 1);
    current_run = session_run_list(filenum, 2);

    Monkey = ProjectRecord{(ProjectRecord.Session==current_session & ProjectRecord.Run==current_run), 'Monkey'};

    if sum(session_run_list(:, 1) == current_session) > 1
        % Currently allows multiple session to be loaded multiple times.
        % Fix & this to prevent mistakes in the future. Currently manually
        % making & sure ot only load the combined session/run file, but
        % this is prone to human error.
        warning ('Loading multiple runs from session %d, but separately', current_session);
    end

    % Extract task timing info
    [iDop, ~, ~, behavior, ~, ~, ~] = ...
        concatenateRuns(reduced_session_run_list, ...
        use_motionCorrection, ...
        'manual_alignment', true, ...
        'project_record_filename', project_record_filename, ...
        'load_all_attempted_trials', false);
    
    % Get task timing for each session
    [fix, memm, mov] = getEpochs(behavior);
    task_timing(filenum,:) = [fix(1) fix(2) memm(2) mov(2)];
    [~, ~, nWindows_session(filenum), nTrials_session(filenum)] = size(iDop);

    partial_filename = sprintf('data_S%dR%d*.mat', current_session, current_run) ;

    cd (fullfile(get_user_data_path('pathType', 'decoding'))) ;

    % Find most recent matching file
    file_listing = dir(partial_filename);
    [~,most_recent_file_index] = max([file_listing.datenum]);
    most_recent_file = file_listing (most_recent_file_index) .name;

    % Load the data
    decoding_results = load (most_recent_file);
    result_horz = decoding_results.result_horz;
    result_vert = decoding_results.result_vert;
    result_combined = decoding_results.result_combined;
    cp_horz = decoding_results.cp_horz;
    cp_vert = decoding_results.cp_vert;
    cp_combined = decoding_results.cp_combined;

    % Fix a few variables
    % Old version of files saved out too many timepoints for angular_error
    % and associated pvalue. Check and fix the saving of each variable if necessary
    for i = 1:length(result_combined)
        angular_error = result_combined{i}.angular_error;
        angular_error_pvalue = result_combined{i}.angular_error_pvalue;
        
        if numel(angular_error)>1
            result_combined{i}.angular_error = angular_error(i);
            result_combined{i}.angular_error_pvalue = angular_error_pvalue(i);
        end
    end
    

    % Extract relevant contents
    % No need to save out cp_horz, cp_vert, cp_combined, but doing for
    % consistency with the `BCI_CV_direction_multicoder_acrossSessions`
    % script. These 3 variables are not used by downstream functions in
    % other scripts to the best of my knowledge
    % Whitney Griggs, January 25, 2025
    across_session_results{filenum} = {result_horz, result_vert, result_combined, reduced_session_run_list, Monkey, cp_horz, cp_vert, cp_combined};
    catch

        warning('Unable to process S%dR%d', current_session, current_run);
    end
end

cd(home_dir);

%% Save file

root_savepath = get_user_data_path('pathType', 'output');
filename = sprintf('AcrossSessionResults_8dir_allMonkeys_dateGenerated_%s.mat', datetime('today', 'Format', 'yyyyMMdd'));
save(fullfile(root_savepath, 'decoding', filename), 'across_session_results', 'bad_sessions', 'task_timing', 'nWindows_session', 'nTrials_session', 'detrend', 'use_motionCorrection');


