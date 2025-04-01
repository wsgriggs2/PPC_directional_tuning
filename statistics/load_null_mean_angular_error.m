function null_mean_angular_error = load_null_mean_angular_error(n_targets, num_repetitions, n_trials)
%% Helper function to load null mean angular error
% Written July 22, 2022
% Author Whitney Griggs

try
    null_distribution = load(sprintf('angular_error_null_distribution_%dtgts_%dreps.mat', n_targets, num_repetitions), 'null_angular_error', 'null_mean_angular_error');
    null_mean_angular_error = null_distribution.null_mean_angular_error(:, 1:n_trials);
catch
    choice = questdlg('This combo of nClasses and num_repetitions did not exist. Do you want to create now?',...
        'Create new null distribution?','yes','no','yes');
    if strcmp(choice,'yes')
        [~, null_mean_angular_error] = generate_null_distribution_angular_error(n_targets, num_repetitions);
    else
        error('The distribution does not exist currently. Aborting.\n');
    end
    null_mean_angular_error = null_mean_angular_error(:, 1:n_trials);
end
              