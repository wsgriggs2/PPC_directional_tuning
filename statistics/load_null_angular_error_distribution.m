function [null_angular_error, null_mean_angular_error] = load_null_angular_error_distribution(nClasses, num_repetitions)
%% Helper function to load the null distributions from file. 
% if the file does not exist, then ask if the user wants to create it.

% Written by Whitney on July 7, 2022
try
    load(sprintf('angular_error_null_distribution_%dtgts_%dreps.mat', nClasses, num_repetitions));
catch
    choice = questdlg('This combo of nClasses and num_repetitions did not exist. Do you want to create now?',...
        'Create new null distribution?','yes','no','yes');
    if strcmp(choice,'yes')
        [null_angular_error, null_mean_angular_error] = generate_null_distribution_angular_error(nClasses, num_repetitions);
    else
        error('The distribution does not exist currently. Aborting.\n');
    end
end
end