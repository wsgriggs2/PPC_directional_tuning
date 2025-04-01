function glm_mask = load_glm_mask(session_run, p_threshold)
% Helper function to load the most recent GLM data and return the
% thresholded mask

start_dir = pwd;

% Logic for filenaming
% Complete filename up until the date generated
partial_filename = sprintf('glm_data_%sgen*.mat', sprintf('S%dR%d_', session_run'));

cd(fullfile(get_user_data_path('pathType', 'glm')));

% Find most recent matching file
file_listing = dir(partial_filename);
[~,most_recent_file_index] = max([file_listing.datenum]);
most_recent_file = file_listing(most_recent_file_index).name;

% Load this most recent file
fprintf('Loading %s \n', most_recent_file);
load(most_recent_file, 'pvalues_corrected');
cd(start_dir);

% To handle older files where I did not separately store F and other useful
% variables
if ~exist('pvalues_corrected', 'var')
    load(most_recent_file, 'glm_data');
    pvalues_corrected = glm_data.p_corrected;
end

% Rename for convenience
glm_pvalues = pvalues_corrected;

% Define GLM mask
glm_mask = ~(glm_pvalues>p_threshold);
end