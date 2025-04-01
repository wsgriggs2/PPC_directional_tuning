
%% Load the Project Record and select individual sessions to analyze

sessions_runs_to_load = specify_sessions_of_interest('project_record_filename', 'ProjectRecord_paper.json');

loadDopplerDataContinuous(whos, 'multiple', true, ...
    'SessionRunList', sessions_runs_to_load,...
    'mc_method', 'normcorre', ...
    'project_record_filename', 'ProjectRecord_paper.json');


[nDepth, nWidth, nTimepoints] = size(data.dop);

dop = data.dop;
behavior = data.behavior;
timestamps = data.timestamps(:, 1);
angiogram = data.angiogram;


%% Confirm that smoothing increases spatial autocorrelation
% Define parameters for the Gaussian smoothing kernel

FWHM = 0;

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
dop = preProcess(dop, 'spatialFilter', {filter_type, filter_size, sigma}, 'zscore', false);


%% Estimate correlation for different sizes of voxel patches
radii_to_test = 0:1:25;
%Do this where I only look at voxels that are between X-1 and X away, not X
% and less away, i.e. approximately a single line of voxels as opposed to a filled circle
% of voxels. This better measures the spatial autocorrelation at differrent
% distances. Otherwise, there are lots of voxels that are closer that
% contribute to the mean correlation.
shape = 'ring';



num_radii_to_test = length(radii_to_test);

% Create place to store the images
voxel_correlations = NaN(nDepth, nWidth, num_radii_to_test);

for radius_ind = 1:num_radii_to_test
    radius = radii_to_test(radius_ind);
    
    voxel_correlations(:, :, radius_ind) = calculate_voxel_patch_spatial_correlation(dop, radius, ...
        'shape', shape);
    fprintf('Calculated radius %d\n', radius);
    
end



% save relevant variables for later possible analysis
save_filename = sprintf('spatial_autocorrelation_FWHM_%d_%s_S%dR%d_dg%s.mat', FWHM, shape, Session, Run, datetime('today', format='yyyyMMdd'));
root_savepath = get_user_data_path('pathType', 'output');
[save_file, file_path] = uiputfile(fullfile(root_savepath, 'spatial_autocorrelation', save_filename));

if save_file ~= 0
    save(fullfile(file_path, save_file), 'voxel_correlations', 'sessionRunList', ...
        'radii_to_test', 'shape', 'FWHM');
end

%% Look at each map just generated

if ~exist('num_radii_to_test', 'var')
    num_radii_to_test = length(radii_to_test);
end


figure;
tld = tiledlayout('flow');
for radius_ind = 1:num_radii_to_test
    nexttile()
    imagesc(voxel_correlations(:, :, radius_ind));%, ...

    colorbar;

end

%% Average spatial correlation across all voxels
% Calculate mean and SEM
[mean_correlation, std_correlation, sem_correlation] = deal(NaN(num_radii_to_test, 1));
for radius_ind = 1:num_radii_to_test
    current_frame = voxel_correlations(:, :, radius_ind);
    mean_correlation(radius_ind) = mean(current_frame, 'all');
    std_correlation(radius_ind) =  std(current_frame, [], 'all');
    sem_correlation(radius_ind) = std_correlation(radius_ind)/sqrt(numel(current_frame));
end

figure;
 shadedErrorBar(radii_to_test, mean_correlation, std_correlation);
%errorbar(radii_to_test, mean_correlation, std_correlation);
ylabel('Mean correlation');
xlabel('Voxel radius');
title(sprintf('S%dR%d spatial autocorrelation\nFWHM - %0.1f', sessionRunList(:, 1), sessionRunList(:, 2), FWHM));
    