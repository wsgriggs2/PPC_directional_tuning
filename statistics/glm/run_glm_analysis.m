function data = run_glm_analysis(dop, behavior, dop_timestamps, varargin)
%% GLM Analysis for 8 target directions
% Will compute F-contrast that finds voxels that respond differently to the 8
% directions.
%
% Written by Whitney Griggs


%% Variable inputs
p = inputParser;
p.addOptional('FWHM',1)
p.addOptional('highpass_freq', 1/128);
p.addOptional('verbose', false);
p.addOptional('hrf', generate_boynton_hrf(16, 'tau', 0.7, 'delta', 2, 'n', 3, 'Fs', 1000));
p.addOptional('scaling_method', 'grand_mean_scaling');

p.parse(varargin{:});
inputs = p.Results;

%% Preprocessing

% Define parameters for the Gaussian smoothing kernel
FWHM = inputs.FWHM; % In voxels; Common values in fMRI is 1-3 voxels
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
Dop_blur = preProcess(dop, 'spatialFilter', {filter_type, filter_size, sigma}, 'zscore', false);

if inputs.verbose
    figure;
    imagesc(makeAngiogram(Dop_blur));
end

% Apply high-pass filter to remove drift
[iDop_highpass, d] = highpass_doppler(Dop_blur, inputs.highpass_freq, 1);

% Scale doppler
iDop_scaled = scale_doppler(iDop_highpass, ...
    'scaling_type', inputs.scaling_method);


%% Extract behavioral things of interest
% Extract information about target position
targetPos = vertcat(behavior.targetPos);

if iscell(targetPos)
    targetPos = cell2mat(targetPos);
end

% Convert from cartesian coordinates to polar coordinates
[angle,~] = cart2pol(targetPos(:,1),targetPos(:,2)); % convert to polar coordinates
angle(angle<0) = angle(angle<0)+2*pi; %Convert to have only positive angles.

% Find unique directions seen.
UniqueAngles = unique(angle);

% create an index showing which trials are for which directions.
TargetPosInd = zeros(length(targetPos),1);
for angleInd = 1:length(UniqueAngles)
    TargetPosInd(angle==UniqueAngles(angleInd)) = angleInd;
end
TargetPosInd_onehot = bsxfun(@eq, TargetPosInd(:), 1:max(TargetPosInd));

%% Create regressor names and labels
boxcar_regressor_names = {'memory period-dir1',    'memory period-dir2', ...
    'memory period-dir3',    'memory period-dir4', ...
    'memory period-dir5',    'memory period-dir6', ...
    'memory period-dir7',    'memory period-dir8', ...
    'sacc+hold period-dir1',      'sacc+hold period-dir2',...
    'sacc+hold period-dir3',      'sacc+hold period-dir4',...
    'sacc+hold period-dir5',      'sacc+hold period-dir6',...
    'sacc+hold period-dir7',      'sacc+hold period-dir8',...
    'juice delivery',       'fixation period'};
boxcar_behavior_start_labels = {'cue',              'cue',   ...
    'cue',              'cue',...
    'cue',              'cue',...
    'cue',              'cue',...
    'target_hold',   'target_hold',...
    'target_hold',   'target_hold',...
    'target_hold',   'target_hold',...
    'target_hold',   'target_hold',...
    'reward',         'initialfixation'};
boxcar_behavior_end_labels = {'target_acquire', 'target_acquire',   ...
    'target_acquire', 'target_acquire',...
    'target_acquire', 'target_acquire',...
    'target_acquire', 'target_acquire',...
    'reward',        'reward',...
    'reward',        'reward',...
    'reward',        'reward',...
    'reward',        'reward',...
    'iti',            'cue'};

boxcar_modulation = cell(1, length(boxcar_regressor_names));

impulse_regressor_names = [];
impulse_behavior_start_labels = [];
impulse_modulation = [];


regressor_types = cell(1, size(boxcar_regressor_names, 2) + size(impulse_regressor_names, 2));
regressor_types(1:size(boxcar_regressor_names, 2)) = {'boxcar'};
regressor_types(size(boxcar_regressor_names, 2)+1:end) = {'impulse'};

regressor_names = [boxcar_regressor_names impulse_regressor_names];
regressor_start_labels = [boxcar_behavior_start_labels impulse_behavior_start_labels];
regressor_end_labels = [boxcar_behavior_end_labels cell(1, size(impulse_regressor_names, 2))];
regressor_modulation = [boxcar_modulation impulse_modulation];

%% Indicate which events are which conditions
impulse_event_ind = cell(size(impulse_regressor_names));
boxcar_event_ind = cell(size(boxcar_regressor_names));
for i = 1:8
    % Cue+Memory period
    boxcar_event_ind{i} = TargetPosInd_onehot(:, i);
    
    % Sacc+Target hold period
    boxcar_event_ind{i+8} = TargetPosInd_onehot(:, i);
end

regressor_event_ind = [boxcar_event_ind impulse_event_ind];

%% create table with specified regressors
regressor_table = table(regressor_types', regressor_names', ...
    regressor_start_labels', regressor_end_labels', regressor_modulation', ...
    regressor_event_ind',...
    'VariableNames', {'type', 'name', 'start_label', 'end_label', ...
    'modulation', 'event_ind'});


%% Generate the regressors
% Generate the non-convolved regressors
regressor_info = get_glm_regressor_info(behavior, regressor_table);

% Convolve regressors with hrf
convolved_regressors = create_glm_regressors(dop_timestamps, regressor_info, inputs.hrf.shape, ...
    'verbose', inputs.verbose);

% Add constant regressor
convolved_regressors = [convolved_regressors ones(size(convolved_regressors, 1), 1)];


%% F-test for difference in activation
% Tests whether any significant difference between different
% directions in memory period
i = 1:7;
j = 1:7;
n = size(convolved_regressors, 2);
Fcontrast = zeros(7, n);
Fcontrast(i, j) = eye(7);
Fcontrast(i, 8) = -1;
if inputs.verbose
    figure; imagesc(Fcontrast); colorbar;
    title('Regressor models to be tested in F-contrast');
end

contrast_red = eye(size(convolved_regressors, 2)) - Fcontrast'*pinv(Fcontrast');
regressors_of_interest_reduced = convolved_regressors * contrast_red;
cols_to_remove = [];
for col = 1:size(regressors_of_interest_reduced, 2)
    if all(regressors_of_interest_reduced(:, col) == 0)
        cols_to_remove = [cols_to_remove, col];
    end
end
regressors_of_interest_reduced(:, cols_to_remove) = [];

% Perform F-test
[F, p, resid, df_model, df_error] = F_test_full_vs_red(iDop_scaled, convolved_regressors, regressors_of_interest_reduced, size(Fcontrast, 1));

% Apply MC correction
pvalues_corrected = apply_MCC(p, ...
    'method', 'FDR', ...
    'verbose', false);


%% F-test for any activation
% Testing which voxels respond during memory period (any activation, not
% difference in activation between groups)
i = 1:8;
j = 1:8;
n = size(convolved_regressors, 2);
Fcontrast = zeros(8, n);
Fcontrast(i, j) = eye(8);
if inputs.verbose
    figure; imagesc(Fcontrast); colorbar;
    title('Regressor models to be tested in F-contrast');
end



contrast_red = eye(size(convolved_regressors, 2)) - Fcontrast'*pinv(Fcontrast');
regressors_of_interest_reduced = convolved_regressors * contrast_red;
cols_to_remove = [];
for col = 1:size(regressors_of_interest_reduced, 2)
    if all(regressors_of_interest_reduced(:, col) == 0)
        cols_to_remove = [cols_to_remove, col];
    end
end
regressors_of_interest_reduced(:, cols_to_remove) = [];

% This is the same thing as above, just easier to write.
convolved_regressors = convolved_regressors(:, [1:8 end]);
regressors_of_interest_reduced = convolved_regressors(:, end);


[F_overall, p_overall, ~, ~, ~] = F_test_full_vs_red(iDop_scaled, convolved_regressors, regressors_of_interest_reduced, size(Fcontrast, 1));


p_overall_corrected = apply_MCC(p_overall, ...
    'method', 'FDR', ...
    'verbose', false);

%% Threshold pvalues for Fcontrast and display
if inputs.verbose
    pvalue_threshold = 1e-2;
    F_thresholded = F;
    F_thresholded(pvalues_corrected>pvalue_threshold) = NaN;
    
    F_overall_thresholded = F_overall;
    F_overall_thresholded(p_overall_corrected>pvalue_threshold) = NaN;
    
    pixelsize = 0.1;
    X_img_mm = pixelsize/2 + (0:size(dop,2)-1)*pixelsize; % So that the center of the image is at true halfway point, not shifted.
    Z_img_mm = pixelsize/2 + (0:size(dop,1)-1)*pixelsize;
    %custom colormap for Positive only activation stats
    lateralizationMapComplete = load('LateralizationP.mat');
    PositiveOnlyColormap = lateralizationMapComplete.LateralizationP;
    PositiveOnlyColormap = PositiveOnlyColormap(1:floor(length(PositiveOnlyColormap(:,1))/2),:);
    
    figure;
    plotDuplexImage(X_img_mm, Z_img_mm, F_thresholded, makeAngiogram(dop),...
        'colormap2use', PositiveOnlyColormap, 'nonlinear_bg', 2, ...
        'showColorbar', true, ...
        'colorbarTitle', 'F-score');
    title(sprintf('GLM F-contrast - Difference in directional activity \nqvalue<%0.10f', pvalue_threshold));
    
    figure;
    cmap = plotDuplexImage(X_img_mm, Z_img_mm, F_overall_thresholded, makeAngiogram(dop),...
        'colormap2use', PositiveOnlyColormap, 'nonlinear_bg', 2, ...
        'showColorbar', true, ...
        'colorbarTitle', 'F-score');
    title(sprintf('GLM F-contrast - Any directional response \nqvalue<%0.10f', pvalue_threshold));
end

data.optional_inputs = inputs;
data.F = F;
data.p_unc = p;
data.p_corrected = pvalues_corrected;
data.dop_scaled = iDop_scaled;
data.regressors_convolved = convolved_regressors;
data.regressors = regressor_info;
data.Fcontrast = Fcontrast;
data.F_overall = F_overall;
data.p_overall_unc = p_overall;
data.p_overall_corr = p_overall_corrected;


