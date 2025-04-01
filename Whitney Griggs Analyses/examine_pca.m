%% For directional tuning analysis where I am looking at PCA
% How many components typically ekpt?
% Any clear relationship to trial variables?
% January 15, 2025 WG

close all

startDir = pwd;
%% hard coded settings
% choose K for K-fold validation
K = 10;
% preprocessing
timeGain = false;
zScore = false;
detrend = true;

% Define parameters for the Gaussian smoothing kernel
FWHM = 1; % In voxels; Common values in fMRI is 1-3 voxels
sigma = FWHM/(sqrt(8*log(2)));
filter_size = 2*ceil(2*sigma)+1; % This is the default filter size used in `imgaussfilt`.
filter_params = {'gaussian', filter_size, sigma};

diskFilter = 2;

% set parameters for computing the mask (if requested)
baselinePeriod = -4:-1;
memPeriod = 2:5;    % relative to memory cue
pThresh = 0.01;     % p value threshold (doesn't matter)
e = 6;             % pixel resolution of map
saveMe = false;

motion_correction_method = 'normcorre';


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


% Apply spatial filter to entire session.
%iDopP = preProcess(iDopP, 'spatialFilter', filter_params, 'zscore', zScore);

if ~exist('fix','var') || ~exist('mem','var') || ~exist('mov','var')
    [fix, memm, mov] = getEpochs(behavior);
end
if exist('mem','var') && ~exist('memm','var')
    memm = mem;
end

[yPix, xPix, nWindows, nTrials] = size(iDopP);

trialTime = fix(1):1/coreParams.framerate:fix(1)+nWindows/coreParams.framerate-1/coreParams.framerate;


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
memTime = ceil(abs(memm(end))*coreParams.framerate);

% Initialize variable
numComponentsToKeep = NaN(nWindows, K, 2); % 2 for horz and vert
pb = ProgressBar(nWindows*K);

disp('Time with for loop');
tic
variance_to_keep = 95;
for time_window = 1:nWindows
    %for time_window = fixTime+5
    % this is the meat of this whole loop where we continue to pile on data
    % (concatenate it, really) as we progress through time. This is on the
    % underlying assumption that the classifier is going to use data from
    % whatever epoch of the task we're currently in to classify data from
    % that same epoch.
    % use all of the data during the fixation period...
    if time_window <= fixTime
        testTimeI = 1:time_window;
        % or use all the data during the memory delay and movement
    else
        testTimeI = fixTime:time_window;
    end
    %testTimeI = i;
    time_window
    
    % prepare the training data for this time window
    data = flattenDoppler2D(iDop_crop, testTimeI);  % nImages x (nPixels*yPixels)
    
    N = size(data,1);           % number of trials
    
    % k-fold & leave one out cross validation
    
    % cross validate here
    indices = crossvalind('Kfold', N, K);
    % for each k-fold
    for k_loop = 1:K
        % create indices for training set & test set
        test = (indices == k_loop);
        train = ~test;
        nan_ind = isnan(train_labels);
        
        % Accounting for NaN labels in case we are not using all possible
        % classes of data.
        horz_train = train & ~nan_ind(:,1);
        vert_train = train & ~nan_ind(:,2);
        
        % Fit and apply z-scoring to train data, then apply to test.
        [data(train, :), mu, sigma] = zscore(data(train, :));
        data(test, :) = (data(test, :) - mu) ./ sigma;
        
        %If nans, then define as 0. NaN values will break downstream
        %functions.
        data(isnan(data)) = 0;
        
        % pca transform
        [numComponentsToKeep(time_window, k_loop, 1), coeffs_horz, scores_horz] = pca_transform(data(horz_train,:), variance_to_keep);
        [numComponentsToKeep(time_window, k_loop, 2), coeffs_vert, scores_vert] = pca_transform(data(vert_train,:), variance_to_keep);
        
        % Iterate counter
        pb.count();
    end
end
delete(pb);

%% Preview top 20 PCA components for last fold of K-fold validation with last timeWindow block

n_timepoints = length(testTimeI);
n_components_to_examine = 10;
%Horizontal and vertical PCA  the same since same trials used for the PCA,
% so only need 1 figure
figure('Units', 'Normalized', 'OuterPosition', [0 0 1 1], 'Renderer', 'Painters');
tld = tiledlayout(n_components_to_examine, n_timepoints, 'TileSpacing', 'none', 'Padding', 'tight');
colormap magma

for timeWindow = 1:n_timepoints
    for i = 1:n_components_to_examine
        nexttile(tld, timeWindow+((i-1)*n_timepoints));
        
        
        data_val_pc = coeffs_horz(:,i);
        data_value_pc_reshape = reshape(data_val_pc, [yPix, xPix, n_timepoints]);
        data_value_pc_reshape_singleTimepoint = data_value_pc_reshape(:, :, timeWindow);
        
        imagesc(data_value_pc_reshape_singleTimepoint)
        axis image
        if timeWindow == 1
            ylabel(sprintf('PC component %d', i));
            yticks([]);
        else
            ax1 = gca;
            ax1.YAxis.Visible = 'off'; % remove y-axis
            
        end
        if i == n_components_to_examine
            xlabel(sprintf('Time: %0.2f', trialTime(testTimeI(timeWindow))));
            xticks([]);
        else
            ax1 = gca;
            ax1.XAxis.Visible = 'off'; % remove y-axis
            xlabel('');
        end
    end
end

title(tld, 'PCA');

if saveMe
    filepath = fullfile(get_user_data_path('pathType', 'output'));
    filename = sprintf('PCAdecoder_first20components_S%dR%d_dg_%s.svg', sessionRunList', datetime('now', 'format', 'yyyyMMdd'));
    saveas(gcf, fullfile(filepath, filename));
end
%% Average across k_loop
numComponentsToKeep_mean = squeeze(mean(numComponentsToKeep, 2));
numComponentsToKeep_std = squeeze(std(numComponentsToKeep, [], 2));

numComponentsToKeep_min = squeeze(min(numComponentsToKeep, [], 2));
numComponentsToKeep_max = squeeze(max(numComponentsToKeep, [], 2));

%% Visualize the number of components

range_above = abs(numComponentsToKeep_max(:, 1) - numComponentsToKeep_mean(:, 1));
range_below = abs(numComponentsToKeep_min(:, 1) - numComponentsToKeep_mean(:, 1));


figure;
shadedErrorBar(trialTime, numComponentsToKeep_mean(:, 1), [range_below'; range_above']);

xlabel('Trial time (s)');
ylabel('# components');
title(sprintf('PCA keeping %d%% variance', variance_to_keep));
legend('Range', 'Mean');

if saveMe
    filepath = fullfile(get_user_data_path('pathType', 'output'));
    filename = sprintf('PCAdecoder_NumComponentsToExplain95Variance_S%dR%d_dg_%s.svg', sessionRunList', datetime('now', 'format', 'yyyyMMdd'));
    saveas(gcf, fullfile(filepath, filename));
end

%% Reshape last coeffs into ntimepoints x nPixels
coeffs_horz_reshape = reshape(coeffs_horz, [nWindows-fixTime+1, yPix*xPix, size(coeffs_horz, 2)]);

%% For Reviewer 4, they requested an additional analysis
% Also do the top PCs correspond to direction tuning? A quick way to look
% at this would be to correlate each voxel’s directional tuning strength
% to their weight onto each component. It would be interesting to know if
% saccade direction explained a large amount of variance across all
% components or specific components.

% Plan is to load GLM dataset with tuning strength alraedy calculated, then
% correlate pixel by pixel with components weights for first 20 components

% Load data
[loadfile, loadpath] = uigetfile(fullfile(get_user_data_path('pathType', 'glm'), sprintf('S%dR%d_8tgt_analysis_mc*dg*.mat', sessionRunList')));

processed_data = load(fullfile(loadpath, loadfile));
pvalue_mask = processed_data.p_cor_registered > 1e-3 | isnan(processed_data.p_cor_registered);

centroidTheta = processed_data.centroidTheta;
centroidTheta_plot = centroidTheta;
centroidTheta_plot(pvalue_mask) = NaN;

centroidRho = processed_data.centroidRho;
centroidRho_plot = centroidRho;
centroidRho_plot(pvalue_mask) = NaN;
nan_ind = isnan(centroidRho_plot);



for j = 1:size(coeffs_horz_reshape, 1)
    % Correlate rho (strength) with different components
    PCA_strength_fig = figure('Units', 'normalized', 'OuterPosition', [0 0 1 1], 'Renderer', 'painter');
    tld = tiledlayout('flow');
    nan_ind = isnan(centroidRho);
    for i = 1:n_components_to_examine
        nexttile(i);
        pc_val = coeffs_horz_reshape(j, :, i);
        corr_value = corrcoef(centroidRho(~nan_ind)', pc_val(~nan_ind)');
        scatter(centroidRho(:), pc_val, 'AlphaData',0.1*ones(size(pc_val)));
        legend(sprintf('Pearson R - %0.2f', corr_value(1, 2)));
        xlabel('Tuning strength');
        ylabel('PC value');
        title(sprintf('PC%d', i));
    end
    title(tld, sprintf('pixelwise correlation between PC weights and tuning strength - timepoint %d', j));
    
    % Correlate F score (statistical significance) with different components
    PCA_fscore_fig = figure('Units', 'normalized', 'OuterPosition', [0 0 1 1], 'Renderer', 'painter');
    tld = tiledlayout('flow');
    for i = 1:n_components_to_examine
        nexttile(i);
        pc_val = coeffs_horz_reshape(j, :, i);
        corr_value = corrcoef(processed_data.F_registered(:)', pc_val');
        scatter(processed_data.F_registered(:), pc_val, 'AlphaData',0.1*ones(size(pc_val)));
        legend(sprintf('Pearson R - %0.2f', corr_value(1, 2)));
        xlabel('F-score');
        ylabel('PC value');
        title(sprintf('PC%d', i));
    end
    title(tld, sprintf('pixelwise correlation between PC weights and F-score - timepoint %d', j));
end

if saveMe
    filepath = fullfile(get_user_data_path('pathType', 'output'));
    filename = sprintf('PCAdecoder_vTuningStrength_timepoint4afterfix_S%dR%d_dg_%s.svg', sessionRunList', datetime('now', 'format', 'yyyyMMdd'));
    saveas(PCA_strength_fig, fullfile(filepath, filename));
end

if saveMe
    filepath = fullfile(get_user_data_path('pathType', 'output'));
    filename = sprintf('PCAdecoder_vFscore_timepoint4afterfix_S%dR%d_dg_%s.svg', sessionRunList', datetime('now', 'format', 'yyyyMMdd'));
    saveas(PCA_fscore_fig, fullfile(filepath, filename));
end

%% Internal functions

function [numComponentsToKeep, pcaCoefficients, data_transform] = pca_transform(data, variance_to_keep)
[pcaCoefficients, pcaScores, ~, ~, explained, ~] = pca(data);
explainedVarianceToKeepAsFraction = variance_to_keep/100;

numComponentsToKeep = find(cumsum(explained)/sum(explained) >= explainedVarianceToKeepAsFraction, 1);
pcaCoefficients = pcaCoefficients(:,1:numComponentsToKeep);
data_transform = pcaScores(:,1:numComponentsToKeep);   % (PCA transformed trainData)

end