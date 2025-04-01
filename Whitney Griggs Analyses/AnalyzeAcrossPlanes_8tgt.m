%% Script to load in multiple planes of data and analyze them
% 1) Load in multiple planes of data
% 2) Create a 4D (or really 5D) dataset.
% 3) Examine how signals propagate through the spatial dimensions compared
%    to time

%% Initialize workspace
close all; clear; clc;


%% Hard code things
% Specify time periods of interest
baselinePeriod = [-1 1];
window_halfsize = 2;
diskFilter = 0;

%% Load colormaps

PositiveNegativeColormap = load('PercentChangeColormap.mat');
PositiveNegativeColormap = PositiveNegativeColormap.PercentChangeColormap;

% Load circular colormap
circle_colormap = load('circular_colormap.mat');
circle_colormap = circle_colormap.colorwheel;


%% Specify which data we want to load
project_record_filename = 'ProjectRecord_paper.json';
ProjectRecord = load_json_as_table(project_record_filename);

session_run_list = specify_sessions_of_interest('project_record_filename', project_record_filename);


% Preallocate space
sessioninfo_indx = NaN(size(session_run_list,1),1);
for i = 1:size(session_run_list,1)
    % set current session/run combo & find row in the ProjectRecord table
    session = session_run_list(i,1);
    run = session_run_list(i,2);
    sessioninfo_indx(i) = intersect(find(ProjectRecord.Session==session),...
        find(ProjectRecord.Run==run));
end

% Check if more than one monkey is used
monkeyPass = isscalar(unique(ProjectRecord.Monkey(sessioninfo_indx)));

if monkeyPass
    Monkey = ProjectRecord.Monkey(sessioninfo_indx(1));
else
    error("No data selected or trying to load more than 1 monkey's data at a time");
end

%% Load in data
% Load in doppler data, where each cell is one coronal slot
loadDopplerData4D(session_run_list, ...
    'project_record_filename', project_record_filename);

% Find all unique slot numbers to examine
ap_planes = unique(session_run_list(:, 3));


%% Preprocess
num_planes = length(ap_planes);
iDopP = cell(length(iDop), 1);

for plane = 1:num_planes
    ap_plane = ap_planes(plane);
    fprintf('Processing ap plane %d \n', ap_plane);
    iDopP{plane} =  preProcess(iDop{plane}, 'timeGain', false, 'spatialFilter', {'disk', diskFilter}, 'zScore', false);
    
end

% convert all to double
for plane = 1:num_planes
    ap_plane = ap_planes(plane);
    iDopP{plane} = double(iDopP{plane});
end

%% Load in GLM masks for each coronal plane
start_dir = pwd;

[F_cell, pvalue_F_unc_cell, pvalue_F_singleSessionCorrected_cell, glm_angiogram_cell] = deal(cell(num_planes, 1));

for plane = 1:num_planes
    ap_plane = ap_planes(plane);
    
    session_load_ind = session_run_list(:, 3) == ap_plane;
    sessionRunList_reduced = session_run_list(session_load_ind, 1:2);
    
    % Logic for filenaming
    % Complete filename up until the date generated
    % Filename will always be stored in ascending order by session, then
    % run, so sort on those features.
    partial_filename = sprintf('glm_data_%sgen*.mat', sprintf('S%dR%d_', sortrows(sessionRunList_reduced, [1, 2])'));
    
    cd(fullfile(get_user_data_path('pathType', 'glm')));
    
    [filename, filepath] = uigetfile(partial_filename, 'Select the GLM data file that matches this session.');
    
    if filepath
        % Can modify glm_data file to save F, pvalues into separate
        % variable for faster access.
        glm_results = load(fullfile(filepath, filename), 'glm_data', 'angiogram', 'coreParams');
        
        % Extract the relevant variables
        F_cell{plane} = glm_results.glm_data.F;
        pvalue_F_unc_cell{plane} = glm_results.glm_data.p_unc;
        pvalue_F_singleSessionCorrected_cell{plane} = glm_results.glm_data.p_corrected;
        glm_angiogram_cell{plane} = glm_results.angiogram;
    else
        error('No GLM data file was chosen. Please try again or create the necessary file using a GLM analysis script.');
    end
end


cd(start_dir);

clear glm_data % This is massive and unnecessary, so we can clear it from the workspace

%% Ideally the GLM mask and the loaded data already match.
% Let's make sure they match
% Align the GLM mask to match the loaded data.
for plane = 1:num_planes
    ap_plane = ap_planes(plane);
    glm_angiogram = glm_angiogram_cell{plane};
    F = F_cell{plane};
    p_unc = pvalue_F_unc_cell{plane};
    p_cor = pvalue_F_singleSessionCorrected_cell{plane};
    
    dop_angiogram = angiogram{plane};
    
    % Use manual alignment to align glm with the dop data
    alignment_tform = align_angiogram(dop_angiogram, glm_angiogram);
    %Convert transform into matlab format
    % Alignment_tform is passed back from ManualImageAlignment
    tform = affine2d(alignment_tform);
    % Align image to match previously loaded data
    glm_angiogram_registered = imwarp(...
        glm_angiogram, ...
        tform, ...
        'OutputView', imref2d(size(dop_angiogram)), ...
        'FillValues', NaN);
    F_registered = imwarp(...
        F, ...
        tform, ...
        'OutputView', imref2d(size(dop_angiogram)), ...
        'FillValues', NaN);
    p_unc_registered = imwarp(...
        p_unc, ...
        tform, ...
        'OutputView', imref2d(size(dop_angiogram)), ...
        'FillValues', NaN);
    p_cor_registered = imwarp(...
        p_cor, ...
        tform, ...
        'OutputView', imref2d(size(dop_angiogram)), ...
        'FillValues', NaN);
    
    % Save back into cell
    F_cell{plane} = F_registered;
    pvalue_F_unc_cell{plane} = p_unc_registered;
    pvalue_F_singleSessionCorrected_cell{plane} = p_cor_registered;
    glm_angiogram_cell{plane} = glm_angiogram_registered;
end

%% Average over each direction
shuffle_direction_labels = false;
anova_effectsize_measure = 'RMSSE';

[iDop_cohenD_8dir_cell, ...
    iDop_percentchange_8dir_cell, iDop_pc_conditions_cell, ...
    pvalue_anova_time_cell, fscore_anova_time_cell, ...
    effectsize_anova_time_cell] = deal(cell(length(ap_planes), 1));
[qvalueFDR_F_cell] = deal(cell(length(ap_planes), 1));

task_timing = nan(length(ap_planes), 4);

for plane = 1:length(ap_planes)
    ap_plane = ap_planes(plane);
    slot_behavior =  behavior{plane};
    
    % Limit analysis to successful trials only
    success = cell2mat({slot_behavior.success});
    slot_behavior(~success) = [];
    
    % Extract timing information about fixation, memory/cue, and movement
    % periods.
    [fix_start_end, memory_end, movement_end] = getEpochs(slot_behavior);
    task_timing(plane, :) = [fix_start_end(1) fix_start_end(2) memory_end(2) movement_end(2)];
    
    % Extract information about targets
    target = cell2mat({slot_behavior.target});
    
    % Extract information about target position
    targetPos = vertcat(slot_behavior.targetPos);
    
    if iscell(targetPos)
        targetPos = cell2mat(targetPos);
    end
    
    % Convert from cartesian coordinates to polar coordinates
    [angle,distance] = cart2pol(targetPos(:,1),targetPos(:,2)); % convert to polar coordinates
    angle(angle<0) = angle(angle<0)+2*pi; %Convert to have only positive angles.
    
    % Find unique directions seen.
    UniqueAngles = unique(angle);
    
    % create an index showing which trials are for which directions.
    % Does not handle distances
    TargetPosInd = zeros(length(target),1);
    for angleInd = 1:length(UniqueAngles)
        TargetPosInd(angle==UniqueAngles(angleInd)) = angleInd;
    end
    
    
    % Assign the direction labels
    conditions = NaN(length(slot_behavior), length(unique(TargetPosInd)));
    % Create an index for each target direction
    for dir = 1:length(unique(TargetPosInd))
        conditions(:, dir) = TargetPosInd == dir; %Left
    end
    
    % Convert into percent change
    % Calculate indices for baseline and memory period
    % Attempting to find the fUS indices that are closest to the bseline
    % period. If using a single point as reference, find the closest. If using
    % multiple points, find all points within desired range.
    if baselinePeriod(1) == baselinePeriod(2)
        baseline_indices = round(baselinePeriod - fix_start_end(1));
    else
        baseline_indices = [ceil(baselinePeriod(1) - fix_start_end(1)) floor(baselinePeriod(2) - fix_start_end(1))];
    end
    
    % Base memStat off length of memory period
    % Round to the nearest integer.
    memStat = [memory_end(2)-window_halfsize memory_end(2)+window_halfsize];
    
    if memStat(1) == memStat(2)
        memory_indices = round(memStat - fix_start_end(1));
    else
        memory_indices = [ceil(memStat(1) - fix_start_end(1)) floor(memStat(2) - fix_start_end(1))];
    end
    
    iDop_pc = convert_to_percent_change(iDopP{plane}, baseline_indices);
    
    % Get trial average for all the conditions
    iDop_pc_conditions = NaN(yPixels(plane), xPixels(plane), nWindows(plane), size(conditions, 2));
    
    for condition = 1:size(conditions, 2)
        trial_ind = logical(conditions(:, condition));
        iDop_pc_conditions(:, :, :, condition) = squeeze(mean(iDop_pc(:, :, :, trial_ind), 4));
    end
    
    iDop_pc_conditions_cell{plane} = iDop_pc_conditions;
    
    
    [iDop_cohenD_8dir, iDop_percentchange_8dir] = deal(NaN(yPixels(plane), xPixels(plane), nWindows(plane), size(conditions, 2)));
    [effectsize_anova_time] = deal(NaN(yPixels(plane), xPixels(plane), nWindows(plane)));
    
    fprintf('### AP Plane: %d ### \n', ap_plane);
    nTimepoints = size(iDop_pc_conditions, 3);
    nConditions = size(conditions, 2);
    zPix = size(iDop_pc_conditions, 1);
    xPix = size(iDop_pc_conditions, 2);
    
    baselineIDop = squeeze(mean(iDop_pc(:, :, baseline_indices, :), 3));
    
    % Moving parfor  down a level, because the code hangs.
    % Depends on size of data loaded. Might be able to move parfor back up
    % with a more pwerful computer with more RAM.
    for condition = 1:nConditions
        fprintf('## Condition: %d ## \n', condition);
        pb = ProgressBar(zPix * xPix * nTimepoints);
        trial_ind = logical(conditions(:, condition));
        parfor timepoint = 1:nTimepoints
            for z = 1:zPix
                for x = 1:xPix
                    compare1 = squeeze(iDop_pc(z, x, timepoint, trial_ind));
                    compare2 = squeeze(baselineIDop(z, x, trial_ind));
                    
                    % Different measures of effect size
                    iDop_cohenD_8dir(z, x, timepoint, condition) = calculate_cohen_d(compare1, compare2);
                    iDop_percentchange_8dir(z, x, timepoint, condition) = calculate_mean_percentchange(compare1, compare2);
                    
                    pb.count;
                end
            end
            
        end
        delete(pb);
        
    end
        
    iDop_cohenD_8dir_cell{plane} = iDop_cohenD_8dir;
    iDop_percentchange_8dir_cell{plane} = iDop_percentchange_8dir;
    
    % Run 1-way anova at each timepoint
    pb = ProgressBar(zPix * xPix * nTimepoints);
    
    [fscore_anova_time, pvalue_anova_time] = deal(NaN(zPix, xPix, nTimepoints));
    
    parfor timepoint = 1:nTimepoints
        for z = 1:zPix
            for x = 1:xPix
                patch_activity = squeeze(iDop_pc(z, x, timepoint, :));
                
                [pval, tbl, ~] = anova1(patch_activity, TargetPosInd, 'off');
                pvalue_anova_time(z, x, timepoint) = pval;
                
                % Convert cell into table.
                anova_table = cell2table(tbl(2:end, :), 'VariableNames', tbl(1, :));
                
                fscore_anova_time(z, x, timepoint) = cell2mat(anova_table.F(1));
                
                % Calculate effect size
                switch anova_effectsize_measure
                    case 'RMSSE'
                        effectsize_anova_time(z, x, timepoint) = calculate_RMS_standardized_effect(patch_activity, TargetPosInd);
                    case 'eta_squared'
                        % eta_squared = SS_between/SS_total
                        effectsize_anova_time(z, x, timepoint) = anova_table.SS(1)/anova_table.SS(3);
                end
                pb.count;
            end
        end
    end
    delete(pb);
    
    pvalue_anova_time_cell{plane} = pvalue_anova_time;
    fscore_anova_time_cell{plane} = fscore_anova_time;
    effectsize_anova_time_cell{plane} = effectsize_anova_time;
end

%% Apply FDR correction
% Collect all pvalues across planes
[pvalue_1D, pvalue_anova_time_1D] = deal([]);
[number_elements_per_slot, number_elements_per_slot_anova_time] = deal(NaN(length(ap_planes), 1));
for plane = 1:length(ap_planes)
    pvalue_temp = pvalue_F_unc_cell{plane};
    pvalue_temp = reshape(pvalue_temp, [], 1);
    number_elements_per_slot(plane) = numel(pvalue_temp);
    pvalue_1D = [pvalue_1D; pvalue_temp];
    
    pvalue_anova_time_temp = pvalue_anova_time_cell{plane};
    pvalue_anova_time_temp = reshape(pvalue_anova_time_temp, [], 1);
    number_elements_per_slot_anova_time(plane) = numel(pvalue_anova_time_temp);
    pvalue_anova_time_1D = [pvalue_anova_time_1D; pvalue_anova_time_temp];
end

disp('Applying FDR correction now');
% Apply FDR correction for planes, conditions, timepoints, and across voxels
[Q] = mafdr(pvalue_1D, 'BHFDR', true);
pvalues_corrected = Q;

[Q] = mafdr(pvalue_anova_time_1D, 'BHFDR', true);
pvalues_anova_time_corrected = Q;

% Reassemble the planes
for plane = 1:length(ap_planes)
    pvalues_for_slot = pvalues_corrected(sum(number_elements_per_slot(1:plane-1))+1:sum(number_elements_per_slot(1:plane-1))+number_elements_per_slot(plane));
    pvalues_anova_time_for_slot = pvalues_anova_time_corrected(sum(number_elements_per_slot_anova_time(1:plane-1))+1:sum(number_elements_per_slot_anova_time(1:plane-1))+number_elements_per_slot_anova_time(plane));
    
    image_size = size(pvalue_F_unc_cell{plane});
    iDop_pvalue = reshape(pvalues_for_slot, image_size(1), image_size(2));
    qvalueFDR_F_cell{plane} = iDop_pvalue;
    
    image_size_anova_time = size(pvalue_anova_time_cell{plane});
    pvalue_anova_time = reshape(pvalues_anova_time_for_slot, image_size_anova_time(1), image_size_anova_time(2), image_size_anova_time(3));
    qvalueFDR_anova_time_cell{plane} = pvalue_anova_time;
end
disp('Done applying FDR correction');


%%  Normalized Center of mass calculation
% Here I scale the magnitudes by the largest amplitude signal (positive or
% negative)

% Temporarily turning off warnings to reduce spurious output to the command
% window
id1 = 'MATLAB:polyshape:repairedBySimplify';
warning('off',id1);
id2 = 'MATLAB:polyshape:boundary3Points';
warning('off',id2);

% Parallelized for speed improvement.
min_windows = min(nWindows(nWindows > 0));
pb = ProgressBar(length(ap_planes)*min_windows*numel(angiogram{end}));
gcp; pctRunOnAll warning('off');
parfor plane = 1:length(ap_planes)
    ap_plane = ap_planes(plane);
    iDop_effectsize_8dir = iDop_cohenD_8dir_cell{plane};
    
    
    [iDop_tuningcurve_theta, iDop_tuningcurve_rho] = deal(NaN(...
        size(iDop_effectsize_8dir, 1), ...
        size(iDop_effectsize_8dir, 2), ...
        size(iDop_effectsize_8dir, 3)));
    
    circular_statistics = cell(min_windows, 1);
    
    for timenum = 1:min_windows
        theta = repmat(UniqueAngles', yPixels(plane) * xPixels(plane), 1);
        rho = reshape(squeeze(iDop_effectsize_8dir(:, :, timenum, :)), [], 8);
        rho_scaled = rho./repmat(max(abs(rho), [], 2), 1, 8);
        
        [x, y] = pol2cart(theta, rho_scaled);
        
        %Find centroid for the tuning curve at each pixel
        [centroidTheta, centroidRho] = deal(NaN(size(x, 1), 1));
        for pixel = 1:size(x,1)
            polyIn = polyshape([x(pixel,:)',y(pixel,:)']);
            [centroidX, centroidY] = centroid(polyIn);
            [centroidTheta(pixel), centroidRho(pixel)] = cart2pol(centroidX,centroidY);
            pb.count;
        end
        
        
        %Reshape the centroids back to the x,z dimensions of iDopP
        centroidTheta = reshape(centroidTheta, yPixels(plane), xPixels(plane));
        centroidRho = reshape(centroidRho, yPixels(plane), xPixels(plane));
        
        % Convert all centroidTheta values to be between [- pi pi]
        greaterThanPIind = centroidTheta > pi;
        lessThanNegPIind = centroidTheta < - pi;
        
        centroidTheta(greaterThanPIind) = centroidTheta(greaterThanPIind) - 2*pi;
        centroidTheta(lessThanNegPIind) = centroidTheta(lessThanNegPIind) + 2*pi;
        
        % Store in iDop_tuningcurve
        iDop_tuningcurve_theta(:, :, timenum) = centroidTheta;
        iDop_tuningcurve_rho(:, :, timenum) = centroidRho;
        
        % Calculate other circular statistics
        circular_statistics{timenum} = get_circular_statistics(theta, rho);
    end
    
    iDop_tuningcurve_theta_cell{plane} = iDop_tuningcurve_theta;
    iDop_tuningcurve_rho_cell{plane} = iDop_tuningcurve_rho;
    iDop_tuningcurve_circ_stats{plane} = circular_statistics;
end
delete(pb);

% Turn back on the warnings
warning('on',id1);
warning('on',id2);


%% Save processed data for later loading
save_workspace = true;

if save_workspace
    data_save_root = get_user_data_path('pathType', 'across session analyses');
    save_fullfile = fullfile(data_save_root, sprintf('%s_8dir_tuningCurve_data_dateGenerated_%s.mat', Monkey, datetime('today', format='yyyyMMdd')));
    save(save_fullfile, 'ap_planes', 'session_run_list',...
        'iDop_tuningcurve_theta_cell', 'iDop_tuningcurve_rho_cell', ...
        'UF', 'angiogram', 'circle_colormap',...
        'Monkey', 'iDop_cohenD_8dir_cell', ...
        'iDop_percentchange_8dir_cell',...
        'qvalueFDR_F_cell', 'glm_angiogram_cell', 'F_cell', ...
        'anova_effectsize_measure',...
        'effectsize_anova_time_cell', ...
        'nWindows', 'task_timing', 'qvalueFDR_anova_time_cell', ...
        'fscore_anova_time_cell', 'pvalue_anova_time_cell', ...
        'pvalue_F_unc_cell', 'iDop_tuningcurve_circ_stats')
end
