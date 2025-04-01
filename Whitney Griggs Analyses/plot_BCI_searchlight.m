%% Plots for BCI_CV_direction_multicoder_searchlight

clear


%% hard coded settings
saveMe = true;

motion_correction_method = 'normcorre';
pvalue_threshold = 1e-3;




%% Load previously created decoding data
% Specify session of interest
project_record_filename = 'ProjectRecord_paper.json';
session_run_list = specify_sessions_of_interest('project_record_filename', project_record_filename, ...
    'single_v_multiple', 'single');

if ~exist('result_combined', 'var')
    [loadfile, loadpath] = uigetfile(fullfile(get_user_data_path('pathType', 'decoding'), sprintf('searchlight_data_S%dR%d_mc*.mat', session_run_list)));
    load_fullfile = fullfile(loadpath, loadfile);
    
    % Load data
    sl_data = load(load_fullfile);
end

%% Load Doppler data to extract behavioral information
% This is overkill and I am sure there is a better way than this, but
% hacking it together right now for paper review. February 1, 2025 WG
loadDopplerData(whos, 'multiple', false,...
    'mc_method', motion_correction_method, ...
    'project_record_filename', project_record_filename, ...
    'SessionRunList', session_run_list);

% Extract behavioral info
[fix, memory_end, mov_end] = getEpochs(behavior);

%% Extract specific variables of interest
Session = session_run_list(1);
Run = session_run_list(2);
sl_result_combined = sl_data.sl_result_combined;
radii_to_test = sl_data.radii_to_test;
timepoints_to_test = sl_data.timepoints_to_test;

%% Extract `empiric_mean_angular_error` from `result_combined`
% This is unnecessary if running entire script, but added so that I can
% load previously saved data instead.
[sl_pvalue_angular_error, sl_empiric_mean_angular_error, ...
    sl_accuracy, sl_pvalue_pc] = deal(NaN(size(sl_result_combined)));
for i = 1:size(sl_result_combined, 1)
    for j = 1:size(sl_result_combined, 2)
        for k = 1:size(sl_result_combined, 3)
            for m = 1:size(sl_result_combined, 4)
                try
                    sl_empiric_mean_angular_error(i, j, k, m) = sl_result_combined{i, j, k, m}.angular_error;
                    sl_accuracy(i, j, k, m) = sl_result_combined{i, j, k, m}.percentCorrect;
                    sl_pvalue_pc(i, j, k, m) = sl_result_combined{i, j, k, m}.p;
                    sl_pvalue_angular_error(i, j, k, m) = sl_result_combined{i, j, k, m}.angular_error_pvalue;
                catch
                    sl_empiric_mean_angular_error(i, j, k, m) = NaN;
                    sl_accuracy(i, j, k, m) = NaN;
                    sl_pvalue_pc(i, j, k, m) = NaN;
                    sl_pvalue_angular_error(i, j, k, m) = NaN;
                end
            end
        end
    end
end

%% Specify which thing we want to plot
metric_to_use = 'angular_error'; % 'angular_error' or 'percent_correct'
switch metric_to_use
    case 'angular_error'
        sl_metric = sl_empiric_mean_angular_error * 180/pi; %in degrees now
        sl_pvalue = sl_pvalue_angular_error;
        %colorbar_limits = [floor(min(sl_metric,[], 'all')) 90];
        colorbar_limits = [30 90];
        colormap2use = flipud(magma);
        colorbar_title = 'Angular error (deg)';
    case 'percent_correct'
        sl_metric = sl_accuracy;
        sl_pvalue = sl_pvalue_pc;
        colorbar_limits = [0 ceil(max(sl_metric,[], 'all'))];
        colormap2use = magma;
        colorbar_title = 'Accuracy (% correct)';
end


%% Visualize searchlight analysis

% Select method for FDR correction - planewise or across all timepoints
FDR_method = 'all'; % 'planewise' or 'all'
FDR_debug = false;
if FDR_debug
    FDR_method = 'planewise'; 
end

n_timepoints = size(sl_result_combined, 3);

trialTime = fix(1):1/coreParams.framerate:fix(1)+nWindows/coreParams.framerate-1/coreParams.framerate;


pixelsize = 0.1;
X_img_mm = pixelsize/2 + (0:size(angiogram,2)-1)*pixelsize; % So that the center of the image is at true halfway point, not shifted.
Z_img_mm = pixelsize/2 + (0:size(angiogram,1)-1)*pixelsize + UF.Depth(1);




circle_center = [ 2 max(Z_img_mm)-2];
for k = 1:length(radii_to_test)
    if strcmp(FDR_method, 'all')
        pvalues_1D = reshape(sl_pvalue(:, :, :, k),[],1);
        nan_ind = isnan(pvalues_1D);
        [Q] = mafdr(pvalues_1D(~nan_ind), 'BHFDR', true);
        pvalues_corrected = NaN(size(pvalues_1D));
        pvalues_corrected(~nan_ind) = Q;
        matrix_size = [size(sl_pvalue_pc, 1) size(sl_pvalue_pc, 2) size(sl_pvalue_pc, 3)];
        sl_pvalue_combined_corrected_all = reshape(pvalues_corrected, matrix_size);
    end


    radius = radii_to_test(k);
    % Plot searchlight analysis results
    fig_searchlight = figure('units','normalized','outerposition',[0 0 1 1]);
    tld = tiledlayout('flow');
        

    if FDR_debug & exist('sl_pvalue_combined_corrected_all', 'var')
        fig_debug = figure('units','normalized','outerposition',[0 0 1 1]);
        tld_debug = tiledlayout('flow');
    end

    for timepoint = 1:n_timepoints
        trial_timepoint = timepoints_to_test(timepoint);

        when_in_trial = trialTime(trial_timepoint);
        if when_in_trial < 0
            task_state = 'fixation';
        elseif when_in_trial == 0
            task_state = 'cue';
        elseif when_in_trial > 0 && when_in_trial < memory_end(2)
            task_state = 'memory';
        elseif when_in_trial > memory_end(2) && when_in_trial < mov_end(2)
            task_state = 'movement';
        elseif when_in_trial > mov_end(2)
            task_state = 'ITI';
        end

        if strcmp(FDR_method, 'planewise')
            % Apply FDR correction for the p-values
            pvalues_1D = reshape(sl_pvalue(:, :, timepoint, k),[],1);
            nan_ind = isnan(pvalues_1D);
            [Q] = mafdr(pvalues_1D(~nan_ind), 'BHFDR', true);
            pvalues_corrected = NaN(size(pvalues_1D));
            pvalues_corrected(~nan_ind) = Q;
            matrix_size = [size(sl_pvalue_pc, 1) size(sl_pvalue_pc, 2)];
            sl_pvalue_combined_corrected = reshape(pvalues_corrected, matrix_size);
        elseif strcmp(FDR_method, 'all')
            sl_pvalue_combined_corrected = sl_pvalue_combined_corrected_all(:, :, timepoint);
        end
    
        % Overlay performance on anatomy
        pvalue_mask = sl_pvalue_combined_corrected < pvalue_threshold;

        if FDR_debug & exist('sl_pvalue_combined_corrected_all', 'var')
            % Check to see difference between the 'all' and 'planewise' FDR
            % methods
                pvalue_mask_from_FDRall = 2 * double(sl_pvalue_combined_corrected_all(:, :, timepoint) < pvalue_threshold);
                difference_in_significant_voxels = pvalue_mask_from_FDRall + double(pvalue_mask);
                difference_in_significant_voxels(difference_in_significant_voxels==0) = NaN;
            

                %nexttile(tld_debug, timepoint);
                figure;
                plotDuplexImage(X_img_mm, Z_img_mm, difference_in_significant_voxels, angiogram,...
                    'colormap2use',distinguishable_colors(3, 'k'),...
                    'nonlinear_bg',2, ...
                    'AutoColorBarLimits',[1 3], ...
                    'showColorbar', true, ...
                    'ColorbarTitle', '1=planewise, 2=all, 3=both');
        end
    
        performance_masked = squeeze(sl_metric(:, :, timepoint, k));
        performance_masked(~pvalue_mask) = NaN;
    
        nexttile(tld, timepoint);
        if timepoint ==n_timepoints
            show_colorbar = true;
        else
            show_colorbar = false;
        end
        cmap = plotDuplexImage(X_img_mm, Z_img_mm, performance_masked, angiogram,...
            'colormap2use',colormap2use, ...
            'nonlinear_bg',2, ...
            'AutoColorBarLimits',colorbar_limits, ...
            'showColorbar', show_colorbar, ...
            'ColorbarTitle', colorbar_title);
        axis off;
        title(sprintf('%0.2f sec - %s', when_in_trial, task_state));
    
        hold on;
        circle_handle = viscircles(circle_center,radius*pixelsize, 'EnhanceVisibility', false, ...
            'Color', 'w');
        hold off;
    end
    title(tld, sprintf('Searchlight analysis S%dR%d (radius = %0.2f)', Session, Run, radius))
end

save_to_file = true;
if save_to_file && saveMe

    suggested_fname = sprintf('searchlight_results_S%dR%d_gen%s.svg', Session, Run, datetime('now', 'Format', 'yyyyMMdd'));

    [filename, filepath] = uiputfile(fullfile(get_user_data_path('pathType', 'decoding'), suggested_fname));
    if any([filename filepath])
        saveas(fig_searchlight, fullfile(filepath, filename));
    end
end

