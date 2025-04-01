%% Compute Sorenson Dice coefficient between searchlight and GLM analyses
% Whitney Griggs December 17, 2024

clear; clc; close all;
%% User settings
debug = true;
saveMe = false;

% Specify the pvalue to threshold maps at
pvalue_threshold = 1e-3;
percent_of_voxels_to_keep = 0.10; % Keep 10% most significant voxels

% For searchlight
display_type = 'pvalue_threshold'; % 'bottom_voxels' or 'top_voxels' Or 'pvalue_threshold'
display_accuracy_or_angularError = 'angular_error'; % 'angular_error' or 'accuracy'
radii_to_test = 2;

% Motion correction method
mc_method = 'normcorre';

if mc_method
    mc_method_string = mc_method;
else
    mc_method_string = 'false';
end

%% Iterate over each requested session comparing the Sorenson-Dice coefficient
% Select which sessions we want to iterate over
project_record_filename = 'ProjectRecord_paper.json';
ProjectRecord = load_json_as_table(project_record_filename);
session_run_list = specify_sessions_of_interest(...
    'project_record_filename', project_record_filename);


%% Loop over each session
keep_session = true(size(session_run_list, 1), 1);

num_sessions = size(session_run_list, 1);
dice_similarity = NaN(num_sessions, 1);
for session_num = 1:num_sessions
    try
        % Load in GLM masks for each coronal plane
        start_dir = pwd;

        % Logic for filenaming
        % Complete filename up until the date generated
        partial_filename = sprintf('glm_data_%sgen*.mat', sprintf('S%dR%d_', session_run_list(session_num, :)'));

        cd(fullfile(get_user_data_path('pathType', 'glm')));

        % Find most recent matching file
        file_listing = dir(partial_filename);
        [~,most_recent_file_index] = max([file_listing.datenum]);
        most_recent_file = file_listing(most_recent_file_index).name;

        % Load this most recent file
        fprintf('Loading %s: file %d/%d\n', most_recent_file, session_num, num_sessions);
        load(most_recent_file, 'F', 'pvalues_corrected', ...
            'angiogram');

        % Rename for convenience
        glm_angiogram = angiogram;
        glm_pvalues = pvalues_corrected;
        glm_F = F;

        % Define GLM mask
        glm_mask = ~(glm_pvalues>pvalue_threshold);

        % Return to original directory
        cd(start_dir);

        % Define axes labels for eventual plots
        pixelsize = 0.1;
        X_img_mm = pixelsize/2 + (0:size(F,2)-1)*pixelsize; % So that the center of the image is at true halfway point, not shifted.
        Z_img_mm = pixelsize/2 + (0:size(F,1)-1)*pixelsize;


        if debug
            % Threshold pvalues for Fcontrast and display
            F_thresholded = F;
            F_thresholded(~glm_mask) = NaN;


            figure;
            transparency_map = isnan(F_thresholded);

            plotDuplexImage(X_img_mm, Z_img_mm, F_thresholded, glm_angiogram,...
                'colormap2use', flipud(magma), 'nonlinear_bg', 2, ...
                'showColorbar', true, ...
                'colorbarTitle', 'F-score',...
                'darkestShade', 0.25,...
                'transparencyMap', ~transparency_map);
            title(sprintf('%s \nqvalue<%0.10f', sprintf('S%dR%d ', session_run_list(session_num, :)'), pvalue_threshold));
        end

        %% Load decoding dataset
        partial_filename = sprintf('data_S%dR%d*.mat', session_run_list(session_num, :)');

        cd(fullfile(get_user_data_path('pathType', 'decoding')));

        % Find most recent matching file
        file_listing = dir(partial_filename);
        [~,most_recent_file_index] = max([file_listing.datenum]);
        most_recent_file = file_listing(most_recent_file_index).name;

        % Load the data
        decoding_results = load(most_recent_file);

        % Extract `empiric_mean_angular_error` from `result_combined`
        [sl_pvalue_angular_error, sl_empiric_mean_angular_error, ...
            sl_accuracy, sl_pvalue] = deal(NaN(size(decoding_results.sl_result_combined)));
        for i = 1:size(decoding_results.sl_result_combined, 1)
            for j = 1:size(decoding_results.sl_result_combined, 2)
                for k = 1:size(decoding_results.sl_result_combined, 3)
                    try
                        sl_empiric_mean_angular_error(i, j, k) = decoding_results.sl_result_combined{i, j, k}.angular_error;
                        sl_accuracy(i, j, k) = decoding_results.sl_result_combined{i, j, k}.percentCorrect;
                        sl_pvalue(i, j, k) = decoding_results.sl_result_combined{i, j, k}.p;
                        sl_pvalue_angular_error(i, j, k) = decoding_results.sl_result_combined{i, j, k}.angular_error_pvalue;
                    catch
                        sl_empiric_mean_angular_error(i, j, k) = NaN;
                        sl_accuracy(i, j, k) = NaN;
                        sl_pvalue(i, j, k) = NaN;
                        sl_pvalue_angular_error(i, j, k) = NaN;
                    end
                end
            end
        end

        radii_index = find(radii_to_test ==2);

        % Apply FDR correction for the p-values
        pvalues_1D = reshape(sl_pvalue,[],1);
        noNaN_ind = isnan(pvalues_1D);
        [Q] = mafdr(pvalues_1D(~noNaN_ind), 'BHFDR', true);
        pvalues_corrected = NaN(size(pvalues_1D));
        pvalues_corrected(~noNaN_ind) = Q;
        matrix_size = size(sl_pvalue);
        sl_pvalue_corrected = reshape(pvalues_corrected, matrix_size);

        % Apply FDR correction for the angular error p-values
        pvalues_1D = reshape(sl_pvalue_angular_error,[],1);
        noNaN_ind = isnan(pvalues_1D);
        [Q] = mafdr(pvalues_1D(~noNaN_ind), 'BHFDR', true);
        pvalues_corrected = NaN(size(pvalues_1D));
        pvalues_corrected(~noNaN_ind) = Q;
        matrix_size = size(sl_pvalue_angular_error);
        sl_pvalue_angular_error_corrected = reshape(pvalues_corrected, matrix_size);


        %% Pick appropriate SL dataset
        colormap_to_use = inferno;
        switch display_accuracy_or_angularError
            case 'angular_error'
                pvalue_to_use = sl_pvalue_angular_error_corrected;
                performance_metric_to_use = sl_empiric_mean_angular_error * 180/pi;
                colorbar_title = 'Angular error (deg)';
                unit_measure = ' deg';
                colorscale = flipud(colormap_to_use);

                colorbar_max = 90;
                colorbar_min = 30;
            case 'accuracy'
                pvalue_to_use = sl_pvalue_corrected;
                performance_metric_to_use = sl_accuracy;
                colorbar_title = 'Performance (% correct)';
                unit_measure = '%';
                colorscale = colormap_to_use;

                colorbar_max = 40;
                colorbar_min = 0;
        end


        % Overlay performance on anatomy
        switch display_type
            case 'pvalue_threshold'
                pvalue_mask = pvalue_to_use(:, :, radii_index) < pvalue_threshold;
                performance_masked = squeeze(performance_metric_to_use(:, :, radii_index));
                performance_masked(~pvalue_mask) = NaN;
            case 'top_voxels'
                % Find top X% of voxels
                performance = performance_metric_to_use(:, :, radii_index);
                pvalue = pvalue_to_use(:,:, radii_index);
                top_quantile_threshold = quantile(performance(:), 1-percent_of_voxels_to_keep);
                performance_mask = performance >= top_quantile_threshold;
                performance_masked = performance;
                performance_masked(~performance_mask) = NaN;
                fprintf('keeping top %d%% of voxels\n', percent_of_voxels_to_keep*100);

                % Report the p-value associated with the respective quantile thresholds
                pvalue_of_quantile_threshold = max(pvalue(performance == top_quantile_threshold),[], 'all');
                fprintf('pvalue_threshold - %d\n', pvalue_of_quantile_threshold);
            case 'bottom_voxels'
                % Find top X% of voxels
                performance = performance_metric_to_use(:, :, radii_index);
                pvalue = pvalue_to_use(:,:, radii_index);
                bottom_quantile_threshold = quantile(performance(:), percent_of_voxels_to_keep);
                performance_mask = performance <= bottom_quantile_threshold;
                performance_masked = performance;
                performance_masked(~performance_mask) = NaN;
                title_string = sprintf('keeping top %d%% of voxels\n', percent_of_voxels_to_keep*100);

                % Report the p-value associated with the respective quantile thresholds
                pvalue_of_quantile_threshold = max(pvalue(performance == bottom_quantile_threshold),[], 'all');
                fprintf('pvalue_threshold - %d\n', pvalue_of_quantile_threshold);
        end

        % Update colorbar limits
        colorbar_max = ceil(max([colorbar_max; performance_masked(:)], [], 'all'));
        colorbar_min = floor(min([colorbar_min; performance_masked(:)], [], 'all'));

        if debug
            fig_searchlight = figure;
            circle_center = [ 2 max(Z_img_mm)-2];
            radius = radii_to_test;

            brightened_image = angiogram./max(angiogram, [], 'all')*1.4;
            brightened_image(brightened_image>1) = 0.99;

            brightened_image = angiogram;

            cmap = plotDuplexImage(X_img_mm, Z_img_mm, performance_masked, brightened_image,...
                'colormap2use', colorscale, 'nonlinear_bg',2, ...
                'AutoColorBarLimits', [colorbar_min colorbar_max], ...
                'showColorbar', true, ...
                'ColorbarTitle', colorbar_title);
            title(sprintf('radius = %0.2f \nS%dR%d', radius, session_run_list(session_num, :)'));

            hold on;
            circle_handle = viscircles(circle_center,radius*pixelsize, 'EnhanceVisibility', false, ...
                'Color', 'w');
        end

        % Create mask from searchlight
        searchlight_map = performance_masked;
        searchlight_mask = ~isnan(searchlight_map);

        %% Check to see how much of image has significant GLM or SL voxels
        eliminate_sessions_with_low_activation = true;
        if eliminate_sessions_with_low_activation
            amt_significant_voxel_threshold = 0.0001;

            percent_voxels_glm_significant = sum(glm_mask(:))/numel(glm_mask);
            percent_voxels_sl_significant = sum(searchlight_mask(:))/numel(searchlight_mask);

            if percent_voxels_glm_significant < amt_significant_voxel_threshold || percent_voxels_sl_significant < amt_significant_voxel_threshold
                keep_session(session_num) = false;
                fprintf('Not enough significant voxels in GLM or searchlight analysis - excluding S%dR%d\n', session_run_list(session_num, :)');
            end
        end

        %% Compute Sorenson Dice coefficient
        % Colorblind palette
        color1 = [100 143 255]./255;
        color2 = [255 176 0]./255;
        color3 = [220 38 127]./255;

        % compute sorenson dice coefficient
        dice_similarity(session_num) = dice(searchlight_mask, glm_mask);

        % Create conjunction mask where each overlay is a different color
        conjunction_mask = nan(size(searchlight_mask));
        conjunction_mask(searchlight_mask) = 1;
        conjunction_mask(glm_mask) = 2;
        conjunction_mask(searchlight_mask & glm_mask) = 3;
        fig_conjunction_map = figure;


        plotDuplexImage(X_img_mm, Z_img_mm, conjunction_mask, angiogram,...
            'colormap2use', [color1; color2; color3], 'nonlinear_bg', 2, ...
            'showColorbar', false, ...
            'colorbarTitle', 'conjunction',...
            'darkestShade', 0.25, ...
            'AutoColorBarLimits', [1 3]);
        title(sprintf('S%dR%d \nSorenson Dice Coefficient: %0.10f', ...
            session_run_list(session_num, :)', dice_similarity(session_num)));

        % Create specified legend
        hold on;
        c1 = scatter(nan, nan, 'filled', 'MarkerFaceColor', color1);
        c2 = scatter(nan, nan, 'filled', 'MarkerFaceColor', color2);
        c3 = scatter(nan, nan, 'filled', 'MarkerFaceColor', color3);
        hold off
        legend([c1, c2, c3], {'PCA+LDA', 'GLM', 'Both'}, 'Location', 'southwest');

        save_to_file = false;
        if save_to_file && saveMe

            suggested_fname = sprintf('dice_sorenson_conjunction_map_multicoder_S%dR%d_gen%s.svg', session_run_list(session_num, :)', datetime('now', 'Format', 'yyyyMMdd'));
            full_filename = fullfile(get_user_data_path('pathType', 'output'), 'dice sorenson', suggested_fname);
            saveas(fig_conjunction_map, full_filename);
        end

        %% Correlated GLM and searchlight values
        % plot GLM vs searchlight voxel values

        % Plot f-score against searchlight metric
        glm_values = glm_F;
        sl_values = performance_metric_to_use(:, :, radii_index);

        % Using the GLM threshold for both for convenience of explanation
        glm_pvalue_ind = glm_pvalues<=pvalue_threshold; %GLM p-values
        sl_pvalue_ind = pvalue_to_use(:, :, radii_index)<=pvalue_threshold; %Searchlight pvalues

        plot_category = ones(size(glm_values)); % neither significant. Default.
        plot_category(glm_pvalue_ind & sl_pvalue_ind) = 4; % Both significant
        plot_category(glm_pvalue_ind & ~sl_pvalue_ind) = 3; % GLM significant, but not SL
        plot_category(~glm_pvalue_ind & sl_pvalue_ind) = 2; % SL significant, but not GLM

        colors_to_use = [216 27 96;
            30 136 229;
            255 193 7;
            0 77 64]./255;

        plot_colors = colors_to_use(plot_category, :);

        switch display_accuracy_or_angularError
            case 'angular_error'
                axis_label = 'Searchlight - angular error (deg)';
            case 'accuracy'
                axis_label = 'Searchlight - percent correct (%)';
        end

        plot_marker_size = 36;
        plot_this = true;
        if plot_this
            if session_num == 1
                fig_correlation = figure('Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
                tld_correlation = tiledlayout('flow');
                title(tld_correlation, sprintf('GLM vs searchlight nsig<=%d', pvalue_threshold));
            end
            % Calculate Pearson correlation
            mdl = fitlm(log(glm_values(:)), sl_values(:));

            nexttile(tld_correlation);
            hold on;
            possible_legend_entries = {'both - ns', 'GLM - ns, SL - s', 'GLM - s, SL - ns', 'both - s'};
            actual_legend_entries = {};
            % Plot each type of point
            for i = 1:4
                plot_ind = plot_category == i;
                if any(plot_ind, 'all')
                    actual_legend_entries{end+1} = possible_legend_entries{i};
                end
                scatter(log(glm_values(plot_ind)), sl_values(plot_ind),plot_marker_size, colors_to_use(i, :), 'filled');
            end
            hold off;

            % Custom limits for paper figure.
            ylim([40 110]);
            xlim([-3.5 5.5]);
            axis square

            % Put R^2 value on plot
            ylimits = ylim;
            xlimits = xlim;
            text(xlimits(1)+1, ylimits(1)+5, sprintf('R^2=%0.2f', mdl.Rsquared.Adjusted));

            % Add title and figure labels
            title(sprintf('S%dR%d', session_run_list(session_num, :)'));
            xlabel('ln(GLM F-score)');
            ylabel(axis_label);
            legend(actual_legend_entries);
        end
    catch
        warning('Unable to analyze session #%d - S%dR%d', ...
            session_num, ...
            session_run_list(session_num, 1),...
            session_run_list(session_num, 2));
    end
end

%% Compute and display summary statistics
median_color = 'r';
marker_color = 'k';
jitter_width = 0.1;

% Hacky way to get monkey names for each row in session_run_list
monkey_names = {'L', 'P'};
monkey_yaxis = [1, 2];
monkey_ind = NaN(num_sessions, 1);
for i = 1:num_sessions
    project_record_indx = ismember( ProjectRecord.Session,session_run_list(i,1)) & ismember(ProjectRecord.Run,session_run_list(i,2));
    monkey_name = ProjectRecord.Monkey{project_record_indx};
    monkey_ind(i) = find(contains(monkey_names, monkey_name));
end

% For each monkey compute statistics and plot
figure_summary = figure;
nan_ind = isnan(dice_similarity);
dice_similarity_noNaN = dice_similarity(~nan_ind & keep_session);
monkey_ind_noNaN = monkey_ind(~nan_ind & keep_session);

% If no monkey 1 or 2, then define a fake point to get it to plot properly.
% Then set ylim to exclude this fake point.
if ~any(monkey_ind_noNaN==1)
    dice_similarity_noNaN(end+1) = -0.1;
    monkey_ind_noNaN(end+1) = 1;
elseif ~any(monkey_ind_noNaN==2)
    dice_similarity_noNaN(end+1) = -0.1;
    monkey_ind_noNaN(end+1) = 2;
end
boxplot(dice_similarity_noNaN, monkey_ind_noNaN);
xlimits = xlim;
xlim_range = diff(xlimits);

for i = 1:2
    dice_similarity_cut = dice_similarity_noNaN(monkey_ind_noNaN == i);
    num_sessions_noNaN = length(dice_similarity_cut);

    quartile_values = prctile(dice_similarity_cut, [25, 50, 75]);
    similarity_mean = mean(dice_similarity_cut);
    similarity_ste = std(dice_similarity_cut)/sqrt(num_sessions_noNaN);
    num_excluded_sessions = sum(~keep_session & monkey_ind == i);

    fprintf('-----SUMMARY STATISTICS (%d sessions)------\n', num_sessions_noNaN);
    fprintf('Monkey - %s\n', monkey_names{i});
    fprintf('[1st quartile, median, 3rd quartile]: [%.02f, %.02f, %.02f]\n', quartile_values');
    fprintf('Mean +/- STE: %0.2f +/- %0.2f\n', similarity_mean, similarity_ste);
    fprintf('Excluded sessions: %d\n', num_excluded_sessions);
end

hold on;
swarmchart(monkey_ind(~nan_ind), dice_similarity(~nan_ind), marker_color, 'filled', 'XJitterWidth', jitter_width);
xticks([1, 2]);
xticklabels(monkey_names);
title('Similarity between searchlight and GLM SPMs');
ylabel('Dice-Sorenson Similarity');
xlabel('Monkey');
ylim([0 1]);
xlim([0.5 2.5]);

save_to_file = true;
if save_to_file && saveMe

    if any(monkey_ind==1) && any(monkey_ind==2)
        suggested_fname = sprintf('dice_sorenson_summaryFigure_allMonkeys_gen%s.svg', datetime('now', 'Format', 'yyyyMMdd'));
    elseif any(monkey_ind==1)
        suggested_fname = sprintf('dice_sorenson_summaryFigure_L_gen%s.svg', datetime('now', 'Format', 'yyyyMMdd'));
    elseif any(monkey_ind==2)
        suggested_fname = sprintf('dice_sorenson_summaryFigure_P_gen%s.svg', datetime('now', 'Format', 'yyyyMMdd'));
    end

    [filename, filepath] = uiputfile(fullfile(get_user_data_path('pathType', 'output'), 'dice sorenson', suggested_fname));
    if any([filename filepath])
        saveas(figure_summary, fullfile(filepath, filename));
    end
end



