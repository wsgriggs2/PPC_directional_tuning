%% Script to plot results of `BCI_CV_direction_multicoder_acrossSessions`

%% Load previously created decoding data
clear; close all; clc;

[loadfile, loadpath] = uigetfile(fullfile(get_user_data_path('pathType', 'decoding'), 'AcrossSessionResults_8dir_*_dateGenerated*.mat'));
load_fullfile = fullfile(loadpath, loadfile);

% Load the entire previous workspace
load(load_fullfile);


%% Loading in saved file and creating necessary variables

% This code block assumes that a specific file was loaded with specific
% sessions/runs in that file. Needs to be updated to handle new files
% Whitney Griggs 01/25/2025
% Still true as of 2025/03/13 - I know it works for old file with old
% Project Record, but need to test for new file.
project_record_filename = 'ProjectRecord_paper.json';
ProjectRecord = load_json_as_table(project_record_filename);


num_of_sessions = length(across_session_results);

bad_sessions = [];

% Remove RT data. It is twice as fast and acquired with different task
% parameters
RT_sessions = 39:49;
nonRT_sessions = 1:38;
across_session_results = across_session_results(nonRT_sessions);
nTrials_session = nTrials_session(nonRT_sessions);
nWindows_session = nWindows_session(nonRT_sessions);
task_timing = task_timing(nonRT_sessions, :);

num_of_sessions = length(across_session_results);


%% Extract monkey and categorize
% P = 1
% L = 2
monkey_names = {'P', 'L'};
monkey_num = NaN(num_of_sessions, 1);
for session_num = 1:num_of_sessions
    current_session_results = across_session_results{session_num};
    current_monkey = current_session_results{5};
    switch current_monkey
        case 'P'
            monkey_num(session_num) = 1;
        case 'L'
            monkey_num(session_num) = 2;
    end
end


%% Plotting

task_timing_copy = task_timing;
monkey_num_copy = monkey_num;

num_of_sessions = length(across_session_results);


% Preallocate
[across_session_percentCorrect, across_session_pvalue, across_session_t, ...
    across_session_angular_error, across_session_angular_error_pvalue] = deal(NaN(num_of_sessions, 10));

for session_num = 1:num_of_sessions
    current_session_results = across_session_results{session_num};
    combined_results = current_session_results{3};
    sessionnum_list(session_num) = unique(current_session_results{4}(:, 1));
    if ~isempty(bad_sessions) && any(ismember(session_num, bad_sessions(:, 1)))
        warning('not using this session');
    else
        n_timepoints = length(combined_results);
        for timepoint = 1:n_timepoints
            current_timepoint_results = combined_results{timepoint};
            across_session_percentCorrect(session_num, timepoint) = current_timepoint_results.percentCorrect;
            across_session_pvalue(session_num, timepoint) = current_timepoint_results.p;
            across_session_t(session_num, timepoint) = current_timepoint_results.t;
            across_session_angular_error(session_num, timepoint) = current_timepoint_results.angular_error;
            across_session_angular_error_pvalue(session_num, timepoint) = current_timepoint_results.angular_error_pvalue;
        end
    end
end

% Remove bad rows
if ~isempty(bad_sessions)
    across_session_percentCorrect(bad_sessions(:, 1), :) = [];
    across_session_pvalue(bad_sessions(:, 1), :) = [];
    across_session_t(bad_sessions(:, 1), :) = [];
    across_session_angular_error(bad_sessions(:, 1), :) = [];
    across_session_angular_error_pvalue(bad_sessions(:, 1), :) = [];

    task_timing_copy(bad_sessions(:, 1), :) = [];
    monkey_num_copy(bad_sessions(:, 1), :) = [];
    sessionnum_list(bad_sessions(:, 1)) = [];
end

monkey_names = {'P', 'L'};

[monkey_indx,~] = listdlg('PromptString','Select monkey(s) to analyze:',...
    'SelectionMode','multiple',...
    'ListString',monkey_names);
monkey_to_use = monkey_names(monkey_indx);
monkey_ind = ismember(monkey_num_copy, find(ismember(monkey_names, monkey_to_use)));
monkey_string = '';
for i = 1:length(monkey_to_use)
    monkey_string = [monkey_string ' ' monkey_to_use{i}];
end
use_session_ind = monkey_ind;


outlier_rows = [];
outlier_ind = true(length(use_session_ind), 1);
outlier_ind(outlier_rows) = false;
use_session_ind = use_session_ind & outlier_ind;

% Find last overlapping timepoint for plotting the time performance
[row, col] = find(across_session_percentCorrect(use_session_ind, :) == 0);
last_overlapping_timepoint = min(col)-1;
%last_overlapping_timepoint = size(across_session_percentCorrect, 2);

% Average across sessions
mean_accuracy_across_time = mean(across_session_percentCorrect(use_session_ind, :), 1);
ste_accuracy_across_time = std(across_session_percentCorrect(use_session_ind, :), 0, 1)/sqrt(nnz(use_session_ind));

mean_angular_error_across_time = mean(across_session_angular_error(use_session_ind, :), 1);
ste_angular_error_across_time = std(across_session_angular_error(use_session_ind, :), 0, 1)/sqrt(nnz(use_session_ind));

%% plot percent correct
sessions_to_use = sessionnum_list(use_session_ind);

%Retrieve plane for each session (but only retrieve first one)
project_record_ind = ismember(ProjectRecord{:, 'Session'}, sessions_to_use);
plane_nums = round(ProjectRecord{project_record_ind, 'ap_plane'});

% Base line color off of plane #
[unique_plane_numbers, ind] = unique(plane_nums);
num_unique_planes = length(unique_plane_numbers);
colors_to_use = distinguishable_colors(num_unique_planes);
colors_for_each_line = colors_to_use(plane_nums - min(plane_nums) + 1, :);

% Create legend entries for each color
legend_entries = cell(num_unique_planes, 1);
for i = 1:num_unique_planes
    legend_entries{i} = sprintf('Plane %d', unique_plane_numbers(i));
end


% Average task onsets
fix_start_average = mean(task_timing(:, 1));
cue_average = 0;
mem_end_average = mean(task_timing(:, 3));
mov_end_average = mean(task_timing(:, 4));

trialTime = round(fix_start_average:1:fix_start_average+length(mean_accuracy_across_time)-1);


% Plot result
figure;
subplot(1, 2, 1);
% Plotting each session
colororder(colors_for_each_line);
plt_handle = plot(trialTime(1:last_overlapping_timepoint), ...
    across_session_percentCorrect(use_session_ind, 1:last_overlapping_timepoint)');

% Hack to add transparency to each line
% Tried to add to the colororder itself, but only accepts mx3 matrices
transparency_level = 0.3;
for i = 1:length(plt_handle)
    plt_handle(i).Color(4) = transparency_level;
end

% Plotting average performance
shadedErrorBar(trialTime(1:last_overlapping_timepoint), ...
    mean_accuracy_across_time(1:last_overlapping_timepoint), ...
    ste_accuracy_across_time(1:last_overlapping_timepoint));

% handle x & y axes
axis([trialTime(1) trialTime(last_overlapping_timepoint) 0 100])
set(gca,'Visible','off')
axes('Position',get(gca,'Position'),...
    'XAxisLocation','bottom',...
    'YAxisLocation','left',...
    'Color','none',...
    'XColor','k','YColor','k',...
    'LineWidth',2,...
    'TickDir','out');
axis([trialTime(1) trialTime(last_overlapping_timepoint) 0 100])


% add labels
yticks(0:20:100)

xticks([0 mem_end_average])
xticklabels({'Cue','Go'})
ylabel('Accuracy (%)')

% add reference lines
hold on; plot([trialTime(1) trialTime(last_overlapping_timepoint)],[100/8 100/8],'--k','LineWidth',2)
text(trialTime(last_overlapping_timepoint), 100/8 + 5, 'chance', 'HorizontalAlignment', 'right', 'FontSize', 15);

% shade fixation & memory period
hFix = patch([trialTime(1) trialTime(1) 0 0], [0 100 100 0], ...
    'k', 'facecolor','b','edgecolor','none');
hFix.FaceAlpha = 0.2;
hMem = patch([0 0 mem_end_average mem_end_average], [0 100 100 0], ...
    'k', 'facecolor','g','edgecolor','none');
hMem.FaceAlpha = 0.2;

% add scale bar
plot([fix_start_average+1 fix_start_average+2], [30 30], 'k', 'LineWidth',5)
text(fix_start_average+1, 35, '1 sec', 'HorizontalAlignment', 'left', 'FontSize', 15)

% Add title
title(sprintf('Percent correct\n%s', monkey_string));

hold off;
legend(plt_handle(ind), legend_entries);


% plot angular error
subplot(1, 2, 2)
colororder(colors_for_each_line);
plt_handle = plot(trialTime(1:last_overlapping_timepoint), ...
    180/pi*across_session_angular_error(use_session_ind, 1:last_overlapping_timepoint)');
% Hack to add transparency to each line
% Tried to add to the colororder itself, but only accepts mx3 matrices
transparency_level = 0.3;
for i = 1:length(plt_handle)
    plt_handle(i).Color(4) = transparency_level;
end

shadedErrorBar(trialTime(1:last_overlapping_timepoint), ...
    180/pi*mean_angular_error_across_time(1:last_overlapping_timepoint), ...
    180/pi*ste_angular_error_across_time(1:last_overlapping_timepoint));


% handle x & y axes
axis([trialTime(1) trialTime(last_overlapping_timepoint) 0 100])
set(gca,'Visible','off')
axes('Position',get(gca,'Position'),...
    'XAxisLocation','bottom',...
    'YAxisLocation','left',...
    'Color','none',...
    'XColor','k','YColor','k',...
    'LineWidth',2,...
    'TickDir','out');
axis([trialTime(1) trialTime(last_overlapping_timepoint) 0 100])


% add labels
yticks(0:30:100)

xticks([0 mem_end_average])
xticklabels({'Cue','Go'})
ylabel('MAAE (deg)')

% add reference lines
hold on; plot([trialTime(1) trialTime(last_overlapping_timepoint)],[90 90],'--k','LineWidth',2)
text(trialTime(last_overlapping_timepoint), 90 + 5, 'chance', 'HorizontalAlignment', 'right', 'FontSize', 15);

% shade fixation & memory period
hFix = patch([trialTime(1) trialTime(1) 0 0], [0 100 100 0], ...
    'k', 'facecolor','b','edgecolor','none');
hFix.FaceAlpha = 0.2;
hMem = patch([0 0 mem_end_average mem_end_average], [0 100 100 0], ...
    'k', 'facecolor','g','edgecolor','none');
hMem.FaceAlpha = 0.2;

% add scale bar
plot([fix_start_average+1 fix_start_average+2], [30 30], 'k', 'LineWidth',5)
text(fix_start_average+1, 35, '1 sec', 'HorizontalAlignment', 'left', 'FontSize', 15)

% Add title
title(sprintf('Angular error\n%s', monkey_string));
hold off;

legend(plt_handle(ind), legend_entries, 'Location', 'best');


%% Look at relationship to session number to see if there is any temporal patttern across days/months

late_memory_indices = [floor(-1*fix_start_average+mem_end_average/2), ceil(-1*fix_start_average+mem_end_average)];
mov_indices = [floor(-1*fix_start_average+mem_end_average), floor(-1*fix_start_average+mem_end_average)+5];

peak_performance_per_session = max(across_session_percentCorrect,[], 2);
mean_performance_lateMemory_per_session = mean(across_session_percentCorrect(:, late_memory_indices), 2);
mean_performance_Movement_per_session = mean(across_session_percentCorrect(:, mov_indices), 2);

x = sessionnum_list(use_session_ind);
y{1} = peak_performance_per_session(use_session_ind)/100;
y{2} = mean_performance_lateMemory_per_session(use_session_ind)/100;
y{3} = mean_performance_Movement_per_session(use_session_ind)/100;

y_labels = {'Peak', 'Mean-lateMemory', 'Mean-earlyMov'};


figure;
tld = tiledlayout(1, 3);
for i = 1:size(sessions_to_use, 1)
    nexttile();

    plot(x, y{i}, '.', 'MarkerSize', 20);
    title(sprintf('%s decoder performance - \n%s', y_labels{i}, monkey_string));
    xlabel('Session Number');
    ylabel('Decoder performance (% correct decode)');
    ylimits = ylim;
    xlimits = xlim;
    ylim([0.1 ylimits(2)]);
    hold on;
    plot([xlimits(1) xlimits(2)]', [1/8 1/8], 'r--');
    hold off;
    xlim(xlimits);

    % Add linear fit to plot
    c = polyfit(x, y{i}, 1);
    mdl = fitlm(x,y{i})

    % Evaluate fit equation using polyval
    y_est = polyval(c,x);
    % Add trend line to plot
    hold on
    plot(x,y_est,'k','LineWidth',2)
    hold off

    legend('Data', 'Chance decode', sprintf('Linear fit - y=%0.2fx + %0.2f \n (r^2 - %0.2f)', c, mdl.Rsquared.Adjusted), ...
        'Location', 'best');
end

%% Look at relationship to what plane it came from
sessions_to_use = sessionnum_list(use_session_ind);

%Retrieve plane for each session
project_record_ind = ismember(ProjectRecord{:, 'Session'}, sessions_to_use);
plane_nums = round(ProjectRecord{project_record_ind, 'ap_plane'});

figh = figure;
tld = tiledlayout(1, 3);
for i = 1:length(y)
    figure(figh);
    nexttile();

    x = plane_nums;

    plot(x, y{i}, '.', 'MarkerSize', 20);
    title(sprintf('%s decoder performance - \n%s', y_labels{i}, monkey_string));
    xlabel('Plane');
    ylabel('Decoder performance (% correct decode)');
    ylimits = ylim;
    xlimits = xlim;
    ylim([0.1 ylimits(2)]);
    hold on;
    plot([xlimits(1) xlimits(2)]', [1/8 1/8], 'r--');
    hold off;
    xlim(xlimits);

    % Find average for each plane
    % Just use peak performance

    unique_planes = unique(plane_nums);
    num_unique_planes = length(unique_planes);
    [mean_performance, sem_performance] = deal(NaN(num_unique_planes, 1));
    for plane_num = 1:num_unique_planes
        peak_performance = y{i};
        mean_performance(plane_num) = mean(peak_performance(plane_nums==unique_planes(plane_num)));
        sem_performance(plane_num) = std(peak_performance(plane_nums==unique_planes(plane_num)))/sqrt(nnz(plane_nums==unique_planes(plane_num)));
    end


    % Plot the mean +/- SEM
    hold on;
    errorbar(unique_planes, mean_performance, sem_performance, '.');

    hold off;
    legend('Data', 'Chance decode', 'Mean +/- SEM', ...
        'Location', 'best');

    % 1-way ANOVA to see if difference between planes
    [p, tbl, stats] = anova1(y{i}, plane_nums);
    multcompare(stats);
end

%% How many sessions were significant?
pvalue_threshold = 0.01;
across_session_pvalue_withNaN = across_session_pvalue;
across_session_pvalue_withNaN(across_session_percentCorrect == 0) = NaN;

across_session_angular_pvalue_withNaN = across_session_angular_error_pvalue;
across_session_angular_pvalue_withNaN(across_session_angular_error == 0) = NaN;

best_pvalue_per_session = min(across_session_angular_pvalue_withNaN, [], 2);
num_significant_sessions = nnz(best_pvalue_per_session(use_session_ind)<pvalue_threshold);
num_total_sessions = nnz(use_session_ind);
fprintf('%d/%d significant sessions\n', num_significant_sessions, num_total_sessions);


%% More specific analysis for specific paper figure
% Get best angular error for each session
across_session_angular_error_withNaN = across_session_angular_error;
across_session_angular_error_withNaN(across_session_angular_error_withNaN == 0) = NaN;
[best_performance_per_session, time_ind] = min(across_session_angular_error_withNaN,[], 2, 'omitnan');
time_ind_linear = sub2ind(size(across_session_angular_pvalue_withNaN), [1:length(across_session_angular_pvalue_withNaN)]', time_ind);


sessions_to_use = sessionnum_list(use_session_ind);
selected_best_performance = best_performance_per_session(use_session_ind) * 180/pi;
pvalues_for_selected_times = across_session_angular_pvalue_withNaN(time_ind_linear);
pvalues_for_selected_times_and_sessions = pvalues_for_selected_times(use_session_ind);
sig_or_not_ind = pvalues_for_selected_times_and_sessions<pvalue_threshold;

%Retrieve plane for each session (but only retrieve first one)
project_record_ind = ismember(ProjectRecord{:, 'Session'}, sessions_to_use);
plane_nums = round(ProjectRecord{project_record_ind, 'ap_plane'});


% Set up options for significant/not significant plot
two_colors = distinguishable_colors(2);
two_labels = [sprintf('p>%0.2f', pvalue_threshold), sprintf('p<=%0.2f', pvalue_threshold)];


figh = figure;
hold on;
if any(~sig_or_not_ind)
    for i = 0:1
        what_to_plot_ind = sig_or_not_ind == i;
        x = plane_nums;
        plot(plane_nums(what_to_plot_ind), selected_best_performance(what_to_plot_ind), ...
            '.', 'MarkerSize', 20, ...
            'Color', two_colors(i+1, :));
    end
    legend(sprintf('p>%0.2f', pvalue_threshold),...
        sprintf('p\\leq%0.2f', pvalue_threshold), ...
        'Chance decode', 'Mean +/- SEM', ...
        'Location', 'best');
else
    x = plane_nums;
    plot(plane_nums, selected_best_performance, ...
        '.', 'MarkerSize', 20);
    legend(sprintf('p\\leq%0.2f', pvalue_threshold), ...
        'Chance decode', 'Mean +/- SEM', ...
        'Location', 'best');
end
title(sprintf('Peak decoder performance - \n%s', monkey_string));
xlabel('Plane');
ylabel('Decoder performance (MAAE (deg))');
ylimits = ylim;
xlimits = xlim;
plot([xlimits(1) xlimits(2)]', [90 90], 'k--');
hold off;
xlim(xlimits);



% Find average for each plane
% Just use peak performance

unique_planes = unique(plane_nums);
num_unique_planes = length(unique_planes);
[mean_performance, sem_performance] = deal(NaN(num_unique_planes, 1));
for plane_num = 1:num_unique_planes
    peak_performance = selected_best_performance;
    mean_performance(plane_num) = mean(peak_performance(plane_nums==unique_planes(plane_num)));
    sem_performance(plane_num) = std(peak_performance(plane_nums==unique_planes(plane_num)))/sqrt(nnz(plane_nums==unique_planes(plane_num)));
end


% Plot the mean +/- SEM

if any(~sig_or_not_ind)
    legend(sprintf('p>%0.2f', pvalue_threshold),...
        sprintf('p\\leq%0.2f', pvalue_threshold), ...
        'Chance decode', 'Mean +/- SEM', ...
        'Location', 'best');
else
    legend(sprintf('p\\leq%0.2f', pvalue_threshold), ...
        'Chance decode', 'Mean +/- SEM', ...
        'Location', 'best');
end


% 1-way ANOVA to see if difference between planes
[p, tbl, stats] = anova1(selected_best_performance, plane_nums);
multcompare(stats);

%% More specific analysis for specific paper figure
% Get best percent correct for each session
across_session_percentCorrect_withNaN = across_session_percentCorrect;
across_session_percentCorrect_withNaN(across_session_percentCorrect_withNaN == 0) = NaN;
[best_performance_per_session, time_ind] = max(across_session_percentCorrect_withNaN,[], 2, 'omitnan');
time_ind_linear = sub2ind(size(across_session_pvalue_withNaN), [1:length(across_session_pvalue_withNaN)]', time_ind);


sessions_to_use = sessionnum_list(use_session_ind);
selected_best_performance = best_performance_per_session(use_session_ind);
pvalues_for_selected_times = across_session_pvalue_withNaN(time_ind_linear);
pvalues_for_selected_times_and_sessions = pvalues_for_selected_times(use_session_ind);
sig_or_not_ind = pvalues_for_selected_times_and_sessions<pvalue_threshold;

%Retrieve plane for each session (but only retrieve first one)
project_record_ind = ismember(ProjectRecord{:, 'Session'}, sessions_to_use);
plane_nums = round(ProjectRecord{project_record_ind, 'ap_plane'});


% Set up options for significant/not significant plot
two_colors = distinguishable_colors(2);
two_labels = [sprintf('p>%0.2f', pvalue_threshold), sprintf('p<=%0.2f', pvalue_threshold)];


figh = figure;
hold on;
if any(~sig_or_not_ind)
    for i = 0:1
        what_to_plot_ind = sig_or_not_ind == i;
        x = plane_nums;
        plot(plane_nums(what_to_plot_ind), selected_best_performance(what_to_plot_ind), ...
            '.', 'MarkerSize', 20, ...
            'Color', two_colors(i+1, :));
    end
    legend(sprintf('p>%0.2f', pvalue_threshold),...
        sprintf('p\\leq%0.2f', pvalue_threshold), ...
        'Chance decode', 'Mean +/- SEM', ...
        'Location', 'best');
else
    x = plane_nums;
    plot(plane_nums, selected_best_performance, ...
        '.', 'MarkerSize', 20);
    legend(sprintf('p\\leq%0.2f', pvalue_threshold), ...
        'Chance decode', 'Mean +/- SEM', ...
        'Location', 'best');
end
title(sprintf('Peak decoder performance - \n%s', monkey_string));
xlabel('Plane');
ylabel('Decoder performance (percent correct)');
ylim([0 100]);
xlimits = xlim;
plot([xlimits(1) xlimits(2)]', 100*[1/8 1/8], 'k--');
hold off;
xlim(xlimits);



% Find average for each plane
% Just use peak performance

unique_planes = unique(plane_nums);
num_unique_planes = length(unique_planes);
[mean_performance, sem_performance] = deal(NaN(num_unique_planes, 1));
for plane_num = 1:num_unique_planes
    peak_performance = selected_best_performance;
    mean_performance(plane_num) = mean(peak_performance(plane_nums==unique_planes(plane_num)));
    sem_performance(plane_num) = std(peak_performance(plane_nums==unique_planes(plane_num)))/sqrt(nnz(plane_nums==unique_planes(plane_num)));
end


% Plot the mean +/- SEM
if any(~sig_or_not_ind)
    legend(sprintf('p>%0.2f', pvalue_threshold),...
        sprintf('p\\leq%0.2f', pvalue_threshold), ...
        'Chance decode', 'Mean +/- SEM', ...
        'Location', 'best');
else
    legend(sprintf('p\\leq%0.2f', pvalue_threshold), ...
        'Chance decode', 'Mean +/- SEM', ...
        'Location', 'best');
end


% 1-way ANOVA to see if difference between planes
[p, tbl, stats] = anova1(selected_best_performance, plane_nums);
multcompare(stats);