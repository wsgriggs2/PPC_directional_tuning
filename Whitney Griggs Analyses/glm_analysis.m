%% Script to perform GLM for single collection of sessions/runs
% Will run a single GLM analysis for whatever sessions/runs are selected,
% even if they are not from same anatomical plane.

% Author: Whitney Griggs
% Date: 2022

%% Initialize
clear; clc; close all;

%% Load colormaps to use
PositiveNegativeColormap = load('PercentChangeColormap.mat');
PositiveNegativeColormap = PositiveNegativeColormap.PercentChangeColormap;

project_record_filename = 'ProjectRecord_paper.json';

%% Load the Project Record and select individual sessions to analyze

sessions_runs_to_load = specify_sessions_of_interest('project_record_filename', project_record_filename);

%% Analyze the selected sessions

% Load in reverse order to keep larger field of view
sessions_runs_to_load = flipud(sessions_runs_to_load);

loadDopplerDataContinuous(whos, 'multiple', true, ...
    'SessionRunList', sessions_runs_to_load,...
    'mc_method', 'normcorre', ...
    'project_record_filename', 'ProjectRecord_paper.json');
[nDepth, nWidth, nTimepoints] = size(data.dop);

dop = data.dop;
behavior = data.behavior;
timestamps = data.timestamps(:, 1);
angiogram = data.angiogram;
%% Run the GLM model

% Create the desired HRF
n_secs = 16;
tau = 1; % Time constant of the gamma function
delta = 1; % A pure delay.
n = 3; % Phase delay of the gamma function
Fs = 1000; % Desired sampling rate resolution
hrf = generate_boynton_hrf(n_secs, 'tau', tau, 'delta', delta, 'n', n, 'Fs', Fs);

% Scaling method
scaling_method = 'voxelwise_mean';


FWHM_size = 0;
fprintf('Running GLM analysis...');
glm_data = run_glm_analysis(dop, behavior, timestamps, ...
    'verbose', false, ...
    'hrf', hrf, ...
    'scaling_method', scaling_method, ...
    'FWHM', FWHM_size);
fprintf('Finished \n');


%% Extract some contents of glm_data
F = glm_data.F;
p = glm_data.p_unc;
pvalues_corrected = glm_data.p_corrected;
analysis_parameters = glm_data.optional_inputs;

F_overall = glm_data.F_overall;
p_overall_unc = glm_data.p_overall_unc;
p_overall_corr = glm_data.p_overall_corr;

%% Save GLM features to file
% Sort in ascending order by session
sessionRunList = sortrows(sessionRunList, [1 2]);

% Logic for filenaming
filename = sprintf('glm_data_%sgen_%s.mat', sprintf('S%dR%d_', sessionRunList'), datetime('today', format='yyyyMMdd'));


savepath = fullfile(get_user_data_path('pathType', 'glm'), filename);

coreParams = data.coreParams;

save(savepath, 'glm_data', 'F', 'p', 'pvalues_corrected', 'sessionRunList',...
    'coreParams', 'analysis_parameters', 'angiogram', ...
    'F_overall', 'p_overall_unc', 'p_overall_corr', '-v7.3');


