%% Script to plot the GLM results from `glm_analysis`
% Author - Whitney Griggs. October 24, 2022

%% Set up some initial values
pvalue_threshold = 1e-5;

%% Find all possible GLM files
start_dir = pwd;

% Get the possible filenames
cd(fullfile(get_user_data_path('pathType', 'glm')));
glm_files = dir('glm_data_*.mat');

glm_filenames = {glm_files.name};

%% Select which GLM files to load

% get Project for which this is a dataset
[indx,~] = listdlg('PromptString','Which GLM results do you want to plot?',...
    'ListString',glm_filenames,...
    'SelectionMode','multiple',...
    'Name','GLM results Selection', ...
    'ListSize', [600 300]);
glm_filename_reduced = glm_filenames(indx);


%% For each file, do the following
for file = 1:size(glm_filename_reduced, 2)
    %% Load data
    load(glm_filename_reduced{file}, 'F', 'p', 'pvalues_corrected', 'angiogram', 'sessionRunList');
    
    %% Threshold pvalues for Fcontrast and display
    F_thresholded = F;
    F_thresholded(pvalues_corrected>pvalue_threshold) = NaN;
    
    pixelsize = 0.1;
    X_img_mm = pixelsize/2 + (0:size(F,2)-1)*pixelsize; % So that the center of the image is at true halfway point, not shifted.
    Z_img_mm = pixelsize/2 + (0:size(F,1)-1)*pixelsize;
    %custom colormap for Positive only activation stats
    lateralizationMapComplete = load('LateralizationP.mat');
    PositiveOnlyColormap = lateralizationMapComplete.LateralizationP;
    PositiveOnlyColormap = PositiveOnlyColormap(1:floor(length(PositiveOnlyColormap(:,1))/2),:);
    
    figure;
    transparency_map = isnan(F_thresholded);
    
    plotDuplexImage(X_img_mm, Z_img_mm, F_thresholded, angiogram,...
        'colormap2use', flipud(magma), 'nonlinear_bg', 2, ...
        'showColorbar', true, ...
        'colorbarTitle', 'F-score',...
        'darkestShade', 0.25,...
        'transparencyMap', ~transparency_map);
    title(sprintf('[%s] \nqvalue<%0.10f', sprintf('S%dR%d ', sessionRunList'), pvalue_threshold));
end
% Change back to initial directory
cd(start_dir)




