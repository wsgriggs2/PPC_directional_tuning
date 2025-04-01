%% Main script for running analysis in directional tuning paper


%% Figure 2
% Generate angiograms
create_session_angiograms;

%% Figure 3 + Supplemental Figure 2A-D
% GLM analysis + ROI analysis for 2 example sessions
glm_analysis;

% Plot glm results
plot_glm_results;

%% Figure 4 + Supplemental Figure 1 E-H + Supplemental Figure 3 
% Single session decoding
% Searchlight (single timepoint)
BCI_CV_direction_multicoder_singleSession;

% Dice Sorenson similarity coefficient
% Correlation between GLM weight and SL accuracy for example sessions
calculate_sorenson_dice_coefficient_across_sessions

%% Figure 5 + Supplemental Figure 6 + Supplemental Figure 7
% Perform pairwise decoding across sessions
BCI_CV_direction_multicoder_pairwise_performance;

% Visualize results of pairwise decoding across sessions
% Supplemental Figure 8 - Similarity of images between decoding sessions
plot_pairwise_decoding_results;

%% Figure 6
% Statistical response across sessions
AnalyzeAcrossPlanes_8tgt

% Plot results
PlotAcrossSessionAnalysis

%% Figure 7
% Decoding across all sessions
BCI_CV_direction_multicoder_acrossSessions;

% Visualize results
plot_across_session_offline_decoding_results;

%% Supplemental Figure 1 and Supplemental Figure 9 and 10
% Pixelwise statistical measures for example sessions
% Pixelwise statistical measures for each coronal strength
Analyze_8tgt_withGLM

%% Supplemental Figure 4
% Searchlight results across decoding time windows
BCI_CV_direction_multicoder_singleSession_searchlight

%% Supplemental Figure 5
% Projection of LDA weights onto vascular images
plot_PCA_LDA_weights

%% Supplemental Figure 11
% Power Doppler data displays spatial autocorrelation
estimate_spatial_autocorrelation

%% Supplemental Figure 12
% # of components to capture 95% of variance in PCA for each example
% session
examine_pca