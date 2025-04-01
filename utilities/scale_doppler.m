function [Dop_scaled, scaling_transform] = scale_doppler(Dop, varargin)
%% Function to scale Doppler data.
% Two possible methods
% 'Voxelwise_mean' - Used by AFNI
% This function uses the same method as AFNI for scaling of data to a
% common scale. See Chen, Taylor, and Cox 2017 under Units section.
% Scales to mean of the voxel.
% 'Grand_mean_scaling' - Used by SPM and FSL
% Normalize entire image to have mean of 1. 

% Input parser
p = inputParser;
addOptional(p, 'scaling_type', 'voxelwise_mean');
addOptional(p, 'baseline_index', []);
parse(p, varargin{:});
result = p.Results;

if ndims(Dop) == 3
    % Z by X by timepoints
    [n_depth, n_width, n_timepoints] = size(Dop);
    
    % flatten the data into a voxel-column
    dataOut = permute(Dop, [2, 1, 3]);% [n_width, n_depth, n_timepoints);
    dataOut = reshape(dataOut, [n_depth*n_width, n_timepoints]);
    % perform the scaling
    switch result.scaling_type
        case 'voxelwise_mean'
            if isempty(result.baseline_index)
                voxel_means = mean(dataOut, 2);
            else
                voxel_means = mean(dataOut(:, result.baseline_index), 2);
            end
            voxel_means_mat = repmat(voxel_means, 1, size(dataOut, 2));
            dataOut = 100*dataOut./voxel_means_mat;
            scaling_transform = voxel_means;
        case 'grand_mean_scaling'
            grand_mean = mean(dataOut, 'all', 'omitnan');
            dataOut = dataOut./grand_mean;
            scaling_transform = grand_mean;
        otherwise 
            error('Not an accepted method for scaling');
    end
    
    % reshape back to original data dims
    dataOut = reshape(dataOut, [n_width, n_depth, n_timepoints]);
    Dop_scaled = permute(dataOut, [2, 1, 3]);
    
elseif ndims(Dop) == 4
    % Z by X by nWindows by nTrials
    [n_depth, n_width, n_windows, n_trials] = size(Dop);
    
    % flatten the data into a voxel-column for zscoring
    dataOut = permute(Dop, [2, 1, 4, 3]);% [n_width, n_depth, n_trials, n_windows]);
    dataOut = reshape(dataOut, [n_depth*n_width, n_windows*n_trials]);
    % perform the scaling
    switch result.scaling_type
        case 'voxelwise_mean'
            if isempty(result.baseline_index)
                voxel_means = mean(dataOut, 2);
            else
                
                % If the baseline_index is 2D, then flatten
                % Assume that rows are trials and columns are timepoints in
                % a trial
                
                if ~isvector(result.baseline_index)
                    result.baseline_index = reshape(result.baseline_index, [], 1);
                end
                
                voxel_means = mean(dataOut(:, result.baseline_index));
            end
            voxel_means_mat = repmat(voxel_means, 1, size(dataOut, 2));
            dataOut = 100*dataOut./voxel_means_mat;
            scaling_transform = voxel_means;
        case 'grand_mean_scaling'
            grand_mean = mean(dataOut, [], 'all');
            dataOut = dataOut./grand_mean;
            scaling_transform = grand_mean;
        otherwise 
            error('Not a supported method for scaling');
    end
    voxel_means = mean(dataOut, 2);
    voxel_means_mat = repmat(voxel_means, 1, size(dataOut, 2));
    dataOut = (dataOut - voxel_means_mat)./voxel_means_mat;
    % reshape back to original data dims
    dataOut = reshape(dataOut, [n_width, n_depth, n_trials, n_windows]);
    Dop_scaled = permute(dataOut, [2, 1, 4, 3]);
end