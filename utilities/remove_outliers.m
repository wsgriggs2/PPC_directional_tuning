function cleaned_data_4D = remove_outliers(data_4D, borderwidth, varargin)
%% Function to remove outliers from percentchange data. 
% Removes only outliers within the specified border region. Most of the
% artifacts seem to be a result of motion correction and edge effects
% likely resulting from interpolation.
%
% Inputs:
%   data_4D - Percentchange version of iDop. 
%   borderwidth - Width of border from which to eliminate outliers
%   varargin - Optional inputs
%       `detection_threshold` - The detection threshold factor for
%                               `isoutlier` function.
%       'dilation_amount' - How much extra around the border to check for
%                           outliers.


% Parse optional inputs
p = inputParser;
addOptional(p, 'detection_threshold', 3);
addOptional(p, 'dilation_amount', 4);
parse(p, varargin{:});
result = p.Results;


% Get dimensions of the data
[yPix, xPix, nTimepoints, nTrials] = size(data_4D);

% Extract std of the data
std_data_4D = squeeze(std(data_4D, 0, 4));
outlier_ind = isoutlier(std_data_4D, 2, ...
    'ThresholdFactor', result.detection_threshold) ...
    | isoutlier(std_data_4D, 1, ...
    'ThresholdFactor', result.detection_threshold);

% Compress across 3rd dimension (time) to eliminate any pixels that are
% outliers in the final image.
outlier_ind = any(outlier_ind, 3);

% Add extra voxels within the border region
outlier_ind([[1:borderwidth], [yPix-borderwidth:yPix]], :) = true;
outlier_ind(:,[[1:borderwidth], [xPix-borderwidth:xPix]]) = true;

% Create mask of the border of the image
dilated_borderwidth = borderwidth + result.dilation_amount;
border_mask = false(size(data_4D, 1), size(data_4D, 2));
border_mask([[1:dilated_borderwidth], [yPix-dilated_borderwidth:yPix]], :) = true;
border_mask(:,[[1:dilated_borderwidth], [xPix-dilated_borderwidth:xPix]]) = true;


% Combine the border mask and the outlier masks. This only eliminates the
% outliers that lie in the border range
combined_mask = outlier_ind & border_mask;
combined_mask = repmat(combined_mask, 1, 1, nTimepoints, nTrials);

cleaned_data_4D = data_4D;
cleaned_data_4D(combined_mask) = NaN;
