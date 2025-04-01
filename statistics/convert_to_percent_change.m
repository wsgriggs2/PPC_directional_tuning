function [DopOut_percent_change] = convert_to_percent_change(DopIn, baseline_indices, varargin)

% Creates a new iDop structure that is percent change from the specified
% baseline for each pixel in DopIn.

% Inputs:
%   * DopIn: 4D Doppler data (depth x width x time x trials)
%   * baseline_indices: [start end]; The time window indices used to calculate 
%     baseline. Percent change is calculated relative to the mean activity across all trials during
%     this time window. This is not aligned to any particular time. Time
%     alignment should happen outside of this function and this function
%     should be passed the relevant iDop indices.
%   * varargin: None currently, but allows for future expansion.
% Outputs:
%   * DopOut_percent_change: 4D doppler data converted to percent change.

% Initially written by Whitney Griggs on November 19, 2020

% Input parser
p = inputParser;
validDopMat = @(x) length(size(x)) == 4;
addRequired(p, 'DopIn',  validDopMat);
addRequired(p, 'baseline_indices');
parse(p, DopIn, baseline_indices, varargin{:});
result = p.Results;

% Convert input parser output into simpler variable names
DopIn = result.DopIn;
baseline_indices = result.baseline_indices;

% Compute percent change
Dop_baseline = mean(DopIn(:, :, baseline_indices, :), 3);
Dop_baseline = repmat(Dop_baseline, 1, 1, size(DopIn, 3), 1);

% Handle case where Dop_baseline equals zero for certain voxels. This
% occurs due to motion correction. If there are parts of the image missing,
% it is filled with zeros. We should change these to be NaNs instead.

Dop_baseline(Dop_baseline == 0) = NaN;

DopOut_percent_change = (DopIn - Dop_baseline)./ Dop_baseline * 100;

end