function [Dop_highpass, d] = highpass_doppler(Dop, fpass, fs, varargin)
%% Function to apply highpass filter on voxelwise basis



%% Variable inputs
p = inputParser;
p.addOptional('impulse_type', 'iir')
p.addOptional('steepness', 0.85);

p.parse(varargin{:});
inputs = p.Results;
%% Perform filtering

if ndims(Dop) == 3
    % Z by X by timepoints
    [n_depth, n_width, n_timepoints] = size(Dop);
    
    % flatten the data into a voxel-column for zscoring
    dataOut = permute(Dop, [2, 1, 3]);% [n_width, n_depth, n_timepoints);
    dataOut = reshape(dataOut, [n_depth*n_width, n_timepoints]);
    % perform the highpass filtering
    [dataOut, d] = highpass(dataOut', fpass, fs, ...
        'ImpulseResponse', inputs.impulse_type, ...
        'Steepness', inputs.steepness);
    dataOut = dataOut';
    % reshape back to original data dims
    dataOut = reshape(dataOut, [n_width, n_depth, n_timepoints]);
    Dop_highpass = permute(dataOut, [2, 1, 3]);
    
elseif ndims(Dop) == 4
    % Z by X by nWindows by nTrials
    [n_depth, n_width, n_windows, n_trials] = size(Dop);
    
    % flatten the data into a voxel-column for zscoring
    dataOut = permute(Dop, [2, 1, 4, 3]);% [n_width, n_depth, n_trials, n_windows]);
    dataOut = reshape(dataOut, [n_depth*n_width, n_windows*n_trials]);
    % perform the highpass filtering
    [dataOut, d] = highpass(dataOut', fpass, fs, ...
        'ImpulseResponse', inputs.impulse_type, ...
        'Steepness', inputs.steepness);
    dataOut = dataOut';
    % reshape back to original data dims
    dataOut = reshape(dataOut, [n_width, n_depth, n_trials, n_windows]);
    Dop_highpass = permute(dataOut, [2, 1, 4, 3]);
end
end