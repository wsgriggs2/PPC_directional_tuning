function effect_size = calculate_voxel_patch_effect_size(...
    Dop1, Dop2, patch_size, varargin)

% Returns the p-value and t-score for each non-overlapping voxel patch
% comparing between Dop1 and Dop2.
% 
% Inputs:
%   * Dop1: 4D Doppler structure; 
%   * Dop2: 4D Doppler structure; must have same depth and width as Dop1
%     patch_size: Size of the patch to be used to calculate the statistics
% Optional inputs:
%   * 'effect_size_measure' - Which effect size measure to calculate. Cohen
%                             d, percent change, [insert future options]
% Outputs:
%   * effect_size
% 

p = inputParser;
validDopMat = @(x) length(size(x)) == 4;
validPatchSize = @(x) isnumeric(x) && isscalar(x) && (x > 0);
addRequired(p, 'Dop1',  validDopMat);
addRequired(p, 'Dop2', validDopMat);
addRequired(p, 'patch_size', validPatchSize);
addOptional(p, 'verbose', true);
addOptional(p, 'effect_size_measure', 'Cohen d');
parse(p, Dop1, Dop2, patch_size, varargin{:});
result = p.Results;

Dop1 = result.Dop1;
Dop2 = result.Dop2;
patch_size = result.patch_size;

% Extract size dimensions from Dop1 and Dop2
[nDepth_Dop1, nWidth_Dop1, nWindows_Dop1, nTrials_Dop1] = size(Dop1);
[nDepth_Dop2, nWidth_Dop2, nWindows_Dop2, nTrials_Dop2] = size(Dop2);

% Check to see if Dop1 is same depth and width as Dop2
if nDepth_Dop1 == nDepth_Dop2 && nWidth_Dop1 == nWidth_Dop2
    nDepth = nDepth_Dop1;
    nWidth = nWidth_Dop1;
else
    error('The depth and width of Dop1 and Dop2 are not the same. These cannot be compared');
end

% Allocate space for the output images
effect_size = NaN(nDepth, nWidth);

% Determine tiling of the patches
if mod(nDepth, patch_size) == 0 && mod(nWidth, patch_size) == 0
    Zmat = (patch_size + 1) / 2 : patch_size : nDepth;
    Xmat = (patch_size + 1) / 2 : patch_size : nWidth;
elseif mod(nDepth, patch_size) == 0
    Zmat = (patch_size + 1) / 2 : patch_size : nDepth;
    Xmat = get_optimal_tiling(patch_size, nWidth);
elseif mod(nWidth, patch_size) == 0
    Zmat = get_optimal_tiling(patch_size, nDepth);
    Xmat = (patch_size + 1) / 2 : patch_size : nWidth;
else
    Zmat = get_optimal_tiling(patch_size, nDepth);
    Xmat = get_optimal_tiling(patch_size, nWidth);
end

if result.verbose
    h = waitbar(0,'Please wait ... computing stat mask');
end

for Z = Zmat
    if result.verbose
        waitbar(Z / length(Zmat))
    end
    for X = Xmat
        iWidth = X - (patch_size - 1)/2 : X + (patch_size - 1)/2;
        iDepth = Z - (patch_size - 1)/2 : Z + (patch_size - 1)/2;
        
        patch1 = squeeze(mean(mean(Dop1(iDepth, iWidth, :, :), 1), 2));
        patch2 = squeeze(mean(mean(Dop2(iDepth, iWidth, :, :), 1), 2));
        
        mean_patch1_activity = mean(patch1, 1);
        mean_patch2_activity = mean(patch2, 1);
        
        switch result.effect_size_measure
            case 'Cohen d'
                % Calculate Cohen's d (measure of effect size) for each voxel patch
                effect_size(iDepth, iWidth) = calculate_cohen_d(mean_patch1_activity, mean_patch2_activity);
            case 'percent change'
                effect_size(iDepth, iWidth) = calculate_mean_percentchange(mean_patch1_activity, mean_patch2_activity);
            case 'median difference'
                effect_size(iDepth, iWidth) = median(mean_patch1_activity) - median(mean_patch2_activity);
        end
    end
end

if result.verbose
    close(h)
    fprintf('Computed new statistics. patch_size=%d \n', patch_size);
end
end

