function patch_centers = get_optimal_tiling(patch_size, image_dimension)
% Generate the optimal tiling in 1D. If mod(image_dimension, patch_size) is
% even, then leave equal space on both sides. If this remainder is odd,
% then leave an extra pixel spacing)relative to the bottom or right spacing
% at the top or left of the image.
% 
% Inputs:
%   patch_size - The scalar representing patch size. Must be an int.
%   image_dimension - The 1D image size. Either depth or width.
% Outputs:
%   patch_centers - The center of the patch.

% Get remainder size
unusable_space = mod(image_dimension, patch_size);
spannable_distance = image_dimension - unusable_space;

% Determine optimal spacing
if mod(unusable_space,2)
    % if remainder is odd
   patch_centers = 1 + (unusable_space - 1) / 2 + (patch_size + 1) / 2 : patch_size : spannable_distance + 1;
else
    % if remainder is even
   patch_centers = unusable_space / 2 + (patch_size + 1) / 2 : patch_size : spannable_distance;
   
end
end