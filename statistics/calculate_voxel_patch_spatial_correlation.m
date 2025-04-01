function voxel_correlation = calculate_voxel_patch_spatial_correlation(...
    Dop, radius, varargin)

%% Return the spatial correlation with neighboring pixels
% This is 100% not the most efficient implementation of this. The most
% efficient implementation would store or cache the pairwise correlation
% between every voxel pair, then pull the relevant pairwise correlation
% when we need it. Here, we keep calculating the same correlations over and
% over and over. This cache/storing should actually be done externally to
% this function.

p = inputParser;
validDopMat = @(x) length(size(x)) == 3;
validPatchSize = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
addRequired(p, 'Dop',  validDopMat);
addRequired(p, 'radius', validPatchSize);
addOptional(p, 'verbose', true);
addOptional(p, 'shape', 'circle');
parse(p, Dop, radius, varargin{:});
result = p.Results;

Dop = result.Dop;
radius = result.radius;
verbose = result.verbose;
shape = result.shape;


% Start the parallel pool
p = gcp;


% Extract size dimensions from Dop1 and Dop2
[nDepth, nWidth, nTimepoints] = size(Dop);
if verbose
    % Create progress bar
    ppm = ProgressBar(nDepth*nWidth);
end
if radius > 0
    % Allocate space for the output images
    voxel_correlation = deal(NaN(nDepth, nWidth));
    
    [columnsInImage, rowsInImage] = meshgrid(1:nWidth, 1:nDepth);
    
    Dop_2D = reshape(Dop, [], size(Dop, 3));
    
    parfor j = 1:nDepth
        for i = 1:nWidth
            switch shape
                case 'circle'
                    ROI_pixels = (rowsInImage - j).^2 ...
                        + (columnsInImage - i).^2 <= radius.^2;


                    mask_ind = false(nDepth, nWidth);
                    mask_ind(ROI_pixels) = true;
                    mask_ind(i, j) = false; % Remove center pixel
                case 'ring'
                    ROI_pixels = (rowsInImage - j).^2 ...
                        + (columnsInImage - i).^2 <= radius.^2;


                    mask_ind = false(nDepth, nWidth);
                    mask_ind(ROI_pixels) = true;
                    
                    % Find perimeter voxels for the desired radius
                    % This almost works, but it keeps the edges of the
                    % image, and I could not think of any quick and elegant
                    % solution that would eliminate the edge pixels that 
                    % were less than the desired radius.
%                     perim_pixels = bwperim(mask_ind);
%                     perim_pixels(j, i) = false;
%                     mask_ind = perim_pixels;
                    
                    % Hack utnil I come up with a better solution for this
                    ROI_pixels_to_exclude = (rowsInImage - j).^2 ...
                        + (columnsInImage - i).^2 < (radius-1).^2;
                    mask_ind(ROI_pixels_to_exclude) = false;
            end
            
            mask_ind_1D = reshape(mask_ind, 1, []);
            
            
            center_pixel = squeeze(Dop(j, i, :));
            surrounding_pixels = squeeze(Dop_2D(mask_ind_1D, :));
            
            
            correlations = corr(center_pixel, surrounding_pixels');
            if any(isnan(correlations), 'all')
                warning('NaN value found. Ignoring it.');
            end
            voxel_correlation(j, i) = mean(correlations, 'omitnan');
            
            if verbose
                ppm.count();
            end
        end
    end
    
else
    voxel_correlation = ones(nDepth, nWidth);
end

if verbose
    delete(ppm);
end
end

