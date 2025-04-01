function dopOut = flattenDoppler2D(dopIn, epochOfInterest, varargin)
% flattenDoppler2D is designed to reshape doppler data into arrays that are
% typically used in the BCI training algorithms.
%
% example use:
% dopOut = flattenDoppler2D(dopIn, *epochOfInterest=1:2*,
% 'verbose', false, 'repBoost', false, 'threeD', true)
% 
% vargout
% dopOut: flattened 2D doppler. If repBoost is on this is...
% (nWindowsUsed*nTrials x nPixels).
% If repboost is off (default) this is...
% (nTrials x nPixels*nWindowsUsed). 
% If threeD is turned on (by option only) this is...
% (nPixels x nWindows x nTrials)
% 
% vargin
% dopIn: 4D doppler data (yPix x xPix x nWindows x nTrials)
% can also pass 3D flattened data, e.g. (yPix*xPix) x nWindows x nTrials
%
% epochOfInterest: time indices to flatten over, e.g. [5 6 7]
%
%
% optional vargin
%
% see vargout for specifics
%
% 'threeD' (bool)
% see vargout for specifics
% 
% verbose: boolean


%% handling varargin
p = inputParser;
p.addOptional('verbose',false,@islogical)
p.addOptional('threeD',false,@islogical)
p.parse(varargin{:});
sets = p.Results;

if ~exist('epochOfInterest','var')
    epochOfInterest = 1:size(dopIn,3);
    if sets.verbose, warning('Epoch of interest not given. Passing all data'), end
end

if ndims(dopIn)==3
    dopIn = reshape(dopIn,[1 size(dopIn)]);
end

%% stack subsequent windows side by side: 4D->3D images of (yPix, xPix*nWindows, nTrials)
% save data dimension sizes for later use
[yPix, xPix, nWindows, nTrials] = size(dopIn); %#ok<ASGLU>
epochLength = length(epochOfInterest);

% if its just a single doppler measurement then squeeze that one
if epochLength==1 
    % this leaves a (yPix, xPix, nTrials) array
    dop3D = squeeze(dopIn(:,:,epochOfInterest,:));
    if size(dopIn,1)==1
        dop3D = reshape(dop3D,[1 size(dop3D)]);
    end
elseif sets.threeD
    % preserves the 4D structure, just selects the epoch
    dop4D = dopIn(:,:,epochOfInterest,:);    
else
    dopTmp = dopIn(:,:,epochOfInterest,:);       
    % this stacks the images in the epoch side by side
    dop3D = reshape(dopTmp, [yPix, epochLength*xPix, nTrials]);
    %dop3D = squeeze(mean(dopIn(:,:,epochOfInterest,:),3)); this just
    %averages across the epoch of interest (depracated)
end
clear dopTmp

%% plot the mean figures
if sets.verbose && ~sets.threeD && ndims(dopIn)==4
    map = mean(dop3D,3);
    %         map = map-min(map(:))/(max(map(:))-min(map(:)));
    map = (map-min(map(:)))/(max(map(:))-min(map(:)));
    map = nthroot(map,4);
    imagesc(map)
    title(['Mean across ' num2str(size(dop3D,3)) ' images'])
    
    % print size of images
    fprintf('%i images found. Image Size: %i x %i\n', ...
        size(dop3D,3), size(dop3D,1), size(dop3D,2))
    input('Press Enter to continue...')
    close all
end

%% reduce dimension
if sets.threeD
    % resize the images into 3D 
    nWindowsKept = size(dop4D,3);
    dopOut = reshape(dop4D, [yPix*xPix, nWindowsKept, nTrials]);    
else
    % this is a bit hacky at this point - putting extra spatial dim back
    if ndims(dopIn)==3
        dop3D = reshape(dop3D,[1 size(dop3D)]);
    end
    % resize the images into 2D
    nImages = size(dop3D,3);
    nPixels = size(dop3D,1)*size(dop3D,2);
    dopOut = zeros(nImages,nPixels); 
    % dopOut is size nImages x (xPixels*yPixels)
    for i = 1:size(dop3D,3)
        currentImage = dop3D(:,:,i);
        dopOut(i,:) = reshape(currentImage,1,nPixels);    
    end
end


end