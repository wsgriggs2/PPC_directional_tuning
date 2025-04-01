function dop  = detrendDrift(dop, baselineWindow)
% dopOut = detrendDrift(dopIn, baselineWindow)
%
% detrendDrift taks in dopIn, which is a 4-D array of dimensions
% yPixels x xPixels x nWindows x nTrials
%
% This function computes baseline within the window defined by argin
% baselineWindow, e.g. baselineWindow = 1:3, and computes a linear
% trendline over the course of the run length (time). It then subtracts
% that linear fit from actual data and passes that abck as dopOut

%% 
[yPixels, xPixels, nWindows, nTrials] = size(dop);

% restrict data to the baseline period only
baselineData = squeeze(mean(dop(:,:,baselineWindow,:),3));

%% get the r squared and percent change maps, voxel by voxel (linear fits)
% [rsq, percChange] = deal(zeros(yPixels,xPixels));

% thePool = gcp;
H = ProgressBar(yPixels, task='Measuring drift. Please wait: ');
for y = 1:yPixels   
    for x = 1:xPixels        
        
        % restrict data to the current voxel only
        voxelData = squeeze(baselineData(y,x,:));       
        
        % compute the line fit
        coeffs = polyfit((1:nTrials)', voxelData, 1);        
        yFit = polyval(coeffs, 1:nTrials);
%         percChange(y,x) = (yFit(end)-yFit(1))/yFit(1)*100;        
        
        % get the r-squared for the fit at this voxel
%         yresid = voxelData'-yFit;    % compute the residual values
%         SSresid = sum(yresid.^2);   % residual sum of squares 
%         % Total sum of squares = variance multiplied by number of observations minus 1:
%         SStotal = (nTrials-1) * var(voxelData);
%         rsq(y,x) = 1 - SSresid/SStotal;  % compute R2 
        
        % compute dopOut
        voxelOut = voxelData-yFit'+yFit(1);               
        dop(y,x,:,:) = repmat(voxelOut,[ 1 nWindows ])';
        
        % Plot the fitted line (debugging only)
%         clf; hold on;
%         plot(voxelData)
%         plot(1:nTrials, yFit, 'r-', 'LineWidth', 3);
%         plot(1:nTrials, voxelOut(1,:), 'g', 'LineWidth', 2);
%         xlabel('successful trial'); ylabel('baseline power'); pause(0.01);
         
    end
    H.count();
end