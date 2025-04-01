function [dopOut, coreParams,template] = correctMotion(dopIn, coreParams, verbose, template)
% this function uses rigid body transforms to correct errant motion in the
% fUS dopppler sequence. It does this on a single-trial basis, i.e. the
% data passed to this function is of the size X x Y x nWindows x nTrials.
% nWindows is the number of fUS activation images in a single trial, e.g.
% ~20 images. 
% 
% inputs:
%
% dopIn: xPixels x yPixels x nWindowsInEpoch x nTrials
% coreParams: the coreparams structure used in createDoppler (and saved in
%             S##R##.mat doppler data files). note: coreParams should
%             contain coreParams.method of a valid choice: 'imregister',
%             'rigid', or 'normcorre'. Otherwise a default is chosen.
%
% output:
%
% dopOut: array xPixels x yPixels x nWindowsInEpoch x nTrials  (corrected)
% coreParams: will now contain the reference frame and trial 

%% check the registration method
if ~isfield(coreParams,'method')
    warning('defaulting to normcorre method!')
    method = 'normcorre';
else 
    method = coreParams.method;
end

if ~strcmp(method,'imregister') && ~strcmp(method,'rigid') && ~strcmp(method,'normcorre')    
end

if ~exist('verbose','var')
    verbose = false;
end

if ~exist('template','var')
    useTemplate = false;
else
    useTemplate = true;
end

%% select reference frame and window
[~, ~, nWindows, nTrials] = size(dopIn);
% refTrial: the trial number to pull the reference (stationary) frame from
% refFrame: the frame/window within the trial  
if strcmp(method,'imregister') || strcmp(method,'rigid') || strcmp(method, 'imregdeform')
    refTrial = ceil(nTrials/4);
    refFrame = 1;
    fixedFrame = zscore(dopIn(:,:,refFrame,refTrial), 0, 'all');
elseif strcmp(method,'normcorre')
    refTrial = 0; 
    refFrame = 0;
else
    error('invalid method chosen.')
end

fprintf('Reference: trial %i/%i, frame %i/%i\n', refTrial, size(dopIn,4), refFrame, size(dopIn,3))
% save out updated coreParams structure
coreParams.motionCorrection_refTrial = refTrial;
coreParams.motionCorrection_refFrame = refFrame;

%% choose the reference frame & initialize data
dopOut = zeros(size(dopIn));

%% imregister method
if strcmp(method,'imregister')
    
    % See comment below about interpolation method. 'Linear' is the default.
    interpolation_method = 'linear'; 
    
    [optimizer, metric] = imregconfig('monomodal');
    
    optimizer.MaximumIterations = 300;
    
    % Set up a matrix to monitor for failed registrations
    failed_registrations = false(size(dopOut, 3), size(dopOut, 4));
    
    % for each frame in the sequence
    ppm = ProgressBar(nTrials*nWindows, 'showWorkerProgress', true);
    parfor trial = 1:nTrials
        for window = 1:nWindows
            
            % select the current frame to work on
            movingFrame = dopIn(:,:,window,trial);
            
            % register the image
            % Although using imregtform followed by imwarp is slower
            % (instead of using imregister), it allows us to specify which
            % interpolation method is used. It was hard to figure out what
            % method was actually used by imregister, but I think it was
            % bilinear interpolation.
            lastwarn('') % Clear last warning message
            img_trans = imregtform(zscore(movingFrame, 0, 'all'), fixedFrame, 'rigid', optimizer, metric);
            dopOut(:,:,window,trial) = imwarp(movingFrame, img_trans, ...
              interpolation_method,...
              'FillValues', 0,...
              'OutputView', imref2d(size(fixedFrame)));
            %temp_frame = imregister(movingFrame, fixedFrame, 'rigid', optimizer, metric);
            %temp_frame(temp_frame == 0) = NaN;
            %dopOut(:,:,window,trial) = temp_frame;
            
            % Check if a warning was thrown
            [warnMsg, ~] = lastwarn;
            if ~isempty(warnMsg)
                failed_registrations(window, trial) = true;
                dopOut(:, :, window, trial) = movingFrame;
            end
            
            % update the waitbar
            ppm.count();
        end
    end
    delete(ppm);
    
    % Fix the failed registrations by taking the closest successful
    % transform
    
    failed_registrations_1D = failed_registrations(:);
    
    for timepoint = 1:length(failed_registrations_1D)
       if failed_registrations(timepoint)
          %Iteratively find closest good timepoint. This ignores that there
          %could be a fair amount of time between trials because we don't use
          %continuous Doppler processing, but rather pull out trials.
          
            % Find index of the misaligned frame
          [window_to_align, trial_to_align] = ind2sub(size(failed_registrations), timepoint);
          
          
          found_replacement_timepoint = false;
          counter = 0;
          while ~found_replacement_timepoint
              counter = counter + 1;
             if (timepoint-counter > 0) && ...
                     ~failed_registrations_1D(timepoint - counter)
                 found_replacement_timepoint = true;
                 % Convert back into trial and window coordinates
                 [window_to_use, trial_to_use] = ind2sub(size(failed_registrations), timepoint-counter);
             elseif (timepoint+counter <= length(failed_registrations_1D)) && ...
                     ~failed_registrations_1D(timepoint + counter)
                 found_replacement_timepoint = true;
                 % Convert back into trial and window coordinates
                 [window_to_use, trial_to_use] = ind2sub(size(failed_registrations), timepoint+counter);
             end
          end
          

          % Find the image transform for that new frame and apply to the
          % misaligned frame
          referenceFrame = dopIn(:, :, window_to_use, trial_to_use);
          movingFrame = dopIn(:, :, window_to_align, trial_to_align);
          img_trans = imregtform(referenceFrame, fixedFrame, 'rigid', optimizer, metric);
          dopOut(:, :, window_to_align, trial_to_align) = imwarp(movingFrame, img_trans, ...
              interpolation_method,...
              'FillValues', 0,...
              'OutputView', imref2d(size(fixedFrame)));
          
          
       end
    end
    
    
end
%% imregdeform method
if strcmp(method,'imregdeform')
    
    
    
    % Set up a matrix to monitor for failed registrations
    failed_registrations = false(size(dopOut, 3), size(dopOut, 4));
    
    % for each frame in the sequence
    ppm = ProgressBar(nTrials*nWindows, 'showWorkerProgress', true);
    parfor trial = 1:nTrials
        for window = 1:nWindows
            
            % select the current frame to work on
            movingFrame = dopIn(:,:,window,trial);
            
            % register the image
            % Although using imregtform followed by imwarp is slower
            % (instead of using imregister), it allows us to specify which
            % interpolation method is used. It was hard to figure out what
            % method was actually used by imregister, but I think it was
            % bilinear interpolation.
            lastwarn('') % Clear last warning message
            [~,img_trans] = imregdeform(zscore(movingFrame, 0, 'all'), fixedFrame, 'DisplayProgress', false);
            dopOut(:,:,window,trial) = img_trans;

            
            % Check if a warning was thrown
            [warnMsg, ~] = lastwarn;
            if ~isempty(warnMsg)
                failed_registrations(window, trial) = true;
                dopOut(:, :, window, trial) = movingFrame;
            end
            
            % update the waitbar
            ppm.count();
        end
    end
    delete(ppm);
    
    % Fix the failed registrations by taking the closest successful
    % transform
    
    failed_registrations_1D = failed_registrations(:);
    
    for timepoint = 1:length(failed_registrations_1D)
       if failed_registrations(timepoint)
          %Iteratively find closest good timepoint. This ignores that there
          %could be a fair amount of time between trials because we don't use
          %continuous Doppler processing, but rather pull out trials.
          
            % Find index of the misaligned frame
          [window_to_align, trial_to_align] = ind2sub(size(failed_registrations), timepoint);
          
          
          found_replacement_timepoint = false;
          counter = 0;
          while ~found_replacement_timepoint
              counter = counter + 1;
             if (timepoint-counter > 0) && ...
                     ~failed_registrations_1D(timepoint - counter)
                 found_replacement_timepoint = true;
                 % Convert back into trial and window coordinates
                 [window_to_use, trial_to_use] = ind2sub(size(failed_registrations), timepoint-counter);
             elseif (timepoint+counter <= length(failed_registrations_1D)) && ...
                     ~failed_registrations_1D(timepoint + counter)
                 found_replacement_timepoint = true;
                 % Convert back into trial and window coordinates
                 [window_to_use, trial_to_use] = ind2sub(size(failed_registrations), timepoint+counter);
             end
          end
          

          % Find the image transform for that new frame and apply to the
          % misaligned frame
          referenceFrame = dopIn(:, :, window_to_use, trial_to_use);
          movingFrame = dopIn(:, :, window_to_align, trial_to_align);
          
          displacement_field = imregdeform(zscore(referenceFrame, 0, 'all'), fixedFrame, 'DisplayProgress', false);
          dopOut(:, :, window_to_align, trial_to_align) = imwarp(movingFrame, displacement_field, ...
              interpolation_method,...
              'FillValues', 0,...
              'OutputView', imref2d(size(fixedFrame)));
          
          
       end
    end
    fprintf('Failed registration frames - %d\n', length(failed_registrations_1D));
    
end

%% normcorre method
if strcmp(method,'normcorre')
    disp('Correcting motion using normcorre')
    if useTemplate
        [dopOut,template] = normcorre_doppler(dopIn, verbose,template);
    else
        [dopOut,template] = normcorre_doppler(dopIn, verbose);
    end
end

end % function

