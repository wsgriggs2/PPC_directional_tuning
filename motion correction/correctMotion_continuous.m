function [dopOut, coreParams,template] = correctMotion_continuous(dopIn, coreParams, verbose, template)
% this function uses rigid body transforms to correct errant motion in the
% fUS dopppler sequence. It does this on a single-trial basis, i.e. the
% data passed to this function is of the size X x Y x nTimepoints.

%
% inputs:
%
% dopIn: xPixels x yPixels x nTimepoints
% coreParams: the coreparams structure used in createDoppler (and saved in
%             S##R##.mat doppler data files). note: coreParams should
%             contain coreParams.method of a valid choice: 'imregister',
%             'rigid', or 'normcorre'. Otherwise a default is chosen.
%
% output:
%
% dopOut: array xPixels x yPixels x nTimepoints  (corrected)
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
[~, ~, nTimepoints] = size(dopIn);
% refTrial: the trial number to pull the reference (stationary) frame from
% refFrame: the frame/window within the trial
if strcmp(method,'imregister') || strcmp(method,'rigid') || strcmp(method, 'imregdeform')
    refFrame = ceil(nTimepoints/4);
    fixedFrame = zscore(dopIn(:,:,refFrame), 0, 'all');
elseif strcmp(method,'normcorre')
    refTrial = 0;
    refFrame = 0;
else
    error('invalid method chosen.')
end

fprintf('Reference: timepoint %i/%i\n', refFrame, size(dopIn,3))
% save out updated coreParams structure
coreParams.motionCorrection_refFrame = refFrame;

%% choose the reference frame & initialize data
dopOut = zeros(size(dopIn));

%% imregister method
if strcmp(method,'imregister')
        
    [optimizer, metric] = imregconfig('monomodal');
    
    %     optimizer.InitialRadius = 0.009; % only for use in multimodal
    %     optimizer.Epsilon = 1.5e-4; % only for use in multimodal
    %     optimizer.GrowthFactor = 1.01; % only for use in multimodal
    optimizer.MaximumIterations = 300;
    
    % Set up a matrix to monitor for failed registrations
    failed_registrations = false(size(dopOut, 3));
    
    % for each frame in the sequence
    ppm = ProgressBar(nTimepoints, 'showWorkerProgress', true);
    parfor timepoint = 1:nTimepoints
        
        % select the current frame to work on
        movingFrame = dopIn(:,:,timepoint);
        
        lastwarn('') % Clear last warning message
        
        % Attempt image registration.
        temp_frame = imregister(movingFrame, fixedFrame, 'rigid', optimizer, metric);
        %temp_frame(temp_frame == 0) = NaN;
        dopOut(:,:,timepoint) = temp_frame;
        
        % Check if a warning was thrown
        [warnMsg, ~] = lastwarn;
        if ~isempty(warnMsg)
            failed_registrations(timepoint) = true;
            dopOut(:, :, timepoint) = movingFrame;
        end
        
        % update the waitbar
        ppm.count();
    end
    delete(ppm);
    
    % Fix the failed registrations by taking the closest successful
    % transform
    
    failed_registrations_1D = failed_registrations(:);
    
    for timepoint = 1:length(failed_registrations_1D)
        if failed_registrations(timepoint)
            %Iteratively find closest good timepoint.
            
            found_replacement_timepoint = false;
            counter = 0;
            while ~found_replacement_timepoint
                counter = counter + 1;
                if (timepoint-counter > 0) && ...
                        ~failed_registrations_1D(timepoint - counter)
                    found_replacement_timepoint = true;
                    counter = -counter;
                elseif (timepoint+counter <= length(failed_registrations_1D)) && ...
                        ~failed_registrations_1D(timepoint + counter)
                    found_replacement_timepoint = true;
                end
            end
            
            
            % Find the image transform for that new frame and apply to the
            % misaligned frame
            referenceFrame = dopIn(:, :, timepoint + counter);
            movingFrame = dopIn(:, :, timepoint);
            img_trans = imregtform(referenceFrame, fixedFrame, 'rigid', optimizer, metric);
            temp_frame = imwarp(movingFrame, img_trans, ...
                'OutputView', imref2d(size(fixedFrame)));
            %temp_frame(temp_frame == 0) = NaN;
            dopOut(:,:,timepoint) = temp_frame;
            
            
        end
        
        
        
    end
end
%% imregdeform method
if strcmp(method,'imregdeform')
    
    % Set up a matrix to monitor for failed registrations
    failed_registrations = false(size(dopOut, 3));
    
    % for each frame in the sequence
    ppm = ProgressBar(nTimepoints, 'showWorkerProgress', true);
    parfor timepoint = 1:nTimepoints
        
        % select the current frame to work on
        movingFrame = dopIn(:,:,timepoint);
        
        lastwarn('') % Clear last warning message
        
        % Attempt image registration.
        [~, temp_frame] = imregdeform(zscore(movingFrame, 0, 'all'), fixedFrame, 'DisplayProgress', verbose);
        %temp_frame(temp_frame == 0) = NaN;
        dopOut(:,:,timepoint) = temp_frame;
        
        % Check if a warning was thrown
        [warnMsg, ~] = lastwarn;
        if ~isempty(warnMsg)
            failed_registrations(timepoint) = true;
            dopOut(:, :, timepoint) = movingFrame;
        end
        
        % update the waitbar
        ppm.count();
    end
    delete(ppm);
    
    % Fix the failed registrations by taking the closest successful
    % transform
    
    failed_registrations_1D = failed_registrations(:);
    
    for timepoint = 1:length(failed_registrations_1D)
        if failed_registrations(timepoint)
            %Iteratively find closest good timepoint.
            
            found_replacement_timepoint = false;
            counter = 0;
            while ~found_replacement_timepoint
                counter = counter + 1;
                if (timepoint-counter > 0) && ...
                        ~failed_registrations_1D(timepoint - counter)
                    found_replacement_timepoint = true;
                    counter = -counter;
                elseif (timepoint+counter <= length(failed_registrations_1D)) && ...
                        ~failed_registrations_1D(timepoint + counter)
                    found_replacement_timepoint = true;
                end
            end
            
            
            % Find the image transform for that new frame and apply to the
            % misaligned frame
            referenceFrame = zscore(dopIn(:, :, timepoint + counter), 0, 'all');
            movingFrame = dopIn(:, :, timepoint);
            img_trans = imregdeform(zscore(referenceFrame, 0, 'all'), fixedFrame, 'DisplayProgress', verbose);
            temp_frame = imwarp(movingFrame, img_trans, ...
                'OutputView', imref2d(size(fixedFrame)));
            %temp_frame(temp_frame == 0) = NaN;
            dopOut(:,:,timepoint) = temp_frame;
            
            
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

