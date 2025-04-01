function [dopOut, behaviorOut, coreParams, angiogram, UFout]...
    = concatenateRuns(sessionRunList, registered, varargin)
% [dopOut, behaviorOut, coreParams] = concatenateRuns(sessionRunList)
%
% session/run list is of the format:
% session run
% [  1     1  ]
% [  1     2  ]
% [  2     1  ]
% [    ...    ]
%
% registered is a boolean whether you want to request motion corrected data

% When aligning motion-corrected data, it first loads in the
% motion-corrected data for each session, then it uses manual alignment to
% optimize the inter-session alignment. Due to slightly different
% positioning and contrast between sessions, the manual alignment is the
% best option Whitney has found yet.

p = inputParser;
addOptional(p, 'manual_alignment',  false, @islogical);
addOptional(p, 'mc_method', 'normcorre');
addOptional(p, 'project_record_filename', '');
addOptional(p, 'load_all_attempted_trials', false);
parse(p, varargin{:});
inputs = p.Results;

if isempty(inputs.project_record_filename)
    inputs.project_record_filename = select_project_record_filename;
end


%% check that data type is consistent across the runs requested
ProjectRecord = load_json_as_table(inputs.project_record_filename);

% Preallocate space
sessioninfo_indx = NaN(size(sessionRunList,1),1);

for i = 1:size(sessionRunList,1)
    % set current session/run combo & find row in the ProjectRecord table
    session = sessionRunList(i,1);
    run = sessionRunList(i,2);
    sessioninfo_indx(i) = intersect(find(ProjectRecord.Session==session),...
        find(ProjectRecord.Run==run));
end

% Check if more than one monkey is used
monkeyPass = isscalar(unique(ProjectRecord.Monkey(sessioninfo_indx)));

% Check if all are same task
taskPass = isscalar(unique(ProjectRecord.Task(sessioninfo_indx)));

% Check what recording system used
system_used = ProjectRecord.RecordingSystem(sessioninfo_indx(1));
systemPass = all(strcmp(ProjectRecord.RecordingSystem(sessioninfo_indx), ProjectRecord.RecordingSystem(sessioninfo_indx(1))));

% Pause if the sessions are deemed incompatible
if ~(monkeyPass && taskPass && systemPass)
    warning('The selected sessions may not be compatible!')
    input('Press any key to continue. (ctrl+c to quit)')
end

%% load the data and concatenate
for i = 1:size(sessionRunList,1)
    % create filename string for this session/run
    session = sessionRunList(i,1);
    run = sessionRunList(i,2);
    dopplerFileName = strcat('doppler_S', num2str(session), '_R', num2str(run), '.mat');

    if inputs.load_all_attempted_trials
       dopplerFileName = strrep(dopplerFileName, '.mat', '_allTrials.mat'); 
    end

    % Load motion corrected data
    if registered
        if ~inputs.manual_alignment
            if size(sessionRunList,1)==1
                % Only load motion corrected data if not concatenated.
                dopplerFileName = strrep(dopplerFileName,'.mat',sprintf('+%s.mat', inputs.mc_method));
            end
        else
            dopplerFileName = strrep(dopplerFileName,'.mat',sprintf('+%s.mat', inputs.mc_method));
        end
    end

    % load data
    try
        fprintf('Loading ''%s'' from %s... \n',dopplerFileName,get_user_data_path('pathType','doppler'))
        load(fullfile(get_user_data_path('pathType','doppler'), dopplerFileName), 'iDop', 'behavior', 'coreParams','UF','angiogram')
        if registered
            data_registered = true;
        end

    catch
        if registered
            warning('Motion corrected data not found. Attempting to load unregistered data.')
            dopplerFileName = strrep(dopplerFileName,sprintf('+%s', inputs.mc_method),'');

            fprintf('Loading ''%s'' from %s... \n',dopplerFileName,get_user_data_path('pathType','doppler'))
            load(fullfile(get_user_data_path('pathType','doppler'),dopplerFileName), 'iDop', 'behavior', 'coreParams','UF');
            data_registered = false;

        else
            error('Cannot load specified file.');
        end
    end

    if registered && ~data_registered
        % If the request is to motion correct and the loaded data is not
        % motion corrected, then motion correct the data.
        coreParams.motionCorrection = true;
        coreParams.method = inputs.mc_method;
        [iDop, coreParams] = correctMotion(iDop, coreParams);
    end

    if i == 1
        % start by initializing to first instance

        dopOut = iDop;
        behaviorOut = behavior;
        UFout = UF;

    else % we are on the second run or later, so concatenate them
        % make sure that we have the same trial length (sometimes, if a run
        % has an epoch length that jitters above or below a full second,
        % this can get rounded up/down based on jitter between runs.

        % First align
        if inputs.manual_alignment

            % For now, will use entire DopOut to make angiogram, but can make
            % improvements as needed. For example, using last 30-40 frames to
            % make angiogram to align to. This would perhaps account for slow
            % drift in anatomy better.
            dopOut_angiogram = makeAngiogram(dopOut);
            iDop_angiogram = makeAngiogram(iDop);

            % Use manual alignment to register between sessions/runs
            alignment_tform = align_angiogram(dopOut_angiogram, iDop_angiogram);
            %Convert transform into matlab format
            % Alignment_tform is passed back from ManualImageAlignment
            tform = affine2d(alignment_tform);
            % Align image to match previously loaded data
            iDop = imwarp(...
                iDop, ...
                tform, ...
                'OutputView', imref2d(size(dopOut_angiogram)));

        end

        % Can either increase size of smaller array to match size of larger
        % arry or decrease size. Tried increasing the size, but it caused
        % lots of downstream problems. Leaving the code in case it is
        % useful in the future, but using the minimum dimensions for now.
        %
        % Note as of 12/5/22, this may be obsolete if the image alignment
        % is used since that should standardize imagesize already.
        use_min_dimensions = true;
        if ~all([size(dopOut,1)==size(iDop,1),size(dopOut,2)==size(iDop,2),size(dopOut,3)==size(iDop,3)])
            if use_min_dimensions
                % Use minimum spatial dimensions and minimum time dimension
                minDims = min(size(dopOut), size(iDop));
                dopOut = dopOut(1:minDims(1), 1:minDims(2), 1:minDims(3), :);
                iDop = iDop(1:minDims(1), 1:minDims(2), 1:minDims(3), :);
            else
                % Use max spatial dimensions and minimum time dimension
                % Assuming that width will be consistently 128 elements (or
                % interpolated version of 128 elements).
                minDims = min(size(dopOut), size(iDop));
                if size(iDop, 1) > size(dopOut, 1)
                    padding = zeros(size(iDop, 1) - size(dopOut, 1), size(dopOut, 2), size(dopOut, 3), size(dopOut, 4));
                    dopOut = [dopOut; padding;];
                elseif size(iDop, 1) < size(dopOut, 1)
                    padding = zeros(size(dopOut, 1) - size(iDop, 1), size(iDop, 2), size(iDop, 3));
                    dopOut = [iDop; padding;];
                end

                dopOut = dopOut(:, :, 1:minDims(3), :);
                iDop = iDop(:, :, 1:minDims(3), :);
            end
        end

        % Concatenate across 4th dimension (trials)
        dopOut = cat(4, dopOut, iDop);


        % update UF
        if ~isfield(UF,'Depth')
            warning('UF missing depth field')
        elseif ~isequal(UF.Depth,UFout.Depth)
            warning('Different depths between the sessions. Using smaller depth');
        end
        % update behavior
        try
            % Handle structs with different number of fields
            behaviorOut_fields = fieldnames(behaviorOut);
            behavior_fields = fieldnames(behavior);
            if length(behaviorOut_fields) >  length(behavior_fields)
                warning('BehaviorOut has more fields than BehaviorIn. Make sure the new fields added make sense.');
                new_fields = behaviorOut_fields(~ismember(behaviorOut_fields, behavior_fields));
                empty_cells = cell(length(behavior), 1);
                for field = 1:length(new_fields)
                    [behavior.(new_fields{field})] = empty_cells{:};
                end
            elseif length(behaviorOut_fields) <  length(behavior_fields)
                warning('BehaviorOut has more fields than BehaviorIn. Make sure the new fields added make sense.');
                behaviorOut_fields = fieldnames(behaviorOut);
                behavior_fields = fieldnames(behavior);
                new_fields = behavior_fields(~ismember(behavior_fields, behaviorOut_fields));
                empty_cells = cell(length(behaviorOut), 1);
                for field = 1:length(new_fields)
                    [behaviorOut.(new_fields{field})] = empty_cells{:};
                end
            end
            behaviorOut = [behaviorOut behavior];
            fprintf('Successful trials concatenated from %i to %i\n', ...
                size(behaviorOut,2)-size(behavior,2)+1, size(behaviorOut,2))
        catch
            disp('Could not concatenate behavior files');
        end % Try catch loop for behavior
    end % If over session number
end % For loop over runs

switch system_used
    case 'prototypeRT'
        coreParams.framerate = 2;
    case 'offlineVantage'
        coreParams.framerate = 1;
end

% save out a pre-image-registration angiogram
if ~exist('angiogram','var')
    angiogram = makeAngiogram(dopOut);
end

% Motion correct entire concatenated dataset. Needs to be done at
% the very end to account for shifts in position between sessions.
if registered && size(sessionRunList,1)>1
    coreParams.motionCorrection = true;
    coreParams.method = inputs.mc_method;
    [dopOut, coreParams] = correctMotion(dopOut, coreParams);
elseif size(sessionRunList, 1) > 1
    warning('You chose to use multiple sessions without alignment.');
end


fprintf('done.\n')


end %Main function end