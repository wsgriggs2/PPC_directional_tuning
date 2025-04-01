function [data]...
    = concatenateContinuousRuns(sessionRunList, registered, varargin)
% data = concatenateRuns(sessionRunList)
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
addOptional(p, 'mc_method', 'imregister');
addOptional(p, 'state_names', []);
addOptional(p, 'project_record_filename', '');
parse(p, varargin{:});
inputs = p.Results;

if isempty(inputs.project_record_filename)
    inputs.project_record_filename = select_project_record_filename;
end

%% check that data type is consistent across the runs requested
ProjectRecord = load_json_as_table(inputs.project_record_filename);

% Preallocate space
ProjectRecord_indx = NaN(size(sessionRunList,1),1);

for i = 1:size(sessionRunList,1)
    % set current session/run combo & find row in the ProjectRecord table
    session = sessionRunList(i,1);
    run = sessionRunList(i,2);
    ProjectRecord_indx(i) = intersect(find(ProjectRecord.Session==session),...
        find(ProjectRecord.Run==run));
end

% Check if more than one monkey is used
monkeyPass = length(unique(ProjectRecord.Monkey(ProjectRecord_indx)))==1;

% Check if all are same task
taskPass = length(unique(ProjectRecord.Task(ProjectRecord_indx)))==1;

% Check that same recording system was used for all sessions
% Check that same recording system was used for all sessions
system_used = ProjectRecord.RecordingSystem(ProjectRecord_indx(1));
systemPass = all(strcmp(ProjectRecord.RecordingSystem(ProjectRecord_indx), ProjectRecord.RecordingSystem(ProjectRecord_indx(1))));

% Pause if the sessions are deemed incompatible
if ~(monkeyPass && taskPass && systemPass)
    fprintf('The selected sessions may not be compatible!\n')
    input('Press any key to continue. (ctrl+c to quit)')
end

%% load the data and concatenate 
for i =1:size(sessionRunList,1) 
    % create filename string for this session/run
    session = sessionRunList(i,1);
    run = sessionRunList(i,2);
    dopplerFileName = strcat('dopplerContinuous_S', num2str(session), '_R', num2str(run), '.mat');
    
    
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
        load(fullfile(get_user_data_path('pathType','doppler'), dopplerFileName));
        if registered
            data_registered = true;
        end
    catch
        if registered
            warning('Motion corrected data not found. Attempting to load unregistered data.')
            dopplerFileName = strrep(dopplerFileName,sprintf('+%s', inputs.mc_method),'');
                fprintf('Loading ''%s'' from %s... \n',dopplerFileName,get_user_data_path('pathType','doppler'))
                load(fullfile(get_user_data_path('pathType','doppler'),dopplerFileName));
                data_registered = false;
        else
            error('Cannot load specified file.');
        end
    end
    
    % This may be a confusing comment in the future, but this statement
    % should work even when using automatic alignment since registered will
    % be true, but the file loaded will be the non-motion corrected files.
    % If those files fail to load, the code should break since there are
    % missing files.
    if registered && ~data_registered
        % If the request is to motion correct and the loaded data is not
        % motion corrected, then motion correct the data.
        coreParams.motionCorrection = true;
        coreParams.method = inputs.mc_method;
        [iDop, coreParams] = correctMotion(iDop, coreParams);
    end
    
    % Add field to struct to represent session number and run number
    [behavior.session] = deal(session);
    [behavior.run] = deal(run);
    
    % Add fields to timestamps to represent session number and run number
    timestamps(:, 2) = session;
    timestamps(:, 3) = run;
    
    % Adjust timing in behavior, timestamps 
    minTime = min([parseTime(behavior(1), 'trialstart'), timestamps(:, 1)'], [], 'all'); 
    timestamps(:, 1) = timestamps(:, 1) - minTime;
    possible_statenames = fieldnames(behavior);
    
    % Note - As of August 31, 2022, this code can be simplified. By just
    % defining the possible statenames using similar logic to what follows,
    % we can reduce the redundant code (namely `adjust_behavior_time`) and
    % improve clarity of this code. 
    if isempty(inputs.state_names)
        try
            % Select all valid statenames
            [indx,tf] = listdlg('PromptString','Select which fields are the state names:',...
            'SelectionMode','multiple',...
            'ListString',possible_statenames);
            state_names = possible_statenames(indx);
            behavior_adjusted = adjust_behavior_time(behavior, possible_statenames(indx), -minTime);
        catch
            while tf
                warning('You (likely) choose incorrect state names. Try again');
                % Select all valid statenames
                [indx,tf] = listdlg('PromptString','Select which fields are the state names:',...
                'SelectionMode','multiple',...
                'ListString',possible_statenames);
                if tf
                    behavior_adjusted = adjust_behavior_time(behavior, possible_statenames(indx), minTime);
                    state_names = possible_statenames(indx);
                end
            end
            if ~tf
                error('No state names chosen');
            end
        end
    else
        try
            behavior_adjusted = adjust_behavior_time(behavior, inputs.state_names, -minTime);
            state_names = inputs.state_names;
        catch
            error('The state names you specified as a variable input are likely invalid. Please try again.');
        end
    end
        
    
    if i == 1 
        % start by initializing to first instance
        
        dopOut = iDop;
        behaviorOut = behavior_adjusted;
        behaviorOut_statenames = state_names;
        UFout = UF;
        timestamps_out = timestamps;
 
        
    else % we are on the second run or later, so concatenate them
        % make sure that we have the same trial length (sometimes, if a run
        % has an epoch length that jitters above or below a full second,
        % this can get rounded up/down based on jitter between runs.
        
        % Find the last timestamp in `behaviorOut`
        last_behavior_entry = behaviorOut(end);
        for k = length(behaviorOut_statenames):-1:1
            possible_last_timestamp = last_behavior_entry.(behaviorOut_statenames{k});
            if ~isempty(possible_last_timestamp)
                last_state_in_final_trial = behaviorOut_statenames{k};
                break;
            end
        end
        
       
        
        % Adjust timing in behavior, timestamps
        % Use last state selected as the final state. Normally this is iti,
        % but trying to avoid hard-coding this for future extensibility.
        maxTime = max([parseTime(behaviorOut(end), last_state_in_final_trial), timestamps_out(:, 1)'], [], 'all');
        timestamps(:, 1) = timestamps(:, 1) + maxTime + 1; % Shift by 1 additional second to avoid overlap of any timestamps
        behavior_adjusted = adjust_behavior_time(behavior_adjusted, state_names, maxTime);
        
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
        
        
        if ~all([size(dopOut,1)==size(iDop,1),size(dopOut,2)==size(iDop,2)])
            % Use minimum spatial dimensions and minimum time dimension
            minDims = min(size(dopOut), size(iDop));
            dopOut = dopOut(1:minDims(1), 1:minDims(2), :);
            iDop = iDop(1:minDims(1), 1:minDims(2), :);
        end
        
        % Concatenate in time
        dopOut =  cat(3, dopOut, iDop);
        
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
            behavior_fields = fieldnames(behavior_adjusted);
            if length(behaviorOut_fields) >  length(behavior_fields)
                warning('BehaviorOut has more fields than BehaviorIn. Make sure the new fields added make sense.');
                new_fields = behaviorOut_fields(~ismember(behaviorOut_fields, behavior_fields));
                empty_cells = cell(length(behavior_adjusted), 1);
                for field = 1:length(new_fields)
                    [behavior_adjusted.(new_fields{field})] = empty_cells{:};
                end
            elseif length(behaviorOut_fields) <  length(behavior_fields)
                warning('BehaviorOut has more fields than BehaviorIn. Make sure the new fields added make sense.');
                behaviorOut_fields = fieldnames(behaviorOut);
                behavior_fields = fieldnames(behavior_adjusted);
                new_fields = behavior_fields(~ismember(behavior_fields, behaviorOut_fields));
                empty_cells = cell(length(behaviorOut), 1);
                for field = 1:length(new_fields)
                    [behaviorOut.(new_fields{field})] = empty_cells{:};
                end
            end
            behaviorOut = [behaviorOut behavior_adjusted];
            fprintf('Successful trials concatenated from %i to %i\n', ...
                size(behaviorOut,2)-size(behavior_adjusted,2)+1, size(behaviorOut,2))
        catch
            disp('Could not concatenate behavior files');
        end % Try catch loop for behavior
        
        % Concatenate the timestamps
        timestamps_out = [timestamps_out; timestamps];
        
    end % If over session number
end % For loop over runs



if inputs.manual_alignment
    % save out a post-image-registration angiogram
    if ~exist('angiogram','var')
        angiogram = makeAngiogram(dopOut);
    end
else
    % save out a pre-image-registration angiogram
    if ~exist('angiogram','var')
        angiogram = makeAngiogram(dopOut);
    end
    
    % Motion correct entire concatenated dataset. Needs to be done at
    % the very end to account for shifts in position between sessions.
    if registered && size(sessionRunList,1)>1
        coreParams.motionCorrection = true;
        coreParams.method = inputs.mc_method;
        [dopOut, coreParams] = correctMotion_continuous(dopOut, coreParams);
    elseif size(sessionRunList, 1) > 1
        warning('You chose to use multiple sessions without alignment.');
    end
end

switch system_used
    case 'prototypeRT'
        coreParams.framerate = 2;
    case 'offlineVantage'
        coreParams.framerate = 1;
end

fprintf('done.\n')

% Choose what information to save out
data.dop = dopOut;
data.coreParams = coreParams;
data.angiogram = angiogram;
data.behavior = behaviorOut;
data.UF = UFout;
data.timestamps = timestamps_out;

end %Main function end