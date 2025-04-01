function data_struct = loadDopplerData4D(session_run_list, varargin)
% load multiple runs of data in for different slot numbers and organize as
% a 4D structure nested within a cell array
%
% Using a cell array avoids the following issues:
% 1. Different imaging plane size
% 2. Different length of trials
% 3. Different number of trials
%
% Eventually, I will need to figure out a way to combine into a single 5D
% structure (slot, z, x, timing, trial).
%
% session_run_list - n X 2 matrix; [session run]. Each row is a
%   single session.
%
% Currently lacks ability to concatenate multiple sessions from same
% coronal plane.


%% Input parser
p = inputParser;
addOptional(p, 'project_record_filename', 'ProjectRecord_paper.json');
parse(p, varargin{:});
inputs = p.Results;

%%

% load Session Info
ProjectRecord = load_json_as_table(inputs.project_record_filename);

% Preallocate space
indx = NaN(size(session_run_list,1),1);

for i = 1:size(session_run_list,1)
    % set current session/run combo & find row in the ProjectRecord table
    session = session_run_list(i,1);
    run = session_run_list(i,2);
    indx(i) = intersect(find(ProjectRecord.Session==session),...
        find(ProjectRecord.Run==run));
end

% Add in slot number for each session
try
    session_ap_pos = [ProjectRecord.ap_plane{indx}];
    if any(isnan(session_ap_pos))
        error('Not all of these specified sessions have an approximate AP Plane. Reverting to using slot #');
    end
    session_run_list(:, 3) = session_ap_pos;
catch
    session_slots_cell = ProjectRecord.Slot(indx);
    session_slots = NaN(length(session_slots_cell), 1);
    for i = 1:length(session_slots_cell)
        if strfind(session_slots_cell{i}, 'cor')
            session_slots(i) = str2num(session_slots_cell{i}(5:end));
        else
            error('This format is not supported currently.');
        end
    end
    session_run_list(:, 3) = session_slots;
end


% Find all unique slot numbers to examine
slotnumbers = unique(session_run_list(:, 3));
num_slotnumbers = length(slotnumbers);

motionCorrectResponse = questdlg('Would you like to load motion corrected data?','select_registration','yes','no','no');
motionCorrectResponse = strcmp(motionCorrectResponse,'yes');

%% get the sessions & runs desired

% create the list of available sessions/runs
Session = ProjectRecord.Session;
Run = ProjectRecord.Run;

% Preallocate variables
[iDop, behavior, coreParams, angiogram, UF] = deal(cell(num_slotnumbers, 1));

% Flip up/down to get larger field of view whenever possible
% Load the last session first because it is probably the largest field
% of view.
session_run_list = sortrows(session_run_list, [1 2], {'descend', 'ascend'});

for i = 1:length(slotnumbers)
    session_load_ind = session_run_list(:, 3) == slotnumbers(i);
    if nnz(session_load_ind) == 1
        % If only one session to load for a given slot
        use_manual_alignment = false;
        
    elseif nnz(session_load_ind) > 1
        % Concatenate  multiple sessions from same slot number.
        % Choose if we want manual or automatic alignment
        alignmentResponse = questdlg('Would you like to use manual alignment across sessions/runs?','Use Manual Alignment?','yes','no','no');
        use_manual_alignment = strcmp(alignmentResponse,'yes');
    end
    
    % Load the data
    [iDop{i}, behavior{i}, coreParams{i},angiogram{i}, UF{i}] = concatenateRuns(...
        session_run_list(session_load_ind, 1:2), ...
        motionCorrectResponse,...
        'manual_alignment', use_manual_alignment, ...
        'project_record_filename', inputs.project_record_filename);
        
   
    % save out a pre-image-registration angiogram
    if ~exist('angiogram','var')  || ~isequal(size(angiogram), size(iDop{i}))
        angiogram{i} = nthroot(squeeze(mean(mean(iDop{i},3),4)),3);
        lowerCutoff = quantile(angiogram{i}(:),0.01);
        angiogram{i}(angiogram{i}<lowerCutoff) = lowerCutoff;
    end
    
    % just some useful values to put in the workspace
    [yPixels(i), xPixels(i), nWindows(i), nTrials(i)] = size(iDop{i});
    
    % image registration (if it didn't load properly)
    if motionCorrectResponse && ~coreParams{i}.motionCorrection
        warning('The data loaded wasn''t motion corrected. Attempting to register now.')
        coreParams{i}.method = 'normcorre';
        disp(['Image registration method: ' coreParams{i}.method])
        iDop{i} = correctMotion(iDop{i}, coreParams);
        coreParams{i}.motionCorrection = true;
    end
    
    % fill in any missing data
    if ~isfield(coreParams{i},'method')
        coreParams{i}.method = 'N/A';
    end
end
%% put variables into workspace
%This can be done by passing output arguments, but to maintain backwards
%compatibility and to keep output line clean, I am using putvar (a MATLAB
%community add-on).

%Turn off this warning message.
warning('off','PUTVAR:overwrite')

%Transfers function variables to workspace.
putvar(Session, Run, yPixels,xPixels,nWindows,nTrials,...
    angiogram, iDop, behavior, coreParams, ...
    ProjectRecord, UF, session_run_list);

data_struct.Session = Session;
data_struct.Run = Run;
data_struct.yPixels = yPixels;
data_struct.xPixels = xPixels;
data_struct.nWindows = nWindows;
data_struct.nTrials = nTrials;
data_struct.angiogram = angiogram;
data_struct.iDop = iDop;
data_struct.behavior = behavior;
data_struct.coreParams = coreParams;
data_struct.ProjectRecord = ProjectRecord;
data_struct.UF = UF;
data_struct.session_run_list = session_run_list;
end

