function data_struct = loadDopplerDataContinuous(allVariablesInWorkspace, singleORmultiple, suppressScreenOutput, varargin)
% loadDopplerDataContinuous -   Load continuous Doppler data from session and
%                               associated variables
%
% INPUTS:
%   allVariablesInWorkspace     struct; If nonempty, checks workspace for existing
%                               variables and avoids overwriting these
%   single_or_multiple:         string; {'single', 'multiple'}; Load one or
%                               multiple sessions?
%   suppressScreenOutput        Do you want to hide command line output?
%   varargin: 
%     SessionRunList:           (1 x 2) double; [session run]' If specified, 
%                               skips selection GUI. Should match session/run
%                               numbers from ProjectRecord JSON.
%     mc_method:                What method do you want to use for motion
%                               correction?

%     project_record_filename:  string: Name of the ProjectRecord*.json
%                               file.



%      
% OUTPUTS:
%   struct_variables:           Struct; The specified loaded data. If
%                               multiple sessions were specified, they are 
%                               all concatenated rather than kept separate.
%
% Authors: Sumner Norman and Whitney Griggs
% Date: 2018

%% Set some defaults

use_default_state_names = true;
% If `state_names` is empty, then will load a GUI later on to let user 
% select appropriate state names
if use_default_state_names
    state_names = {'trialstart', 'initialfixation', 'fixationhold', 'cue', ...
        'memory', 'target_acquire', 'target_hold', 'reward', 'iti'}; 
else
    state_names = [];
end

%% Parse inputs
if ~exist('singleORmultiple','var')
    singleORmultiple = 'multiple';
end
if ~exist('suppressScreenOutput','var')
    suppressScreenOutput = false;
end

%% Input parser
p = inputParser;
addOptional(p, 'SessionRunList', []);
addOptional(p, 'mc_method', NaN); %normcorre, imregister, or NaN
% If user does not specify which project record JSON file to use, then ask
% them.
addOptional(p, 'project_record_filename', '');
parse(p, varargin{:});
inputs = p.Results;

if isempty(inputs.project_record_filename)
    inputs.project_record_filename = select_project_record_filename;
end


%% if the data is already in the workspace, skip loading and return
VariablesReqd = {'data','Session','Run'};
if exist('allVariablesInWorkspace','var')
    varsExist = ismember(VariablesReqd,{allVariablesInWorkspace.name});
    if all(varsExist)
        if ~isempty(inputs.SessionRunList)
            choice = questdlg('Would you like to (re)load the data?',...
                'Data (reload)?', 'Yes','No','Yes');
            if ~strcmp(choice, 'Yes')
                disp('doppler data already loaded. continuing')
                if nargout > 0
                    error('This function cannot pass back the desired data. Fix your code.');
                end
                return
            end
        else
            disp('doppler data already loaded. continuing')
            if nargout > 0
                error('This function cannot pass back the desired data. Fix your code.');
            end
            return
        end
    end
end

% load Session Info
ProjectRecord = load_json_as_table(inputs.project_record_filename);
if ~suppressScreenOutput
    disp(ProjectRecord)
end

%% get the sessions & runs desired
if isempty(inputs.SessionRunList)
    % create the list of available sessions/runs
    runStrings = cell(size(ProjectRecord,1),1);

    for i = 1:size(ProjectRecord,1)
        Session = ProjectRecord.Session(i);
        Run = ProjectRecord.Run(i);
        Slot = ProjectRecord.Slot(i);
        Monkey = ProjectRecord.Monkey(i);
        nTargets = ProjectRecord.nTargets(i);
        Project = ProjectRecord.Project(i);
        
        runStrings{i} = ['Monkey ' char(Monkey) ', Session ' num2str(Session) ', Run ' num2str(Run) ', Slot ' num2str(Slot) ', Project ' char(Project) ', nTargets ' num2str(nTargets)];

    end
    
    % use a GUI to select which ones you want
    [indx,tf] = listdlg('PromptString','Select Session/Run(s) to load:',...
        'SelectionMode',singleORmultiple,...
        'ListString',runStrings, ...
        'ListSize', [600 300]);
    
    sessionRunList = [ProjectRecord.Session(indx) ProjectRecord.Run(indx)];
else
    sessionRunList = inputs.SessionRunList;
    
    tf = true;
end


if size(sessionRunList, 1) > 1
    % Choose if we want manual or automatic alignment
    alignmentResponse = questdlg('Would you like to use manual alignment across sessions/runs?','Use Manual Alignment?','yes','no','no');
    use_manual_alignment = strcmp(alignmentResponse,'yes');
else
    use_manual_alignment = false;
end

if tf
    %% use concatenateRuns to load data
    if isnan(inputs.mc_method)
        % choose if you want motion corrected data
        motionCorrectResponse = questdlg('Would you like to load motion corrected data?','select_registration','yes','no','no');
        motionCorrectResponse = strcmp(motionCorrectResponse,'yes');
        mc_method = 'normcorre';
    else
        motionCorrectResponse = any(inputs.mc_method);
        mc_method = inputs.mc_method;
    end
    

    % Load the last session first because it is probably the largest field
    % of view.
    sessionRunList = sortrows(sessionRunList, [1 2], {'descend', 'ascend'});
    
    data = concatenateContinuousRuns(...
        sessionRunList, ...
        motionCorrectResponse, ...
        'manual_alignment', use_manual_alignment, ...
        'mc_method', mc_method, ...
        'state_names', state_names, ...
        'project_record_filename', inputs.project_record_filename);

    %% save out a pre-image-registration angiogram
    if ~isfield(data, 'angiogram') || ~isequal(size(data.angiogram), size(data.dop))
        data.angiogram = makeAngiogram(data.dop);
    end
      
    %% image registration (if it didn't load properly)
    if motionCorrectResponse && ~data.coreParams.motionCorrection
        warning('The data loaded wasn''t motion corrected. Attempting to register now.')
        data.coreParams.method = 'normcorre';
        disp(['Image registration method: ' data.coreParams.method])
        data.dop = correctMotion_continuous(data.dop, data.coreParams);
        data.coreParams.motionCorrection = true;
    end
    
    % fill in any missing data
    if ~isfield(data.coreParams,'method')
        data.coreParams.method = 'N/A';
    end
    
    % set session/run (run=0 for display purposes if we have multiple runs)
    Session = sessionRunList(end,1);
    Run = sessionRunList(end,2);
    if size(sessionRunList,1)>1
        Run = 0;
    end
    
    %% clear some memory
    clear moreFlag fileName allVariables
    
    %% put variables into workspace
    %This can be done by passing output arguments, but to maintain backwards
    %compatibility and to keep output line clean, I am using putvar (a MATLAB
    %community add-on).
    
    %Turn off this warning message.
    warning('off','PUTVAR:overwrite')
    
    %Transfers function variables to workspace.
    putvar(Session, Run, data, sessionRunList, ProjectRecord);
end
putvar(tf);

data_struct.Session = Session;
data_struct.Run = Run;
data_struct.data = data;
data_struct.sessionRunlist = sessionRunList;
data_struct.ProjectRecord = ProjectRecord;
