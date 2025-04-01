function struct_variables = loadDopplerData(allVariablesInWorkspace, singleORmultiple, suppressScreenOutput, varargin)
% loadDopplerData -             Loads task-aligned doppler data and associated
%                               variables
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
%     load_all_attempted_trials:For the real-time BMI data, do you want to
%                               load all trials regardless of successful 
%                               decoding?


%      
% OUTPUTS:
%   struct_variables:           Struct; The specified loaded data. If
%                               multiple sessions were specified, they are 
%                               all concatenated rather than kept separate.
%
% Authors: Sumner Norman and Whitney Griggs
% Date: 2018

%% Parse inputs
if ~exist('singleORmultiple','var')
    singleORmultiple = 'multiple';
end
if ~exist('suppressScreenOutput','var')
    suppressScreenOutput = true;
end

p = inputParser;
addOptional(p, 'SessionRunList', []);
addOptional(p, 'mc_method', NaN);

% If user does not specify which project record JSON file to use, then ask
% them.
addOptional(p, 'project_record_filename', '');
addOptional(p, 'load_all_attempted_trials', false);
parse(p, varargin{:});
inputs = p.Results;

if isempty(inputs.project_record_filename)
    inputs.project_record_filename = select_project_record_filename;
end

%% if the data is already in the workspace, skip loading and return
VariablesReqd = {'iDop','Session','Run','behavior'};
if exist('allVariablesInWorkspace','var') && ~isempty(allVariablesInWorkspace)
    varsExist = ismember(VariablesReqd,{allVariablesInWorkspace.name});
    if all(varsExist)
        disp('doppler data already loaded. continuing')
        return
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
    % choose if you want motion corrected data
    if isnan(inputs.mc_method)
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

    [iDop, behavior, coreParams,angiogram, UF] = ...
        concatenateRuns(...
        sessionRunList, ...
        motionCorrectResponse, ...
        'manual_alignment', use_manual_alignment, ...
        'mc_method', mc_method, ...
        'project_record_filename', inputs.project_record_filename, ...
        'load_all_attempted_trials', inputs.load_all_attempted_trials);

    %% save out a pre-image-registration angiogram
    if ~exist('angiogram','var') || ~isequal(size(angiogram), size(iDop))
        angiogram = nthroot(squeeze(mean(mean(iDop,3),4)),3);
        lowerCutoff = quantile(angiogram(:),0.01);
        angiogram(angiogram<lowerCutoff) = lowerCutoff;
    end

    %% useful values to put in the workspace
    [yPixels, xPixels, nWindows, nTrials] = size(iDop);

    %% image registration (if it didn't load properly)
    if motionCorrectResponse && ~coreParams.motionCorrection
        warning('The data loaded wasn''t motion corrected. Attempting to register now.')
        coreParams.method = 'normcorre';
        disp(['Image registration method: ' coreParams.method])
        iDop = correctMotion(iDop, coreParams);
        coreParams.motionCorrection = true;
    end

    % fill in any missing data
    if ~isfield(coreParams,'method')
        coreParams.method = 'N/A';
    end

    % set session/run (run=0 for display purposes if we have multiple runs)
    Session = sessionRunList(end,1);
    Run = sessionRunList(end,2);
    if size(sessionRunList,1)>1, Run = 0; end

    %% clear some memory
    clear moreFlag fileName allVariables

    %% put variables into workspace
    %This can be done by passing output arguments, but to maintain backwards
    %compatibility and to keep output line clean, I am using putvar (a MATLAB
    %community add-on).

    %Turn off this warning message.
    warning('off','PUTVAR:overwrite')

    %Transfers function variables to workspace.
    putvar(Session, Run, yPixels, xPixels, nWindows, nTrials,...
        angiogram, iDop, behavior, coreParams, sessionRunList, ...
        ProjectRecord, UF);

    % Pass back as struct as well for cleanliness
    struct_variables.Session = Session;
    struct_variables.Run = Run;
    struct_variables.yPixels = yPixels;
    struct_variables.xPixels = xPixels;
    struct_variables.nWindows = nWindows;
    struct_variables.nTrials = nTrials;
    struct_variables.angiogram = angiogram;
    struct_variables.iDop = iDop;
    struct_variables.behavior = behavior;
    struct_variables.coreParams = coreParams;
    struct_variables.sessionRunList = sessionRunList;
    struct_variables.ProjectRecord = ProjectRecord;
    struct_variables.UF = UF;
end
putvar(tf);
