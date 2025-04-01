function session_run_list = specify_sessions_of_interest(varargin)
%% Helper function to specify which sessions/runs we want to load
% Loads specified ProjectRecord JSON file, gives user an option to select
% which sessions/runs to load, and then returns `session_run_list`
% 
% Author: Whitney Griggs
% Date: March 18, 2025

%% Input parser
p = inputParser;
addOptional(p, 'project_record_filename', 'ProjectRecord_paper.json');
addOptional(p, 'single_v_multiple', 'multiple');
parse(p, varargin{:});
inputs = p.Results;

if isempty(inputs.project_record_filename)
    inputs.project_record_filename = select_project_record_filename;
end

%% Load Project Record
ProjectRecord = load_json_as_table(inputs.project_record_filename);


%% get the sessions & runs desired

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
    'SelectionMode',inputs.single_v_multiple,...
    'ListString',runStrings, ...
    'ListSize', [600 300]);

session_run_list = [ProjectRecord.Session(indx) ProjectRecord.Run(indx)];

end