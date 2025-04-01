%% Plot specified vasculature background images - meso angiograms

clear
close all
pixelsize = 0.1;

project_record_filename = 'ProjectRecord_paper.json';
ProjectRecord = load_json_as_table(project_record_filename);

sessionRunList = specify_sessions_of_interest('project_record_filename', project_record_filename);


% choose if you want motion corrected data
motionCorrectResponse = questdlg('Would you like to load motion corrected data?','select_registration','yes','no','no');
motionCorrectResponse = strcmp(motionCorrectResponse,'yes');

combine_sessions = questdlg('Would you like to load combine runs within sessions?','combine_session_runs?','yes','no','no');
combine_sessions = strcmp(combine_sessions,'yes');

if combine_sessions
    [session_nums, ia, ic] = unique(sessionRunList(:, 1));
    
else
    ic = 1:size(sessionRunList, 1);
    session_nums = sessionRunList(:, 1);
end

%% use concatenateRuns to load data
figure;
tld = tiledlayout('flow');

for session = 1:length(session_nums)
    [~, ~, ~, ~, ~, angiogram, UF] = ...
        concatenateRuns(...
        sessionRunList(ic == session, :), ...
        motionCorrectResponse, ...
        'manual_alignment', true, ...
        'project_record_filename', project_record_filename);
    
    if length(sessionRunList(ic == session, 1)) > 1
        ProjectRecord_ind = find(ismember(ProjectRecord.Session, sessionRunList(ic == session, 1)) & ismember(ProjectRecord.Run, sessionRunList(ic == session, 2)), 1);
        Session = ProjectRecord.Session(ProjectRecord_ind);
        Run = 0; % If multiple runs, then use zero to specify this
    else
        ProjectRecord_ind = ProjectRecord.Session == sessionRunList(ic == session, 1) & ProjectRecord.Run == sessionRunList(ic == session, 2);
        Session = ProjectRecord.Session(ProjectRecord_ind);
        Run = ProjectRecord.Run(ProjectRecord_ind);
    end
    
    
    % Scale angiogram to better mesovasculature
    angiogram_scaled = nthroot(angiogram, 5);
    
    X_img_mm = pixelsize/2 + (0:size(angiogram,2)-1)*pixelsize;
    try
        Z_img_mm = pixelsize/2 + (0:size(angiogram,1)-1)*pixelsize + UF.Depth(1);
    catch
        warning('Missing info S%dR%d', Session, Run);
        Z_img_mm = pixelsize/2 + (0:size(angiogram,1)-1)*pixelsize;
    end
    
    nexttile();
    
   figure_title = sprintf('S%dR%d', Session, Run);
    
    plot_angiogram(angiogram_scaled, X_img_mm, Z_img_mm, ...
        'title', figure_title,...
        'colormap', 'inferno', ...
        'show_colormap', false)
    
    drawnow;
end

tld.TileSpacing = 'tight';
tld.Padding = 'tight';