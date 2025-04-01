%% Script to generate angiograms for specific session

%% Specify which sessions we want to generate angiograms for
project_record_filename = 'ProjectRecord_paper.json';
session_run_list = specify_sessions_of_interest('project_record_filename', project_record_filename);
n_sessions = size(session_run_list, 1);

%% Retrieve and plot angiograms within a single figure

% Open figure and create tiledlayout
figure;
tld = tiledlayout('flow');

% For each session, retrieve and plot angiogram as one tile in figure
for session = 1:n_sessions
    % Load data
    loadDopplerData([], 'multiple', true, ...
        'sessionRunList', session_run_list(session, :), ...
        'project_record_filename', 'ProjectRecord_paper.json', ...
        'mc_method', 'normcorre');

    % Define image size
    pixelsize = 0.1;
    X_img_mm = pixelsize/2 + (0:size(iDop,2)-1)*pixelsize; 
    Z_img_mm = pixelsize/2 + (0:size(iDop,1)-1)*pixelsize + UF.Depth(1);
    
    % Plot angiogram
    nexttile(session);
    plot_angiogram(angiogram, X_img_mm, Z_img_mm, ...
        'title', sprintf('S%dR%d', session_run_list(session, :)));

end