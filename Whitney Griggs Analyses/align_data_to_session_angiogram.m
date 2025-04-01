function aligned_angiogram = align_data_to_session_angiogram(angiogram, session_run_list)
%% Align the angiogram to a specified subsession
% This is similar to `align_anatomical_data.m`, but uses one of the
% sessions' angiogram instead of the anatomical images acquired on a
% separate date. It ensures a better match especially when there are only a
% few sessions that are being aligned into a composite functional image.

%% Load some things
processed_data_path = get_user_data_path('pathType', 'doppler');


%% Load and align each angiogram

% Find all the unique AP planes
ap_planes = unique(session_run_list(:, 3));
num_planes = length(ap_planes);

% Pre-allocate space
aligned_angiogram = cell(num_planes, 1);

% Iterate over each AP plane
for plane = 1:num_planes
   plane_ind = session_run_list(:, 3) == ap_planes(plane);
   
   session_run_list_reduced = session_run_list(plane_ind, :);
   
   if size(session_run_list_reduced, 1) > 1
       % GUI to select which session to use for the angiogram
       runStrings = cell(nnz(plane_ind), 1);
       for i = 1:size(session_run_list_reduced,1)
            session = session_run_list_reduced(i, 1);
            run = session_run_list_reduced(i, 2);
            runStrings{i} = ['Session ' num2str(session) ', Run ' num2str(run)];
       end

        % use a GUI to select which ones you want
        [indx,tf] = listdlg('PromptString','Select Session/Run(s) to load:',...
            'SelectionMode','single',...
            'ListString',runStrings, ...
            'ListSize', [600 300]);
   
   
   % Load the specified angiogram
   session = session_run_list_reduced(indx, 1);
   run = session_run_list_reduced(indx, 2);
   filename = strcat('doppler_S', num2str(session), '_R', num2str(run), '+normcorre.mat');
   session_angiogram = load(fullfile(processed_data_path, filename), 'angiogram');
   session_angiogram = session_angiogram.angiogram;
   
   alignment_tform = align_angiogram(angiogram{plane}, session_angiogram);
    %Convert transform into matlab format
    % Alignment_tform is passed back from ManualImageAlignment
    tform = affine2d(alignment_tform);
    % Align image to match previously loaded data
    aligned_angiogram{plane} = imwarp(...
        session_angiogram, ...
        tform, ...
        'OutputView', imref2d(size(angiogram{plane})));
   else
         aligned_angiogram{plane} = angiogram{plane};
   end 
end
