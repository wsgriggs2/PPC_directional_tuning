%% Preview Power Doppler data

loadDopplerData(whos, 'multiple', false, ...
    'mc_method', 'normcorre', ...
    'project_record_filename', 'ProjectRecord_paper.json');

preview_iDop_video(iDop, session_run_list, coreParams);