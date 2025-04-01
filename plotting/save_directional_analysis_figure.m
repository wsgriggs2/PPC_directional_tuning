function save_directional_analysis_figure(fig, filename_suffix, monkey_name)
    base_filename = sprintf('across_session_analysis_%s_', monkey_name);
    base_filepath = fullfile(get_user_data_path('pathType', 'across session analyses'));

    full_filename = fullfile(base_filepath, sprintf('%s%s_dg%s.svg', ...
        base_filename, ...
        filename_suffix, ...
        datetime('today', format='yyyyMMdd')));
    saveas(fig, full_filename, 'svg');
end