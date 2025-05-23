function json_table = load_json_as_table(filename)

%% Function to load in json session or project records as a table
% Inputs:
%   filename - String indicating filename of the json record to be loaded.
% Outputs:
%   json_table - Table; one row per session with metadata as columns.
%
% Author: Whitney Griggs
% Date: March 18, 2025

json_file = fileread(filename);
json_info = jsondecode(json_file);

if size(json_info, 1) > 1
    json_table = struct2table(json_info);
else
    json_table = struct2table(json_info, 'AsArray', true);
end

column_names = json_table.Properties.VariableNames;

% convert all char arrays to strings (except for the Notes entries)
for column = 1:width(json_table)
    if iscellstr(json_table{:, column}) && ...
            ~strcmp(json_table(:, column).Properties.VariableNames, 'Notes')   
        json_table.(column_names{column}) = string(json_table{:, column}); 
    end    
end