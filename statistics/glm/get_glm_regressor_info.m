function regressor_table = get_glm_regressor_info(behavior, regressors, varargin)
%% Helper function to generate the start time, duration, and modulation for specified regressors
% Inputs
%   behavior - behavior struct
%   regressors - table with at least 5 columns; type, name, start_label,
%                end_label, event_ind;
% Outputs
%   regressor_table - table; name, start_time, duration, and modulation


% Notes
% Supports "event_indices" if the start and end labels are insufficient
%       to identify conditions. Useful if for example a subset of the
%       trials are for a specific condition but has identical start/end
%       labels to another condition.

% Written by Whitney Griggs September 14, 2022


%% Variable inputs
p = inputParser;
p.addOptional('iti_regressor', true);

p.parse(varargin{:});
inputs = p.Results;

%% Create GLM regressors

% To be modeled as impulse functions
impulse_regressor_ind = strcmp(regressors.type, 'impulse');
num_impulse_regressors = nnz(impulse_regressor_ind);
if num_impulse_regressors
    impulse_regressor_names = regressors.name(impulse_regressor_ind);
    impulse_behavior_labels = regressors.start_label(impulse_regressor_ind);
    impulse_event_ind = regressors.event_ind(impulse_regressor_ind);
    impulse_modulation = regressors.modulation(impulse_regressor_ind);
else
    impulse_regressor_names = [];
    impulse_behavior_labels = [];
    impulse_event_ind = [];
    impulse_modulation = [];
end


% To be modeled as boxcar functions
boxcar_regressor_ind = strcmp(regressors.type, 'boxcar');
num_boxcar_regressors = nnz(boxcar_regressor_ind);
if num_boxcar_regressors
    boxcar_regressor_names = regressors.name(boxcar_regressor_ind);
    boxcar_behavior_start_labels = regressors.start_label(boxcar_regressor_ind);
    boxcar_behavior_end_labels = regressors.end_label(boxcar_regressor_ind);
    boxcar_event_ind = regressors.event_ind(boxcar_regressor_ind);
    boxcar_modulation = regressors.modulation(boxcar_regressor_ind);
else
    boxcar_regressor_names = [];
    boxcar_behavior_start_labels = [];
    boxcar_behavior_end_labels = [];
    boxcar_event_ind = [];
    boxcar_modulation = [];
end

% handle impulse functions first
if num_impulse_regressors
    [impulse_event_timing, impulse_event_duration, impulse_event_modulation] = deal(cell(length(impulse_behavior_labels), 1));
    for event = 1:num_impulse_regressors
        impulse_event_timing{event} = parseTime(behavior, impulse_behavior_labels{event});
        impulse_event_duration{event} = [];
        impulse_event_modulation{event} = impulse_modulation{event};
    end
    
    % Correct for any provided event indices
    for regressor = 1:num_impulse_regressors
        if ~isempty(impulse_event_ind{regressor})
            impulse_event_timing{regressor} = impulse_event_timing{regressor}(impulse_event_ind{regressor});
            impulse_event_modulation{regressor} = impulse_modulation{regressor}(impulse_event_ind{regressor});
        end
    end
else
    [impulse_event_timing, impulse_event_duration, impulse_event_modulation] = deal([]);
end


% Handle boxcar functions
if num_boxcar_regressors
    [boxcar_event_timing, boxcar_event_duration, boxcar_event_modulation] = deal(cell(length(boxcar_behavior_start_labels), 1));
    
    for event = 1:num_boxcar_regressors
        boxcar_event_timing{event} = parseTime(behavior, boxcar_behavior_start_labels{event});
        boxcar_event_duration{event} = parseTime(behavior, boxcar_behavior_end_labels{event}) - parseTime(behavior, boxcar_behavior_start_labels{event});
        boxcar_event_modulation{event} = boxcar_modulation{event};
    end
    
    % Correct for any provided event indices
    for regressor = 1:num_boxcar_regressors
        if ~isempty(boxcar_event_ind{regressor})
            boxcar_event_timing{regressor} = boxcar_event_timing{regressor}(boxcar_event_ind{regressor});
            boxcar_event_duration{regressor} =  boxcar_event_duration{regressor}(boxcar_event_ind{regressor});
            if ~isempty(boxcar_event_modulation{regressor})
                boxcar_event_modulation{regressor} = boxcar_event_modulation{regressor}(boxcar_event_ind{regressor});
            end
        end
    end
else
    [boxcar_event_timing, boxcar_event_duration, boxcar_event_modulation] = deal([]);
end


combined_event_labels = [impulse_regressor_names; boxcar_regressor_names];
combined_event_timing = [impulse_event_timing; boxcar_event_timing];
combined_event_duration = [impulse_event_duration; boxcar_event_duration];
combined_event_modulation = [impulse_event_modulation; boxcar_event_modulation];

% Add iti period if we want
% This has to be processed differently because it spans two trials
if inputs.iti_regressor
    iti_regressor_name = {'ITI'};
    iti_start_label = 'iti';
    iti_end_label = 'trialstart';
    
    iti_event_timing = parseTime(behavior, iti_start_label);
    % Discard the very last ITI because we won't always have a trialstart
    % for it.
    iti_event_timing = iti_event_timing(1:end-1);
    iti_end_timing = parseTime(behavior, iti_end_label);
    
    iti_event_duration = iti_end_timing(2:end) - iti_event_timing;
    
    combined_event_labels = [combined_event_labels; iti_regressor_name];
    combined_event_timing = [combined_event_timing; {iti_event_timing}];
    combined_event_duration = [combined_event_duration; {iti_event_duration}];
    combined_event_modulation = [combined_event_modulation; cell(1, 1)];
end

% Remove all NaNs from the regressors
nRegressors = length(combined_event_labels);
for regressor = 1:nRegressors
    
    % Pull out info about current regressor
    regressor_timing = combined_event_timing{regressor};
    regressor_duration = combined_event_duration{regressor};
    regressor_modulation = combined_event_modulation{regressor};
    
    % Find NaNs
    nan_ind = isnan(regressor_timing);
    if ~isempty(regressor_duration)
        nan_ind = nan_ind | isnan(regressor_duration);
    end
    if ~isempty(regressor_modulation)
        nan_ind = nan_ind | isnan(regressor_modulation);
    end
    
    % Remove NaNs
    combined_event_timing{regressor} = regressor_timing(~nan_ind);
    if ~isempty(regressor_duration)
        combined_event_duration{regressor} = regressor_duration(~nan_ind);
    end
    if ~isempty(regressor_modulation)
        combined_event_modulation{regressor} = regressor_modulation(~nan_ind);
    end
end

regressor_table = table(combined_event_labels, ...
    combined_event_timing, ...
    combined_event_duration, ...
    combined_event_modulation, ...
    'VariableNames', {'regressorName', 'timestamps', 'duration', 'modulation'});