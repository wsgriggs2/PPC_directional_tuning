function [empiric_mean_angular_error, null_mean_angular_error, pvalue, ste_angular_error] = calculate_angular_error(predicted_class, actual_class, varargin)
%% Generates final mean angular error and associated statistics
% Inputs:
%   predicted_class             - (nx1); The predicted class labels
%   actual_class                - (nx1); The actual class labels
%   varargin
%       num_repetitions         - (scalar); How many replicates for the null
%                                 distribution
%       ignore_middle_points    - (bool); How should we handle the middle points
%                                 for the 8-tgt case?
%       null_mean_angular_error - (nxnum_repetitions); What is the mean angular error of
%                                 the null distribution? This should be
%                                 the mean angular error for
%                                 each trial. It will then grab the
%                                 appropriate final trial to use.
% Outputs:
%   empiric_mean_angular_error  - Actual mean angular error seen in
%                                 experiment
%   null_mean_angular_error     - Null distribution for angular error
%                                 assuming a uniform distribution.
%   pvalue                      - Empiric p-value. What percentage of the
%                                 null distribution exceeds the empiric
%                                 mean angular error.
%
% Written by Whitney Griggs July 1, 2022


%% Handling varargin
p = inputParser;
p.addOptional('num_repetitions',1000);
p.addOptional('ignore_middle_points', false);
p.addOptional('null_mean_angular_error', []);
p.parse(varargin{:});
inputs = p.Results;


%% Clean up data depending on nClasses
% Find how many classes there are.
nClasses = length(unique(actual_class));

switch nClasses
    case 8
        % Find all the middle predictions
        middle_ind = predicted_class==5;
        
        if inputs.ignore_middle_points
            % Remove middle predictions
            direction_predicted = predicted_class(~middle_ind);
            direction_predicted(direction_predicted > 5) = direction_predicted(direction_predicted > 5) - 1;
            
            % Remove middle predictions
            direction_actual = actual_class(~middle_ind);
            direction_actual(direction_actual > 5) = direction_actual(direction_actual > 5) - 1;
        else
            % Assign middle predictions to the worst possible angular error
            % (pi).
            direction_predicted = predicted_class;
            direction_actual = actual_class;
            
            % Map back to 1:8 instead of 1:9
            direction_predicted(direction_predicted>5) = direction_predicted(direction_predicted > 5) - 1;
            direction_actual(direction_actual > 5) = direction_actual(direction_actual > 5) - 1;
            
            % Assign the middle ind to actual+pi for prediction. This is
            % worst case angular error.
            direction_predicted(middle_ind) = mod(direction_actual(middle_ind) + 4, 8);
            direction_predicted(direction_predicted == 0) = 8;
        end
        
        % Map from predicted number to direction
        angle_lookup_table = 0:pi/4:7*pi/4;
        angle_lookup_table = angle_lookup_table([6 7 8 5 1 4 3 2]);
        
    case 2
        direction_predicted = predicted_class;
        direction_actual = actual_class;
        
        % Map from predicted number to direction
        angle_lookup_table = [0 pi];
        
end

%% Load null distribution
nTrials = length(direction_actual);
if isempty(inputs.null_mean_angular_error)
    
    [null_angular_error, null_mean_angular_error] = load_null_angular_error_distribution(nClasses, inputs.num_repetitions);
else
    % if passed as an optional argument
    null_mean_angular_error =  inputs.null_mean_angular_error;
    if size(null_mean_angular_error, 1) ~= inputs.num_repetitions
       warning('Mismatch between specified and actual number of repetitions');
       inputs.num_repetitions = size(null_mean_angular_error, 1);
    end
end

% Should already be loaded, we only need the first nTrials
null_mean_angular_error = null_mean_angular_error(:, nTrials);


% Calculate angular error for each trial
empiric_angular_error = abs(angdiff(angle_lookup_table(direction_predicted), angle_lookup_table(direction_actual)));
empiric_mean_angular_error = mean(empiric_angular_error);
ste_angular_error = std(empiric_angular_error)/sqrt(length(empiric_angular_error));

pvalue = 1 - nnz(null_mean_angular_error > empiric_mean_angular_error)/inputs.num_repetitions;

end



