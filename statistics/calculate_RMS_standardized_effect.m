function psi = calculate_RMS_standardized_effect(data, condition_labels)

% Calculate psi, root-mean-square standardized effect
% This is based upon Wikipedia implementation.

% Inputs:
%   * data - 1D vector of all the observations
%   * condition_labels - Condition label for each observation
% Outputs:
%   * psi - The root-mean-square standardized effect


% Calculate group-specific quantities
unique_labels = unique(condition_labels(~isnan(condition_labels)));
n_groups = length(unique_labels);
if n_groups > 0
    [group_means, squared_error, group_size] = deal(NaN(n_groups, 1));
    for group = 1:n_groups
        group_data = data(condition_labels==unique_labels(group));
        group_data(isnan(group_data)) = [];
        group_size(group) = numel(group_data);
        group_means(group) = mean(group_data);
        squared_error(group) = sum((group_data - group_means(group)).^2);
    end

    % Calculate overall variables
    population_mean = mean(data, 'omitnan');
    MSE = sum(squared_error)/(sum(group_size)-n_groups);
    squared_difference = sum((group_means(group) - population_mean)^2);

    % Calculate RMSSE
    psi = sqrt(1/(n_groups-1)*squared_difference/MSE);
else
    psi = NaN;
    warning('No unique groups found. Returning NaN');
end