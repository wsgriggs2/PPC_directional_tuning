function d = calculate_cohen_d(x1, x2)

% Calculate Cohen's D for effect size of the difference between two means
% Inputs
%   * x1, x2 - Vectors to be compared.
%
% Outputs
%   * d - Effect size measure

% Calculate basic quantities
n1 = nnz(~isnan(x1));
n2 = nnz(~isnan(x2));
mean_x1 = mean(x1, 'omitnan');
mean_x2 = mean(x2, 'omitnan');
var_x1 = var(x1, 0, 'omitnan');
var_x2 = var(x2, 0, 'omitnan');


% Calculate pooled standard deviation (assuming two independent samples)
pooled_std = sqrt(((n1-1)*var_x1 + (n2-1)*var_x2)/(n1+n2-2));

%Calculate Cohen's d
d = (mean_x1 - mean_x2)/pooled_std;