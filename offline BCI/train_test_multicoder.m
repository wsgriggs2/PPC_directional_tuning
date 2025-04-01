function [cp_horz, cp_vert, cp_combined, custom_counting_matrix] = train_test_multicoder(data, labels, train_ind, test_ind, varargin)
%% No cross-validation here. Train and test on differnet portions of data for multicoder architecture

%% handling varargin
p = inputParser;
p.addOptional('verbose',false,@islogical)
p.addOptional('validationMethod','kFold')
p.addOptional('K',10);
p.addOptional('N',175);
p.addOptional('classificationMethod','CPCA+LDA')
p.addOptional('testData',NaN)
p.addParameter('m', 1)
p.addParameter('variance_to_keep', 95);
p.parse(varargin{:});
sets = p.Results;

%%

if any(train_ind & test_ind)
    warning('There is a data leak. Testing on some of same samples as the training.');
end

train_data = data(train_ind, :);
test_data = data(test_ind, :);

% Z-score the data
[train_data, mu, sigma] = zscore(train_data);
test_data = (test_data - train_data) ./ sigma;

% Create the class performance variables
cp_horz = classperf(labels(:,1));     % initializing class performance var
cp_vert = classperf(labels(:,2));     % initializing class performance var

% Find how many unique classes
unique_horz_classes = unique(labels(~isnan(labels(:, 1))));
unique_vert_classes = unique(labels(~isnan(labels(:, 2))));
possible_combined_labels = length(unique_horz_classes) * length(unique_vert_classes);
num_unique_vert_classes = length(unique_vert_classes);

% Initialize custom counting matrix
custom_counting_matrix = zeros(possible_combined_labels);

% Create the combined labels and performance variable
labels_combined = labels(:,1) + length(unique_vert_classes) * (labels(:,2)-1);
if all(~isnan(labels_combined))
    cp_combined = classperf(labels_combined);
end

%% k-fold & leave one out cross validation

class_horz = classifyDoppler(train_data, ...
    labels(train_ind,1), ...
    test_data, ...
    'method', sets.classificationMethod, ...
    'm', sets.m, ...
    'variance_to_keep', sets.variance_to_keep);
class_vert = classifyDoppler(train_data, ...
    labels(train_ind,1), ...
    test_data, ...
    'method', sets.classificationMethod, ...
    'm', sets.m, ...
    'variance_to_keep', sets.variance_to_keep);

if all(~isnan(labels_combined))
    % Fix with dynamically assigned # of classes instead of 3
    class_combined = class_horz + num_unique_vert_classes * (class_vert-1);
end
classperf(cp_horz, class_horz, test_ind);
classperf(cp_vert, class_vert, test_ind);
if all(~isnan(labels_combined))
    classperf(cp_combined, class_combined, test_ind);
    
    custom_counting_matrix = update_counting_matrix(custom_counting_matrix, class_combined, labels_combined(test_ind));
end

percentCorrect_combined = cp_combined.CorrectRate*100;
nCorrect_combined = sum(diag(cp_combined.CountingMatrix));
nCounted_combined = sum(cp_combined.CountingMatrix(:));
chance_combined = 1/length(cp_combined.ClassLabels);
p_combined = binomialTest(nCorrect_combined, nCounted_combined, chance_combined, 'one');
