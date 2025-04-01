function [cp_horz, p_horz, cp_vert, p_vert, cp_combined, p_combined, countingmatrix_custom] = crossValidate_multicoder(data, labels, varargin)
%
% [cp, p] = crossValidate(data, labels, varargin)
%
% crossValidate will load data from a set of session/run combinations
% (using the concatenateRuns.m function) and then cross validate those
% results using your method of choice
%
%
% vargin:
%   data =  nImages x (nPixels*yPixels)
%   labels = 2D array; x, y columns
%       For horizontal and vertical axes, use the following labels
%       1 = Negative (down or left)
%        2 = At center (along vertical or horizontal axis)
%        3 = Positive (up or right)
%
% optional varargin:
% 'verbose': boolean
%
% 'classificationMethod': string. Options are:
%   'CPCA+LDA'
%   'PCA+LDA'
%
% 'validationMethod': string. Options are:
%   'leaveOneOut'
%   'kFold' (note: defaults to 10-fold unless you pass 'K')
%
% 'K': integer, e.g. '10' for 10-fold validation
%
% 'testData' = nImages x (nPixels*yPixels)
% note: this MUST be the same size as 'data' and must use the same labels!
% this optional data varargin was created with the intention of passing
% data from a different time window for creating dynamic training/decodes.
%
% 'm' - Positive integer; determines subspace dimensionality for cPCA. Max
% is (number of classes - 1). Default is 1.
%
% varargout:
% cp, 'class performance', created using classperf
% note that, among other things, cp contains:
% confusion matrix: 'countingMatrix', e.g. [1 0 ; 0 1]
% correct rate, e.g. 0.95
%
% p: p value (based on binomial test), e.g. 0.05

%% handling varargin
p = inputParser;
p.addOptional('verbose',false,@islogical);
p.addOptional('validationMethod','kFold');
p.addOptional('K',10);
p.addOptional('classificationMethod','CPCA+LDA');
p.addOptional('testData',NaN);
p.addParameter('m', 1);
p.addParameter('variance_to_keep', 95);
p.parse(varargin{:});
sets = p.Results;

%% normalize to z score & get useful vars
display_progress = sets.verbose;
sets.testData = zscore(sets.testData);
N = size(data, 1);

cp_horz = classperf(labels(~isnan(labels(:,1)),1));     % initializing class performance var
cp_vert = classperf(labels(~isnan(labels(:,2)),2));     % initializing class performance var

unique_horz_classes = unique(labels(~isnan(labels(:, 1))));
unique_vert_classes = unique(labels(~isnan(labels(:, 2))));

possible_combined_labels = length(unique_horz_classes) * length(unique_vert_classes);

countingmatrix_custom = zeros(possible_combined_labels);

labels_combined = labels(:,1) + length(unique(labels(~isnan(labels(:,1)),1))) * (labels(:,2)-1);

if all(~isnan(labels_combined))
    cp_combined = classperf(labels_combined);
end

%% k-fold & leave one out cross validation
    
    % leave one out is just a special case of K-Fold wehre K=N;
    if strcmp(sets.validationMethod,'leaveOneOut')
        sets.K = N;
        display_progress  = true;
    end
    
    % cross validate here
    indices = crossvalind('Kfold', N, sets.K);
    % for each k-fold
    for i = 1:sets.K
        % create indices for training set & test set
        test = (indices == i);
        train = ~test;
        nan_ind = isnan(labels);
        
        % Accounting for NaN labels in case we are not using all possible
        % classes of data.
        horz_train = train & ~nan_ind(:,1);
        vert_train = train & ~nan_ind(:,2);
        horz_test = test & ~nan_ind(:,1);
        vert_test = test & ~nan_ind(:,2);

        % Fit and apply z-scoring to train data, then apply to test.
        [data(train, :), mu, sigma] = zscore(data(train, :));
        data(test, :) = (data(test, :) - mu) ./ sigma;
        
        %If nans, then define as 0. NaN values will break downstream
        %functions.
        data(isnan(data)) = 0;
        
        % classify
        if ~isnan(sets.testData)
            class_horz = classifyDoppler(data(horz_train,:), ...
                labels(horz_train,1), ...
                sets.testData(horz_test,:), ...
                'method', sets.classificationMethod, ...
                'm', sets.m, ...
                'variance_to_keep', sets.variance_to_keep);
            class_vert = classifyDoppler(data(vert_train,:), ...
                labels(vert_train,2) , ...
                sets.testData(vert_test,:), ...
                'method', sets.classificationMethod, ...
                'm', sets.m, ...
                'variance_to_keep', sets.variance_to_keep);
        else
            class_horz = classifyDoppler(data(horz_train,:), ...
                labels(horz_train,1), ...
                data(horz_test,:), ...
                'method', sets.classificationMethod, ...
                'm', sets.m, ...
                'variance_to_keep', sets.variance_to_keep);
            class_vert = classifyDoppler(data(vert_train,:), ...
                labels(vert_train,2), ...
                data(vert_test,:), ...
                'method', sets.classificationMethod, ...
                'm', sets.m, ...
                'variance_to_keep', sets.variance_to_keep);
        end
        if all(~isnan(labels_combined))
            class_combined = class_horz + length(unique(labels(~isnan(labels(:,1)),1))) * (class_vert-1);
        end
        classperf(cp_horz, class_horz, horz_test(~isnan(labels(:,1))));
        classperf(cp_vert, class_vert, vert_test(~isnan(labels(:,2))));
        if all(~isnan(labels_combined))
            classperf(cp_combined, class_combined, test);
            
            countingmatrix_custom = update_counting_matrix(countingmatrix_custom, class_combined, labels_combined(test));
        end
        if display_progress
            if i == 1
                f = waitbar(0, 'Leave-one-out analysis - Please wait...');
            else
                waitbar(i/sets.K, f);
            end
        end
    end
    if display_progress
        close(f);
    end


%% calculate classification accuracy measures
% May be slight mismatch between percentCorrect and nCorrect/nCounted due
% to `inconclusive` samples. MATLAB ignores them. We keep them in nCounted.
percentCorrect_horz = cp_horz.CorrectRate*100;
nCorrect_horz = sum(diag(cp_horz.CountingMatrix));
nCounted_horz = sum(cp_horz.CountingMatrix(:));
chance_horz = 1/length(cp_horz.ClassLabels);
p_horz = binomialTest(nCorrect_horz, nCounted_horz, chance_horz, 'one');

percentCorrect_vert = cp_vert.CorrectRate*100;
nCorrect_vert = sum(diag(cp_vert.CountingMatrix));
nCounted_vert = sum(cp_vert.CountingMatrix(:));
chance_vert = 1/length(cp_vert.ClassLabels);
p_vert = binomialTest(nCorrect_vert, nCounted_vert, chance_vert, 'one');

percentCorrect_combined = cp_combined.CorrectRate*100;
nCorrect_combined = sum(diag(cp_combined.CountingMatrix));
nCounted_combined = sum(cp_combined.CountingMatrix(:));
num_classes_combined = (length(cp_vert.ClassLabels)*length(cp_horz.ClassLabels));
chance_combined = 1/num_classes_combined;

% For directional tuning paper, using modified conservative estimate of
% chance level due to class imbalance, i.e., ideal chance is 11.11% but due
% to no middle category, chance is actually closer to 12.5%
chance_combined = 1/length(cp_combined.ClassLabels);

p_combined = binomialTest(nCorrect_combined, nCounted_combined, chance_combined, 'one');

%% display measures if verbose is on
if sets.verbose
    % classification accuracy (%)
    fprintf('\nHorizontal Classification Accuracy: \n%i / %i trials correctly classified (%2.2f%% correct)\t',...
        nCorrect_horz, nCounted_horz, percentCorrect_horz)
    % p-value
    if p_horz<0.0001
        fprintf('(binomial test: p < 0.001)\n')
    else
        fprintf('(binomial test: p = %1.3f)\n', p_horz)
    end
    
    % classification accuracy (%)
    fprintf('\nVertical Classification Accuracy: \n%i / %i trials correctly classified (%2.2f%% correct)\t',...
        nCorrect_vert, nCounted_vert, percentCorrect_vert)
    % p-value
    if p_vert<0.0001
        fprintf('(binomial test: p < 0.001)\n')
    else
        fprintf('(binomial test: p = %1.3f)\n', p_vert)
    end
    
    % classification accuracy (%)
    fprintf('\nCombined Classification Accuracy: \n%i / %i trials correctly classified (%2.2f%% correct)\t',...
        nCorrect_combined, nCounted_combined, percentCorrect_combined)
    % p-value
    if p_vert<0.0001
        fprintf('(binomial test: p < 0.001)\n')
    else
        fprintf('(binomial test: p = %1.3f)\n', p_combined)
    end
end

end
