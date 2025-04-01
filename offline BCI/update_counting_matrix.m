function updated_matrix = update_counting_matrix(counting_matrix, class_combined, test_labels)
% Hack for countingmatrix in classperformance
% In case of multicoder, some classes are never seen, and not sure how to
% get classperformance to handle this elegantly without labelling them as
% inconclusive

countingmatrix_temp = zeros(size(counting_matrix));

for i = 1:9 % Prediction
    for j = 1:9 % Truth
        countingmatrix_temp(i, j) = nnz(class_combined == i & test_labels == j);
    end
end

updated_matrix = counting_matrix + countingmatrix_temp;