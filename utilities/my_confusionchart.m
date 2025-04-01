function plot_handle = my_confusionchart(confusionmat, varargin)

%% Create a confusion matrix according to WHitney's specifications



%% Parse variable input
p = inputParser;
addOptional(p, 'climits',  [0 100]);
addOptional(p, 'Normalization', 'row-normalized');
addOptional(p, 'colormap', inferno);
addOptional(p, 'labels', '');
addOptional(p, 'superimpose_text', 'count');
addOptional(p, 'FontSize', 16);
addOptional(p, 'FontColor', 'k');
addOptional(p, 'title', '');
parse(p, varargin{:});
result = p.Results;


%% Extract some variables
% How many classes are there?
nClasses_true = size(confusionmat, 1);
nClasses_predicted = size(confusionmat, 2);

%% Compute desired confusion matrix
confusionmat_norm = confusionmat;
switch result.Normalization
    case 'row-normalized'
        row_total = sum(confusionmat, 2);
        confusionmat_norm = confusionmat ./ row_total *100;
        colorbar_string = '% correct (row-normalized)';
    case 'absolute'
        result.climits = [0 max(confusionmat, [], 'all')];
                colorbar_string = 'count';
    case 'column-normalized'
        column_total = sum(confusionmat, 1);
        confusionmat_norm = confusionmat ./ column_total *100;
                colorbar_string = '% correct (column-normalized)';
    case 'total-normalized'
        confusionmat_norm = confusionmat./sum(confusionmat, 'all')*100;
        colorbar_string = '% correct (total-normalized)';
end

%% Create confusion matrix
plot_handle = imagesc(confusionmat_norm);
colormap(result.colormap);
caxis(result.climits);
cb = colorbar;
cb.Label.String = colorbar_string;
cb.Color = result.FontColor;

xticks(1:nClasses_predicted);
yticks(1:nClasses_true);
ylabel('True class');
xlabel('Predicted class');

%% Customize confusion matrix
if ~isempty(result.labels)
   tickLabels = result.labels;
   xticklabels(tickLabels(1:nClasses_predicted));
   yticklabels(tickLabels(1:nClasses_true));
end

if result.superimpose_text
    for i = 1:nClasses_true
        for j = 1:nClasses_predicted
            count = confusionmat(i, j); 
            row_total = sum(confusionmat(i, :));
            switch result.superimpose_text
                case 'count'
                    label = strcat(num2str(count), '/', num2str(row_total));
                case 'percentage'
                    label = sprintf('%0.1f%%', confusionmat_norm(i, j));
            end
            if confusionmat_norm(i, j) < 50
                text_color = 'w';
            else
                text_color = 'k';
            end
            text(j,i,label, ...
                'HorizontalAlignment', 'center', ...
                'Color', text_color)
        end
    end   
end

title(result.title);
set(gca, 'FontSize',result.FontSize);
set(gca, 'XColor', result.FontColor);
set(gca, 'YColor', result.FontColor);
axis square;