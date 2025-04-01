function plt = imagesc_confusion_matrix(counting_matrix, varargin)
%% Function to create a simple confusion matrix
% X-axis is predicted value
% Y-axis is true value


p = inputParser;
p.addOptional('label_strings', []);
p.addOptional('row_normalize', true);
p.addOptional('title', []);
p.addOptional('FontColor', 'k');
p.addOptional('superimpose_text', '');
p.addOptional('colorbar_limits', [0 100]);
p.addOptional('colormap', inferno);
p.addOptional('show_colorbar', true);
p.addOptional('tiledlayout', []);
p.parse(varargin{:});
inputs = p.Results;

if ~isempty(inputs.tiledlayout)
    nexttile(inputs.tiledlayout);
end

n_true_classes = size(counting_matrix, 1);
n_pred_classes = size(counting_matrix, 2);

% row normalize
if inputs.row_normalize
    counting_matrix_plot = counting_matrix;
    counting_matrix_plot = counting_matrix_plot./sum(counting_matrix_plot, 2);
    counting_matrix_plot(isnan(counting_matrix_plot)) = 0;
else
    counting_matrix_plot = counting_matrix;
end

plt = imagesc(counting_matrix_plot*100);
yticks(1:n_true_classes);
xticks(1:n_pred_classes);
if inputs.show_colorbar
    clb = colorbar;
    clb.Label.String = 'Percent (row-normalized)';
    clb.Color = inputs.FontColor;
end
colormap(gca, inputs.colormap);
caxis(inputs.colorbar_limits);
if ~isempty(inputs.label_strings)
    yticklabels(inputs.label_strings(1:n_true_classes));
    xticklabels(inputs.label_strings(1:n_pred_classes));
end
xlabel('Predicted class');
ylabel('True class');
title(inputs.title, 'Color', inputs.FontColor);
set(gca, 'XColor', inputs.FontColor);
set(gca, 'YColor', inputs.FontColor);
text_color = inputs.FontColor;

if inputs.superimpose_text
    for i = 1:n_true_classes
        for j = 1:n_pred_classes
            count = counting_matrix(i, j); 
            row_total = sum(counting_matrix(i, :));
            switch inputs.superimpose_text
                case 'count'
                    label = strcat(num2str(count), '/', num2str(row_total));
                case 'percentage'
                    label = sprintf('%0.1f%%', 100*counting_matrix_plot(i, j));
                    if counting_matrix_plot(i, j) < .5
                        text_color = 'w';
                    else
                        text_color = 'k';
                    end
            end
            
            text(j,i,label, ...
                'HorizontalAlignment', 'center', ...
                'Color', text_color)
        end
    end   
end

axis image;

