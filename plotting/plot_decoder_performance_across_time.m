function plot_decoder_performance_across_time(trial_time, decoder_performance, pvalues, chance_level, mov_start_time, varargin)
%% Function that creates plot showing decoder performance across trial time
% Inputs:
%   trial_time          - (n x 1) Vector; Time for each entry in `decoder_performance`
%   decoder_performance - (n x 1) Vector; Percent correct (or other metric) at each
%                         timepoint
%   pvalues             - (n x 1) Vector; Probability of outcome for each
%                         timepoint
%   mov_start_time      - Scalar; Indicates when the movement start period
%                         began
%   chance_level        - Scalar; What is the expected chance performance
%                         rate?
%   varargin
%       p_threshold     - Double; What threshold do you want to consider
%                         statistically significant?
%       title           - Plot title
%
% Extracted into function form by Whitney on 2022/02/10


p = inputParser;
p.addOptional('p_threshold', 0.05);
p.addOptional('title', 'Decoder Performance across trial time');
p.addOptional('tiledlayout', []);
p.addOptional('ylim', [0 100]);
p.addOptional('color_line_by_pvalue', false);
p.addOptional('show_colorbar', false);
p.parse(varargin{:});
inputs = p.Results;

ylimits = inputs.ylim;

if ~isempty(inputs.tiledlayout)
   nexttile(inputs.tiledlayout);
end

% shade fixation & memory period
hFix = patch([trial_time(1) trial_time(1) 0 0], [0 100 100 0], ...
    'k', 'facecolor','b','edgecolor','none');
hFix.FaceAlpha = 0.2;
hMem = patch([0 0 mov_start_time mov_start_time], [0 100 100 0], ...
    'k', 'facecolor','g','edgecolor','none');
hMem.FaceAlpha = 0.2;

hold on
if inputs.color_line_by_pvalue
    % When you use this option, the resulting .svg will be an image instead
    % of line graphics.
    % Can fix this by using figure('Renderer', 'Painters')
    x = trial_time;
    y = decoder_performance;
        
   color_gradient = log10(pvalues);
   color_gradient(color_gradient<-6 | isinf(color_gradient)) = -5; % This is the limit of our bootstrap for MAAE
   surface([x;x], [y;y], [color_gradient;color_gradient], ...
       'facecol', 'no', ...
       'edgecol', 'interp', ...
       'linew', 4);
   
   % Define colormap
   plasma_colormap = plasma(256);
   cut_plasma = plasma_colormap(1:200, :);
   colormap(gca, flipud(cut_plasma));
   ax = gca;
   ax.CLim = [-5 0];
   
    if inputs.show_colorbar       
        cb = colorbar('Location', 'east');
        cb.Label.String = 'log(pvalue)';

    end
else
    plot(trial_time, decoder_performance, '-k', ...
        'LineWidth',4);
end
hold off

% handle x & y axes
% This extra `axis` command is needed otherwise it doesn't
% rescale properly otherwise.
axis([trial_time(1) trial_time(end) ylimits(1) ylimits(2)]);
set(gca,'Visible','off')
axes('Position',get(gca,'Position'),...
    'XAxisLocation','bottom',...
    'YAxisLocation','left',...
    'Color','none',...
    'XColor','k','YColor','k',...
    'LineWidth',2,...
    'TickDir','out');
axis([trial_time(1) trial_time(end) ylimits(1) ylimits(2)]);

% add labels
yticks(0:25:ylimits(2))
xticks([0 mov_start_time])
xticklabels({'Cue','Go'})
ylabel('Accuracy (%)')

% add reference lines
hold on; plot([trial_time(1) trial_time(end)],[chance_level chance_level],'--k','LineWidth',2)
text(trial_time(end), chance_level + 5, 'chance', 'HorizontalAlignment', 'right', 'FontSize', 15);

% add scale bar
plot([mov_start_time+3 mov_start_time+4], [20 20], 'k', 'LineWidth',5)
text(mov_start_time+3, 25, '1 sec', 'HorizontalAlignment', 'left', 'FontSize', 15)

% plot p value & labels
if ~inputs.color_line_by_pvalue && any(pvalues<inputs.p_threshold/length(trial_time))
    plot(trial_time(pvalues<inputs.p_threshold/length(trial_time)), ylimits(2),'*k','MarkerSize',10)
end

% title
title(inputs.title);

% set bigger font size
set(gca,'FontSize', 15,'FontName','Times');


end