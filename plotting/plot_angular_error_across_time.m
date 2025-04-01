function plot_angular_error_across_time(trial_time, angular_error, pvalues, chance_level, mov_start_time, varargin)
%% Function that creates plot showing decoder performance across trial time
% Inputs:
%   trial_time          - (n x 1) Vector; Time for each entry in `decoder_performance`
%   Angular Error       - (n x 1) Vector; Mean angular error (rad)
%   pvalues             - (n x 1) Vector; Probability of outcome for each
%                         timepoint
%   mov_start_time      - Scalar; Indicates when the movement start period
%                         began
%   chance_level        - Scalar; What is the expected chance performance
%                         rate? (In radians)
%   varargin
%       p_threshold     - Double; What threshold do you want to consider
%                         statistically significant?
%       title           - Plot title
%       errorbars       - Plot errorbars if specified (in radians)
%
% Extracted into function form by Whitney on 2022/02/10
% Added errorbar option WG 2025/01/28


p = inputParser;
p.addOptional('p_threshold', 0.05);
p.addOptional('title', 'Decoder Performance across trial time');
p.addOptional('tiledlayout', []);
p.addOptional('ylim', [0 pi]);
p.addOptional('color_line_by_pvalue', false);
p.addOptional('show_colorbar', false);
p.addOptional('errorbars', []);
p.parse(varargin{:});
inputs = p.Results;

if ~isempty(inputs.tiledlayout)
    nexttile(inputs.tiledlayout);
end


%% Convert rad to deg
angular_error = angular_error*360/(2*pi);
chance_level = chance_level*360/(2*pi);
ylimits = inputs.ylim*180/pi;

if ~isempty(inputs.errorbars)
    inputs.errorbars = inputs.errorbars * 180/pi;
end
%% Plot performance
% shade fixation & memory period
hFix = patch([trial_time(1) trial_time(1) 0 0], [0 180 180 0], ...
    'k', 'facecolor','b','edgecolor','none');
hFix.FaceAlpha = 0.2;
hMem = patch([0 0 mov_start_time mov_start_time], [0 180 180 0], ...
    'k', 'facecolor','g','edgecolor','none');
hMem.FaceAlpha = 0.2;

if inputs.color_line_by_pvalue
    % When you use this option, the resulting .svg will be an image instead
    % of line graphics.
    % Can fix this by using figure('Renderer', 'Painters')
    x = trial_time;
    y = angular_error;

    color_gradient = log10(pvalues);
    color_gradient(color_gradient<-6 | isinf(color_gradient)) = -5; % This is the limit of our bootstrap.

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
    if ~isempty(inputs.errorbars)
        hold on
        %  surface([x;x], [y;y], [color_gradient;color_gradient], ...
        % 'facecol', 'no', ...
        % 'edgecol', 'interp', ...
        % 'linew', 4);
        shadedErrorBar(x,y,inputs.errorbars, ...
            'patchSaturation', 0.1);
        % plot(x, y+inputs.errorbars, 'k');
        % plot(x, y-inputs.errorbars, 'k');
        hold off
    end
else
    if isempty(inputs.errorbars)

        plot(trial_time, angular_error, '-k', ...
            'LineWidth',4);
    else
        shadedErrorBar(trial_time, angular_error, inputs.errorbars);
    end
end

% handle x & y axes
axis([trial_time(1) trial_time(end) ylimits(1) ylimits(2)])
set(gca,'Visible','off')
axes('Position',get(gca,'Position'),...
    'XAxisLocation','bottom',...
    'YAxisLocation','left',...
    'Color','none',...
    'XColor','k','YColor','k',...
    'LineWidth',2,...
    'TickDir','out');
axis([trial_time(1) trial_time(end) ylimits(1) ylimits(2)])

% add labels
yticks(0:45:ylimits(2))
xticks([0 mov_start_time])
xticklabels({'Cue','Go'})
ylabel('Angular error (deg)')

% add reference lines
hold on; plot([trial_time(1) trial_time(end)],[chance_level chance_level],'--k','LineWidth',2)
text(trial_time(end), chance_level + 5, 'chance', 'HorizontalAlignment', 'right', 'FontSize', 15);

% add scale bar
plot([-3 -2], [22 22], 'k', 'LineWidth',5)
text(-3, 32, '1 sec', 'HorizontalAlignment', 'left', 'FontSize', 15)

% plot p value & labels
if ~inputs.color_line_by_pvalue && any(pvalues<inputs.p_threshold/length(trial_time))
    plot(trial_time(pvalues<inputs.p_threshold/length(trial_time)), ylimits(2),'*k','MarkerSize',10)
end

% title
title(inputs.title);

% set bigger font size
set(gca,'FontSize', 15,'FontName','Times')


end