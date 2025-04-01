function plot_angiogram(backgroundImage,X_img_mm, Z_img_mm, varargin)
%% Helper function to plot the angiogram with correct data aspect ratio
%
% Author: Whitney Griggs
% Date: March 18, 2025
 
%% Parse variable inputs
p = inputParser;
p.addOptional('title', 'Average Doppler Image');
p.addOptional('colormap', 'inferno');
p.addOptional('show_colormap', false);
p.addOptional('colorbar_title', '');
p.addOptional('data_aspect_ratio', [1 1 1]);
p.parse(varargin{:});
inputs = p.Results;

%% Plot the angiogram
imagesc(X_img_mm, Z_img_mm, backgroundImage);
daspect(inputs.data_aspect_ratio);
axis image;
xlabel('mm');
ylabel('mm');
ax = gca;
colormap(inputs.colormap);
title(inputs.title, 'interpreter', 'none');

if inputs.show_colormap
    cb = colorbar;
    if ~isempty(inputs.colorbar_title)
       cb.Label.String = inputs.colorbar_title;
    end
end