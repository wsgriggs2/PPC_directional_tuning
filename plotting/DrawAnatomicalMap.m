function DrawAnatomicalMap(backgroundImage,X_img_mm, Z_img_mm)

%Easy way to plot anatomical image with millimeter measurements shown.


Ibg = flipud(backgroundImage); % background image
    
    
    % Match Foreground image scale to Background image scale
    Ibg = (Ibg - min(Ibg(:)))./(max(Ibg(:))-min(Ibg(:))); % normalise background data [0 1]
    Ibg = nthroot(Ibg,4);           % nonlinear scale between [0 1]
    Ibg = 2*Ibg-1;                  % scale to [-1 1]
    
    
    
    % Set figure background to black, removes tick labels, sets ticks to figure frame color
    framecolor = [204/256 204/256 204/256]; % color frame matlab figure (RGB)
    
    ca = gca;
    ca.Color = 'k';
    %ca.XColor = framecolor;
    %ca.YColor = framecolor;
    
    % Background display 60 dB range (transparency scale)
    set(ca,'ALim',[min(Ibg(:)) max(Ibg(:))]);
    set(ca,'CLim',[0.0001 0.01]); % display 60 dB range (CData)
    
    % Define image color mask (m x n x 3 array of RGB values)
    Mask = ones(size(Ibg,1),size(Ibg,2),3); % White color mask
    
    % Display Background with image object syntax (property name/property value pairs)
    % CData == 3D array of values specifying the color of each pixel of the image
    % AlphaData == m-by-n matrix of non-NaNvalues specifying the transparency of each pixel
    % AlphaDataMapping scaled == Scaling the elements of AlphaData torange between the min and max values of the axes ALim property
    imagesc('CData',Mask, 'AlphaData',Ibg,'AlphaDataMapping','scaled'); % here we defined the Background as a transparent image
    axis image
    
    % Display Foreground
    % CDataMapping == scales the values according to the values of the axes CLim property
    
    % fix the ticks
    [~, newTicksI] = find(mod(X_img_mm,1)==0);
    xticks(newTicksI);
    xticklabels(X_img_mm(newTicksI()))
    xlabel('mm');
    [~, newTicksI] = find(mod(Z_img_mm,1)==0);
    yticks(newTicksI);
    yticklabels(fliplr(Z_img_mm(newTicksI())))
    ylabel('mm');
    ax = gca;
    ax.FontSize = 18;
    ca.XColor = [0 0 0];
    ca.YColor = [0 0 0];
    
end