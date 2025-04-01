function alignment_tform = align_angiogram(angiogram1, angiogram2)

angiogram1_zscore = zscore(angiogram1, 0, 'all');
angiogram2_zscore = zscore(angiogram2, 0, 'all');

% Align angiogram2 to angiogram1
app = image_alignment_app(angiogram1_zscore, angiogram2_zscore);
pauseCount=0;

%Mandatory wait period because other the code
%keeps running forward without waiting for user input
%from the manual image alignment GUI.
while isempty(app.alignment_tform)
    pause(1);
    pauseCount=pauseCount+1;
    if pauseCount == 24*60*60*7 % Wait up to 7 day
        error('Alignment_tform was not created within 7 days. Stopping code');
        break;
    end
end
alignment_tform = app.alignment_tform;
close(app.UIFigure);
clear app;
