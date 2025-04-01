function [aligned_iDop, tform] = align_iDop(iDop1, iDop2)
% Align iDop2 to iDop1

dop1_angiogram = makeAngiogram(iDop1);
dop2_angiogram = makeAngiogram(iDop2);

% Use manual alignment to register between sessions/runs
alignment_tform = align_angiogram(zscore(dop1_angiogram), zscore(dop2_angiogram));
%Convert transform into matlab format
% Alignment_tform is passed back from ManualImageAlignment
tform = affine2d(alignment_tform);
% Align image to match previously loaded data
aligned_iDop = imwarp(...
    iDop2, ...
    tform, ...
    'OutputView', imref2d(size(dop1_angiogram)));
            