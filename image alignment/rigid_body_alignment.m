function [rigidbody_tform, image2_tform, vertshift, horzshift, rotation, image2_transformed] ...
    = rigid_body_alignment(image1, image2)
% rigid_body_alignment   Rigid body alignment of two neurovascular images
%
% INPUTS:
%   image1:              2D matrix; Neurovascular map #1
%   image2:              2D matrix; Neurovascular map #2
%   
% OUTPUTS:
%   rigidbody_tform:     (3 x 3) double matrix; Transform matrix
%   image2_tform:        MATLAB format for transform matrix
%   vertshift:           Pixel shift - vertical
%   horzshift:           Pixel shift - horizontal
%   rotation:            Angular rotation
%   image2_transformed   Image2 with transforms applied


[optimizer,metric] = imregconfig('monomodal'); %Initialize the image alingment for monomodal image
rigidbody_tform = imregtform(image2, image1, 'rigid', optimizer,metric); %Find best alignment using image intensity
image2_tform = rigidbody_tform.T;
[vertshift, horzshift, rotation] = extractTFormInfo(image2_tform); %Extract info from transform matrix
image2_transformed = imwarp(image2, rigidbody_tform,'OutputView',imref2d(size(image1))); %Apply transform matrix to image

end


%% Local function
function [vertshift,horzshift,rotation_in_degrees] = extractTFormInfo(tform)
    %Extract info from MATLAB affine2d transform matrix
    vertshift = -tform(3,2);
    horzshift = tform(3,1);
    rotation_in_degrees = asind(tform(1,2));
end