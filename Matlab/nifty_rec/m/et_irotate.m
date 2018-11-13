function out_image = et_irotate(in_image, rotation, center, use_gpu, background)
%ET_ROTATE
%    Rotation of 2D and 3D images. GPU accelerated. It performs the inverse rotation of et_rotate()
%    when given the same input parameters. 
%
%Description
%    OUT_IMAGE = ET_IROTATE(IN_IMAGE, ROTATION, CENTER, USE_GPU, BACKGROUND)
%
%    IN_IMAGE is a 2D or 3D matrix.
%
%    ROTATION specifies rotation in radians along x,y,z axis.
%
%    CENTER is the center of rotation (xc,yc,zc). The unit measure is pixels, but it does not 
%    necessarily need to be integer.
%
%    USE_GPU is optional and it enables GPU acceleration if a compatible GPU 
%    device is installed in the system. By default use_gpu is set to 0 (disabled).
%
%    BACKGROUND is the value the background is set to when performing rotation.
%    It defaults to 0.
%
%GPU acceleration
%    If a CUDA compatible Grahpics Processing Unit (GPU) is installed, 
%    the projection algorithm can take advantage of it. Set use_gpu parameter
%    to 1 to enable GPU acceleration. If GPU acceleration is not available, 
%    the value of the parameter is uninfluential.
%
%Algorithm notes
%    Rotation is performed by rotating back the support of the target image and resampling 
%    the source image by trilinear interpolation. 
%
%
%Example
%   N = 128;
%   use_gpu = 1;
%   background = 0;
%   in_image = ones(N,N,N);
%   rotation = [0.0, pi/8, pi/4];
%   center = [N/2, N/2, N/2];
%   rot_image = et_rotate(in_image,rotation,center,use_gpu,background);
%   out_image = et_irotate(rot_image,rotation,center,use_gpu,background);
%
%See also
%   ET_ROTATE
%
% 
%Stefano Pedemonte
%Copyright 2009-2013 CMIC-UCL.
%Gower Street, London, UK


if not(exist('use_gpu'))
    use_gpu = 0;
end
    
if not(exist('background'))
    background = 0;
end

out_image = et_irotate_mex(in_image, rotation, center, use_gpu, background);

