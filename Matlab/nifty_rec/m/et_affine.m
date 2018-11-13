function out_image = et_affine(in_image, affine, use_gpu, background)
%ET_AFFINE
%    Apply Affine Transformation to a 2D or 3D image. GPU accelerated.
%
%Description
%    OUT_IMAGE = ET_AFFINE(IN_IMAGE, AFFINE, USE_GPU, BACKGROUND)
%
%    IN_IMAGE is a 2D or 3D matrix.
%
%    AFFINE specifies a [4x4] Affine Transformation Matrix.
%
%    USE_GPU is optional and it enables GPU acceleration if a compatible GPU 
%    device is installed in the system. By default use_gpu is set to 0 (disabled).
%
%    BACKGROUND is the value the background is set to when performing the affine transform.
%    It defaults to 0.
%
%GPU acceleration
%    If a CUDA compatible Grahpics Processing Unit (GPU) is installed, 
%    the projection algorithm can take advantage of it. Set use_gpu parameter
%    to 1 to enable GPU acceleration. If GPU acceleration is not available, 
%    the value of the parameter is uninfluential.
%
%Algorithm notes
%    Affine is performed by transforming the support of the target image and resampling 
%    the source image by trilinear interpolation. 
%
%Example
%   N = 128;
%   use_gpu = 1;
%   background = 0;
%   in_image = ones(N,N,N);
%   affine = eye(4);
%   out_image = et_affine(in_image,affine,use_gpu,background);
%
%See also
%   ET_ROTATE, ET_CONVOLVE, ET_PROJECT, ET_BACKPROJECT,
%   ET_SET_GPU, ET_LIST_GPUS,
%
% 
%Stefano Pedemonte
%Copyright 2009-2010 CMIC-UCL.
%Gower Street, London, UK


if not(exist('use_gpu'))
    use_gpu = 0;
end
    
if not(exist('background'))
    background = 0;
end

out_image = et_affine_mex(in_image, affine, use_gpu, background);


