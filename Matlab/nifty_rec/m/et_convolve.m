function out_image = et_convolve(in_image, kernels, use_gpu)
%ET_CONVOLVE
%    Convolve multiple 1D arrays or 2D images with multiple kernels. FFT based. GPU accelerated.
%
%Description
%    OUT_IMAGE = ET_CONVOLVE(IN_IMAGE, KERNELS, USE_GPU)
%
%    IN_IMAGE is a 2D or 3D matrix.
%
%    KERNELS is a 2D or 3D matrix.
%
%    USE_GPU is optional and it enables GPU acceleration if a compatible GPU 
%    device is installed in the system. By default use_gpu is set to 0 (disabled).
%
%GPU acceleration
%    If a CUDA compatible Grahpics Processing Unit (GPU) is installed, 
%    the projection algorithm can take advantage of it. Set use_gpu parameter
%    to 1 to enable GPU acceleration. If GPU acceleration is not available, 
%    the value of the parameter is uninfluential.
%
%Algorithm notes
%    FFT based convolution.
%
%
%Example
%   N = 128;
%   use_gpu = 1;
%   in_image = ones(N,N,N);
%   kernels = ones(3,3,N);
%   out_image = et_convolve(in_image,kernels,use_gpu);
%
%See also
%   ET_AFFINE, ET_ROTATE, ET_PROJECT, ET_BACKPROJECT
%
% 
%Stefano Pedemonte
%Copyright 2009-2010 CMIC-UCL.
%Gower Street, London, UK


if not(exist('use_gpu'))
    use_gpu = 0;
end

out_image = et_convolve_mex(in_image, kernels, use_gpu);

