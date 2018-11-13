%% et_convolve
% Convolve multiple 1D arrays or 2D images with multiple kernels. FFT based. GPU accelerated.
%
%% Description
% |out_image = et_convolve(in_image, kernels, use_gpu)|
%
% |in_image| is a 2D or 3D matrix.
%
% |kernels| is a 2D or 3D matrix.
%
% |use_gpu| is optional and it enables GPU acceleration if a compatible GPU 
% device is installed in the system. By default use_gpu is set to 0 (disabled).
%
%% GPU acceleration
% If a CUDA compatible Grahpics Processing Unit (GPU) is installed, 
% the projection algorithm can take advantage of it. Set use_gpu parameter
% to 1 to enable GPU acceleration. If GPU acceleration is not available, 
% the value of the parameter is uninfluential.
%
%% Algorithm notes
% FFT based convolution.
%
%
%% Example
N = 128;
use_gpu = 1;
in_image = ones(N,N,N);
kernels = ones(3,3,N);
out_image = et_convolve(in_image,kernels,use_gpu);

%% See also
% <et_rotate_help.html |et_rotate|>, <et_affine_help.html |et_affine|>,
% <et_project_help.html |et_project|>, <et_backproject_help.html |et_backproject|>,
% <et_set_gpu_help.html |et_set_gpu|>, <et_list_gpus_help.html |et_list_gpus|>,
%
%
%% Source
% Stefano Pedemonte
%
% Copyright 2009-2010 CMIC-UCL.
% Gower Street, London, UK

