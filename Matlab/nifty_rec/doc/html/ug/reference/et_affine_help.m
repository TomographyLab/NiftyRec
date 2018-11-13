%% et_affine
% Apply Affine Transformation to a 2D or 3D image. GPU accelerated.
%
%% Description
% |out_image = et_affine(in_image, affine, use_gpu, background)|
%
% in_image is a 2D or 3D matrix.
%
% |affine| specifies a |[4x4]| Affine Transformation Matrix.
%
% |use_gpu| is optional and it enables GPU acceleration if a compatible GPU 
% device is installed in the system. By default use_gpu is set to |0| (disabled).
%
% |background| is the value the background is set to when performing the affine transform.
% It defaults to |0|.
%
%% GPU acceleration
% If a CUDA compatible Grahpics Processing Unit (GPU) is installed, 
% the projection algorithm can take advantage of it. Set use_gpu parameter
% to 1 to enable GPU acceleration. If GPU acceleration is not available, 
% the value of the parameter is uninfluential.
%
%% Algorithm notes
% Affine is performed by transforming the support of the target image and resampling 
% the source image by trilinear interpolation. 
%
%% Example
N = 128;
use_gpu = 1;
background = 0;
in_image = ones(N,N,N);
affine = eye(4);
out_image = et_affine(in_image,affine,use_gpu,background);

%% See also
% <et_rotate_help.html |et_rotate|>, <et_convolve_help.html |et_convolve|>, 
% <et_project_help.html |et_project|>, <et_backproject_help.html |et_backproject|>,
% <et_set_gpu_help.html |et_set_gpu|>, <et_list_gpus_help.html |et_list_gpus|>,
% 
%
%% Source
% Stefano Pedemonte
%
% Copyright 2009-2010 CMIC-UCL.
% Gower Street, London, UK

