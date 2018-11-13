%% et_rotate
% Rotation of 2D and 3D images. GPU accelerated.
%
%% Description
% |out_image = et_rotate(in_image, rotation, center, use_gpu, background)
%
% |in_image| is a 2D or 3D matrix.
%
% |rotation| specifies rotation in radians along |(x,y,z)| axis.
%
% |center| is the center of rotation |(xc,yc,zc)|. The unit measure is pixels, but it does not 
% necessarily need to be integer.
%
% |use_gpu| is optional and it enables GPU acceleration if a compatible GPU 
% device is installed in the system. By default use_gpu is set to 0 (disabled).
%
% |background| is the value the background is set to when performing rotation.
% It defaults to |0|.
%
%% GPU acceleration
% If a CUDA compatible Grahpics Processing Unit (GPU) is installed, 
% the projection algorithm can take advantage of it. Set use_gpu parameter
% to |1| to enable GPU acceleration. If GPU acceleration is not available, 
% the value of the parameter is uninfluential.
%
%% Algorithm notes
% Rotation is performed by rotating back the support of the target image and resampling 
% the source image by trilinear interpolation. 
%
%
%% Example
N = 128;
use_gpu = 1;
background = 0;
in_image = ones(N,N,N);
rotation = [0.0:pi/8:pi/4];
center = [N/2, N/2, N/2];
out_image = et_rotate(in_image,rotation,center,use_gpu,background);

%% See also
% <et_affine_help.html |et_affine|>,
% <et_project_help.html |et_project|>, <et_backproject_help.html |et_backproject|>,
% <et_set_gpu_help.html |et_set_gpu|>, <et_list_gpus_help.html |et_list_gpus|>,
%
% 
%% Source
% Stefano Pedemonte
%
% Copyright 2009-2010 CMIC-UCL.
% Gower Street, London, UK

