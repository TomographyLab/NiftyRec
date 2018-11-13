%% et_project
% Projection function for Emission Tomographic reconstruction
%
%% Description
% Function for projection of activity into detector space.
%
% |sinogram = et_project(activity,cameras,psf,use_gpu,background)|
%
% |activity| is a 2D or 3D matrix of activity.
%
% |cameras| specifies camera positions and it can be expressed in two forms: 
% a matrix of size |[n,3]| representing angular position of each camera 
% with respect of x,y,z axis; or a column vector of length |n| where for each 
% camera, only rotation along z axis is specified. This is the most common 
% situation in PET and SPECT where the gantry rotates along a single axis.
%
% |psf| is a Depth-Dependent Point Spread Function.
%
% |use_gpu| is optional and it enables GPU acceleration if a compatible GPU 
% device is installed in the system. By default use_gpu is set to 0 (disabled).
%
% |background| is the value the background is set to when performing rotation.
% It defaults to 0.
%
%% GPU acceleration
% If a CUDA compatible Grahpics Processing Unit (GPU) is installed, 
% the projection algorithm can take advantage of it. Set use_gpu parameter
% to 1 to enable GPU acceleration. If GPU acceleration is not available, 
% the value of the parameter is uninfluential.
%
%% Algorithm notes
% Rotation based projection algorithm with trilinear interpolation.
% Depth-Dependent Point Spread Function is applyed in the frequency domain.
%
%% Reference
% Pedemonte, Bousse, Erlandsson, Modat, Arridge, Hutton, Ourselin, 
% "GPU Accelerated Rotation-Based Emission Tomography Reconstruction", NSS/MIC 2010
%
%% Example
N = 128;
use_gpu = 1;
activity = ones(128,128,128);
PSF = ones(3,3,128);
cameras = [0:pi/100:pi]';
sinogram = et_project(activity,cameras,PSF,use_gpu);

%% See also
% <et_backproject_help.html |et_backproject|>, <et_mlem_reconstruct_help.html |et_mlem_reconstruct|>,
% <et_list_gpus_help.html |et_list_gpus|>, <et_set_gpu_help.html |et_set_gpu|>,

%% Source
% Stefano Pedemonte
%
% Copyright 2009-2010 CMIC, UCL
% Gower Street, London, UK


