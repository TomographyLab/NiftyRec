%% et_backproject
% Backprojection function for Emission Tomographic reconstruction
%
%% Description
% Function for backprojection of activity from detector space.
%
% |image = et_backproject(sinogram,cameras,psf,use_gpu)|
%
% |sinogram| is a 2D or 3D sinogram.
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
n_cameras = 100;
use_gpu = 1;
sinogram = ones(N,N,n_cameras);
PSF = ones(3,3,n_cameras);
cameras = [0:pi/n_cameras:pi-pi/n_cameras]';
%image = et_backproject(sinogram,cameras,PSF,use_gpu);

%% See also
% <et_project_help.html |et_project|>, <et_mlem_reconstruct_help.html |et_mlem_reconstruct|>,
% <et_list_gpus_help.html |et_list_gpus|>, <et_set_gpu_help.html |et_set_gpu|>,

%% Source
% Stefano Pedemonte
%
% Copyright 2009-2010 CMIC, UCL
% Gower Street, London, UK

