function [activity_new, projection, normalization, update] = et_mapem_step(activity_old, normalization, sinogram, cameras, attenuation, psf, beta, gradient_prior, GPU, background, background_attenuation, epsilon)
%ET_MAPEM_STEP
%    Step of Maximum A Posteriori iterative reconstsruction algorithm for Emission Tomography
%
%Description
%    This function computes an estimate of the activity, given the previous estimate and the gradient 
%    of the prior distribution.
%
%    [NEW_ACTIVITY, PROJECTION, NORMALIZATION, UPDATE] = ET_MAPEM_STEP(ACTIVITY, NORM, SINO, CAMERAS, ATTENUATION, PSF, BETA, GRAD_PRIOR, USE_GPU, BACKGROUND, EPSILON)
%
%    ATIVITY is a 2D or 3D matrix of activity, typically estimated in the previous MAPEM step
%
%    NORM specifies a normalization volume, See examples.
%
%    SINO is a 2D or 3D sinogram.
%
%    CAMERAS specifies camera positions and it can be expressed in two forms: 
%    a matrix of size [n,3] representing angular position of each camera 
%    with respect of x,y,z axis; or a column vector of length n where for each 
%    camera, only rotation along z axis is specified. This is the most common 
%    situation in PET and SPECT where the gantry rotates along a single axis.
%
%    ATTENUATION is the attenuation map (adimensional). Refer to the programming manual. 
%    If it's set to a scalar then attenuation is not applied (faster).
%
%    PSF is a Depth-Dependent Point Spread Function. If it's set to a
%    scalar then the PSF i considered ideal (ray-casting). 
%
%    BETA parameter for sensitivity of the prior term
%
%    GRAD_PRIOR gradient of the prior probability distribution
%
%    USE_GPU is optional and it enables GPU acceleration if a compatible GPU 
%    device is installed in the system. By default use_gpu is set to 0 (disabled).
%
%    BACKGROUND is the value the background is set to when performing rotation.
%    It defaults to 0.
%
%    BACKGROUND_ATTENUATION is the value the attenuation background is set
%    to when performing rotation. 
%    It defaults to 0.
%
%    EPSILON is a small value that is added to the projection in order to
%    avoid division by zero. It defaults to 1e-6.
%
%GPU acceleration
%    If a CUDA compatible Grahpics Processing Unit (GPU) is installed, 
%    the projection algorithm can take advantage of it. Set use_gpu parameter
%    to 1 to enable GPU acceleration. If GPU acceleration is not available, 
%    the value of the parameter is uninfluential.
%
%Algorithm notes
%    Rotation based projection algorithm with trilinear interpolation.
%    Depth-Dependent Point Spread Function is applyed in the frequency domain.
%
%Reference
%    Pedemonte, Bousse, Erlandsson, Modat, Arridge, Hutton, Ourselin, 
%    "GPU Accelerated Rotation-Based Emission Tomography Reconstruction", NSS/MIC 2010
%
%Example
%   N = 128;
%   N_cameras = 120;
%   mlem_steps = 30;
%   USE_GPU = 1;
%   psf = ones(11,11,N);
%   cameras = linspace(0,pi,N_cameras)';
%   mask = et_spherical_phantom(N,N,N,N/2,1,0,(N+1)/2,(N+1)/2,(N+1)/2);
%   phantom = et_spherical_phantom(N,N,N,N/8,100,10,N/4,N/3,N/2) .* mask;
%   attenuation = et_spherical_phantom(N,N,N,N/8,0.002,0.001,N/4,N/3,N/2) .* mask;
%   sinogram = et_poissrnd(et_project(phantom,cameras,attenuation,psf,USE_GPU));
%   norm = et_backproject(ones(N,N,N_cameras),cameras,attenuation,psf,USE_GPU);
%   activity = ones(N,N,N);  %initial activity estimate
%   for i=1:mlem_steps
%       activity = et_mapem_step(activity, norm, sinogram, cameras, attenuation, psf, 0, 0, USE_GPU);
%       imagesc(activity(:,:,N/2)); pause(0.1)
%   end
%
%See also
%   ET_PROJECT, ET_BACKPROJECT, ET_MLEM_DEMO, ET_MAPEM_MRF_DEMO
%
% 
%Stefano Pedemonte
%Copyright 2009-2012 CMIC-UCL
%Gower Street, London, UK

if not(exist('beta','var'))
    beta = 0;
end
if not(exist('gradient_prior','var'))
    gradient_prior = 0;
end
if not(exist('GPU','var'))
    GPU = 0;
end
if not(exist('background','var'))
    background = 0;
end
if not(exist('background_attenuation','var'))
    background_attenuation = 0;
end
if not(exist('espilon','var'))
    epsilon = 1e-8;
end
if not(exist('psf','var'))
    psf = 0;
end
if not(exist('attenuation','var'))
    attenuation = 0;
end

projection = et_project(activity_old, cameras, attenuation, psf, GPU, background, background_attenuation);
projection(projection<0) = 0 ;
update = et_backproject((sinogram+epsilon) ./ (projection+epsilon), cameras, attenuation, psf, GPU, background, background_attenuation);
activity_new = activity_old .* (update+epsilon);
normalization = normalization - beta * gradient_prior;
normalization(normalization<0) = 0;
activity_new = activity_new ./ (normalization+epsilon);

return 

