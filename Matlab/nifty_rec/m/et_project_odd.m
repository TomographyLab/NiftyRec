function sinogram = et_project_odd(activity, cameras, attenuation, psf, use_gpu, background, background_attenuation, truncate_negative_values)
%ET_PROJECT_ODD
%    Projection function for Emission Tomographic reconstruction. Odd
%    volume size
%
%Description
%    Function for projection of activity into detector space.
%
%    SINOGRAM = ET_PROJECT_ODD(ACTIVITY, CAMERAS, ATTENUATION, PSF, USE_GPU, BACKGROUND, BACKGROUND_ATTENUATION, TRUNCATE_NEGATIVE_VALUES)
%
%    ATIVITY is a 2D or 3D matrix of activity.
%
%    CAMERAS specifies camera positions and it can be expressed in two forms: 
%    a matrix of size [n,3] representing angular position of each camera 
%    with respect of x,y,z axis; or a column vector of length n where for each 
%    camera, only rotation along z axis is specified. This is the most common 
%    situation in PET and SPECT where the gantry rotates along a single axis.
%
%    ATTENUATION specifies attenuation coefficients. It must be the same size as ACTIVITY.
%
%    PSF is a Depth-Dependent Point Spread Function. 
%
%    USE_GPU is optional and it enables GPU acceleration if a compatible GPU 
%    device is installed in the system. By default use_gpu is set to 0 (disabled).
%
%    BACKGROUND is the value the background is set to when performing rotation. 
%    It defaults to 0. 
%
%    BACKGROUND_ATTENUATION is the value the attenuation background is set 
%    to when performing rotation. It defaults to 0.
%
%    TRUNCATE_NEGATIVE_VALUES defaults to 1. If 1 then negative results are truncated to 0. 
%    If 0 there is no truncation. 
%
%GPU acceleration
%    If a CUDA compatible Grahpics Processing Unit (GPU) is installed, 
%    the projection algorithm can take advantage of it. Set use_gpu parameter
%    to 1 to enable GPU acceleration. If GPU acceleration is not available, 
%    the value of the parameter is uninfluential.
%
%Algorithm notes
%    Rotation based projection algorithm with trilinear interpolation.
%    FFT-based Depth-Dependent Point Spread Function. 
%
%Reference
%    Pedemonte, Bousse, Erlandsson, Modat, Arridge, Hutton, Ourselin, 
%    "GPU Accelerated Rotation-Based Emission Tomography Reconstruction", NSS/MIC 2010
%
%Example
%   N = 127;
%   use_gpu = 1;
%   activity = ones(N,N,N);
%   attenuation = 0;
%   PSF = ones(7,7,N);
%   cameras = [0:pi/100:pi]';
%   sinogram = et_project_odd(activity,cameras,attenuation,PSF,use_gpu);
%
%See also
%   ET_PROJECT
%
% 
%Stefano Pedemonte
%Copyright 2009-2013 CMIC-UCL
%Gower Street, London, UK


if not(exist('attenuation','var'))
    attenuation = 0;
end

if not(exist('background_attenuation','var'))
    background_attenuation = 0;
end

if not(exist('psf','var'))
    psf = 0;
end

if not(exist('use_gpu','var'))
    use_gpu = 0;
end
    
if not(exist('background','var'))
    background = 0;
end

if not(exist('truncate_negative_values','var'))
    truncate_negative_values = 1;
end

if not(isscalar(activity))
    activity_even = zeros(size(activity,1)+1,size(activity,2),size(activity,3)); 
    activity_even(1:end-1,:,:) = activity; 
else
    activity_even = activity; 
end
if not(isscalar(attenuation))
    attenuation_even = zeros(size(attenuation,1)+1,size(attenuation,2),size(attenuation,3)); 
    attenuation_even(1:end-1,:,:) = attenuation; 
else
    attenuation_even = attenuation; 
end

if not(size(cameras,2)==3)
    fprintf('et_project_odd() expects the parameter "cameras" to be of size [n_cameras,3]\n');
    return;
end
centers = repmat([(size(activity_even,1)-2)/2, (size(activity_even,2)-1)/2, (size(activity_even,3)-2)/2],size(cameras,1),1); 
cameras = [cameras,centers];
sinogram = et_project_mex(activity_even, cameras, attenuation_even, psf, use_gpu, background, background_attenuation, truncate_negative_values);
sinogram = sinogram(1:end-1,:,:); 


