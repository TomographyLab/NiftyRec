function [sinogram,partial] = et_project_partial(activity, cameras, attenuation, psf, use_gpu, rotation_of_partial_integrals, background, background_image, background_attenuation, truncate_negative_values)
%ET_PROJECT
%    Experimental - Projector for Emission Tomographic reconstruction. Returns projection and
%    partial line integrals. 
%
%Description
%    Function for projection of activity into detector space.
%
%    [SINOGRAM,PARTIAL] = ET_PROJECT_PARTIAL(ACTIVITY, CAMERAS, ATTENUATION, PSF, USE_GPU, ROTATION_OF_PARTIAL_INTEGRALS, BACKGROUND, BACKGROUND_IMAGE, BACKGROUND_ATTENUATION, TRUNCATE_NEGATIVE_VALUES)
%
%    SINOGRAM is the projection 
%
%    PARTIAL are the partial line intergrals. It's a volume of size
%    [Nx,Ny,Nz,Ncameras] where Nx,Ny,Nz is the size of the activity and
%    attenuation volumes, Ncameras is the number of camera positions. 
%    Each 3D volume [:,:,:,i] contains the partial line integrals for the
%    i-th camera. 
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
%    ROTATION_OF_PARTIAL_INTEGRALS if enabled the partial integrals are rotated 
%    back to the image frame. Defailts to 0 (disabled). 
%
%    BACKGROUND is the value the background is set to when performing rotation. 
%    It defaults to 0. 
%
%    BACKGROUND_IMAGE Background light pattern. This parameters is intended e.g. for coded aperture imaging. 
%    It's an image of size [Nx,Ny]. If it is set to a scalar value then it is ignored. 
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
%   N = 128;
%   use_gpu = 1;
%   activity = ones(N,N,N);
%   attenuation = zeros(N,N,N);
%   PSF = ones(7,7,N);
%   cameras = [0:pi/100:pi]';
%   [sinogram,partial] = et_project_partial(activity,cameras,attenuation,PSF,use_gpu);
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

if not(exist('psf','var'))
    psf = 0;
end

if not(exist('use_gpu','var'))
    use_gpu = 0;
end
   
if not(exist('rotation_of_partial_integrals','var'))
    rotation_of_partial_integrals = 0;
end
 
if not(exist('background','var'))
    background = 0;
end

if not(exist('background_image','var'))
    background_image = 0;
end

if not(exist('background_attenuation','var'))
    background_attenuation = 0;
end

if not(exist('truncate_negative_values','var'))
    truncate_negative_values = 1;
end

[sinogram,partial] = et_project_partial_mex(activity, cameras, attenuation, psf, use_gpu, rotation_of_partial_integrals, background, background_image, background_attenuation, truncate_negative_values);

