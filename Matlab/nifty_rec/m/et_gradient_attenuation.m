function gradient = et_gradient_attenuation(activity, sinogram, cameras, attenuation, psf, use_gpu, background, background_attenuation, truncate_negative_values)

%ET_GRADIENT_ATTENUATION
%    Computes gradient of the attenuation map. 
%
%Description
%    Computes gradient of the atteunation map. 
%
%    IMAGE = ET_GRADIENT_ATTENUATION(ACTIVITY, SINOGRAM, CAMERAS, ATTENUATION, PSF, USE_GPU, BACKGROUND, BACKGROUND_ATTENUATION, TRUNCATE_NEGATIVE_VALUES)
%
%    ACTIVITY is 2D or 3D activity (estimate)
%
%    SINOGRAM is a 2D or 3D sinogram.
%
%    CAMERAS specifies camera positions and it can be expressed in two forms: 
%    a matrix of size [n,3] representing angular position of each camera 
%    with respect of x,y,z axis; or a column vector of length n where for each 
%    camera, only rotation along z axis is specified. This is the most common 
%    situation in PET and SPECT where the gantry rotates along a single axis.
%
%    ATTENUATION specifies attenuation coefficient map (estimate). 
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
%    If 0 there is no truncati
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
%
%See also
%    ET_PROJECT, ET_MAPEM_STEP, ET_MLEM_DEMO
%    ET_LIST_GPUS, ET_SET_GPU
%
% 
%Stefano Pedemonte
%Copyright 2009-2012 CMIC-UCL.
%Gower Street, London, UK


if not(exist('attenuation'))
    attenuation = 0;
end

if not(exist('psf'))
    psf = 0;
end

if not(exist('use_gpu'))
    use_gpu = 0;
end

if not(exist('background'))
    background = 0;
end

if not(exist('background_attenuation'))
    background_attenuation = 0;
end

if not(exist('truncate_negative_values'))
    truncate_negative_values = 1;
end

gradient = et_gradient_attenuation_mex(activity, sinogram, cameras, attenuation, psf, use_gpu, background, background_attenuation, truncate_negative_values); 

