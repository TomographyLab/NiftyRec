function projection = tt_project_perspective(attenuation, image_origin, detector_origin, detector_shape, psf, use_gpu, background)
%ET_PROJECT
%    Projection function for Transmission Tomography
%
%Description
%
%    PROJECTION = TT_PROJECT_PERSPECTIVE(ATTENUATION, IMAGE_ORIGIN, DETECTOR_ORIGIN, DETECTOR_SHAPE, PSF, USE_GPU, BACKGROUND)
%
%    ATTENUATION is a 3D matrix of attenuation.
%
%    IMAGE_ORIGIN is a matrix of size [n,3] representing the coordinates of the 
%    center of the attenuation image voxel at [1,1]. Refer to the Programming Manual and the exemples. 
%
%    DETECTOR_ORIGIN is a matrix of size [n,3] representing the coordinates of the 
%    detector, which can be different at each of the n projections. 
%
%    DETECTOR_SHAPE is a matrix of size [n,2] representing the size of the dectors. 
%
%    PSF is a Depth-Dependent Point Spread Function.
%
%    USE_GPU is optional and it enables GPU acceleration if a compatible GPU 
%    device is installed in the system. By default use_gpu is set to 0 (disabled).
%
%    BACKGROUND is the value the background is set to when performing rotation.
%    It defaults to 0.
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
%    Refer to the Programming Manual of NiftyRec.
%
%Reference
%    Pedemonte, Bousse, Erlandsson, Modat, Arridge, Hutton, Ourselin, 
%    "GPU Accelerated Rotation-Based Emission Tomography Reconstruction", NSS/MIC 2010
%
%Example
%
%See also
%   TT_BACKPROJECT, ET_PROJECT
%
% 
%Stefano Pedemonte
%Copyright 2009-2012 CMIC-UCL
%Gower Street, London, UK


if not(exist('psf'))
    psf = 0;
end

if not(exist('use_gpu'))
    use_gpu = 0;
end
    
if not(exist('background'))
    background = 0;
end

projection = tt_project_perspective_mex(attenuation, image_origin, detector_origin, detector_shape, psf, use_gpu, background);

