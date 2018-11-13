function [activity_new, projection, normalization, update] = et_osmapem_step(subset_order, activity_old, sinogram, cameras, attenuation, psf, beta, gradient_prior, GPU, background, background_attenuation, epsilon)

%ET_OSMAPEM_STEP
%    Step of Maximum A Posteriori Ordered Subset iterative reconstsruction algorithm for Emission Tomography
%
%Description
%    This function computes an estimate of the activity, given the previous estimate and the gradient 
%    of the prior distribution.
%
%    [NEW_ACTIVITY, PROJECTION, NORMALIZATION, UPDATE] = ET_OSMAPEM_STEP(SUBSET_ORDER, ACTIVITY, SINO, CAMERAS, ATTENUATION, PSF, BETA, GRAD_PRIOR, USE_GPU, BACKGROUND, BACKGOUND_ATTENUATION, EPSILON)
%
%    SUBSET_ORDER is the subset of cameras to be used for each iteration.
%    SUBSET_ORDER=8 means that 1/8 of the cameras are used, SUBSET_ORDER=16
%    means that 1/16 of the cameras are used. If the division is not exact,
%    the result is rounded.
%    
%    ATIVITY is a 2D or 3D matrix of activity, typically estimated in the
%    previous MAPEM step
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
%    PSF is a Depth-Dependent Point Spread Function.
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
%    BACKGROUND_ATTENUATION is the value the attenuation background is set to when performing rotation.
%    It defaults to 0.
%
%    EPSILON is a small value that is added to the projection in order to avoid division by zero.
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
%   n_cameras = 120;
%   mapem_steps = 20;
%   subset_order = 16;
%   USE_GPU = 1;
%   phantom = et_spherical_phantom(N,N,N,N/8,100,0);
%   attenuation = 0;
%   psf = ones(5,5,N);
%   cameras = [0:2*pi/n_cameras:2*pi-2*pi/n_cameras]';
%   sinogram = et_poissrnd(et_project(phantom,cameras,attenuation,psf,USE_GPU));
%   activity = ones(N,N,N);  %initial activity
%   for i=1:mapem_steps
%       fprintf('OSMAPEM Step: %d\n',i);
%       activity = et_osmapem_step(subset_order, activity, sinogram, cameras, attenuation, psf, 0, 0, USE_GPU);
%       imagesc(activity(:,:,N/2)); pause(0.1);
%   end
%
%See also
%   ET_MAPEM_STEP, ET_MLEM_DEMO
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


%Extract random cameras: select randomly the first camera, 
%then select successive camera by drawing from a Gaussian 
%centred 'subset_order' indexes away 
%and procede that way in the same direction until 'N_cameras/subset_order' 
%cameras are selected
N1 = size(activity_old,1);
N2 = size(activity_old,2);
N_cameras = size(cameras,1);
if subset_order<2
    subset_order=2;
    disp 'Warning: Minimun subset order is 2, using subset_order=2'
end
N_cameras_sub = round(N_cameras/subset_order);
cameras_indexes=zeros(1,N_cameras_sub);
cameras_indexes(1)=round(rand()*N_cameras);
for cam = 2:N_cameras_sub
    cameras_indexes(cam)=ceil(normrnd(cameras_indexes(1)+(cam-1)*subset_order,subset_order/2));
    while not(cameras_indexes(cam)>cameras_indexes(cam-1))
        cameras_indexes(cam)=ceil(normrnd(cameras_indexes(1)+(cam-1)*subset_order,subset_order/2));
    end

end
cameras_indexes=rem(cameras_indexes,N_cameras);
cameras_indexes(cameras_indexes<=0)=N_cameras;

%Random without replacement
%N_cameras_sub = round(N_cameras/subset_order);
%cameras_indexes = round(0.5 + (N_cameras-0.5) * rand(1,N_cameras_sub)); 
cameras_sub = cameras(cameras_indexes,:);

%compute sensitivity for the subset
normalization = et_backproject(ones(N1,N2,N_cameras_sub), cameras_sub, attenuation, psf, GPU, background, background_attenuation);

%project and backproject
projection = et_project(activity_old, cameras_sub, attenuation, psf, GPU, background, background_attenuation);
projection(projection<0) = 0 ;
update = et_backproject((sinogram(:,:,cameras_indexes)+epsilon) ./ (projection+epsilon), cameras_sub, attenuation, psf, GPU, background, background_attenuation);
activity_new = activity_old .* (update+epsilon);
normalization = (normalization - beta * gradient_prior);
normalization(normalization<0) = 0;
activity_new = activity_new ./ (normalization+epsilon);

return

