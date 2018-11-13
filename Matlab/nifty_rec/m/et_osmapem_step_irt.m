function [activity_new, projection, normalization, update] = et_osmapem_step_irt(subset_order, activity_old, sinogram, cameras, attenuation, psf, beta, gradient_prior,epsilon)

%ET_OSMAPEM_STEP_IRT
%    Step of Maximum A Posteriori Ordered Subset iterative reconstsruction algorithm for Emission Tomography
%    Identical to 'et_osmapem_step' but based on IRT. 
%
%Description
%    This function computes an estimate of the activity, given the previous estimate and the gradient 
%    of the prior distribution.
%
%    [NEW_ACTIVITY, PROJECTION, NORMALIZATION, UPDATE] = ET_OSMAPEM_STEP_IRT(SUBSET_ORDER, ACTIVITY, SINO, CAMERAS, ATTENUATION, PSF, BETA, GRAD_PRIOR)
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
%Reference
%    IRT reconstruction toolbox, Fessler. 
%
%Example
%   N = 128;
%   n_cameras = 120;
%   mapem_steps = 20;
%   subset_order = 16;
%   phantom = et_spherical_phantom(N,N,N,N/8,100,0);
%   attenuation = 0;
%   psf = ones(5,5,N);
%   cameras = [0:2*pi/n_cameras:2*pi-2*pi/n_cameras]';
%   sinogram = et_poissrnd(et_project(phantom,cameras,attenuation,psf,USE_GPU));
%   activity = ones(N,N,N);  %initial activity
%   for i=1:mapem_steps
%       fprintf('OSMAPEM Step: %d\n',i);
%       activity = et_osmapem_step_irt(subset_order, activity, sinogram, cameras, attenuation, psf, 0, 0);
%       imagesc(activity(:,:,N/2)); pause(0.1);
%   end
%
%See also
%   ET_OSMAPEM_STEP, ET_PROJECT_IRT, ET_BACKPROJECT_IRT
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
if not(exist('psf','var'))
    psf = 0;
end
if not(exist('attenuation','var'))
    attenuation = 0;
end
if not(exist('espilon','var'))
    epsilon = 1e-6;
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
normalization = et_backproject_irt(ones(N1,N2,N_cameras_sub), cameras_sub, attenuation, psf);

%project and backproject
projection = et_project_irt(activity_old, cameras_sub, attenuation, psf);
projection(projection<epsilon) = epsilon ;
update = et_backproject_irt(sinogram(:,:,cameras_indexes) ./ projection, cameras_sub, attenuation, psf);
activity_new = activity_old .* update;
normalization = (normalization - beta * gradient_prior);
normalization(normalization<epsilon) = epsilon;
activity_new = activity_new ./ (normalization);

return

