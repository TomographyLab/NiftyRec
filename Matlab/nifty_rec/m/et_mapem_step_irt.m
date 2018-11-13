function [activity_new, projection, normalization, update] = et_mapem_step_irt(activity_old, normalization, sinogram, cameras, attenuation, psf, beta, gradient_prior, epsilon)
%ET_MAPEM_STEP_IRT
%    Step of Maximum A Posteriori iterative reconstsruction algorithm for Emission Tomography
%    This function is equivalent to 'et_mapem_step' but based on the IRT toolbox. 
%
%Description
%    This function computes an estimate of the activity, given the previous estimate and the gradient 
%    of the prior distribution.
%
%    [NEW_ACTIVITY, PROJECTION, NORMALIZATION, UPDATE] = ET_MAPEM_STEP_IRT(ACTIVITY, NORM, SINO, CAMERAS, ATTENUATION, PSF, BETA, GRAD_PRIOR, EPSILON)
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
%    EPSILON is a small value that is added to the projection in order to
%    avoid division by zero. It defaults to 1e-6.
%
%Reference
%    IRT reconstruction toolbox, Fessler. 
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
%   norm = et_backproject_irt(ones(N,N,N_cameras),cameras,attenuation,psf,USE_GPU);
%   activity = ones(N,N,N);  %initial activity estimate
%   for i=1:mlem_steps
%       activity = et_mapem_step_irt(activity, norm, sinogram, cameras, attenuation, psf, 0, 0, USE_GPU);
%       imagesc(activity(:,:,N/2)); pause(0.1)
%   end
%
%See also
%   ET_MAPEM_STEP, ET_PROJECT_IRT, ET_BACKPROJECT_IRT
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
if not(exist('espilon','var'))
    epsilon = 1e-6;
end
if not(exist('psf','var'))
    psf = 0;
end
if not(exist('attenuation','var'))
    attenuation = 0;
end

projection = et_project_irt(activity_old, cameras, attenuation, psf);
projection(projection<epsilon) = epsilon ;
update = et_backproject_irt(sinogram ./ projection, cameras, attenuation, psf);
activity_new = activity_old .* update;
normalization = normalization - beta * gradient_prior;
normalization(normalization<epsilon) = epsilon;
activity_new = activity_new ./ (normalization);

return 

