function [fisher,fisher_prior] = et_fisher_grid(activity, cameras, grid, attenuation, psf, use_gpu, epsilon, background, background_attenuation)

%ET_FISHER_GRID
%    Approximate the FIsher Information on a grid of points
%
%Description
%
%    FISHER = ET_FISHER_GRID(ACTIVITY, CAMERAS, GRID, ATTENUATION, PSF, USE_GPU, BACKGROUND, EPSILON, BACKGROUND, BACKGROUND_ATTENUATION)
%
%    ATIVITY is a 2D or 3D matrix of activity.
%
%    CAMERAS specifies camera positions and it can be expressed in two forms: 
%    a matrix of size [n,3] representing angular position of each camera 
%    with respect of x,y,z axis; or a column vector of length n where for each 
%    camera, only rotation along z axis is specified. This is the most common 
%    situation in PET and SPECT where the gantry rotates along a single axis.
%
%    GRID is the grid of points where the Fisher information matrix (FIM) is computed. 
%    The grid must have the same size as ACTIVITY. It is filled with 0
%    everywhere except for the voxels where the activity is considered to be uncertain. 
%    The FIM of the non-zero voxels will be computed. 
%    Each non-zero voxel has to contain a unique integer between 1 and the
%    number of non-zero voxels. Thi number will be the index of the FIM for
%    that voxel. Please refer to the example in the inline help. 
%
%    ATTENUATION specifies attenuation coefficients. It must be the same size as ACTIVITY.
%
%    PSF is a Depth-Dependent Point Spread Function. 
%
%    USE_GPU is optional and it enables GPU acceleration if a compatible GPU 
%    device is installed in the system. By default use_gpu is set to 0 (disabled).
%
%    EPSILON is small number. If not specified, it defaults to 1e-12. The
%    elements of the FIM are based on a summation of terms that contain the
%    inverse of the expected counts in each detector bin. The expected
%    counts are obtained internally by projecting ACTIVITY and are positive
%    values, however they might be 0, in which case the inverse diverges. 
%    Expected counts that are smaller then EPSILON are set to EPSILON. 
%
%    BACKGROUND is the value the background is set to when performing rotation. 
%    It defaults to 0. 
%
%    BACKGROUND_ATTENUATION is the value the attenuation background is set 
%    to when performing rotation. It defaults to 0.
%
%GPU acceleration
%    If a CUDA compatible Grahpics Processing Unit (GPU) is installed, 
%    the projection algorithm can take advantage of it. Set use_gpu parameter
%    to 1 to enable GPU acceleration. If GPU acceleration is not available, 
%    the value of the parameter is uninfluential.
%
%
%Example
%   N = 128;
%   use_gpu = 1;
%   activity = ones(N,N,N);
%   attenuation = zeros(N,N,N);
%   PSF = ones(7,7,N);
%   cameras = [0:pi/100:pi]';
%   grid = zeros(N,N,N);
%   counter = 0;
%   for i=4:12:128
%       for j=4:12:128
%           for k=4:12:128
%               counter=counter+1; 
%               grid(i,j,k)=counter;
%   end; end; end;             
%   fisher = et_fisher_grid(activity,cameras,grid,attenuation,PSF,use_gpu);
%   cov = inv(fisher);
%   var = reshape(diag(cov),11,11,11);
%   imagesc(var(:,:,5))
%
%See also
%   ET_FISHER_GRID_INVPROJECTION
%
% 
%Stefano Pedemonte
%Copyright 2009-2012 CMIC-UCL
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
    
if not(exist('epsilon'))
    epsilon = 0;
end

if not(exist('background'))
    background = 0;
end

if not(exist('background_attenuation'))
    background_attenuation = 0;
end

[fisher,fisher_prior] = et_fisher_grid_mex(activity, cameras, grid, attenuation, psf, use_gpu, epsilon, background, background_attenuation);

return 
