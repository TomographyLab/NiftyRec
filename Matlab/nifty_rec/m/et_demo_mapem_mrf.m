

% ET_DEMO_MAPEM_MRF
%     NiftyRec SPECT demo: Maximum a posteriori Expectation Maximisation reconstruction 
%     with regularisation based on quadratic Markov Random Field. 
%
%See also
%   ET_DEMO_MLEM, ET_MAPEM_STEP, ET_DEMO_OSEM
%
% 
%Stefano Pedemonte
%Copyright 2009-2012 CMIC-UCL
%Gower Street, London, UK


%% Parameters
N          = 128;
N_cameras  = 120;
cameras    = linspace(0,2*pi,N_cameras)';
psf        = ones(5,5,N);
N_counts   = 50e6;

beta_mrf   = 10;

iter_mapem = 30;
GPU        = 1;

%% Simulate SPECT scan 
disp('Creating synthetic sinogram..');
mask = et_spherical_phantom(N,N,N,N*0.45,1,0,(N+1)/2,(N+1)/2,(N+1)/2);
phantom = et_spherical_phantom(N,N,N,N/8,100,10,N/4,N/3,N/2) .* mask;
attenuation = et_spherical_phantom(N,N,N,N/8,0.00002,0.00001,N/4,N/3,N/2) .* mask;
ideal_sinogram = et_project(phantom, cameras, attenuation, psf, GPU);
ideal_sinogram = ideal_sinogram/sum(ideal_sinogram(:))*N_counts;
sinogram = et_poissrnd(ideal_sinogram);


%% Reconstruction:

%Create normalization volume
disp('Creating normalization..');
norm = et_backproject(ones(N,N,length(cameras)), cameras, attenuation, psf, GPU);

%Create kernel for Total Variation gradient
k1 = [-1,-1,-1;-1,-1,-1;-1,-1,-1];
k2 = [-1,-1,-1;-1,26,-1;-1,-1,-1];
k3 = [-1,-1,-1;-1,-1,-1;-1,-1,-1];
kernel = zeros(3,3,3);
kernel(:,:,1) = k1; kernel(:,:,2) = k2; kernel(:,:,3) = k3;

%Reconstruction
disp('Reconstructing..');
activity = ones(N,N,N);
for i=1:iter_mapem
    fprintf('\nMAPEM step: %d / %d',i,iter_mapem);
    prior_gradient = - convn(activity,kernel,'same');
    activity = et_mapem_step(activity, norm, sinogram, cameras, attenuation, psf, beta_mrf, prior_gradient, GPU, 0, 0.0001);
    imagesc(activity(:,:,floor(N/2))); colormap gray; axis square; pause(0.2)
end

disp('Done');

if GPU
    et_reset_gpu();
end

