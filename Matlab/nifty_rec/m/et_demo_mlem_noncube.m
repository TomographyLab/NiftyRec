
% ET_DEMO_MLEM_noncube
%     NiftyRec Demo: MLEM SPECT reconstruction with non cubic volume
%
%See also
%     ET_DEMO_MLEM, ET_DEMO_OSEM, ET_DEMO_MAPEM_MRF, ET_MAPEM_STEP
%
% 
%Stefano Pedemonte
%Copyright 2009-2013 CMIC-UCL
%Gower Street, London, UK


%% Parameters
N          = 128;
m          = 64;
N_cameras  = 120;
psf        = ones(5,5,N);
N_counts   = 200e6;
iter_mlem  = 50;
GPU        = 1;

cameras    = zeros(N_cameras,3); 
cameras(:,2)  = linspace(0,2*pi,N_cameras); 


%% Simulate SPECT scan 
disp('Creating synthetic sinogram..');
phantom = et_spherical_phantom(N,m,N,N/8,30,10,N/4,N/3,N/2)+et_spherical_phantom(N,m,N,N/9,15,0,3*N/4,m/2,N/2);
attenuation = et_spherical_phantom(N,m,N,N/8,0.00002,0.00001,N/4,m/3,N/2);
ideal_sinogram = et_project(phantom, cameras, attenuation, psf, GPU);
ideal_sinogram = ideal_sinogram/sum(ideal_sinogram(:))*N_counts;
sinogram = et_poissrnd(ideal_sinogram);


%% Reconstruction:
%Compute normalization volume
disp('Computing normalization..');
norm = et_backproject(ones(N,m,length(cameras)), cameras, attenuation, psf, GPU) ;

%Reconstruction
disp('Reconstructing..');
hFig = figure(); set(hFig, 'Position', get(0,'ScreenSize')/2); 
activity = ones(N,m,N);
for i=1:iter_mlem
    fprintf('\nMLEM step: %d',i);
    activity = et_mapem_step(activity, norm, sinogram, cameras, attenuation, psf, 0, 0, GPU, 0, 0.0001);
    scale = 900;
    subplot(1,3,1); image(scale*flipud(squeeze(activity(:,:,floor(N/2)))')); colormap gray; axis equal tight; 
    subplot(1,3,2); image(scale*flipud(squeeze(activity(:,floor(m/2),:))')); colormap gray; axis equal tight; 
    subplot(1,3,3); image(scale*flipud(squeeze(activity(floor(N/4),:,:))')); colormap gray; axis equal tight; pause(0.1);    
end


%% Cleanup 
if GPU
    et_reset_gpu();
end

disp('MLEM Done.');
