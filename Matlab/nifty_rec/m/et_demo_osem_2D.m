
% ET_DEMO_OSEM_2D
%     NiftyRec Demo: Ordered Subsets Expectation Maximisation (OSEM) SPECT
%     reconstruction - 2D.
%
%See also
%   ET_DEMO_OSEM, ET_DEMO_MLEM, ET_DEMO_MAPEM_MRF, ET_OSMAPEM_STEP
%
% 
%Stefano Pedemonte
%Copyright 2009-2013 CMIC-UCL
%Gower Street, London, UK


%% Parameters
N          = 256;
N_cameras  = 180;
cameras = zeros(N_cameras,3);
cameras(:,2)=(pi/180)*(0:180/N_cameras:180-180/N_cameras);
psf        = ones(1,1,N);
N_counts   = 1e10/128;

iter_osem    = 300;
subset_order = 16;
GPU          = 0;

%% Simulate SPECT scan 
disp('Creating synthetic sinogram..');
load niftyrec_logo; 
phantom = reshape(niftyrec_logo,N,1,N);
attenuation = 0;
ideal_sinogram = et_project(phantom, cameras, attenuation, psf, GPU);
ideal_sinogram = ideal_sinogram/sum(ideal_sinogram(:))*N_counts;
sinogram =  et_poissrnd(ideal_sinogram);

figure(1); subplot(1,2,1); imagesc(squeeze(sinogram)); colormap gray; axis square tight off; title('Sinogram'); 
figure(1); subplot(1,2,2); image(zeros(N,N)); axis square tight off; title(sprintf('OSEM-%d reconstruction',subset_order)); pause(3);

%% Reconstruction:
disp('Reconstructing..');
activity = ones(N,1,N);
for i=1:iter_osem
    fprintf('\nOSEM step: %d',i);
    activity = et_osmapem_step(subset_order, activity, sinogram, cameras, attenuation, psf, 0, 0, GPU, 0, 0.0001);
    figure(1); subplot(1,2,2); imagesc(squeeze(activity(:,1,:))); colormap gray; axis square tight off; title(sprintf('OSEM-%d reconstruction - iteration %d',subset_order,i)); pause(0.1)
end
disp('Done');

if GPU
    et_reset_gpu();
end

