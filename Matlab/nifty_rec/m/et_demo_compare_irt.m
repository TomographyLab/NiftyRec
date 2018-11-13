
% ET_DEMO_COMPARE_IRT
%
%See also
%   ET_DEMO_MLEM
%
% 
%Stefano Pedemonte
%Copyright 2009-2013 CMIC-UCL
%Gower Street, London, UK


%% Parameters
N_x        = 256;
N_y        = 1; 
N_z        = 256;
N_cameras  = 240;
iter_mlem  = 100;
GPU        = 0;

%% Simulate SPECT scan 
cameras = zeros(N_cameras,3);
cameras(:,2)=(pi/180)*(0:360/N_cameras:360-360/N_cameras);

disp('Creating synthetic sinogram..');
mask = et_spherical_phantom(N_x,N_y,N_z,N_x*0.45,1,0,(N_x+1)/2,(N_y+1)/2,(N_z+1)/2);
N_spheres  = 8;
max_radius = 25;
min_radius = 1.5;
spiral_radius = N_x/3; 
spiral_center_x = N_x/2; 
spiral_center_z = N_z/2; 
phantom = zeros(N_x,N_y,N_z);
i=0;
for radius = linspace(min_radius,max_radius,N_spheres)
    theta = 2*pi / N_spheres * i;
    loc_sphere = [cos(theta),sin(theta);-sin(theta),cos(theta)]*spiral_radius*[0.5*sqrt(2),0.5*sqrt(2)]';
    phantom = phantom + et_spherical_phantom(N_x,N_y,N_z,radius,10,0.00001,loc_sphere(1)+spiral_center_x,1,loc_sphere(2)+spiral_center_z); 
    i=i+1;
end
phantom = mask.*phantom; 
attenuation = 0;

%% Reconstruction IRT, NO PSF:
figure; title('Reconstruction IRT, NO PSF');
psf = repmat([0,0,0;0,1,0;0,0,0],[1,1,N_x]);
sinogram = et_project_irt(phantom, cameras, attenuation, psf);
normalisation = et_backproject_irt(ones(N_x,N_y,length(cameras)), cameras, attenuation, psf) ;
disp('Reconstructing with IRT, NO PSF..');
activity_irt_nopsf = ones(N_x,N_y,N_z); 
log_likelihood_irt_nopsf = zeros(1,iter_mlem);
mse_irt_nopsf = zeros(1,iter_mlem); 
for i=1:iter_mlem
    fprintf('\n1/4 MLEM step: %d',i);
    [activity_irt_nopsf,projection] = et_mapem_step_irt(activity_irt_nopsf, normalisation, sinogram, cameras, attenuation, psf, 0, 0);
    log_likelihood_irt_nopsf(i) = et_log_likelihood(sinogram,projection); 
    mse_irt_nopsf(i) = norm(phantom(:)-activity_irt_nopsf(:));
    subplot(2,1,1); imagesc(squeeze(activity_irt_nopsf(:,ceil(N_y/2),:))); colormap gray; axis square; 
    subplot(2,1,2); hold off; plot(squeeze(phantom(N_x/2,ceil(N_y/2),:)),'r'); hold on; plot(squeeze(activity_irt_nopsf(N_x/2,ceil(N_y/2),:)),'b'); pause(0.2);
end
disp('Done');

%% Reconstruction NiftyRec, NO PSF:
figure; title('Reconstruction NiftyRec, NO PSF');
psf=0;
sinogram = et_project(phantom, cameras, attenuation, psf, GPU);
normalisation = et_backproject(ones(N_x,N_y,length(cameras)), cameras, attenuation, psf, GPU) ;
disp('Reconstructing with NiftyRec, NO PSF..');
activity_niftyrec_nopsf = ones(N_x,N_y,N_z);
log_likelihood_niftyrec_nopsf = zeros(1,iter_mlem);
mse_niftyrec_nopsf = zeros(1,iter_mlem); 
for i=1:iter_mlem
    fprintf('\n2/4 MLEM step: %d',i);
    [activity_niftyrec_nopsf,projection] = et_mapem_step(activity_niftyrec_nopsf, normalisation, sinogram, cameras, attenuation, psf, 0, 0);
    log_likelihood_niftyrec_nopsf(i) = et_log_likelihood(sinogram,projection); 
    mse_niftyrec_nopsf(i) = norm(phantom(:)-activity_niftyrec_nopsf(:));
    subplot(2,1,1); imagesc(squeeze(activity_niftyrec_nopsf(:,ceil(N_y/2),:))); colormap gray; axis square; 
    subplot(2,1,2); hold off; plot(squeeze(phantom(N_x/2,ceil(N_y/2),:)),'r'); hold on; plot(squeeze(activity_niftyrec_nopsf(N_x/2,ceil(N_y/2),:)),'b'); pause(0.2);
end
disp('Done');

%% Reconstruction IRT, PSF:
figure; title('Reconstruction IRT, with PSF');
psf = repmat(fspecial('gaussian',3,0.6),[1,1,N_x]);
sinogram = et_project_irt(phantom, cameras, attenuation, psf);
normalisation = et_backproject_irt(ones(N_x,N_y,length(cameras)), cameras, attenuation, psf) ;
disp('Reconstructing with IRT, with PSF..');
activity_irt_withpsf = ones(N_x,N_y,N_z);
log_likelihood_irt_withpsf = zeros(1,iter_mlem);
mse_irt_withpsf = zeros(1,iter_mlem); 
for i=1:iter_mlem
    fprintf('\n3/4 MLEM step: %d',i);
    [activity_irt_withpsf,projection] = et_mapem_step_irt(activity_irt_withpsf, normalisation, sinogram, cameras, attenuation, psf);
    log_likelihood_irt_withpsf(i) = et_log_likelihood(sinogram,projection); 
    mse_irt_withpsf(i) = norm(phantom(:)-activity_irt_withpsf(:));    
    subplot(2,1,1); imagesc(squeeze(activity_irt_withpsf(:,ceil(N_y/2),:))); colormap gray; axis square; 
    subplot(2,1,2); hold off; plot(squeeze(phantom(N_x/2,ceil(N_y/2),:)),'r'); hold on; plot(squeeze(activity_irt_withpsf(N_x/2,ceil(N_y/2),:)),'b'); pause(0.2);
end
disp('Done');

%% Reconstruction NiftyRec, PSF:
figure; title('Reconstruction NiftyRec, with PSF');
psf = repmat(fspecial('gaussian',3,0.6),[1,1,N_x]);
sinogram = et_project(phantom, cameras, attenuation, psf);
normalisation = et_backproject(ones(N_x,N_y,length(cameras)), cameras, attenuation, psf) ;
disp('Reconstructing with NiftyRec, with PSF..');
activity_niftyrec_withpsf = ones(N_x,N_y,N_z);
log_likelihood_niftyrec_withpsf = zeros(1,iter_mlem);
mse_niftyrec_withpsf = zeros(1,iter_mlem); 
for i=1:iter_mlem
    fprintf('\n4/4 MLEM step: %d',i);
    [activity_niftyrec_withpsf,projection] = et_mapem_step(activity_niftyrec_withpsf, normalisation, sinogram, cameras, attenuation, psf);
    log_likelihood_niftyrec_withpsf(i) = et_log_likelihood(sinogram,projection); 
    mse_niftyrec_withpsf(i) = norm(phantom(:)-activity_niftyrec_withpsf(:));     
    subplot(2,1,1); imagesc(squeeze(activity_niftyrec_withpsf(:,ceil(N_y/2),:))); colormap gray; axis square; 
    subplot(2,1,2); hold off; plot(squeeze(phantom(N_x/2,ceil(N_y/2),:)),'r'); hold on; plot(squeeze(activity_niftyrec_withpsf(N_x/2,ceil(N_y/2),:)),'b'); pause(0.2);
end
disp('Done');

%save reconstructions.mat activity_irt_nopsf activity_irt_withpsf activity_niftyrec_nopsf activity_niftyrec_withpsf phantom mask sinogram cameras; 
%save plots_data.mat log_likelihood_irt_nopsf log_likelihood_irt_withpsf log_likelihood_niftyrec_nopsf log_likelihood_niftyrec_withpsf mse_irt_nopsf mse_irt_withpsf mse_niftyrec_nopsf mse_niftyrec_withpsf; 
