

% ET_DEMO_ESTIMATE_ATTENUATION
%     NiftyRec Demo: estimation of the attenuation map from the photon counts with the MLEM algorithm. 
%
%See also
%   ET_DEMO_MLEM, ET_DEMO_OSEM, ET_DEMO_MAPEM_MRF
%
% 
%Stefano Pedemonte
%Copyright 2009-2012 CMIC-UCL
%Gower Street, London, UK


%% Parameters
N_x          = 128;
N_y          = 128;
N_z          = 128;
N_cameras    = 120;
cameras      = zeros(N_cameras,3);
cameras(:,2) =(pi/180)*(0:360/N_cameras:360-360/N_cameras);
psf          = repmat(fspecial('gaussian',3,1),[1,1,N_x])
N_counts     = 200e6;
save_images  = 0;
iter_mlem    = 200;
subset_order = 8; 
GPU          = 1; 

%% Simulate SPECT scan 
disp('Creating synthetic sinogram..');
mask = et_spherical_phantom(N_x,N_z,N_y,N_x*0.45,1,0,(N_x+1)/2,(N_z+1)/2,(N_y+1)/2);

radius_sphere_1   = N_x/6;
xyz_sphere_1      = [7*N_x/12, (N_z+1)/2, 7*N_y/12];
activity_sphere_1 = 1;
radius_sphere_2   = N_x/8;
xyz_sphere_2      = [1*N_x/3, (N_z+1)/2, 1*N_y/3];
attenuation_sphere_2 = 0.02;      

phantom = mask.* (et_spherical_phantom(N_x,N_z,N_y,radius_sphere_1,activity_sphere_1,activity_sphere_1/5,xyz_sphere_1(1),xyz_sphere_1(2),xyz_sphere_1(3)) + ...
                  activity_sphere_1/10); 
attenuation_phantom = mask.* et_spherical_phantom(N_x,N_z,N_y,radius_sphere_2,attenuation_sphere_2,attenuation_sphere_2/10,xyz_sphere_2(1),xyz_sphere_2(2),xyz_sphere_2(3)); 

ideal_sinogram = et_project(phantom, cameras, attenuation_phantom, psf, GPU);
ratio = sum(ideal_sinogram(:))/N_counts;
ideal_sinogram = ideal_sinogram/ ratio;
sinogram = et_poissrnd(ideal_sinogram);

h1=figure();
h2=figure();
figure(h1); imagesc(reshape(phantom(:,floor(N_y/2),:),N_x,N_y)); colormap gray; axis equal tight off; 
if save_images; print('-dpng','-r200', sprintf('./video_spheres/phantom_activity_y.png')); end
figure(h2); imagesc(reshape(attenuation_phantom(:,floor(N_y/2),:),N_x,N_y)); colormap gray; axis equal tight off; pause(0.2); 
if save_images; print('-dpng','-r200', sprintf('./video_spheres/phantom_attenuation_y.png')); end
for i=1:N_cameras
    figure(h1); image(0.15*squeeze(sinogram(:,:,i))'); colormap gray; axis tight equal off; pause(0.1); 
    if save_images; print('-dpng','-r200', sprintf('./video_spheres/sinogram_%d.png',i)); end
end

%% Reconstruction:
disp('Reconstructing..');
attenuation_epsilon = 0.00001;
activity = ones(N_x,N_z,N_y);
attenuation = 0.01*ones(128,128,128);
for i=1:iter_mlem
    fprintf('\nMLEM step: %d',i);
    % estimate activity 
    activity = et_osmapem_step(subset_order, activity, sinogram, cameras, attenuation, psf, 0, 0, GPU, 0, 0.0001); 
    % estimate attenuation 
    l = 1./attenuation; 
    gradient = mask.*et_gradient_attenuation(mask.*activity,sinogram./(et_project(activity,cameras,attenuation,psf,1)+1e-30),cameras,attenuation,0,1);
    gradient_norm = mask.*et_gradient_attenuation(mask.*activity,1+0*sinogram,cameras,attenuation,0,1);
    gradient_norm(gradient_norm<1e-6) = 1e-6;
    a = gradient ./ gradient_norm;
    a(a<1e-6)=1e-6;
    l = l.*a;
    attenuation = mask./l;

    figure(h1); imagesc(reshape(activity(:,floor(N_y/2),:),N_x,N_y)); colormap gray; axis equal tight off; title('Activity');
    if save_images; print('-dpng','-r200', sprintf('./video_spheres/activity_%d.png',i)); end
    figure(h2); imagesc(reshape(attenuation(:,floor(N_y/2),:),N_x,N_y)); colormap gray; axis equal tight off; title('Attenuation'); pause(0.2); 
    if save_images; print('-dpng','-r200', sprintf('./video_spheres/attenuation_%d.png',i)); end
end
disp('Done');

if GPU
    et_reset_gpu();
end

