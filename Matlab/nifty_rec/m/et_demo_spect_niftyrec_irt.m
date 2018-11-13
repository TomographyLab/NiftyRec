
% ET_DEMO_SPECT_NIFTYREC_IRT
%     NiftyRec Demo: Compare NiftyRec and IRT reconstructions. 
%     NiftyRec has a NiftyRec style interface to IRT to enable interoperability and comparison. 
%
%See also
%     ET_PROJECT_IRT, ET_BACKPROJECT_IRT, ET_PROJECT, ET_BACKPROJECT
%
% 
%Stefano Pedemonte
%Copyright 2009-2013 CMIC-UCL
%Gower Street, London, UK


N_x = 64;
N_y = 64;
N_z = 64;
N_projections = 120;
theta = 180;
N_counts = 50e6;
subset_order = 16; 
N_iter = 50;
use_gpu = 1;

cameras = zeros(N_projections,3);
cameras(:,2)=(pi/180)*(0:theta/N_projections:theta-theta/N_projections);

% Point Spread Function ___________________________________________________
dx = 2.46;               % pixel size in y direction zubal
d  = 1.5;                % d is the collimator hole aperture [mm]
l  = 35;                 % l is the septa length [mm]
h  = 100;                % h is the source-to-collimator distance
t  = 0.2;                % t is the septal thickness
intrinsic = 3.8;         % intrinsic is the FWHM of the intrinsic resolution
obj2det = (N_y/2)*dx+50;  %detector orbit radius of rotation
   
[FWHM slope efficiency] = et_collimator_parallelholes(d,l,h,t,intrinsic,obj2det);
psf = make_3s_psfs(N_y, dx , obj2det, intrinsic, slope);

% Phantom _________________________________________________________________
fprintf('Creating synthetic data.. \n');
radius_sphere_1   = N_x/6;
xyz_sphere_1      = [7*N_x/12, (N_z+1)/2, 7*N_y/12];
activity_sphere_1 = 1;
radius_sphere_2   = N_x/8;
xyz_sphere_2      = [1*N_x/3, (N_z+1)/2, 1*N_y/3];
activity_sphere_2 = 2;
phantom = et_spherical_phantom(N_x,N_z,N_y,radius_sphere_1,activity_sphere_1,0,xyz_sphere_1(1),xyz_sphere_1(2),xyz_sphere_1(3)) + ...
          et_spherical_phantom(N_x,N_z,N_y,radius_sphere_2,activity_sphere_2,0,xyz_sphere_2(1),xyz_sphere_2(2),xyz_sphere_2(3)); 
attenuation = 0;

sinogram = et_project(phantom,cameras,attenuation,psf,use_gpu);
sinogram_irt = et_project_irt(phantom,cameras,attenuation,psf);

sinogram = et_poissrnd(sinogram*N_counts/sum(sinogram(:)));
sinogram_irt = et_poissrnd(sinogram_irt*N_counts/sum(sinogram_irt(:)));

% Reconstruct _____________________________________________________________
fprintf('Reconstruction.. \n');
activity = ones(N_x,N_z,N_y);
activity_irt = ones(N_x,N_z,N_y);
norm = et_backproject(ones(N_x,N_z,N_projections), cameras, attenuation, psf, use_gpu) ;
norm_irt = et_backproject_irt(ones(N_x,N_z,N_projections), cameras, attenuation, psf) ;
for iter = 1:N_iter
    fprintf('OSEM step: %d/%d (subset order %d) \n',iter,N_iter,subset_order);
    tic; activity = et_mapem_step(activity, norm, sinogram, cameras, attenuation, psf, 0, 0, use_gpu);
    fprintf('  - NiftyRec iteration: %g \n', toc);
    tic; activity_irt = et_mapem_step_irt(activity_irt, norm_irt, sinogram_irt, cameras, attenuation, psf, 0, 0);    
    fprintf('  - IRT iteration:      %g \n', toc);    
    subplot(1,2,1); imagesc(reshape(activity(:,floor(N_z/2),:),N_x,N_y)); colormap gray; axis equal tight off; title('NiftyRec');
    subplot(1,2,2); imagesc(reshape(activity_irt(:,floor(N_z/2),:),N_x,N_y)); colormap gray; axis equal tight off; title('IRT'); pause(0.1)
end

if use_gpu
    et_reset_gpu();
end


