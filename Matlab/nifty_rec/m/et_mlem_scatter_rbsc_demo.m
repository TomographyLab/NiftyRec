
%ET_MLEM_SCATTER_RBSC_DEMO
%Demo Reconstruction Based Scatter Compensation (RBSC) MLEM


%% Parameters
cameras    = [0:360/120:360-360/120]'*pi/180;
psf        = ones(15,15,N);
N          = 128;
scatter = ones(N,N,length(cameras));

iter_mlem = 30;
GPU        = 1;

disp('Loading sinogram..');
%sinogram   = load_nii('TOMO_I123_EM001_DS.img');
%sinogram   = double(sinogram.img);
%for i=1:120
%    sinogram(:,:,i) = sinogram(:,:,i)';
%end

sinogram = et_poissrnd( abs(et_project(et_spherical_phantom(N,N,N,N/8,100,0,N/4,N/3,N/2), cameras, psf, GPU)));

%% Reconstruction:

%Create normalization volume
disp('Creating normalization..');
norm = et_backproject(ones(N,N,length(cameras)), cameras, psf, GPU) ;
activity = et_reorder_activity_in(ones(N,N,N));

%Reconstruction
disp('Reconstructing..');
for i=1:iter_mlem
    fprintf('\nMLEM step: %d',i);
    proj = et_project(activity, cameras, psf, GPU, 0) + scatter + 0.0001 ;
    update = et_backproject(sinogram ./ proj, cameras, psf, GPU, 0) ;
    activity = activity .* update;
    activity = activity ./ norm;
    a = et_reorder_activity_out(activity);
    imagesc(a(:,:,floor(N/4))); colormap gray; axis square; pause(0.2)
end
activity = et_reorder_activity_out(activity);

disp('Done');

