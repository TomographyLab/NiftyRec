
%% Test et_rotate

%% Parameters
N             = 128;
N_projections = 120;
theta_first   = 0.0;   
theta_last    = 2*pi; 
cameras      = zeros(N_projections,3);
cameras(:,1) = linspace(theta_first,theta_last,N_projections);
psf           = ones(5,5,N);

GPU           = 1;

phantom_type  = 1;  % 0 for 'brain FDG PET'; 1 for 'sphere in uniform background'

%% Make phantom
disp('Creating synthetic sinogram..');
mask = et_spherical_phantom(N,N,N,N*0.45,1,0,(N+1)/2,(N+1)/2,(N+1)/2);
if phantom_type == 0
    phantom = et_load_nifti('activity_128.nii'); 
    phantom = phantom.img .* mask;
else
    phantom = et_spherical_phantom(N,N,N,N/8,30,10,N/4,N/3,N/2) .* mask;
end
attenuation = et_spherical_phantom(N,N,N,N/8,0.00002,0.00001,N/4,N/3,N/2) .* mask;

%% Simulate SPECT scan 
% Project using et_project (fast): 
tic;
projection_1 = et_project(phantom,cameras,attenuation,psf,GPU); 
fprintf('et_project: %2.2f seconds\n',toc); 

% Project using et_rotate and et_project (slow): 
tic;
projection_2 = zeros(N,N,N_projections);
for i = 1:N_projections
    rotated_volume = et_rotate(phantom,cameras(i,:),[(N+1)/2,(N+1)/2,(N+1)/2],GPU); 
    projection_2(:,:,i) = et_project(rotated_volume,[0,0,0],attenuation,psf,GPU);
end
fprintf('et_rotate -> et_project: %2.2f seconds\n',toc); 

figure; 
for i=1:N_projections
    subplot(1,2,1); imagesc(projection_1(:,:,i)); axis square tight; title('et\_project()'); 
    subplot(1,2,2); imagesc(projection_2(:,:,i)); axis square tight; title('et\_rotate() -> et\_project()'); 
    pause(0.1); 
end



