
% ET_DEMO_FBP
%     NiftyRec Demo: FBP SPECT reconstruction 
%
%See also
%   ET_DEMO_MLEM, ET_DEMO_OSEM, ET_DEMO_MAPEM_MRF
%
% 
%Stefano Pedemonte
%Copyright 2009-2013 CMIC-UCL
%Gower Street, London, UK


%% Parameters
N              = 128;
N_projections  = 120;
cameras        = linspace(0,2*pi,N_cameras)';
psf            = ones(5,5,N);
N_counts       = 50e6;

iter_mlem      = 30;
GPU            = 1;

phantom_type   = 3;  % 0 for 'brain FDG PET'; 1 for sphere in uniform background; 2 for spiral

%% Initialise plot
figure(1); subplot(2,3,1); image(zeros(N,N)); colormap gray; axis square off tight; 
figure(1); subplot(2,3,2); image(zeros(N,N)); colormap gray; axis square off tight; 
figure(1); subplot(2,3,3); image(zeros(N,N)); colormap gray; axis square off tight; 
figure(1); subplot(2,3,4); image(zeros(N,N)); colormap gray; axis square off tight; 
figure(1); subplot(2,3,5); image(zeros(N,N)); colormap gray; axis square off tight; 
figure(1); subplot(2,3,6); image(zeros(N,N)); colormap gray; axis square off tight; pause(0.1); 

%% Load or generate phantom 
disp('Creating synthetic sinogram..');
mask = et_spherical_phantom(N,N,N,N*0.45,1,0,(N+1)/2,(N+1)/2,(N+1)/2);
if phantom_type == 0
    phantom_file = 'activity_128.nii'; 
    if not(exist(phantom_file,'file'))
        fprintf('Synthetic data file %s cannot be found. \nPlease make sure that NiftyRec data is installed. \nIf NiftyRec was built and installed from the source code, \nrun CMake, enable the INCLUDE_DATA flag and run "make install".\n',phantom_file); 
        return;
    end
    phantom = et_load_nifti(phantom_file); 
    phantom = phantom.img .* mask;
elseif phantom_type == 1
    phantom = et_spherical_phantom(N,N,N,N/8,60,20,N/4,N/2,N/3) .* mask;
else
    phantom = et_spiral_phantom(N,N,N,60,20,8,2,10,40,2) .* mask;
end
attenuation = 0;
figure(1); subplot(2,3,1); imagesc(squeeze(phantom(:,floor(N/2),:))); colormap gray; axis square off tight;  pause(0.1); 

%% Simulate SPECT scan 
ideal_sinogram = et_project(phantom, cameras, attenuation, psf, GPU);
ideal_sinogram = ideal_sinogram/sum(ideal_sinogram(:))*N_counts;
sinogram = et_poissrnd(ideal_sinogram);

disp('Visualising sinogram..');
for i=1:N_projections  
    figure(1); subplot(2,3,2); image(0.3*squeeze(sinogram(:,:,i))'); colormap gray; axis square off tight;      
    if i<N/2
        figure(1); subplot(2,3,3); image(0.3*squeeze(sinogram(:,i,:))); colormap gray; axis square off tight;  
    end
    pause(0.05); 
end


%% Reconstruction - FBP:
disp('Filtering ..');
H = linspace(0, 2*pi, N/2);
H_ramp = [H H(end-1:-1:2)];
h_ramp = ifftshift( ifft( H_ramp ) ); h_ramp = h_ramp / max(h_ramp); 
h_ramp = h_ramp(N/2-10:N/2+10);
sinogram_filtered = sinogram; 
for cam = 1:N_projections
    for i =1:N
        sinogram_filtered(:,i,cam) = conv(sinogram(:,i,cam),h_ramp,'same'); 
    end
end

disp('Back-projecting ..');
tic; 
activity = et_backproject(sinogram_filtered, cameras, attenuation, psf, GPU); 
time_fbp = toc; 

scale = 0.004; 
figure(1); subplot(2,3,4); image(scale*squeeze(mask(:,floor(N/2),:).*activity(:,floor(N/2),:))); colormap gray; axis square off tight;     
figure(1); subplot(2,3,5); image(scale*squeeze(mask(floor(N/2),:,:).*activity(floor(N/2),:,:))); colormap gray; axis square off tight;
figure(1); subplot(2,3,6); image(scale*squeeze(mask(:,:,floor(N/2)).*activity(:,:,floor(N/2))')); colormap gray; axis square off tight; pause(0.05);     




%% Cleanup 
if GPU
    et_reset_gpu();
end

fprintf('FBP done in %d sec.\n',time_fbp);





