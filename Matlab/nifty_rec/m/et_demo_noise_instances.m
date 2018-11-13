

% ET_DEMO_NOISE_INSTANCES
%     NiftyRec SPECT demo: Reconstruction of multiple instances of nosie (2D) in parallel 
%     on the Graphics Processing Unit. 
%
%See also
%   ET_DEMO_MLEM, ET_DEMO_OSEM, ET_DEMO_FISHER_INFORMATION 
%
% 
%Stefano Pedemonte
%Copyright 2009-2012 CMIC-UCL
%Gower Street, London, UK


%% The following computes 'N_noise_instances' noise instances in parallel with the GPU by stacking them
%% in a 3D image. No PSF. 
%% Stacking the 2D images in a 3D volume allowes us to explit the processing power of the GPU in full, 
%% enabling 10000 MLEM iterations of 10240 noise instances in about one day on NVidia GTX285 with images 64x64. 

%% Parameters
N_noise_instances = 2048; 
N_repetitions     = 5;
N                 = 64;
N_cameras         = 100;
N_counts_128      = 150e6; 
N_counts          = N_counts_128* (N/128)^2 / N;  %(1/N) in order to use only one plane - 2D; 
cameras           = zeros(N_cameras,3);
cameras(:,2)      = (pi/180)*(0:180/N_cameras:180-180/N_cameras);
iter_mlem         = 10000;
GPU               = 1;
which_phantom     = 0; %0,1,2

%% Simulate SPECT scan 
fprintf('Creating synthetic sinogram..\n');
mask = et_spherical_phantom(N,N,N,N*0.45,1,0,(N+1)/2,(N+1)/2,(N+1)/2);
if which_phantom==1
    %1)spheres:
    phantom = mask .* (et_spherical_phantom(N,N,N, N/8  ,30 ,0  ,N/2   ,N/2,N/4      ) + ...
                       et_spherical_phantom(N,N,N, N/12 ,20 ,0  ,N/2   ,N/2,1.1*3*N/4) + ...
                       et_spherical_phantom(N,N,N, N/10 ,40 ,0  ,N/2   ,N/2,1.1*N/2  ) + ... 
                       et_spherical_phantom(N,N,N, N/10 ,3  ,15 ,3*N/4 ,N/2,1.0*N/2  ) ) ;
    phantom = phantom.*et_spherical_phantom(N,N,N, N/10 ,0  ,1  ,N/4   ,N/2,1.0*N/2  );
    phantom = phantom(:,N/2,:);
elseif which_phantom==2
    %2)uniform activity:
    phantom = mask;
    phantom = phantom(:,N/2,:);
elseif which_phantom==3
    %3)ncat phantom
    load ncat_phantom.mat
    phantom = mask.*ncat_phantom;
    phantom = phantom(:,64,:); 
end

attenuation = 0;

phantom_stack = zeros(N,N_noise_instances,N);
for i=1:N_noise_instances
    phantom_stack(:,i,:)=phantom; 
end
ideal_sinogram_stack = et_project(phantom_stack, cameras, attenuation, 0, GPU);
ideal_sinogram_stack = ideal_sinogram_stack/sum(ideal_sinogram_stack(:))*N_counts*N_noise_instances;

time_elapsed_sec=0;
for repetition=1:N_repetitions
    fprintf('================= Repetition: %d/%d ================\n',repetition, N_repetitions); tic;
    sinogram_stack = et_poissrnd(ideal_sinogram_stack);
    fprintf(' Noisy sinogram generated in %3.3f seconds. \n',toc); tic; 
    %% Reconstruction:
    activity_stack = ones(N,N_noise_instances,N);
    norm_stack = et_backproject(ones(N,N_noise_instances,N_cameras),cameras,attenuation,0,GPU);
    fprintf(' Normalisation generated in %3.3f seconds. \n',toc); 
    for iter=1:iter_mlem 
        activity_stack = et_mapem_step(activity_stack, norm_stack, sinogram_stack, cameras, attenuation, 0, 0, 0, GPU, 0, 0.0001); 
        time_left_sec = time_elapsed_sec*((N_repetitions-repetition)*iter_mlem+(iter_mlem-iter)); 
        time_left_h = floor(time_left_sec/60/60);
        time_left_m = floor((time_left_sec-60*60*time_left_h)/60); 
        time_left_s = floor(time_left_sec-60*60*time_left_h-60*time_left_m); 
        fprintf(' - Repetition %d/%d, iteration %d/%d, %3.3f seconds for 1 MLEM iteration, %d:%d:%d left\n',repetition, N_repetitions, iter,iter_mlem, time_elapsed_sec, time_left_h,time_left_m,time_left_s); 
        if (floor(iter/100)*100==iter)
            time_elapsed_sec = toc/100; 
            save(sprintf('./out_64_ncat_2/instance%d_iter%d.mat',repetition,iter),'activity_stack'); 
            tic;
        end
    end 
end

if GPU
    et_reset_gpu();
end

